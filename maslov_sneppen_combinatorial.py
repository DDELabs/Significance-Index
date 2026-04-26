"""
maslov_sneppen_combinatorial.py
================================
Combinatorial node-removal Maslov-Sneppen null model with synergy analysis
for EAS and EUR populations (ovary tissue; extensible to all tissues).

Manuscript
----------
Title : Cross-population gene network topology reveals ancestry-specific
        hub genes in metabolic and neurological tissues
Preprint: https://doi.org/10.21203/rs.3.rs-8610143/v1

Description
-----------
Beyond single-gene hub analysis, this script tests whether *combinations*
of hub genes exert a greater-than-additive (super-additive, i.e. synergistic)
effect on network fragmentation.

For each population/tissue/combination triplet, the script:
  1. Loads the crosstalk edge list and builds the undirected gene network.
  2. Single-node null  — runs the Maslov-Sneppen null for every gene that
     appears in at least one combination (needed to compute synergy).
  3. Combinatorial null — measures the disconnection fraction after
     simultaneously removing all genes in each combination across
     N_PERMUTATIONS rewired networks.
  4. Synergy null — computes
         synergy = f_combo − Σ f_individual
     in both the observed network and every rewired network, yielding
     a synergy z-score and empirical p-value.
  5. Saves every null distribution as a compressed .npz file so that
     figures can be regenerated without re-running.

Gene combinations tested
------------------------
Population  Tissue   Combination
----------  -------  --------------------------------
EAS         Ovary    (TP53, GNAS)
EAS         Ovary    (YWHAE, TP53, GNAS)
EAS         Ovary    (YWHAE, TP53, RUNX1, CTNNB1)

EUR         Ovary    (CDK9, CDK1, ELAVL1)
EUR         Ovary    (XRCC1, NEIL2, CDK9)

Add further tissues by extending POPULATION_TISSUE_MAP below.

Input files (place in working directory)
-----------------------------------------
  crosstalk_edges_EAS_ovary.csv   (columns: gene_A, gene_B)
  crosstalk_edges_EUR_ovary.csv

Outputs (written to results/<POPULATION>/<TISSUE>/)
----------------------------------------------------
  null_combo_<TISSUE>_<GENES>.npz     -- positional null array + f_obs
  null_synergy_<TISSUE>_<GENES>.npz   -- synergy null array + observed synergy
  null_single_<TISSUE>_<GENE>.npz     -- per-gene null (for synergy)
  positional_<POPULATION>_<TISSUE>.csv
  synergy_<POPULATION>_<TISSUE>.csv
  Fig_combinatorial_<POPULATION>_<TISSUE>.pdf

Usage
-----
  python maslov_sneppen_combinatorial.py

  To replot from saved .npz files without re-running:
      from maslov_sneppen_combinatorial import replot_from_npz
      replot_from_npz("EAS", "Ovary")

Dependencies
------------
  networkx >= 2.8
  numpy >= 1.22
  pandas >= 1.4
  matplotlib >= 3.5
  statsmodels >= 0.13

Reproducibility
---------------
  SEED = 42   (controls all random number generation)

License
-------
  MIT License. See repository LICENSE file.
"""

import os
import random
import warnings

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
N_PERMUTATIONS = 1000   # number of rewired null networks
Q              = 10     # rewiring factor (Q * |E| swap attempts)
SEED           = 42     # global random seed
RESULTS_DIR    = "results"

# Map: population → tissue → {edge_file, combinations}
# Each combination is a list of gene names to remove simultaneously.
POPULATION_TISSUE_MAP = {
    "EAS": {
        "Ovary": {
            "edge_file": "crosstalk_edges_EAS_ovary.csv",
            "combinations": [
                ("TP53", "GNAS"),
                ("YWHAE", "TP53", "GNAS"),
                ("YWHAE", "TP53", "RUNX1", "CTNNB1"),
            ],
        },
        # Add further EAS tissues here, e.g.:
        # "Adipose": {
        #     "edge_file": "crosstalk_edges_EAS_adipose.csv",
        #     "combinations": [("APP", "SMAD2")],
        # },
    },
    "EUR": {
        "Ovary": {
            "edge_file": "crosstalk_edges_EUR_ovary.csv",
            "combinations": [
                ("CDK9", "CDK1", "ELAVL1"),
                ("XRCC1", "NEIL2", "CDK9"),
            ],
        },
    },
}
# ─────────────────────────────────────────────────────────────────────────────


# ── DIRECTORY HELPERS ─────────────────────────────────────────────────────────

def combo_tag(combo):
    """Return a filename-safe string for a gene-combination tuple."""
    return "+".join(combo)


def out_dir(population, tissue):
    path = os.path.join(RESULTS_DIR, population, tissue)
    os.makedirs(path, exist_ok=True)
    return path


def combo_npz_path(population, tissue, combo, kind="combo"):
    """
    Return the .npz file path for a given combo and result kind.

    Parameters
    ----------
    kind : str
        'combo'   → positional null
        'synergy' → synergy null
        'single'  → individual-gene null
    """
    tag = combo_tag(combo) if isinstance(combo, tuple) else combo
    fname = f"null_{kind}_{tissue}_{tag}.npz"
    return os.path.join(out_dir(population, tissue), fname)


# ── NETWORK ───────────────────────────────────────────────────────────────────

def load_network(edge_file):
    """Load undirected gene network from a CSV edge list."""
    df = pd.read_csv(edge_file)
    G  = nx.from_pandas_edgelist(df, source="gene_A", target="gene_B")
    G.remove_edges_from(nx.selfloop_edges(G))
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ── DISCONNECTION FRACTIONS ───────────────────────────────────────────────────

def disconnection_fraction(G, node):
    """
    Fraction of remaining nodes NOT in the largest connected component
    after removing *node*.
    """
    if node not in G:
        return np.nan
    H = G.copy()
    H.remove_node(node)
    n = H.number_of_nodes()
    if n == 0:
        return 1.0
    return (n - len(max(nx.connected_components(H), key=len))) / n


def disconnection_fraction_combo(G, nodes):
    """
    Fraction of remaining nodes NOT in the largest connected component
    after simultaneously removing all genes in *nodes*.
    """
    H = G.copy()
    for n in nodes:
        if n in H:
            H.remove_node(n)
    n = H.number_of_nodes()
    if n == 0:
        return 1.0
    return (n - len(max(nx.connected_components(H), key=len))) / n


# ── MASLOV-SNEPPEN REWIRING ────────────────────────────────────────────────────

def maslov_sneppen_rewire(G, Q=10, seed=None):
    """
    Degree-preserving random rewiring (Maslov & Sneppen, Science 2002).

    Performs Q * |E| double-edge swaps while avoiding self-loops and
    multi-edges.  Returns (rewired_graph, n_successful_swaps).
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    G2      = G.copy()
    n_swaps = Q * G2.number_of_edges()
    done    = 0
    attempts = 0

    while done < n_swaps:
        attempts += 1
        if attempts > n_swaps * 20:
            print(f"    [Rewire] Warning: only {done}/{n_swaps} swaps completed")
            break
        edges = list(G2.edges())
        if len(edges) < 2:
            break
        i1, i2 = random.sample(range(len(edges)), 2)
        u1, v1 = edges[i1]
        u2, v2 = edges[i2]
        if len({u1, v1, u2, v2}) < 4:
            continue
        if G2.has_edge(u1, v2) or G2.has_edge(u2, v1):
            continue
        if u1 == v2 or u2 == v1:
            continue
        G2.remove_edge(u1, v1)
        G2.remove_edge(u2, v2)
        G2.add_edge(u1, v2)
        G2.add_edge(u2, v1)
        done += 1

    return G2, done


# ── NULL RUNS ─────────────────────────────────────────────────────────────────

def run_single_node_null(G, genes, n_permutations=1000, Q=10, seed=42):
    """
    Run the Maslov-Sneppen null for each gene in *genes* individually.

    Returns a dict {gene: np.array of null disconnection fractions}.
    """
    np.random.seed(seed)
    seeds = np.random.randint(0, 100_000, size=n_permutations)
    null_single = {}
    for gene in genes:
        if gene not in G:
            print(f"  WARNING: {gene} not in network")
            continue
        vals = []
        for i in range(n_permutations):
            Gr, _ = maslov_sneppen_rewire(G, Q=Q, seed=int(seeds[i]))
            vals.append(disconnection_fraction(Gr, gene))
        null_single[gene] = np.array(vals)
        print(f"  ✓ single null: {gene}")
    return null_single


def run_combinatorial_null(G, combos, n_permutations=1000, Q=10, seed=42):
    """
    Run the Maslov-Sneppen null for each combination in *combos*.

    Returns a dict {combo_tuple: np.array of null disconnection fractions}.
    """
    np.random.seed(seed)
    seeds      = np.random.randint(0, 100_000, size=n_permutations)
    null_combo = {c: [] for c in combos}

    for i in range(n_permutations):
        if (i + 1) % 200 == 0:
            print(f"  Permutation {i+1}/{n_permutations}")
        Gr, _ = maslov_sneppen_rewire(G, Q=Q, seed=int(seeds[i]))
        for combo in combos:
            null_combo[combo].append(disconnection_fraction_combo(Gr, list(combo)))

    return {k: np.array(v) for k, v in null_combo.items()}


# ── SAVE / LOAD NPZ ───────────────────────────────────────────────────────────

def save_combo_null(population, tissue, combo, null_array, f_obs):
    """Save positional null array for a gene combination."""
    path = combo_npz_path(population, tissue, combo, kind="combo")
    np.savez_compressed(path, null_array=null_array, f_obs=np.array([f_obs]))
    print(f"  [Saved positional npz]  {path}")


def save_synergy_null(population, tissue, combo, synergy_null, synergy_obs):
    """Save synergy null array for a gene combination."""
    path = combo_npz_path(population, tissue, combo, kind="synergy")
    np.savez_compressed(path, null_array=synergy_null,
                        synergy_obs=np.array([synergy_obs]))
    print(f"  [Saved synergy npz]     {path}")


def save_single_null(population, tissue, gene, null_array, f_obs_single):
    """Save individual-gene null array (required to reconstruct synergy)."""
    path = combo_npz_path(population, tissue, gene, kind="single")
    np.savez_compressed(path, null_array=null_array,
                        f_obs=np.array([f_obs_single]))
    print(f"  [Saved single npz]      {path}")


def load_combo_null(population, tissue, combo):
    path = combo_npz_path(population, tissue, combo, kind="combo")
    d    = np.load(path)
    return d["null_array"], float(d["f_obs"][0])


def load_synergy_null(population, tissue, combo):
    path = combo_npz_path(population, tissue, combo, kind="synergy")
    d    = np.load(path)
    return d["null_array"], float(d["synergy_obs"][0])


def load_single_null(population, tissue, gene):
    path = combo_npz_path(population, tissue, gene, kind="single")
    d    = np.load(path)
    return d["null_array"], float(d["f_obs"][0])


# ── STATISTICS ────────────────────────────────────────────────────────────────

def compute_combo_stats(G, combo, null_combo_dist):
    """
    Positional statistics for a gene combination.

    Returns a dict with f_obs, mu_null, sigma_null, z_score, p_emp, etc.
    """
    f_obs = disconnection_fraction_combo(G, list(combo))
    mu    = np.mean(null_combo_dist)
    sigma = np.std(null_combo_dist, ddof=1)
    z     = (f_obs - mu) / sigma if sigma > 0 else np.inf
    p     = float(np.mean(null_combo_dist >= f_obs))
    p95   = np.percentile(null_combo_dist, 95)
    return dict(
        combo=combo_tag(combo), combo_size=len(combo),
        f_obs=f_obs, mu_null=mu, sigma_null=sigma,
        z_score=z, p_emp=p, p95_null=p95, n_perm=len(null_combo_dist),
    )


def compute_synergy(G, combo, f_obs_single, null_single, null_combo_dist):
    """
    Synergy statistics for a gene combination.

    synergy_obs = f_obs_combo − Σ f_obs_individual

    The null synergy distribution is derived from the same permutations
    used for the combinatorial null, making the test fully non-parametric.

    Returns a dict including synergy_obs, z_synergy, p_synergy,
    super_additive flag, and the full synergy_null_dist array.
    """
    f_obs_c           = disconnection_fraction_combo(G, list(combo))
    indiv             = np.array([f_obs_single[g] for g in combo])
    syn_obs           = f_obs_c - np.sum(indiv)
    null_indiv_arrays = np.array([null_single[g] for g in combo])
    syn_null          = null_combo_dist - np.sum(null_indiv_arrays, axis=0)
    mu_syn  = np.mean(syn_null)
    sig_syn = np.std(syn_null, ddof=1)
    z_syn   = (syn_obs - mu_syn) / sig_syn if sig_syn > 0 else np.inf
    p_syn   = float(np.mean(syn_null >= syn_obs))
    return dict(
        combo=combo_tag(combo),
        f_obs_combo=f_obs_c, f_obs_sum_indiv=float(np.sum(indiv)),
        synergy_obs=syn_obs, mu_synergy_null=mu_syn,
        sigma_synergy=sig_syn, z_synergy=z_syn, p_synergy=p_syn,
        super_additive=bool(syn_obs > 0),
        synergy_null_dist=syn_null,
    )


# ── PLOTTING ──────────────────────────────────────────────────────────────────

def plot_combo(population, tissue, combos, positional_results,
               synergy_results, null_combo, null_synergy_dists):
    """
    For each gene combination, produce a two-panel figure:
      Left  : positional null distribution (combinatorial disconnection).
      Right : synergy null distribution (super- vs sub-additivity).

    Saves a single multi-panel PDF.
    """
    n_rows = len(combos)
    fig, axes = plt.subplots(n_rows, 2, figsize=(14, 4.5 * n_rows))
    if n_rows == 1:
        axes = [axes]

    fig.suptitle(
        f"Combinatorial node removal null model\n"
        f"Maslov-Sneppen ({N_PERMUTATIONS} permutations, Q={Q})"
        f" · {population} · {tissue}",
        fontsize=13, fontweight="bold",
    )

    for row, combo in enumerate(combos):
        tag         = combo_tag(combo)
        pr          = positional_results[combo]
        sr          = synergy_results[combo]
        nc          = null_combo[combo] * 100
        ns          = null_synergy_dists[combo] * 100
        f_obs_pct   = pr["f_obs"] * 100
        mu_pct      = pr["mu_null"] * 100
        syn_obs_pct = sr["synergy_obs"] * 100

        # Left: positional null
        ax = axes[row][0]
        ax.hist(nc, bins=30, color="#4a90d9", alpha=0.72,
                edgecolor="white", lw=0.5)
        ax.axvline(f_obs_pct, color="#d62728", lw=2, ls="--",
                   label=f"Observed: {f_obs_pct:.1f}%")
        ax.axvline(mu_pct, color="#2ca02c", lw=1.5, ls=":",
                   label=f"Null mean: {mu_pct:.1f}%")
        sig_col = "#d62728" if pr["p_emp"] < 0.05 else "#555"
        ax.text(
            0.97, 0.97,
            f"Z = {pr['z_score']:.2f}\np = {pr['p_emp']:.4f}",
            transform=ax.transAxes, va="top", ha="right", fontsize=10,
            color=sig_col,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.85, edgecolor=sig_col, lw=0.9),
        )
        ax.set_title(
            f"Positional null: {tag}\n"
            f"f_obs={f_obs_pct:.1f}% vs null_mean={mu_pct:.1f}%",
            fontsize=10,
        )
        ax.set_xlabel("Disconnection fraction (%)")
        ax.set_ylabel("Permutation count")
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        # Right: synergy null
        ax = axes[row][1]
        ax.hist(ns, bins=30, color="#e69c24", alpha=0.80,
                edgecolor="white", lw=0.5)
        ax.axvline(syn_obs_pct, color="#d62728", lw=2, ls="--",
                   label=f"Observed synergy: {syn_obs_pct:.1f}%")
        ax.axvline(0, color="gray", lw=1.2, ls="-",
                   label="Zero synergy (additive)")
        label    = "SUPER-additive" if sr["super_additive"] else "SUB-additive"
        sig_col2 = "#d62728" if sr["p_synergy"] < 0.05 else "#555"
        ax.text(
            0.97, 0.97,
            f"{label}\nZ = {sr['z_synergy']:.2f}\np = {sr['p_synergy']:.3f}",
            transform=ax.transAxes, va="top", ha="right", fontsize=10,
            color=sig_col2,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.85, edgecolor=sig_col2, lw=0.9),
        )
        f_sum_pct = sr["f_obs_sum_indiv"] * 100
        ax.set_title(
            f"Synergy null: {tag}\n"
            f"f_obs_combo={f_obs_pct:.1f}% vs sum_indiv={f_sum_pct:.1f}%",
            fontsize=10,
        )
        ax.set_xlabel("Synergy = f_combo − Σf_individual (%)")
        ax.set_ylabel("Permutation count")
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    fig_path = os.path.join(
        out_dir(population, tissue),
        f"Fig_combinatorial_{population}_{tissue}.pdf",
    )
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [Plot saved] {fig_path}")


# ── REPLOT FROM SAVED NPZ ─────────────────────────────────────────────────────

def replot_from_npz(population, tissue):
    """
    Regenerate the figure from saved .npz files without re-running the null model.

    Requires that main() was run at least once with this script so that
    all .npz files exist in results/<population>/<tissue>/.

    Example
    -------
    >>> from maslov_sneppen_combinatorial import replot_from_npz
    >>> replot_from_npz("EAS", "Ovary")
    """
    combos = [
        tuple(c)
        for c in POPULATION_TISSUE_MAP[population][tissue]["combinations"]
    ]
    edge_file = POPULATION_TISSUE_MAP[population][tissue]["edge_file"]
    G = load_network(edge_file)

    all_genes    = sorted(set(g for c in combos for g in c))
    null_single  = {}
    f_obs_single = {}
    for gene in all_genes:
        arr, fobs         = load_single_null(population, tissue, gene)
        null_single[gene] = arr
        f_obs_single[gene] = fobs

    null_combo         = {}
    null_synergy_dists = {}
    positional_results = {}
    synergy_results    = {}

    for combo in combos:
        arr, fobs          = load_combo_null(population, tissue, combo)
        null_combo[combo]  = arr
        syn_arr, syn_obs   = load_synergy_null(population, tissue, combo)
        null_synergy_dists[combo] = syn_arr

        mu = np.mean(arr); sigma = np.std(arr, ddof=1)
        positional_results[combo] = dict(
            f_obs=fobs, mu_null=mu, sigma_null=sigma,
            z_score=(fobs - mu) / sigma if sigma > 0 else np.inf,
            p_emp=float(np.mean(arr >= fobs)),
        )
        f_sum = sum(f_obs_single[g] for g in combo)
        mu_s  = np.mean(syn_arr); sig_s = np.std(syn_arr, ddof=1)
        synergy_results[combo] = dict(
            synergy_obs=syn_obs, f_obs_combo=fobs, f_obs_sum_indiv=f_sum,
            z_synergy=(syn_obs - mu_s) / sig_s if sig_s > 0 else np.inf,
            p_synergy=float(np.mean(syn_arr >= syn_obs)),
            super_additive=bool(syn_obs > 0),
        )

    plot_combo(population, tissue, combos, positional_results,
               synergy_results, null_combo, null_synergy_dists)
    print(
        f"Replot done → "
        f"{out_dir(population, tissue)}/"
        f"Fig_combinatorial_{population}_{tissue}.pdf"
    )


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    for population, tissues in POPULATION_TISSUE_MAP.items():
        for tissue, cfg in tissues.items():
            edge_file = cfg["edge_file"]
            combos    = [tuple(c) for c in cfg["combinations"]]

            print(f"\n{'='*60}")
            print(f"  {population} | {tissue}")
            print(f"{'='*60}")

            if not os.path.exists(edge_file):
                print(f"  [SKIP] {edge_file} not found")
                continue

            G         = load_network(edge_file)
            all_genes = sorted(set(g for c in combos for g in c))

            # ── 1. Single-node nulls (needed for synergy) ─────────────────
            print("\n  [1/3] Single-node nulls...")
            null_single  = run_single_node_null(G, all_genes, N_PERMUTATIONS, Q, SEED)
            f_obs_single = {g: disconnection_fraction(G, g) for g in all_genes}
            for gene in all_genes:
                if gene in null_single:
                    save_single_null(population, tissue, gene,
                                     null_single[gene], f_obs_single[gene])

            # ── 2. Combinatorial null ─────────────────────────────────────
            print("\n  [2/3] Combinatorial null...")
            null_combo = run_combinatorial_null(G, combos, N_PERMUTATIONS, Q, SEED)

            positional_results = {}
            synergy_results    = {}
            null_synergy_dists = {}

            for combo in combos:
                f_obs_c = disconnection_fraction_combo(G, list(combo))
                save_combo_null(population, tissue, combo,
                                null_combo[combo], f_obs_c)

                null_indiv = np.array([null_single[g] for g in combo])
                syn_null   = null_combo[combo] - np.sum(null_indiv, axis=0)
                indiv_sum  = sum(f_obs_single[g] for g in combo)
                syn_obs    = f_obs_c - indiv_sum
                save_synergy_null(population, tissue, combo, syn_null, syn_obs)

                null_synergy_dists[combo] = syn_null
                positional_results[combo] = compute_combo_stats(
                    G, combo, null_combo[combo]
                )
                synergy_results[combo] = compute_synergy(
                    G, combo, f_obs_single, null_single, null_combo[combo]
                )

            # ── 3. Save CSVs and plot ─────────────────────────────────────
            print("\n  [3/3] Saving CSVs and figure...")
            pos_rows = list(positional_results.values())
            for r in pos_rows:
                r.update(population=population, tissue=tissue)
            syn_rows = [
                {k: v for k, v in synergy_results[c].items()
                 if k != "synergy_null_dist"}
                for c in combos
            ]

            pos_df = pd.DataFrame(pos_rows)
            syn_df = pd.DataFrame(syn_rows)

            if len(combos) > 1:
                pos_df["q_BH"] = multipletests(
                    pos_df["p_emp"], method="fdr_bh"
                )[1]
                syn_df["q_BH_synergy"] = multipletests(
                    syn_df["p_synergy"], method="fdr_bh"
                )[1]

            d = out_dir(population, tissue)
            pos_df.to_csv(
                os.path.join(d, f"positional_{population}_{tissue}.csv"),
                index=False,
            )
            syn_df.to_csv(
                os.path.join(d, f"synergy_{population}_{tissue}.csv"),
                index=False,
            )
            plot_combo(
                population, tissue, combos,
                positional_results, synergy_results,
                null_combo, null_synergy_dists,
            )

    print("\nDone.")
    print("\nTo replot any tissue without re-running the null model:")
    print("  from maslov_sneppen_combinatorial import replot_from_npz")
    print('  replot_from_npz("EAS", "Ovary")')


if __name__ == "__main__":
    main()
