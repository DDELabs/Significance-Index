"""
maslov_sneppen_single_node.py
==============================
Single-node Maslov-Sneppen configuration-model null analysis
for EAS and EUR populations across multiple tissues.

Description
-----------
For each tissue/population pair, this script:
  1. Loads the crosstalk edge list (CSV with columns gene_A, gene_B).
  2. Builds an undirected gene network.
  3. Computes the observed disconnection fraction for a query (hub) gene.
  4. Runs a Maslov-Sneppen degree-preserving rewiring null model
     (N_PERMUTATIONS permutations, Q * |E| edge-swap attempts each).
  5. Computes z-score and empirical p-value; applies BH-FDR across tissues.
  6. Saves per-tissue results and a combined summary table.

Input files (place in working directory)
-----------------------------------------
  crosstalk_edges_<POPULATION>_<tissue>.csv   (columns: gene_A, gene_B)

Outputs (written to results/<POPULATION>/<tissue>/)
----------------------------------------------------
  null_<TISSUE>_<GENE>.npz       -- null distribution + observed value
  stats_<TISSUE>_<GENE>.csv      -- per-tissue statistics
  Fig_<TISSUE>_<GENE>.pdf        -- null distribution figure

Combined output (results/)
--------------------------
  EAS_all_tissues.csv
  EUR_all_tissues.csv

Usage
-----
  python maslov_sneppen_single_node.py

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
# Each entry: (population, tissue_label, edge_csv_file, query_gene)
# Comment-out any row to skip that tissue during a run.

TISSUE_MAP = [
    # EAS population
    ("EAS", "Ovary",     "crosstalk_edges_EAS_ovary.csv",     "TP53"),
    ("EAS", "Adipose",   "crosstalk_edges_EAS_adipose.csv",   "APP"),
    ("EAS", "Liver",     "crosstalk_edges_EAS_liver.csv",     "APP"),
    ("EAS", "Muscle",    "crosstalk_edges_EAS_muscle.csv",    "SMAD2"),
    ("EAS", "Pancreas",  "crosstalk_edges_EAS_pancreas.csv",  "SUOX"),
    ("EAS", "Pituitary", "crosstalk_edges_EAS_pituitary.csv", "GNAS"),
    ("EAS", "Cortex",    "crosstalk_edges_EAS_cortex.csv",    "CDK2"),
]

N_PERMUTATIONS = 1000   # number of rewired null networks
Q              = 10     # rewiring factor (Q * |E| swap attempts)
SEED           = 42     # global random seed
RESULTS_DIR    = "results"

# Optional sanity checks (set True to enable)
RUN_REWIRE_SANITY_CHECK = False
RUN_MIXING_DIAGNOSTICS  = False
# ─────────────────────────────────────────────────────────────────────────────


# ── DIRECTORY HELPERS ─────────────────────────────────────────────────────────

def tissue_dir(population, tissue):
    path = os.path.join(RESULTS_DIR, population, tissue)
    os.makedirs(path, exist_ok=True)
    return path


# ── NETWORK ───────────────────────────────────────────────────────────────────

def load_network(edge_file):
    """Load undirected gene network from a CSV edge list."""
    edge_df = pd.read_csv(edge_file)
    G = nx.from_pandas_edgelist(edge_df, source="gene_A", target="gene_B")
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


# ── TOPOLOGY ──────────────────────────────────────────────────────────────────

def compute_topology(G):
    """Return a DataFrame of degree, betweenness, and clustering per node."""
    degree_dict      = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G, normalized=True)
    clustering_dict  = nx.clustering(G)

    node_profile = pd.DataFrame({
        "gene"        : list(degree_dict.keys()),
        "degree"      : list(degree_dict.values()),
        "betweenness" : [betweenness_dict[n] for n in degree_dict],
        "clustering"  : [clustering_dict[n]  for n in degree_dict],
    }).sort_values("degree", ascending=False).reset_index(drop=True)

    node_profile["degree_rank"]      = node_profile["degree"].rank(ascending=False).astype(int)
    node_profile["betweenness_rank"] = node_profile["betweenness"].rank(ascending=False).astype(int)
    return node_profile


# ── DISCONNECTION FRACTION ────────────────────────────────────────────────────

def disconnection_fraction(G, node_to_remove):
    """
    Fraction of remaining nodes NOT in the largest connected component
    after removing *node_to_remove*.
    """
    if node_to_remove not in G:
        return np.nan
    H = G.copy()
    H.remove_node(node_to_remove)
    n_remaining = H.number_of_nodes()
    if n_remaining == 0:
        return 1.0
    largest = len(max(nx.connected_components(H), key=len))
    return (n_remaining - largest) / n_remaining


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

    G_rewired    = G.copy()
    n_edges      = G_rewired.number_of_edges()
    n_swaps      = Q * n_edges
    n_successful = 0
    n_attempts   = 0

    while n_successful < n_swaps:
        n_attempts += 1
        if n_attempts > n_swaps * 20:
            print(f"    [Rewire] Warning: only {n_successful}/{n_swaps} swaps completed")
            break

        edges = list(G_rewired.edges())
        if len(edges) < 2:
            break

        idx1, idx2 = random.sample(range(len(edges)), 2)
        u1, v1 = edges[idx1]
        u2, v2 = edges[idx2]

        if len({u1, v1, u2, v2}) < 4:
            continue
        if G_rewired.has_edge(u1, v2) or G_rewired.has_edge(u2, v1):
            continue
        if u1 == v2 or u2 == v1:
            continue

        G_rewired.remove_edge(u1, v1)
        G_rewired.remove_edge(u2, v2)
        G_rewired.add_edge(u1, v2)
        G_rewired.add_edge(u2, v1)
        n_successful += 1

    return G_rewired, n_successful


# ── NULL MODEL ────────────────────────────────────────────────────────────────

def run_null_model(G, query_gene, n_permutations=1000, Q=10, seed=42):
    """
    Run the Maslov-Sneppen null model for a single query gene.

    Returns an array of disconnection fractions from *n_permutations*
    independently rewired networks.
    """
    null_vals = []
    np.random.seed(seed)
    seeds = np.random.randint(0, 100_000, size=n_permutations)

    for i in range(n_permutations):
        if (i + 1) % 200 == 0:
            print(f"    Permutation {i+1}/{n_permutations}")
        G_rw, _ = maslov_sneppen_rewire(G, Q=Q, seed=int(seeds[i]))
        null_vals.append(disconnection_fraction(G_rw, query_gene))

    return np.array(null_vals)


# ── SAVE / LOAD NULL RESULTS ──────────────────────────────────────────────────

def save_null(population, tissue, gene, null_array, f_obs):
    """Persist null array and observed value as a compressed .npz file."""
    path = os.path.join(tissue_dir(population, tissue),
                        f"null_{tissue}_{gene}.npz")
    np.savez_compressed(path, null_array=null_array, f_obs=np.array([f_obs]))
    print(f"    [Saved] {path}")
    return path


def load_null(population, tissue, gene):
    """Load a previously saved null distribution."""
    path = os.path.join(RESULTS_DIR, population, tissue,
                        f"null_{tissue}_{gene}.npz")
    data = np.load(path)
    return data["null_array"], float(data["f_obs"][0])


# ── STATISTICS ────────────────────────────────────────────────────────────────

def compute_statistics(gene, f_obs, null_array):
    """
    Compute z-score, empirical p-value, and null summary statistics.
    """
    mu      = np.mean(null_array)
    sigma   = np.std(null_array, ddof=1)
    z_score = (f_obs - mu) / sigma if sigma > 0 else (np.inf if f_obs > mu else 0.0)
    p_emp   = float(np.mean(null_array >= f_obs))
    p95     = np.percentile(null_array, 95)
    return dict(
        gene=gene, f_obs=f_obs, mu_null=mu, sigma_null=sigma,
        z_score=z_score, p_emp=p_emp, p95_null=p95,
        n_perm=len(null_array),
    )


# ── PLOTTING ──────────────────────────────────────────────────────────────────

def plot_null_distribution(result, null_array, population, tissue, out_dir):
    """
    Plot the null distribution of disconnection fractions with the
    observed value annotated. Saves a PDF to *out_dir*.
    """
    gene  = result["gene"]
    f_obs = result["f_obs"]
    z, p  = result["z_score"], result["p_emp"]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(null_array * 100, bins=30, color="#4a90d9", alpha=0.72,
            edgecolor="white", linewidth=0.5,
            label=f"Null (n={len(null_array)})")
    ax.axvline(f_obs * 100, color="#d62728", lw=2.0, ls="--",
               label=f"Observed: {f_obs*100:.1f}%")
    ax.axvline(result["mu_null"] * 100, color="#2ca02c", lw=1.5, ls=":",
               label=f"Null mean: {result['mu_null']*100:.1f}%")
    ax.axvline(result["p95_null"] * 100, color="#9467bd", lw=1.2,
               ls="-.", alpha=0.85,
               label=f"Null 95th: {result['p95_null']*100:.1f}%")

    sig_color = "#d62728" if p < 0.05 else "#555555"
    ax.text(
        0.97, 0.97, f"Z = {z:.2f}\np_emp = {p:.4f}",
        transform=ax.transAxes, va="top", ha="right", fontsize=10,
        color=sig_color,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  alpha=0.85, edgecolor=sig_color, lw=0.9),
    )

    ax.set_title(
        f"{gene}  [{population} — {tissue}]\n"
        f"Maslov-Sneppen null  (n={len(null_array)} permutations, Q={Q})",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlabel("Nodes disconnected upon removal (%)", fontsize=10)
    ax.set_ylabel("Count (permutations)", fontsize=10)
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()

    out_fig = os.path.join(out_dir, f"Fig_{tissue}_{gene}.pdf")
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    [Plot]  {out_fig}")


# ── REPLOT FROM SAVED DATA ────────────────────────────────────────────────────

def replot(population, tissue, gene):
    """
    Regenerate a figure from a saved .npz file without re-running the null model.

    Example
    -------
    >>> from maslov_sneppen_single_node import replot
    >>> replot("EAS", "Ovary", "TP53")
    """
    null_array, f_obs = load_null(population, tissue, gene)
    result = compute_statistics(gene, f_obs, null_array)
    out = tissue_dir(population, tissue)
    plot_null_distribution(result, null_array, population, tissue, out)
    print(f"Replot complete → {out}/Fig_{tissue}_{gene}.pdf")


# ── OPTIONAL DIAGNOSTICS ──────────────────────────────────────────────────────

def rewire_sanity_check(G, tissue):
    """Assert degree sequence is preserved and no self-loops are introduced."""
    G_test, _ = maslov_sneppen_rewire(G, Q=Q, seed=SEED)
    deg_orig    = sorted(dict(G.degree()).values(), reverse=True)
    deg_rewired = sorted(dict(G_test.degree()).values(), reverse=True)
    assert deg_orig == deg_rewired, f"[{tissue}] Degree sequence changed!"
    assert not list(nx.selfloop_edges(G_test)), f"[{tissue}] Self-loops present!"
    shared = set(G.edges()) & set(G_test.edges())
    print(
        f"    [Sanity] Degree preserved ✓ | No self-loops ✓ | "
        f"Turnover: {1 - len(shared)/G.number_of_edges():.1%}"
    )


def assess_mixing(G, tissue, Q_values=(10, 50, 100, 200), n_trials=10):
    """Assess edge turnover across different Q values to verify mixing."""
    print(f"    [Mixing] {tissue}")
    n_edges = G.number_of_edges()
    for Qv in Q_values:
        turnovers = [
            1 - len(set(G.edges()) &
                    set(maslov_sneppen_rewire(G, Q=Qv, seed=t)[0].edges())) / n_edges
            for t in range(n_trials)
        ]
        print(
            f"      Q={Qv:4d} | "
            f"turnover = {np.mean(turnovers):.1%} ± {np.std(turnovers):.1%}"
        )


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Collect results grouped by population for the summary tables
    all_stats_by_pop = {}

    for population, tissue, edge_file, query_gene in TISSUE_MAP:
        print(f"\n{'='*60}")
        print(f"  {population} | {tissue} | gene = {query_gene}")
        print(f"{'='*60}")

        if not os.path.exists(edge_file):
            print(f"  [SKIP] File not found: {edge_file}")
            continue

        G = load_network(edge_file)
        print(
            f"  Nodes: {G.number_of_nodes()} | "
            f"Edges: {G.number_of_edges()} | "
            f"Components: {nx.number_connected_components(G)}"
        )

        if query_gene not in G.nodes():
            print(f"  [SKIP] {query_gene} not in network — check node names.")
            continue

        if RUN_REWIRE_SANITY_CHECK:
            rewire_sanity_check(G, tissue)
        if RUN_MIXING_DIAGNOSTICS:
            assess_mixing(G, tissue)

        node_profile = compute_topology(G)

        f_obs = disconnection_fraction(G, query_gene)
        print(f"  f_obs({query_gene}) = {f_obs:.4f}")

        print(f"  Running {N_PERMUTATIONS} permutations...")
        null_array = run_null_model(G, query_gene, N_PERMUTATIONS, Q, SEED)

        out = tissue_dir(population, tissue)
        save_null(population, tissue, query_gene, null_array, f_obs)

        result = compute_statistics(query_gene, f_obs, null_array)
        result["population"] = population
        result["tissue"]     = tissue

        gene_row = node_profile[node_profile["gene"] == query_gene]
        if not gene_row.empty:
            result["degree"]           = int(gene_row["degree"].values[0])
            result["degree_rank"]      = int(gene_row["degree_rank"].values[0])
            result["betweenness"]      = float(gene_row["betweenness"].values[0])
            result["betweenness_rank"] = int(gene_row["betweenness_rank"].values[0])

        stats_df  = pd.DataFrame([result])
        stats_csv = os.path.join(out, f"stats_{tissue}_{query_gene}.csv")
        stats_df.to_csv(stats_csv, index=False)
        print(f"    [Saved] {stats_csv}")

        plot_null_distribution(result, null_array, population, tissue, out)

        all_stats_by_pop.setdefault(population, []).append(result)

    # ── combined summary tables (one per population) ───────────────────────────
    for population, stats_list in all_stats_by_pop.items():
        if not stats_list:
            continue
        summary = pd.DataFrame(stats_list)

        _, q_vals, _, _ = multipletests(summary["p_emp"], alpha=0.05, method="fdr_bh")
        summary["q_BH"]        = q_vals
        summary["significant"] = summary["q_BH"] < 0.05

        cols = [
            "population", "tissue", "gene",
            "degree", "degree_rank",
            "f_obs", "mu_null", "sigma_null",
            "z_score", "p_emp", "q_BH", "significant",
        ]
        display = summary[[c for c in cols if c in summary.columns]].copy()
        for c in ["f_obs", "mu_null", "sigma_null"]:
            if c in display.columns:
                display[c] = (display[c] * 100).round(2)
        display["z_score"] = display["z_score"].round(3)
        display["p_emp"]   = display["p_emp"].round(4)
        display["q_BH"]    = display["q_BH"].round(4)

        out_table = os.path.join(RESULTS_DIR,
                                 f"TableS6_{population}_all_tissues.csv")
        display.to_csv(out_table, index=False)

        print(f"\n{'='*60}")
        print(f"COMBINED RESULTS — {population} all tissues")
        print("=" * 60)
        print(display.to_string(index=False))
        print(f"\n[Saved] {out_table}")

    print("\n✓ All tissues complete.")
    print("\nTo replot any result without re-running:")
    print("  from maslov_sneppen_single_node import replot")
    print('  replot("EAS", "Ovary", "TP53")')


if __name__ == "__main__":
    main()
