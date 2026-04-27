"""
maslov_sneppen_combinatorial.py
================================
Combinatorial node-removal Maslov-Sneppen null model with synergy analysis.

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

Input files (place in working directory)
-----------------------------------------
  crosstalk_edges_<POPULATION>_<tissue>.csv   (columns: gene_A, gene_B)

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
import networkx as nx
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from statsmodels.stats.multitest import multipletests
from scipy import stats
import warnings

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

N_PERMUTATIONS = 1000
Q = 10
SEED = 42
RESULTS_DIR = "combinatorial_results_liver"

# Population-specific combinations
COMBINATIONS = {
    'EAS': {
        'liver': []
    },
    'EUR': {
        'liver': [
            ('CD74', 'HLA-B'),
            ('CD74', 'APP', 'HLA-B'),
            ('CD74', 'APP', 'PA2G4', 'HLA-B'),
        ]
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# NETWORK UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def load_network(edge_file):
    """Load network from CSV file with gene_A, gene_B columns."""
    edge_df = pd.read_csv(edge_file)
    G = nx.from_pandas_edgelist(edge_df, source="gene_A", target="gene_B")
    G.remove_edges_from(nx.selfloop_edges(G))
    print(f"  Loaded network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def disconnection_fraction(G, node_to_remove):
    """
    Fraction of nodes disconnected from largest component after removing one node.
    
    Parameters
    ----------
    G : nx.Graph
    node_to_remove : str
        Node label
    
    Returns
    -------
    float : disconnected fraction [0, 1]
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


def disconnection_fraction_combo(G, nodes_to_remove):
    """
    Remove all nodes in nodes_to_remove simultaneously from G.
    Return fraction of remaining nodes NOT in largest connected component.
    
    Parameters
    ----------
    G : nx.Graph
    nodes_to_remove : list of node labels
    
    Returns
    -------
    float : disconnected fraction [0, 1]
    """
    # Only remove nodes that actually exist
    valid = [n for n in nodes_to_remove if n in G]
    if not valid:
        return 0.0
    
    H = G.copy()
    H.remove_nodes_from(valid)
    n_remaining = H.number_of_nodes()
    if n_remaining == 0:
        return 1.0
    
    components = sorted(nx.connected_components(H), key=len, reverse=True)
    largest = len(components[0])
    return (n_remaining - largest) / n_remaining

# ─────────────────────────────────────────────────────────────────────────────
# REWIRING (MASLOV-SNEPPEN)
# ─────────────────────────────────────────────────────────────────────────────

def maslov_sneppen_rewire(G, Q=10, seed=None):
    """
    Maslov-Sneppen edge-swapping preserves degree distribution.
    
    Parameters
    ----------
    G : nx.Graph
    Q : int
        Number of swaps = Q × number of edges
    seed : int
        Random seed
    
    Returns
    -------
    G_rewired : nx.Graph
    n_successful : int
        Number of successful swaps
    """
    import random
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    G_rewired = G.copy()
    n_edges = G_rewired.number_of_edges()
    n_swaps = Q * n_edges
    n_successful = 0
    n_attempts = 0
    
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
        
        # Check constraints
        if len({u1, v1, u2, v2}) < 4:
            continue
        if G_rewired.has_edge(u1, v2) or G_rewired.has_edge(u2, v1):
            continue
        if u1 == v2 or u2 == v1:
            continue
        
        # Perform swap
        G_rewired.remove_edge(u1, v1)
        G_rewired.remove_edge(u2, v2)
        G_rewired.add_edge(u1, v2)
        G_rewired.add_edge(u2, v1)
        n_successful += 1
    
    return G_rewired, n_successful

# ─────────────────────────────────────────────────────────────────────────────
# SINGLE NODE NULL (needed for synergy analysis)
# ─────────────────────────────────────────────────────────────────────────────

def run_single_node_null(G, query_genes, n_permutations=1000, Q=10, seed=42):
    """
    Run Maslov-Sneppen null for individual genes.
    
    Returns
    -------
    dict : {gene: np.array of null distribution}
    """
    print(f"\n  Running single-node null for {len(query_genes)} genes...")
    null_single = {}
    
    np.random.seed(seed)
    rng_seeds = np.random.randint(0, 100000, size=n_permutations)
    
    for gene in query_genes:
        if gene not in G:
            print(f"    WARNING: {gene} not in network, skipping")
            continue
        
        null_vals = []
        for i in range(n_permutations):
            if (i + 1) % 200 == 0:
                print(f"    {gene}: {i+1}/{n_permutations}")
            
            G_rewired, _ = maslov_sneppen_rewire(G, Q=Q, seed=int(rng_seeds[i]))
            null_vals.append(disconnection_fraction(G_rewired, gene))
        
        null_single[gene] = np.array(null_vals)
        print(f"    ✓ {gene}: f_obs={disconnection_fraction(G, gene):.4f}")
    
    return null_single

# ─────────────────────────────────────────────────────────────────────────────
# COMBINATORIAL NULL
# ─────────────────────────────────────────────────────────────────────────────

def run_combinatorial_null(G_original, query_combinations, n_permutations=1000, Q=10, seed=42):
    """
    For each combination, generate n_permutations rewired networks
    and measure disconnection upon simultaneous removal of all genes
    in the combination.
    
    Returns
    -------
    null_combo : dict
        {combo_tuple : np.array of length n_permutations}
    """
    # Validate
    all_genes = set().union(*[set(c) for c in query_combinations])
    missing = [g for g in all_genes if g not in G_original.nodes()]
    if missing:
        raise ValueError(f"Genes not in network: {missing}")
    
    null_combo = {combo: [] for combo in query_combinations}
    
    np.random.seed(seed)
    rng_seeds = np.random.randint(0, 100000, size=n_permutations)
    
    print(f"\n  Combinatorial null: {len(query_combinations)} combinations × {n_permutations} permutations")
    
    for i in range(n_permutations):
        if (i + 1) % 200 == 0:
            print(f"    Permutation {i+1}/{n_permutations}")
        
        # ONE rewired network per permutation --- shared across all combinations
        # This is the key efficiency: one rewiring, multiple combos tested
        G_rewired, _ = maslov_sneppen_rewire(
            G_original,
            Q=Q,
            seed=int(rng_seeds[i])
        )
        
        for combo in query_combinations:
            f = disconnection_fraction_combo(G_rewired, list(combo))
            null_combo[combo].append(f)
    
    return {k: np.array(v) for k, v in null_combo.items()}

# ─────────────────────────────────────────────────────────────────────────────
# SYNERGY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def compute_synergy(G_original, combo, f_obs_single, null_single, null_combo_dist):
    """
    Compute synergy = f_obs(combo) - sum(f_obs(individual))
    
    Also compute null synergy distribution:
    For each permutation k:
        synergy_null[k] = f_null_combo[k] - sum(f_null_single_i[k])
    
    This tests whether the real combination is MORE synergistic
    than degree-preserved random combinations.
    
    Returns dict with synergy statistics.
    """
    # Observed values
    f_obs_c = disconnection_fraction_combo(G_original, list(combo))
    f_obs_individuals = np.array([f_obs_single[g] for g in combo])
    synergy_obs = f_obs_c - np.sum(f_obs_individuals)
    
    # Null synergy distribution
    # For each permutation, synergy_null[k] = f_null_combo[k] - sum(f_null_single_i[k])
    n_perm = len(null_combo_dist)
    null_individual_arrays = np.array([null_single[g] for g in combo])
    
    # Shape: (k_genes, n_perm) → sum over genes → (n_perm,)
    null_sum_individuals = np.sum(null_individual_arrays, axis=0)
    synergy_null = null_combo_dist - null_sum_individuals
    
    mu_syn = np.mean(synergy_null)
    sig_syn = np.std(synergy_null, ddof=1)
    z_syn = (synergy_obs - mu_syn) / sig_syn if sig_syn > 0 else np.inf
    p_syn = np.mean(synergy_null >= synergy_obs)
    
    return {
        'combo': combo,
        'f_obs_combo': f_obs_c,
        'f_obs_sum_indiv': float(np.sum(f_obs_individuals)),
        'synergy_obs': synergy_obs,
        'mu_synergy_null': mu_syn,
        'sigma_synergy': sig_syn,
        'z_synergy': z_syn,
        'p_synergy': p_syn,
        'super_additive': synergy_obs > 0,
        'synergy_null_dist': synergy_null,
    }

# ─────────────────────────────────────────────────────────────────────────────
# POSITIONAL NULL STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_combo_statistics(G_original, combo, null_combo_dist):
    """
    Z-score and empirical p-value for one combination's positional null.
    
    Tests: is this combination more disruptive than the same combination
    in a topologically randomized (but degree-preserved) network?
    """
    f_obs = disconnection_fraction_combo(G_original, list(combo))
    f_null = null_combo_dist
    
    mu_null = np.mean(f_null)
    sigma_null = np.std(f_null, ddof=1)
    z_score = (f_obs - mu_null) / sigma_null if sigma_null > 0 else np.inf
    p_emp = np.mean(f_null >= f_obs)
    
    return {
        'combo': '+'.join(combo),
        'combo_size': len(combo),
        'genes': list(combo),
        'f_obs': f_obs,
        'mu_null': mu_null,
        'sigma_null': sigma_null,
        'z_score': z_score,
        'p_emp': p_emp,
        'p95_null': np.percentile(f_null, 95),
        'n_perm': len(f_null),
    }

# ─────────────────────────────────────────────────────────────────────────────
# EMPIRICAL BROWN'S METHOD
# ─────────────────────────────────────────────────────────────────────────────

def empirical_browns_method(p_values, null_distributions):
    """
    Combine p-values using Empirical Brown's Method.
    
    Brown's method accounts for correlation between tests by estimating
    the covariance structure from the null distributions.
    
    Parameters
    ----------
    p_values : list or np.array
        Individual p-values to combine
    null_distributions : list of np.arrays
        Each array is the null distribution for one test (same length)
    
    Returns
    -------
    dict with combined_p_value and other statistics
    """
    k = len(p_values)
    
    if k == 1:
        return {
            'combined_p': p_values[0],
            'df_brown': 1,
            'scaling_factor': 1.0,
            'method': 'single_test'
        }
    
    # Convert p-values to chi-square statistics
    # χ² = -2 × ln(p)
    chi2_stats = -2 * np.log(np.array(p_values))
    
    # Under independence: sum of χ² follows χ²(2k) distribution
    # But if tests are correlated, we need to estimate the effective df
    
    # Compute covariance matrix from null distributions
    null_matrix = np.column_stack(null_distributions)  # shape: (n_perm, k)
    cov_matrix = np.cov(null_matrix.T)
    
    # Brown's method scaling factor
    # E[χ²] = 2k under independence
    # Var[χ²] = 4k under independence
    # With correlation: Var[sum of χ²] = 4k + 2×sum(cov_ij)
    
    expected_var_independent = 4 * k
    
    # Estimate variance from covariance matrix
    # Var(sum X_i) = sum(Var(X_i)) + 2×sum(Cov(X_i, X_j))
    total_variance = np.sum(cov_matrix)
    
    # Scaling factor c
    c = total_variance / (2 * k)
    
    # Effective degrees of freedom
    df_effective = 2 * k / c if c > 0 else 2 * k
    
    # Combined test statistic
    T_brown = np.sum(chi2_stats) / c if c > 0 else np.sum(chi2_stats)
    
    # P-value from chi-square distribution
    combined_p = 1 - stats.chi2.cdf(T_brown, df=df_effective)
    
    return {
        'combined_p': combined_p,
        'df_brown': df_effective,
        'scaling_factor': c,
        'test_statistic': T_brown,
        'individual_p_values': p_values,
        'method': 'empirical_browns'
    }

# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────

def plot_combinatorial_results(stats_list, synergy_list, null_combo, population, tissue):
    """
    Two-panel figure per combination:
    Left  — positional null: histogram + observed
    Right — synergy null: histogram + observed synergy
    """
    n = len(stats_list)
    fig, axes = plt.subplots(n, 2, figsize=(12, 4 * n))
    
    if n == 1:
        axes = axes.reshape(1, 2)
    
    for i, (stat, syn) in enumerate(zip(stats_list, synergy_list)):
        combo_label = stat['combo']
        f_null = null_combo[tuple(stat['genes'])]
        syn_null = syn['synergy_null_dist']
        
        # ── Left panel: positional null ──────────────────────────────────────
        ax = axes[i, 0]
        ax.hist(
            f_null * 100, bins=30,
            color='#4a90d9', alpha=0.75, edgecolor='white', linewidth=0.4
        )
        ax.axvline(
            stat['f_obs'] * 100,
            color='#d62728', linewidth=2.2, linestyle='--',
            label=f"Observed: {stat['f_obs']*100:.1f}%"
        )
        ax.axvline(
            stat['mu_null'] * 100,
            color='#2ca02c', linewidth=1.5, linestyle=':',
            label=f"Null mean: {stat['mu_null']*100:.1f}%"
        )
        
        sig_col = '#d62728' if stat['p_emp'] < 0.05 else '#555'
        ax.text(
            0.97, 0.97,
            f"Z = {stat['z_score']:.2f}\np = {stat['p_emp']:.4f}",
            transform=ax.transAxes, va='top', ha='right',
            fontsize=9, color=sig_col,
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                     ec=sig_col, lw=0.8, alpha=0.9)
        )
        
        ax.set_title(
            f"Positional null: {combo_label}\n"
            f"f_obs={stat['f_obs']*100:.1f}% vs null_mean={stat['mu_null']*100:.1f}%",
            fontsize=10, fontweight='bold'
        )
        ax.set_xlabel('Disconnection fraction (%)')
        ax.set_ylabel('Permutation count')
        ax.legend(fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        
        # ── Right panel: synergy null ────────────────────────────────────────
        ax = axes[i, 1]
        ax.hist(
            syn_null * 100, bins=30,
            color='#ff7f0e', alpha=0.75, edgecolor='white', linewidth=0.4
        )
        ax.axvline(
            syn['synergy_obs'] * 100,
            color='#d62728', linewidth=2.2, linestyle='--',
            label=f"Observed synergy: {syn['synergy_obs']*100:.1f}%"
        )
        ax.axvline(0, color='#888', linewidth=1.0, linestyle='-',
                   label='Zero synergy (additive)')
        
        syn_col = '#d62728' if syn['p_synergy'] < 0.05 else '#555'
        direction = 'SUPER-additive' if syn['synergy_obs'] > 0 else 'Sub-additive'
        ax.text(
            0.97, 0.97,
            f"{direction}\nZ = {syn['z_synergy']:.2f}\np = {syn['p_synergy']:.3f}",
            transform=ax.transAxes, va='top', ha='right',
            fontsize=9, color=syn_col,
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                     ec=syn_col, lw=0.8, alpha=0.9)
        )
        
        ax.set_title(
            f"Synergy null: {combo_label}\n"
            f"f_obs_combo={syn['f_obs_combo']*100:.1f}% vs "
            f"sum_indiv={syn['f_obs_sum_indiv']*100:.1f}%",
            fontsize=10, fontweight='bold'
        )
        ax.set_xlabel('Synergy = f_combo − Σf_individual (%)')
        ax.set_ylabel('Permutation count')
        ax.legend(fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
    
    plt.suptitle(
        f"Combinatorial node removal null model\n"
        f"Maslov-Sneppen ({N_PERMUTATIONS} permutations, Q={Q}) · {population} · {tissue}",
        fontsize=12, fontweight='bold', y=1.01
    )
    plt.tight_layout()
    
    fname = f'Fig_combinatorial_{population}_{tissue}.pdf'
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"  [Saved] {fname}")

# ─────────────────────────────────────────────────────────────────────────────
# MASTER RUNNER
# ─────────────────────────────────────────────────────────────────────────────

def run_full_combinatorial_analysis(
    edge_file,
    population='EAS',
    tissue='liver',
    seed=42
):
    """
    Runs the complete combinatorial null for the combinations reported
    in the manuscript.
    
    Parameters
    ----------
    edge_file : str
        Path to CSV file with columns "gene_A", "gene_B"
    population : str
        'EAS' or 'EUR'
    tissue : str
        Tissue name (e.g., 'liver')
    seed : int
        Random seed
    """
    print(f"\n{'='*70}")
    print(f"Combinatorial Analysis: {population} · {tissue}")
    print(f"{'='*70}")
    
    # Load network
    G_original = load_network(edge_file)
    
    # Get combinations for this population/tissue
    query_combinations = COMBINATIONS.get(population, {}).get(tissue, [])
    
    if not query_combinations:
        print(f"  No combinations defined for {population} {tissue}")
        return None
    
    # Filter to combinations whose genes exist in this network
    valid_combos = []
    for combo in query_combinations:
        missing = [g for g in combo if g not in G_original.nodes()]
        if missing:
            print(f"  Skipping {combo}: {missing} not in network")
        else:
            valid_combos.append(tuple(combo))
    
    if not valid_combos:
        print("  No valid combinations for this tissue/population")
        return None
    
    print(f"  Testing {len(valid_combos)} combinations:")
    for combo in valid_combos:
        print(f"    • {'+'.join(combo)}")
    
    # ── Step 1: Run single-node null (needed for synergy) ──────────────────
    all_genes = list(set().union(*valid_combos))
    null_single = run_single_node_null(
        G_original,
        query_genes=all_genes,
        n_permutations=N_PERMUTATIONS,
        Q=Q,
        seed=seed
    )
    
    # Observed single-gene values
    f_obs_single = {
        g: disconnection_fraction(G_original, g)
        for g in all_genes if g in G_original
    }
    
    # ── Step 2: Run positional null (Maslov-Sneppen) ───────────────────────
    null_combo = run_combinatorial_null(
        G_original,
        query_combinations=valid_combos,
        n_permutations=N_PERMUTATIONS,
        Q=Q,
        seed=seed
    )
    
    # ── Step 3: Compute positional statistics ──────────────────────────────
    stats_list = []
    for combo in valid_combos:
        stat = compute_combo_statistics(G_original, combo, null_combo[combo])
        stats_list.append(stat)
    
    stats_df = pd.DataFrame([
        {k: v for k, v in s.items() if k != 'genes'}
        for s in stats_list
    ])
    
    # BH-FDR across all combinations
    if len(stats_df) > 1:
        _, q_vals, _, _ = multipletests(stats_df['p_emp'], method='fdr_bh')
        stats_df['q_BH'] = q_vals
    else:
        stats_df['q_BH'] = stats_df['p_emp']
    
    stats_df['significant_positional'] = stats_df['q_BH'] < 0.05
    stats_df['population'] = population
    stats_df['tissue'] = tissue
    
    print("\n── Positional null results ──")
    print(stats_df[[
        'combo', 'combo_size', 'f_obs', 'mu_null',
        'z_score', 'p_emp', 'q_BH', 'significant_positional'
    ]].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    
    # ── Step 4: Synergy analysis ───────────────────────────────────────────
    synergy_list = []
    for combo in valid_combos:
        # Only compute synergy if all genes have single-node null distributions
        if all(g in null_single for g in combo):
            syn = compute_synergy(
                G_original, combo,
                f_obs_single, null_single, null_combo[combo]
            )
        else:
            # Simplified: just report observed synergy without null
            f_c = disconnection_fraction_combo(G_original, list(combo))
            f_sum = sum(f_obs_single.get(g, 0) for g in combo)
            syn = {
                'combo': combo,
                'f_obs_combo': f_c,
                'f_obs_sum_indiv': f_sum,
                'synergy_obs': f_c - f_sum,
                'z_synergy': np.nan,
                'p_synergy': np.nan,
                'super_additive': (f_c - f_sum) > 0,
                'synergy_null_dist': np.zeros(N_PERMUTATIONS),
            }
        synergy_list.append(syn)
    
    synergy_df = pd.DataFrame([
        {k: v for k, v in s.items()
         if k not in ('synergy_null_dist', 'combo')}
        for s in synergy_list
    ])
    synergy_df.insert(0, 'combo', ['+'.join(c) for c in valid_combos])
    
    print("\n── Synergy analysis ──")
    print(synergy_df[[
        'combo', 'f_obs_combo', 'f_obs_sum_indiv',
        'synergy_obs', 'z_synergy', 'p_synergy', 'super_additive'
    ]].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    
    # ── Step 5: Empirical Brown's Method ───────────────────────────────────
    if len(valid_combos) > 1:
        print("\n── Empirical Brown's Method (combined p-value) ──")
        
        # Combine positional p-values
        p_values_positional = stats_df['p_emp'].values
        null_dists_positional = [null_combo[combo] for combo in valid_combos]
        
        brown_result = empirical_browns_method(p_values_positional, null_dists_positional)
        
        print(f"  Individual p-values: {[f'{p:.4f}' for p in p_values_positional]}")
        print(f"  Combined p-value (Brown's): {brown_result['combined_p']:.6f}")
        print(f"  Effective df: {brown_result['df_brown']:.2f}")
        print(f"  Scaling factor: {brown_result['scaling_factor']:.4f}")
        
        # Add to summary
        stats_df['brown_combined_p'] = brown_result['combined_p']
        stats_df['brown_df'] = brown_result['df_brown']
        
        # Combine synergy p-values (if available)
        synergy_p_values = synergy_df['p_synergy'].dropna().values
        if len(synergy_p_values) > 1:
            synergy_null_dists = [s['synergy_null_dist'] for s in synergy_list 
                                  if not np.isnan(s['p_synergy'])]
            brown_synergy = empirical_browns_method(synergy_p_values, synergy_null_dists)
            print(f"  Synergy combined p-value (Brown's): {brown_synergy['combined_p']:.6f}")
            synergy_df['brown_synergy_combined_p'] = brown_synergy['combined_p']
    
    # ── Step 6: Visualize ───────────────────────────────────────────────────
    plot_combinatorial_results(
        stats_list, synergy_list, null_combo,
        population, tissue
    )
    
    # ── Step 7: Save tables ─────────────────────────────────────────────────
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    positional_file = os.path.join(RESULTS_DIR, f'positional_{population}_{tissue}.csv')
    stats_df.drop(columns=['genes'], errors='ignore').to_csv(positional_file, index=False)
    print(f"  [Saved] {positional_file}")
    
    synergy_file = os.path.join(RESULTS_DIR, f'synergy_{population}_{tissue}.csv')
    synergy_df.to_csv(synergy_file, index=False)
    print(f"  [Saved] {synergy_file}")
    
    return {
        'positional': stats_df,
        'synergy': synergy_df,
        'null_dists': null_combo,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    
    # Example usage:
    
    # ── EAS liver ──────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("RUNNING: EAS liver")
    print("="*70)
    
    results_eas_liver = run_full_combinatorial_analysis(
        edge_file='crosstalk_edges_EAS_liver.csv',
        population='EAS',
        tissue='liver',
        seed=42
    )
    
    # ── EUR liver ──────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("RUNNING: EUR liver")
    print("="*70)
    
    results_eur_liver = run_full_combinatorial_analysis(
        edge_file='crosstalk_edges_EUR_liver.csv',
        population='EUR',
        tissue='liver',
        seed=42
    )
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved in '{RESULTS_DIR}/' directory")
    print("  • Positional null tables: positional_<POPULATION>_<TISSUE>.csv")
    print("  • Synergy tables: synergy_<POPULATION>_<TISSUE>.csv")
    print("  • Figures: Fig_combinatorial_<POPULATION>_<TISSUE>.pdf")