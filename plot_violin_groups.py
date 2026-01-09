import argparse
import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams

# Argument parsing
parser = argparse.ArgumentParser(description="Plot marker genes for PR precursors and Rods per sample")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('markers', help='Path to file with marker gene list (one per line)')
args = parser.parse_args()

myObject = args.myObject
markers = args.markers

# Set scanpy figure directory
sc.settings.figdir = 'figures/'
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300

# Load marker genes
with open(markers, 'r') as f:
    marker_genes = [g.strip() for g in f if g.strip()]
print("Loaded marker genes:", marker_genes)

# Read AnnData
adata = sc.read_h5ad(myObject)

# Filter to only include PR precursors and Rods
target_celltypes = ['PR precursors', 'Rod']
adata_filtered = adata[adata.obs['celltype'].isin(target_celltypes)].copy()

print(f"\nFiltered to: {adata_filtered.shape[0]} cells ({target_celltypes})")

# Filter marker genes to those present
present = [g for g in marker_genes if g in adata_filtered.var_names]
if not present:
    print("Error: No valid marker genes found")
    exit(1)

print(f"\nGenes to plot: {present}")

# Get samples
samples = sorted(adata_filtered.obs['renamed_samples'].unique())
print(f"\nSamples: {samples}")

# Colors for celltypes
celltype_colors = {'PR precursors': '#1f77b4', 'Rod': '#ff7f0e'}

# ===============================
# 1. RIDGEPLOT (all genes in one figure)
# ===============================
print("\nCreating ridgeplot...")

# Create one ridgeplot figure with all genes
n_genes = len(present)
fig, axes = plt.subplots(n_genes, 1, figsize=(10, 3 * n_genes))
if n_genes == 1:
    axes = [axes]

for idx, gene in enumerate(present):
    ax = axes[idx]
    
    # Prepare data for ridge plot
    plot_data = []
    group_labels = []
    
    for sample in samples:
        for celltype in target_celltypes:
            # Create combined label
            label = f"{sample} - {celltype}"
            
            # Filter cells
            mask = (adata_filtered.obs['renamed_samples'] == sample) & \
                   (adata_filtered.obs['celltype'] == celltype)
            
            if mask.sum() > 0:
                expr = adata_filtered[mask, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                elif hasattr(expr, 'A'):
                    expr = expr.A.flatten()
                else:
                    expr = expr.flatten()
                
                # Remove zeros for cleaner plot
                expr = expr[expr > 0]
                
                if len(expr) > 0:
                    plot_data.append(expr)
                    group_labels.append((label, celltype_colors[celltype]))
    
    # Create ridge plot using kde
    from scipy.stats import gaussian_kde
    
    # Calculate y positions
    y_positions = np.arange(len(plot_data))
    
    # Plot each density
    for i, (data, (label, color)) in enumerate(zip(plot_data, group_labels)):
        if len(data) > 1:
            kde = gaussian_kde(data)
            x_grid = np.linspace(0, np.max(data) * 1.1, 200)
            density = kde(x_grid)
            
            # Normalize and offset
            normalized_density = density / density.max() * 0.8
            ax.fill_between(x_grid, i + normalized_density, i, 
                           alpha=0.7, color=color)
            
            # Add label
            ax.text(-0.02, i + 0.4, label, ha='right', va='center', 
                   fontsize=9, transform=ax.get_yaxis_transform())
    
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.set_xlabel('Expression' if idx == n_genes - 1 else '')
    ax.set_title(gene, fontsize=12, pad=10)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')

plt.suptitle('Ridge Plots: PR precursors vs Rod per Sample', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/ridgeplot_{os.path.basename(myObject).replace('.h5ad', '')}.png", 
           bbox_inches='tight')
plt.close()

# ===============================
# 2. VIOLIN PLOT - 2 VIOLINS FOR EACH SAMPLE (all genes in one figure)
# ===============================
print("\nCreating violin plot with 2 violins per sample...")

fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 6))
if len(present) == 1:
    axes = [axes]

for idx, gene in enumerate(present):
    ax = axes[idx]
    
    # Prepare data for violin plot: 2 violins per sample (PR precursors and Rod)
    violin_data = []
    positions = []
    colors = []
    
    # Position tracking
    pos = 1
    sample_offset = 0.3  # Space between the 2 violins within a sample
    
    for sample in samples:
        # Check if we have cells for both celltypes in this sample
        has_data = False
        
        for celltype_idx, celltype in enumerate(target_celltypes):
            # Filter cells for this sample and celltype
            mask = (adata_filtered.obs['renamed_samples'] == sample) & \
                   (adata_filtered.obs['celltype'] == celltype)
            
            if mask.sum() > 0:
                expr = adata_filtered[mask, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                elif hasattr(expr, 'A'):
                    expr = expr.A.flatten()
                else:
                    expr = expr.flatten()
                
                if len(expr) > 0:
                    violin_data.append(expr)
                    
                    # Calculate position: PR precursors on left, Rod on right within sample
                    if celltype == 'PR precursors':
                        position = pos - sample_offset/2
                        color = '#1f77b4'
                    else:  # Rod
                        position = pos + sample_offset/2
                        color = '#ff7f0e'
                    
                    positions.append(position)
                    colors.append(color)
                    has_data = True
        
        # Only increment position if we had data for this sample
        if has_data:
            pos += 1.5  # Space between samples
    
    # Create violin plot
    if violin_data:
        parts = ax.violinplot(violin_data, positions=positions, showmeans=False, showmedians=True)
        
        # Color the violins
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(0.7)
        
        # Color the medians
        for partname in ('cmedians',):
            if partname in parts:
                parts[partname].set_edgecolor('black')
                parts[partname].set_linewidth(1.5)
    
    # Set x-tick labels (sample names in the middle)
    # Calculate middle position for each sample
    x_ticks = []
    x_labels = []
    
    pos = 1
    for sample in samples:
        # Check if this sample has any data
        sample_mask = adata_filtered.obs['renamed_samples'] == sample
        if sample_mask.sum() > 0:
            x_ticks.append(pos)
            x_labels.append(sample)
            pos += 1.5
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    
    # Add vertical lines to separate samples
    for tick in x_ticks:
        ax.axvline(x=tick + 0.75, color='gray', alpha=0.3, linewidth=0.5, linestyle='--')
    
    ax.set_title(gene)
    ax.set_xlabel('')
    ax.set_ylabel('Expression' if idx == 0 else '')
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add legend
    from matplotlib.patches import Patch
    if idx == 0:
        legend_elements = [
            Patch(facecolor='#1f77b4', alpha=0.7, label='PR precursors'),
            Patch(facecolor='#ff7f0e', alpha=0.7, label='Rod')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

plt.suptitle('Violin Plots: PR precursors vs Rod per Sample (2 violins per sample)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/violin_{os.path.basename(myObject).replace('.h5ad', '')}.png", 
           bbox_inches='tight')
plt.close()

# ===============================
# 3. BOX PLOT with jitter (2 boxes per sample, all genes in one figure)
# ===============================
print("\nCreating box plot with jitter...")

fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 6))
if len(present) == 1:
    axes = [axes]

for idx, gene in enumerate(present):
    ax = axes[idx]
    
    # Prepare data for box plot
    box_data = []
    positions = []
    x_labels = []
    colors = []
    
    pos = 1
    sample_offset = 0.15  # Smaller offset for boxes vs violins
    
    for sample in samples:
        has_data = False
        
        for celltype_idx, celltype in enumerate(target_celltypes):
            # Filter cells for this sample and celltype
            mask = (adata_filtered.obs['renamed_samples'] == sample) & \
                   (adata_filtered.obs['celltype'] == celltype)
            
            if mask.sum() > 0:
                expr = adata_filtered[mask, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                elif hasattr(expr, 'A'):
                    expr = expr.A.flatten()
                else:
                    expr = expr.flatten()
                
                if len(expr) > 0:
                    box_data.append(expr)
                    
                    # Calculate position with offset for celltypes
                    if celltype == 'PR precursors':
                        position = pos - sample_offset
                        color = '#1f77b4'
                    else:  # Rod
                        position = pos + sample_offset
                        color = '#ff7f0e'
                    
                    positions.append(position)
                    colors.append(color)
                    has_data = True
                    
                    # Add jitter points
                    jitter_x = np.random.normal(position, 0.01, len(expr))
                    ax.scatter(jitter_x, expr, alpha=0.3, s=10, color=color, edgecolor='none')
        
        # Add sample label
        if has_data:
            x_labels.append((pos, sample))
            pos += 1
    
    # Create box plot
    if box_data:
        bp = ax.boxplot(box_data, positions=positions, widths=0.1, 
                       patch_artist=True, showfliers=False)
        
        # Color boxes
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.7)
        
        # Color medians
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(1.5)
    
    # Set x-tick labels
    if x_labels:
        tick_positions = [pos for pos, _ in x_labels]
        tick_labels = [label for _, label in x_labels]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    
    # Add vertical lines to separate samples
    for tick in tick_positions:
        ax.axvline(x=tick + 0.5, color='gray', alpha=0.3, linewidth=0.5, linestyle='--')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', alpha=0.7, label='PR precursors'),
        Patch(facecolor='#ff7f0e', alpha=0.7, label='Rod')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    ax.set_title(gene)
    ax.set_xlabel('')
    ax.set_ylabel('Expression' if idx == 0 else '')
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle('Box Plots with Jitter: PR precursors vs Rod per Sample', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/boxplot_{os.path.basename(myObject).replace('.h5ad', '')}.png", 
           bbox_inches='tight')
plt.close()

print(f"\nDone! 3 plots saved to figures/ folder:")
print(f"1. ridgeplot_{os.path.basename(myObject).replace('.h5ad', '')}.png")
print(f"2. violin_{os.path.basename(myObject).replace('.h5ad', '')}.png")
print(f"3. boxplot_{os.path.basename(myObject).replace('.h5ad', '')}.png")
