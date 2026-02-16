import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Plot marker genes for specific celltypes per sample")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('markers', help='Path to file with marker gene list (one per line)')
parser.add_argument('celltypes', help='Path to file with celltype list (one per line)')
args = parser.parse_args()

myObject = args.myObject
markers = args.markers
celltypes_file = args.celltypes

# Load celltypes
with open(celltypes_file, 'r') as f:
    target_celltypes = [ct.strip() for ct in f if ct.strip()]
print("Target celltypes:", target_celltypes)

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

# Filter to only include specified celltypes
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

# Colors for samples (use a colormap for samples)
sample_colors = {}
cmap = plt.cm.Set3  # Good for categorical data
for i, sample in enumerate(samples):
    sample_colors[sample] = cmap(i / max(1, len(samples) - 1))

# Create a string with gene names for filenames
gene_string = "_".join(present)
base_filename = os.path.basename(myObject).replace('.h5ad', '')

# ===============================
# 1. RIDGEPLOT - Group by celltype, samples within each celltype
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
    colors = []

    # Process celltypes in order from file
    for celltype in target_celltypes:
        for sample in samples:
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
                    group_labels.append(f"{sample}\n{celltype}")
                    colors.append(sample_colors[sample])

    # Create ridge plot using kde
    from scipy.stats import gaussian_kde

    # Calculate y positions
    y_positions = np.arange(len(plot_data))

    # Plot each density
    for i, (data, label, color) in enumerate(zip(plot_data, group_labels, colors)):
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

    # Add vertical lines to separate celltypes
    # Find separation points between celltypes
    current_celltype = None
    for i, label in enumerate(group_labels):
        celltype = label.split('\n')[1]
        if current_celltype is None:
            current_celltype = celltype
        elif celltype != current_celltype:
            ax.axhline(y=i - 0.5, color='black', alpha=0.5, linestyle='--', linewidth=0.5)
            current_celltype = celltype

plt.suptitle(f'Ridge Plots: Samples within {", ".join(target_celltypes)}', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/{gene_string}_ridgeplot_{base_filename}.png",
           bbox_inches='tight')
plt.close()

# ===============================
# 2. VIOLIN PLOT - Group by celltype, samples within each
# ===============================
print("\nCreating violin plot grouped by celltype...")

fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 6))
if len(present) == 1:
    axes = [axes]

for idx, gene in enumerate(present):
    ax = axes[idx]

    # Prepare data
    violin_data = []
    positions = []
    colors = []
    group_labels = []

    # Position tracking
    start_positions = {}
    current_pos = 1
    
    # Calculate starting positions for each celltype
    for celltype_idx, celltype in enumerate(target_celltypes):
        start_positions[celltype] = current_pos
        # Update for next celltype (leave gap of 1 between celltypes)
        if celltype_idx < len(target_celltypes) - 1:
            current_pos += len(samples) + 1
        else:
            current_pos += len(samples)

    # Process each celltype in order
    for celltype in target_celltypes:
        start_pos = start_positions[celltype]
        
        for i, sample in enumerate(samples):
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
                    positions.append(start_pos + i)
                    colors.append(sample_colors[sample])
                    group_labels.append(sample)

    # Create violin plot
    if violin_data:
        parts = ax.violinplot(violin_data, positions=positions, showmeans=False, showmedians=True)

        # Color the violins by sample
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(0.7)

    # Set x-tick labels and positions
    if len(samples) > 0:
        # Set ticks at sample positions
        sample_ticks = list(positions)
        sample_tick_labels = group_labels

        ax.set_xticks(sample_ticks)
        ax.set_xticklabels(sample_tick_labels, rotation=45, ha='right', fontsize=8)

        # Add group labels for each celltype
        for celltype in target_celltypes:
            start_pos = start_positions[celltype]
            if any(celltype in ct for ct in target_celltypes):
                middle_pos = start_pos + (len(samples) - 1) / 2
                ax.text(middle_pos, -0.15, celltype, ha='center', va='top',
                        transform=ax.get_xaxis_transform(), fontsize=10, fontweight='bold')

        # Add vertical lines to separate groups
        for celltype_idx, celltype in enumerate(target_celltypes[:-1]):
            start_pos = start_positions[celltype]
            separator = start_pos + len(samples) + 0.5
            ax.axvline(x=separator, color='black', alpha=0.5, linestyle='--', linewidth=1)

    ax.set_title(gene)
    ax.set_xlabel('')
    ax.set_ylabel('Expression' if idx == 0 else '')
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle(f'Violin Plots: Samples within {", ".join(target_celltypes)}', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/{gene_string}_violin_{base_filename}.png",
           bbox_inches='tight')
plt.close()

# ===============================
# 3. BOX PLOT with jitter - Group by celltype, samples within each
# ===============================
print("\nCreating box plot with jitter...")

fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 6))
if len(present) == 1:
    axes = [axes]

for idx, gene in enumerate(present):
    ax = axes[idx]

    # Prepare data
    box_data = []
    positions = []
    colors = []
    group_labels = []

    # Position tracking (same as violin plot)
    start_positions = {}
    current_pos = 1
    
    # Calculate starting positions for each celltype
    for celltype_idx, celltype in enumerate(target_celltypes):
        start_positions[celltype] = current_pos
        # Update for next celltype (leave gap of 1 between celltypes)
        if celltype_idx < len(target_celltypes) - 1:
            current_pos += len(samples) + 1
        else:
            current_pos += len(samples)

    # Process each celltype in order
    for celltype in target_celltypes:
        start_pos = start_positions[celltype]
        
        for i, sample in enumerate(samples):
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
                    positions.append(start_pos + i)
                    colors.append(sample_colors[sample])
                    group_labels.append(sample)

                    # Add jitter points
                    jitter_x = np.random.normal(start_pos + i, 0.02, len(expr))
                    ax.scatter(jitter_x, expr, alpha=0.3, s=10, color=sample_colors[sample], edgecolor='none')

    # Create box plot
    if box_data:
        bp = ax.boxplot(box_data, positions=positions, widths=0.6,
                       patch_artist=True, showfliers=False)

        # Color boxes by sample
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.7)

        # Color medians
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(1.5)

    # Set x-tick labels and group labels
    if len(samples) > 0:
        # Set ticks at sample positions
        sample_ticks = list(positions)
        ax.set_xticks(sample_ticks)
        ax.set_xticklabels(group_labels, rotation=45, ha='right', fontsize=8)

        # Add group labels for each celltype
        for celltype in target_celltypes:
            start_pos = start_positions[celltype]
            middle_pos = start_pos + (len(samples) - 1) / 2
            ax.text(middle_pos, -0.15, celltype, ha='center', va='top',
                    transform=ax.get_xaxis_transform(), fontsize=10, fontweight='bold')

        # Add vertical lines to separate groups
        for celltype_idx, celltype in enumerate(target_celltypes[:-1]):
            start_pos = start_positions[celltype]
            separator = start_pos + len(samples) + 0.5
            ax.axvline(x=separator, color='black', alpha=0.5, linestyle='--', linewidth=1)

    ax.set_title(gene)
    ax.set_xlabel('')
    ax.set_ylabel('Expression' if idx == 0 else '')
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle(f'Box Plots with Jitter: Samples within {", ".join(target_celltypes)}', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(f"figures/{gene_string}_boxplot_{base_filename}.png",
           bbox_inches='tight')
plt.close()

print(f"\nDone! 3 plots saved to figures/ folder:")
print(f"1. figures/{gene_string}_ridgeplot_{base_filename}.png")
print(f"2. figures/{gene_string}_violin_{base_filename}.png")
print(f"3. figures/{gene_string}_boxplot_{base_filename}.png")
