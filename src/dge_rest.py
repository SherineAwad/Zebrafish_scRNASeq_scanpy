import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text  # For automatic non-overlapping labels with arrows

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
    parser.add_argument("-c", "--csv", required=True, help="Output CSV file for full DE table")
    parser.add_argument("-g", "--groupby", default=None, help="Column to group by (kept for style)")
    parser.add_argument("-o", "--out", required=True, help="Output H5AD file with DE results")
    parser.add_argument("-n", "--topn", type=int, default=None, help="Number of top genes per group")
    parser.add_argument("-m", "--metric", type=str, default="wilcoxon_score",
                        help="Metric to rank genes by (wilcoxon_score, logfoldchange, pval, pval_adj)")
    parser.add_argument("--sample1", required=True, help="Group to test (Control)")
    parser.add_argument("-k", "--topk", type=int, help="Top K genes to label in volcano plot")

    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if "renamed_samples" not in adata.obs:
        raise KeyError("Column 'renamed_samples' not found in adata.obs")

    # ---------------- SUBSET: CONTROL + REST ----------------
    all_samples = adata.obs["renamed_samples"].unique().tolist()
    rest_samples = [s for s in all_samples if s != args.sample1]

    # create a new column with only two groups: Control vs rest
    adata.obs['Control_vs_rest'] = adata.obs["renamed_samples"].apply(
        lambda x: args.sample1 if x == args.sample1 else "rest"
    )
    sub = adata.copy()
    sub.uns["log1p"] = {"base": None}

    # ---------------- RANK GENES ----------------
    sc.tl.rank_genes_groups(
        sub,
        groupby="Control_vs_rest",
        groups=[args.sample1],
        reference="rest",
        method="wilcoxon",
        n_genes=None
    )

    # ---------------- EXTRACT RESULTS ----------------
    r = sub.uns["rank_genes_groups"]
    g = args.sample1
    rows = []
    for i in range(len(r["names"][g])):
        rows.append({
            "group": f"{args.sample1}_vs_rest",
            "gene": r["names"][g][i],
            "wilcoxon_score": r["scores"][g][i],
            "pval": r["pvals"][g][i],
            "pval_adj": r["pvals_adj"][g][i],
            "logfoldchange": r["logfoldchanges"][g][i]
        })

    df = pd.DataFrame(rows)

    # ---------------- TOP-N ----------------
    if args.topn is not None:
        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not in DE table")
        finite_df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchange"])
        ascending = args.metric in ["pval", "pval_adj"]
        df = finite_df.sort_values(args.metric, ascending=ascending).head(args.topn)

    # ---------------- ADD TOPK ----------------
    if args.topk is not None:
        df["topk"] = args.topk

    # ---------------- SAVE CSV AND H5AD ----------------
    df.to_csv(args.csv, index=False)
    sub.uns["dge_table"] = df
    sub.write_h5ad(args.out)
    print(f"Saved DGE table → {args.csv}")
    print(f"Saved DE-containing object → {args.out}")

    # ---------------- VOLCANO PLOT ----------------
    plt.figure(figsize=(6,6))
    x = r["logfoldchanges"][g]
    y = -np.log10(np.array(r["pvals_adj"][g]) + 1e-300)
    plt.scatter(x, y, s=5, alpha=0.7)
    plt.xlabel("log2 fold change")
    plt.ylabel("-log10(adj pval)")
    plt.title(f"Volcano plot: {args.sample1} vs rest")
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
    plt.axvline(x=0, color='grey', linestyle='--')
    plt.xlim(-20, 30)  # <-- ONLY CHANGE: Set x-axis limits from -20 to 30

    # ---------------- LABEL TOP-K GENES WITH adjustText (all get arrows) ----------------
    if args.topk is not None:
        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not in DE table")
        ascending = args.metric in ["pval", "pval_adj"]
        finite_df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchange"])
        top_genes = finite_df.sort_values(args.metric, ascending=ascending).head(args.topk)

        texts = []
        for _, row in top_genes.iterrows():
            x_val = row["logfoldchange"]
            y_val = -np.log10(row["pval_adj"] + 1e-300)
            # place the text at the same spot as the point, arrow will always be drawn
            texts.append(plt.text(x_val, y_val, row["gene"], fontsize=6))

        # force all arrows to be drawn by setting 'force_text=True' and tweaking expand parameters
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
            expand_points=(2.0,2.0),
            force_text=True
        )

    fig_dir = "figures"
    os.makedirs(fig_dir, exist_ok=True)
    csv_prefix = os.path.splitext(os.path.basename(args.csv))[0]
    fig_path = os.path.join(fig_dir, f"{csv_prefix}_volcano.png")
    plt.savefig(fig_path, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"Saved Volcano plot → {fig_path}")

if __name__ == "__main__":
    main()
