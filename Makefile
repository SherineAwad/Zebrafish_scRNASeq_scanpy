


AC_combined.h5ad:
	python merge_clusters.py -i filtered_corrected_AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad -o AC_combined.h5ad --annotations AC_combined_annotations.txt




figures/cabp5a_cabp5b_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png:
	python plot_per_celltypes.py annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad cabp5a_cabp5b.txt BCgroup.txt


figures/gap43_marcksa_marcksb_stmn2a_stmn2b_klf4_klf17_klf9_atoh7_pou4f2_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png:
	python plot_per_celltypes.py annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad gap43_marcksa.txt RGCgroup.txt





figures/umap_annotationsON.png:
	python annotate_subtypes.py --input AC_combined.h5ad --annotations AC_subtypes_annotations.txt --output AC_subtypes_annotated.h5ad --remove_cluster 45



AC_combined_cell_ratios.csv:
	python stats.py -i AC_subtypes_annotated.h5ad -c combined_leiden


figures/umapAC_combined_colored_by_sample.png:
	python plotAC.py -i AC_subtypes_annotated.h5ad 

figures/umap_AC_subtypes_annotated_vsx2_Harmony.png:
	python plot_markers.py AC_subtypes_annotated.h5ad markers.txt 


figures/AC_subtypes_annotated_dotplot_sorted.png:
	python plot_dotplot.py -i AC_subtypes_annotated.h5ad
