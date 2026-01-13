


AC_combined.h5ad:
	python merge_clusters.py -i filtered_corrected_AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad -o AC_combined.h5ad --annotations AC_combined_annotations.txt


figures/umap_AC_combined_vsx2_Harmony.png:
	python plot_markers.py AC_combined.h5ad markers.txt
