# Zebrafish Project 

We analyse Zebrafish LD, NMDA and control scRNASeq samples using scanpy. 




##  Quality Control Plots

These violin plots visualize cell-level quality metrics **before and after filtering** low-quality cells.

---

###  `violin_QC.png`

This plot shows **raw quality control metrics** (e.g., number of genes per cell, mitochondrial percentage) **before filtering**. It helps ide.pngy thresholds for removing low-quality or dead cells.

![violin_QC.png](figures/violin_QC.png)

---

###  `violin_AfterQC.png`

This plot shows the same QC metrics **after filtering**. It confirms that poor-quality cells were successfully removed based on chosen thresholds.

![violin_AfterQC.png](figures/violin_AfterQC.png)


## UMAP After Batch Effect Removal (Harmony)

### `umap_sHarmonyzebrafishes.png`

This UMAP plot shows zebrafish single-cell data **after batch effect correction** using **Harmony**.

By aligning shared biological structure across different batches or samples, Harmony improves the clustering and visualization. Cells are now grouped primarily by **biological similarity** rather than technical batch differences, enabling more accurate downstream analysis.

![umap_sHarmonyzebrafishes.png](figures/umap_sHarmonyzebrafishes.png)

### UMAP per sample

Below are the UMAP plots split by individual samples after reclustering:

### UMAP per sample

Below are the UMAP plots split by individual samples after reclustering:

<img src="figures/umap_reclustered_LD.png" alt="UMAP LD sample" width="33%"><img src="figures/umap_reclustered_NMDA.png" alt="UMAP NMDA sample" width="33%"><img src="figures/umap_reclustered_Control.png" alt="UMAP Control sample" width="33%">

## Clustering 

### UMAP after clustering (resolution: default =1)  

<img src="figures/originalgeneByGene/umap_clustersNumbers.png" alt="Umap Annotations" width="90%">

### Dot plot for Marker Genes

<img src="figures/originalgeneByGene/dotplot_markerGenes.png" alt="Umap Annotations" width="90%">


### Marker Genes UMAP 

<img src="figures/originalgeneByGene/umap_acta2.png" alt="acta2" width="33%"><img src="figures/originalgeneByGene/umap_apoea.png" alt="apoea" width="33%"><img src="figures/originalgeneByGene/umap_apoeb.png" alt="apoeb" width="33%">
<img src="figures/originalgeneByGene/umap_aqp4.png" alt="aqp4" width="33%"><img src="figures/originalgeneByGene/umap_arr3a.png" alt="arr3a" width="33%"><img src="figures/originalgeneByGene/umap_ascl1a.png" alt="ascl1a" width="33%">
<img src="figures/originalgeneByGene/umap_atoh7.png" alt="atoh7" width="33%"><img src="figures/originalgeneByGene/umap_bhlhe23.png" alt="bhlhe23" width="33%"><img src="figures/originalgeneByGene/umap_cabp5a.png" alt="cabp5a" width="33%">
<img src="figures/originalgeneByGene/umap_cabp5b.png" alt="cabp5b" width="33%"><img src="figures/originalgeneByGene/umap_calb1.png" alt="calb1" width="33%"><img src="figures/originalgeneByGene/umap_calb2a.png" alt="calb2a" width="33%">
<img src="figures/originalgeneByGene/umap_calb2b.png" alt="calb2b" width="33%"><img src="figures/originalgeneByGene/umap_cbln4.png" alt="cbln4" width="33%"><img src="figures/originalgeneByGene/umap_ccnd1.png" alt="ccnd1" width="33%">
<img src="figures/originalgeneByGene/umap_cdh2.png" alt="cdh2" width="33%"><img src="figures/originalgeneByGene/umap_chata.png" alt="chata" width="33%"><img src="figures/originalgeneByGene/umap_cdk1.png" alt="cdk1" width="33%">
<img src="figures/originalgeneByGene/umap_crx.png" alt="crx" width="33%"><img src="figures/originalgeneByGene/umap_csf2rb.png" alt="csf2rb" width="33%"><img src="figures/originalgeneByGene/umap_dla.png" alt="dla" width="33%">
<img src="figures/originalgeneByGene/umap_ebf3a.png" alt="ebf3a" width="33%"><img src="figures/originalgeneByGene/umap_elavl3.png" alt="elavl3" width="33%"><img src="figures/originalgeneByGene/umap_fgf19.png" alt="fgf19" width="33%">
<img src="figures/originalgeneByGene/umap_gad1a.png" alt="gad1a" width="33%"><img src="figures/originalgeneByGene/umap_gad1b.png" alt="gad1b" width="33%"><img src="figures/originalgeneByGene/umap_gad2.png" alt="gad2" width="33%">
<img src="figures/originalgeneByGene/umap_gfap.png" alt="gfap" width="33%"><img src="figures/originalgeneByGene/umap_gli1.png" alt="gli1" width="33%"><img src="figures/originalgeneByGene/umap_gnat2.png" alt="gnat2" width="33%">
<img src="figures/originalgeneByGene/umap_guca1b.png" alt="guca1b" width="33%"><img src="figures/originalgeneByGene/umap_her12.png" alt="her12" width="33%"><img src="figures/originalgeneByGene/umap_her4.2.png" alt="her4.2" width="33%">
<img src="figures/originalgeneByGene/umap_her4.3.png" alt="her4.3" width="33%"><img src="figures/originalgeneByGene/umap_her4.4.png" alt="her4.4" width="33%"><img src="figures/originalgeneByGene/umap_her6.png" alt="her6" width="33%">
<img src="figures/originalgeneByGene/umap_igf2a.png" alt="igf2a" width="33%"><img src="figures/originalgeneByGene/umap_igf2b.png" alt="igf2b" width="33%"><img src="figures/originalgeneByGene/umap_insm1a.png" alt="insm1a" width="33%">
<img src="figures/originalgeneByGene/umap_isl1.png" alt="isl1" width="33%"><img src="figures/originalgeneByGene/umap_isl2a.png" alt="isl2a" width="33%"><img src="figures/originalgeneByGene/umap_isl2b.png" alt="isl2b" width="33%">
<img src="figures/originalgeneByGene/umap_kcnj8.png" alt="kcnj8" width="33%"><img src="figures/originalgeneByGene/umap_lhx1a.png" alt="lhx1a" width="33%"><img src="figures/originalgeneByGene/umap_mbpa.png" alt="mbpa" width="33%">
<img src="figures/originalgeneByGene/umap_mpeg1.1.png" alt="mpeg1.1" width="33%"><img src="figures/originalgeneByGene/umap_neflb.png" alt="neflb" width="33%"><img src="figures/originalgeneByGene/umap_nefla.png" alt="nefla" width="33%">
<img src="figures/originalgeneByGene/umap_nefma.png" alt="nefma" width="33%"><img src="figures/originalgeneByGene/umap_nefmb.png" alt="nefmb" width="33%"><img src="figures/originalgeneByGene/umap_neurod1.png" alt="neurod1" width="33%">
<img src="figures/originalgeneByGene/umap_notch1a.png" alt="notch1a" width="33%"><img src="figures/originalgeneByGene/umap_notch1b.png" alt="notch1b" width="33%"><img src="figures/originalgeneByGene/umap_nr2e3.png" alt="nr2e3" width="33%">
<img src="figures/originalgeneByGene/umap_nrl.png" alt="nrl" width="33%"><img src="figures/originalgeneByGene/umap_ompa.png" alt="ompa" width="33%"><img src="figures/originalgeneByGene/umap_opn1mw1.png" alt="opn1mw1" width="33%">
<img src="figures/originalgeneByGene/umap_opn1mw2.png" alt="opn1mw2" width="33%"><img src="figures/originalgeneByGene/umap_opn1mw3.png" alt="opn1mw3" width="33%"><img src="figures/originalgeneByGene/umap_opn1mw4.png" alt="opn1mw4" width="33%">
<img src="figures/originalgeneByGene/umap_opn1sw1.png" alt="opn1sw1" width="33%"><img src="figures/originalgeneByGene/umap_opn1sw2.png" alt="opn1sw2" width="33%"><img src="figures/originalgeneByGene/umap_pax2a.png" alt="pax2a" width="33%">



## After annotation using marker genes

<img src="figures/umapZebrafishes_annotations.png" alt="Umap Annotations" width="90%">  
<img src="figures/umapZebrafishes_annotationsON.png" alt="Umap Annotations ON" width="90%">

## We need to increase the resolution to 2.5 to seperate the tiny branch of Rod Cells

🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑

## 🚨 Starting from here, we have plots of resolution 2.5 🚨

🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑


### UMAP after clustering 

<img src="figures/umap_clusters.png" alt="Umap Annotations" width="90%">


## After re-Clustering annotation using marker genes

<img src="figures/Zebrafishes_annotations.png" alt="Umap Annotations" width="90%">
<img src="figures/Zebrafishes_annotationsON.png" alt="Umap Annotations ON" width="90%">


## Remove Celltype Cones_MG_MGPC_PostMitotic as likely doublet 

<img src="figures/Zebrafishes_Nannotations.png" alt="Umap Annotations" width="90%">
<img src="figures/Zebrafishes_NannotationsON.png" alt="Umap Annotations ON" width="90%">
 

## Cell ratio after increasing resolution and re-annotations 

| Sample   | MG   | MGPC | PR precursors | Rod  | Cones | BC    | AC    | HC    | RGC  | Microglia_ImmuneCells | RPE  | Melanocyte | Endothelial | Perycites | Oligodenrocyte |
|----------|------|------|---------------|------|-------|-------|-------|-------|------|------------------------|------|-------------|--------------|-----------|----------------|
| Control  | 6.49 | 1.13 | 0.88          | 4.76 | 7.54  | 46.85 | 10.24 | 13.40 | 3.20 | 1.74                   | 1.34 | 0.27        | 0.57         | 1.09      | 0.51           |
| LD       | 43.12| 15.54| 14.44         | 7.90 | 3.32  | 4.86  | 3.05  | 6.17  | 1.40 | 0.14                   | 0.03 | 0.00        | 0.00         | 0.02      | 0.02           |
| NMDA     | 29.89| 14.74| 13.87         | 13.74| 2.75  | 5.50  | 9.23  | 6.94  | 2.96 | 0.22                   | 0.00 | 0.02        | 0.02         | 0.05      | 0.07           |


<img src="figures/Restacked_bar_sample_by_celltype.png" width="90%">

## Umap for Control,LD, and NMDA merged 

<img src="figures/umap_merged_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png" width="90%"> 

## Umap per sample 
<img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Control.png" width="33%"><img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_LD.png" width ="33%"><img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_NMDA.png" width="33%">



### No. of Cells per Sample 

| renamed_samples | Number of cells |
|-----------------|-----------------|
| Control         | 134,237         |
| NMDA            | 15,424          |
| LD              | 11,214          |


### Mini Dotplot 

<img src="figures/dotplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes_markerGenes.png" width="90%">

### Selected Markers UMAP per request 

<img src="figures/umap_annotated_clustered_corrected_prdm1a_Harmony.png" width="33%"><img src="figures/umap_annotated_clustered_corrected_prdm1b_Harmony.png" width="33%"><img src="figures/umap_annotated_clustered_corrected_prkcaa_Harmony.png" width="33%">

<img src="figures/umap_annotated_clustered_corrected_nr2e3_Harmony.png?v=1" width="33%">


### Selected Markers Violin per request

<img src="figures/violin_annotated_clustered_corrected_ascl1a_Harmony.png?v=2" width="33%"><img src="figures/violin_annotated_clustered_corrected_atoh7_Harmony.png?v=2" width="33%"><img src="figures/violin_annotated_clustered_corrected_cdk1_Harmony.png?v=2" width="33%">

<img src="figures/violin_annotated_clustered_corrected_nr2e3_Harmony.png?v=2" width="33%"><img src="figures/violin_annotated_clustered_corrected_nrl_Harmony.png?v=2" width="33%">



### Selected genes Comparisons between groups  


#### elavl3/gad2/slc6a9 in AC
![](figures/elavl3_gad2_slc6a9_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/elavl3_gad2_slc6a9_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


#### cabp5a/cabp5b in BC
![](figures/cabp5a_cabp5b_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/cabp5a_cabp5b_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


#### prdm1a/nr2e3 in Rod and PR precursors 
![](figures/prdm1a_nr2e3_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/prdm1a_nr2e3_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)

#### gnat2/arr3a in Cones
![](figures/gnat2_arr3a_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/gnat2_arr3a_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


#### isl2b/rbpms2b in RGC
![](figures/isl2b_rbpms2b_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/isl2b_rbpms2b_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


#### ompa/gad2 in HC
![](figures/ompa_gad2_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/ompa_gad2_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


#### rho/guca1b in Rod and PR precursors
![](figures/rho_guca1b_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)
![](figures/rho_guca1b_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=3)


### gap43/marcksa/marcksb/stmn2a/stmn2b/klf4/klf17/klf9/atoh7/pou4f2 in RGC
![](figures/gap43_marcksa_marcksb_stmn2a_stmn2b_klf4_klf17_klf9_atoh7_pou4f2_violin_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=1)
![](figures/gap43_marcksa_marcksb_stmn2a_stmn2b_klf4_klf17_klf9_atoh7_pou4f2_boxplot_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png?v=1)

### UMAPs of Major Retinal Cell Subtypes

The following UMAP plots illustrate the clustering of subtypes within each major retinal cell type. These visualizations highlight transcriptional diversity captured by single-cell RNA sequencing.

#### Cones Photoreceptor Subtypes before batch correction 
<img src="figures/umapCones.png" width="70%"> 
<img src="figures/umap_merged_Cones.png" width ="70%">

<img src= "figures/umap_Cones_NMDA.png" width ="33%"><img src= "figures/umap_Cones_LD.png" width ="33%"><img src= "figures/umap_Cones_Control.png" width ="33%">


#### Cones Photoreceptor Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony.png" width ="90%">

<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_LD.png" width ="33%"><img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_NMDA.png" width ="33%"><img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_Control.png" width ="33%">

### QC after batch correction (STRICT PARAMETERS)
<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_leiden.png" width ="90%">
<img src="figures/violin_Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_QC.png" width ="90%">

### Cones Marker Genes after batch correction 

<img src="figures/umap_corrected_Cones_gnat2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_grk7b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_tgfa_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_tbx2a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_opn1sw1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_mpzl2b_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_opn1sw2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_opn1mw4_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_opn1lw1_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_snap25a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_gngt2a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_arr3a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_thrb_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_opn1lw2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_rho_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_nr2e3_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_vsx1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_cabp5a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_ompa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_rlbp1a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_aqp4_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_arr3a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_elavl3_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_isl2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_pou4f2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_mpeg1.1_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_Cones_mbpa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_Cones_rpe65a_sHarmony.png" width="30%" />



#### Amacrine Cell (AC) Subtypes
<img src="figures/umapAC.png" width="70%">
<img src="figures/umap_merged_AC.png" width ="70%">

<img src= "figures/umap_AC_NMDA.png" width ="33%"><img src= "figures/umap_AC_LD.png" width ="33%"><img src= "figures/umap_AC_Control.png" width ="33%">


#### AC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony.png" width ="90%">

<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_LD.png" width ="33%"><img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_NMDA.png" width ="33%"><img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_Control.png" width ="33%">

### QC after batch correction (STRICT PARAMETERS)
<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_leiden.png" width ="90%">
<img src="figures/violin_AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_QC.png" width ="90%">

### AC Marker Genes after batch correction 

<img src="figures/umap_corrected_AC_tfap2a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_chata_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_gad1a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_gad2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_slc6a9_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_th_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_slc18a3a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_chgb_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_bhlhe22_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_bhlhe23_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_fosab_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_tkta_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_hmgb2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_sox2_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_sox4a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_rho_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_nr2e3_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_vsx1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_cabp5a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_ompa_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_gnat2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_rlbp1a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_aqp4_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_arr3a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_elavl3_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_isl2b_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_pou4f2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_mpeg1.1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_AC_mbpa_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_AC_rpe65a_sHarmony.png" width="30%" />



### AC clusters after removing clusters 11, 48, and 49
<img src="figures/corrected_AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_filteredannotationsON.png" width ="90%"> 


### Recombining subclusters 1 and 45 in AC

<img src="figures/umapAC_combined_colored_by_sample.png?v=3" width="70%">
<img src="figures/umapAC_combined_combined_leiden.png?v=3" width="70%">

<div style="display: flex; gap: 10px;">
  <img src="figures/AC_combined_Control_no_leiden.png?v=5" width="30%">
  <img src="figures/AC_combined_LD_no_leiden.png?v=5" width="30%">
  <img src="figures/AC_combined_NMDA_no_leiden.png?v=5" width="30%">
</div>


<div style="display: flex; gap: 10px;">
  <img src="figures/umapAC_combined_Control.png?v=3" width="30%">
  <img src="figures/umapAC_combined_LD.png?v=3" width="30%">
  <img src="figures/umapAC_combined_NMDA.png?v=3" width="30%">
</div>


### Dotplot for subclusters in AC 


![Dotplot](figures/AC_combined_dotplot_sorted.png?v=1)

![Dotplot Full](figures/dotplot__AC_combined_annotated_dotplot_sorted_full.png?v=1) 

### Feature plots for subclusters in AC after recombining 1 and 45 

<img src="figures/umap_AC_combined_apoeb_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_cx52.7_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_isl1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_opn1sw2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc5a7a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_tfap2a_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_aqp4_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_cx52.9_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_isl2b_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_pcp4a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc6a1b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_tgfa_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_arr3a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_eomesa_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_lhx1a_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_rbpms2a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc6a1l_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_th_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_bhlhe22_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_fosab_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_mafaa_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_rbpms2b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc6a5_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_thrb_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_bhlhe23_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gad1a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_meis2a_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_rho_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc6a9_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_tkta_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_cabp5a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gad1b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_mpzl2b_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_rlbp1a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_snap25a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_vamp1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_cabp5b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gad2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_nos1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_robo2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_sox2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_vip_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_calb1_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gfap_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_npy_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_slc17a6a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_sox4a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_vsx1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_calb2a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gnat2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_nr2e3_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_slc17a6b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_sst1.1_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_vsx2_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_calb2b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_gngt2a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_ompa_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_slc17a7a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_sst2_Harmony.png?v=2" width="30%" /> 
<img src="figures/umap_AC_combined_chata_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_grk7b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_opn1lw1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_slc17a8_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_stat3_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_chgb_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_grm6a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_opn1lw2_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc17a9a_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_tac1_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_crhb_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_grm6b_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_opn1mw4_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_slc17a9b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_tbx2a_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_cx52.6_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_hmgb2b_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_opn1sw1_Harmony.png?v=2" width="30%" />
<img src="figures/umap_AC_combined_slc18a3a_Harmony.png?v=2" width="30%" /> <img src="figures/umap_AC_combined_tcf7l2_Harmony.png?v=2" width="30%" />

### Cell counts and cell ratio in AC subset

[[View stats summary for AC here](https://docs.google.com/spreadsheets/d/1nlbB0la8o_QE1IjOtlR64hxlaQ4oTQQTXWKVcIf17dQ/edit?usp=sharing)

#### Retinal Ganglion Cell (RGC) Subtypes
<img src="figures/umapRGC.png" width="70%">
<img src="figures/umap_merged_RGC.png" width ="70%">

<img src="figures/umap_RGC_NMDA.png" width ="33%"><img src="figures/umap_RGC_LD.png" width ="33%"><img src="figures/umap_RGC_Control.png" width ="33%">



#### RGC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony.png" width ="90%">

<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_LD.png" width ="33%"><img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_NMDA.png" width ="33%"><img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_Control.png" width ="33%">

### QC after batch correction (STRICT PARAMETERS)
<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_leiden.png" width ="90%">
<img src="figures/violin_RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_QC.png" width ="90%">


### RGC Marker Genes after batch correction

<img src="figures/umap_corrected_RGC_robo2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_isl2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_rbpms2a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_RGC_rbpms2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_mafaa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_eomesa_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_RGC_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_arr3a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_elavl3_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_RGC_isl2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_pou4f2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_mpeg1.1_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_RGC_mbpa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_RGC_rpe65a_sHarmony.png" width="30%" /><img src="figures/umap_corrected_RGC_tbr1b_sHarmony.png" width="30%" />


#### Horizontal Cell (HC) Subtypes
<img src="figures/umapHC.png" width="70%">
<img src="figures/umap_merged_HC.png" width ="70%">

<img src= "figures/umap_HC_NMDA.png" width ="33%"><img src= "figures/umap_HC_LD.png" width ="33%"><img src= "figures/umap_HC_Control.png" width ="33%">


#### HC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony.png" width ="90%">

<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_LD.png" width ="33%"><img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_NMDA.png" width ="33%"><img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_Control.png" width ="33%">

### QC after batch correction (STRICT PARAMETERS)
<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_leiden.png" width ="90%">
<img src="figures/violin_HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_QC.png" width ="90%">

### HC Marker Genes after batch correction

<img src="figures/umap_corrected_HC_ompa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_cx52.6_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_cx52.9_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_cx52.7_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_pcp4a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_lhx1a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_isl1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_gad2_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_rho_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_nr2e3_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_vsx1_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_cabp5a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_gnat2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_rlbp1a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_aqp4_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_arr3a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_elavl3_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_isl2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_pou4f2_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_HC_mpeg1.1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_mbpa_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_HC_rpe65a_sHarmony.png" width="30%" />

### HC clusters after removing clusters 14, 15, and 16 

<img src="figures/corrected_HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_filteredannotationsON.png" width ="90%">


#### Bipolar Cell (BC) Subtypes
<img src="figures/umapBC.png" width="70%">
<img src="figures/umap_merged_BC.png" width ="70%">

<img src= "figures/umap_BC_NMDA.png" width ="33%"><img src= "figures/umap_BC_LD.png" width ="33%"><img src= "figures/umap_BC_Control.png" width ="33%">


### BC Marker Genes without batch correction 

<img src="figures/umap_BC_annotated_vsx2_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_cabp5a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_cabp5b_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_bhlhe23_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_vamp1_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_vsx1_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_grm6a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_grm6b_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_calb1_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_calb2a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_calb2b_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_rho_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_nr2e3_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_ompa_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_gnat2_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_rlbp1a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_aqp4_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_gad1b_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_arr3a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_elavl3_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_isl2b_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_pou4f2_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_mpeg1.1_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_mbpa_Harmony.png" width="30%" />
<img src="figures/umap_BC_annotated_rpe65a_Harmony.png" width="30%" /> <img src="figures/umap_BC_annotated_prkcaa_Harmony.png" width="30%" />


#### MG Subtypes
<img src="figures/umapMG.png" width="70%">
<img src="figures/umap_merged_MG.png" width ="70%">

<img src= "figures/umap_MG_NMDA.png" width ="33%"><img src= "figures/umap_MG_LD.png" width ="33%"><img src= "figures/umap_MG_Control.png" width ="33%">


#### MG Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony.png" width ="90%">

<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_LD.png" width ="33%"><img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_NMDA.png" width ="33%"><img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_Control.png" width ="33%">


### QC after batch correction (STRICT PARAMETERS)
<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_leiden.png" width ="90%">
<img src="figures/violin_MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_sHarmony_QC.png" width ="90%">

### MG marker Genes after batch correction

<img src="figures/umap_corrected_MG_gfap_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_apoeb_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_aqp4_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_gfap_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_rlbp1a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_stat3_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_sox2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_rho_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_nr2e3_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_vsx1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_cabp5a_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_ompa_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_gnat2_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_gad1b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_arr3a_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_elavl3_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_isl2b_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_pou4f2_sHarmony.png" width="30%" />
<img src="figures/umap_corrected_MG_mpeg1.1_sHarmony.png" width="30%" /> <img src="figures/umap_corrected_MG_mbpa_sHarmony.png" width="30%" />


### MG clusters after removing clusters 11 and 13 

<img src="figures/corrected_MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_filteredannotationsON.png" width ="90%">


## Pearson correlation 

- The script first calculates the average expression of each gene within each cell type, separately for the two conditions you’re comparing. So for every gene, you have an average expression value for each cell type under condition 1 and condition 2.

- Once these averages are calculated, the script looks across all cell types to see how highly each gene is expressed in at least one cell type. It takes the maximum of these average expressions for each gene. This gives a single value per gene representing its highest mean expression in any cell type.

- Then, the script filters out genes whose maximum mean expression is below the threshold you set. In other words, genes that are not expressed strongly in any cell type are removed from the analysis. Only genes that have at least one cell type where their mean expression exceeds the threshold are kept.

- This way, the correlation analysis later is done only on genes that are meaningfully expressed, avoiding noise from genes that are essentially silent across all cell types.




```bash
python pearson_heatmap.py annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad renamed_samples Control LD control_vs_ld_heatmap.png --min_mean_expr n
``` 
| Cell Type | Control | LD |
|------------|:-------:|:--:|
| AC | ✅ | ✅ |
| BC | ✅ | ✅ |
| Cones | ✅ | ✅ |
| HC | ✅ | ✅ |
| MG | ✅ | ✅ |
| MGPC | ✅ | ✅ |
| Microglia_ImmuneCells | ✅ | ✅ |
| Oligodenrocyte | ✅ | ✅ |
| PR precursors | ✅ | ✅ |
| Perycites | ✅ | ✅ |
| RGC | ✅ | ✅ |
| RPE | ✅ | ✅ |
| Rod | ✅ | ✅ |
| **Melanocyte** | ✅ | ❌ |
| **Endothelial** | ✅ | ❌ |

## Control vs LD 


##### > 0.05 
![pearson correlation Heatmap](control_vs_ld_heatmap_0.05.png?v=3)


##### > 0.1 
![pearson correlation Heatmap](control_vs_ld_heatmap_0.1.png?v=3)

##### > 0.2 
![pearson correlation Heatmap](control_vs_ld_heatmap_0.2.png?v=3)


##### > 0.5 
![pearson correlation Heatmap](control_vs_ld_heatmap_0.5.png?v=3)



## Control vs NMDA




##### > 0.05 

![control vs nmda 0.05](control_vs_nmda_heatmap_0.05.png?v=1)


##### > 0.1 

![control vs nmda 0.1](control_vs_nmda_heatmap_0.1.png?v=1)


##### > 0.2 

![control vs nmda 0.2](control_vs_nmda_heatmap_0.2.png?v=1)

#####  > 0.5 

![control vs nmdaI 0.5](control_vs_nmda_heatmap_0.5.png?v=1)


## LD vs NMDA


##### > 0.05 

![ld vs nmda 0.05](ld_vs_nmda_heatmap_0.05.png?v=1)


##### > 0.1 

![ld vs nmda 0.1](ld_vs_nmda_heatmap_0.1.png?v=1)


##### > 0.2

![ld vs nmda 0.2](ld_vs_nmda_heatmap_0.2.png?v=1)
 

##### > 0.5 

![ld vs nmda 0.5](ld_vs_nmda_heatmap_0.5.png?v=1)



# Differential gene expression: per celltype 

```
 sc.tl.rank_genes_groups(
        adata_for_dge,
        groupby=groupby, #groupby can be celltype or leiden ..etc  
        method="wilcoxon",
        reference="rest",
        use_raw=True
    )
```
#### groupby = celltype 

[Download all DGE here](https://docs.google.com/spreadsheets/d/1yGYvD2I4IOZQUterO7G_7qOr6lOub_5T3ggRZiGe8nY/edit?usp=sharing)


### Heatmap by Wilcoxon score 

![Zebrafish DGE heatmap](figures/zebrafish_dge_heatmap.png?v=3)

## Differential gene expression per subset 

#### AC: groubpy = leiden 

[AC DGE Top 50 Wilcoxon score](https://docs.google.com/spreadsheets/d/1e0EvUsbGpTTTzivpPWfO1h9YnlFL4DihGgNPsW5IpiU/edit?usp=sharing)

[AC DGE Top 50 logFC](https://docs.google.com/spreadsheets/d/1Zu66WvspVaGxdWB2jrI-2t2_f3AgjaE7euTH59H3pl4/edit?usp=sharing)

[AC DGE Top 50 adj Pvalue](https://docs.google.com/spreadsheets/d/1UFqpXxLPI6HPbHYcpuci8zIg9Y7f6t8OmN4ezhrx75o/edit?usp=sharing)

##### Warning: very big file
[Download AC DGE here](https://drive.google.com/file/d/1915_eRknU12MrrVRrWtWEOee-vEhoXS9/view?usp=sharing)

  
### Heatmap by Wilcoxon score
![AC heatmap](figures/AC_dge_heatmap.png?v=4)


# :rotating_light: Debugging section

### For testing purposes: Control vs Control 

![control vs control](control_vs_control_heatmap_0.0.png?v=4) 


# 🐍 For Debugging Issues

### 💾 Command & Output

```
# Name                     Version          Build            Channel
scanpy                     1.9.3            pypi_0           pypi
```

```bash
python print_h5ad.py annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad 
Loading AnnData object from: annotated_clustered_corrected_doubletRemoved_Zebrafishes.h5ad

==================================================
BASIC AnnData INFORMATION
==================================================
Overall shape: (160875, 24597) (cells: 160875, genes: 24597)
Observation names (cells): ['AAACCAAAGCTTAACG-1', 'AAACCAAAGTACGCCG-1', 'AAACCATTCGCTCATT-1', 'AAACCCGCACCTCACG-1', 'AAACCCGCATCGGGAC-1']...
Variable names (genes): ['fgfr1op2', 'si:dkey-21h14.12', 'si:dkey-285e18.2', 'znf1114', 'si:dkey-21h14.10']...

=== adata.X (Primary Expression Matrix) ===
Shape: (160875, 24597) (cells: 160875, genes: 24597)
Type: dense
Data type: float32
Min: -6.042159
Max: 381.777557
Mean: -0.001693
Std: 0.995661
Non-zero values: 3957042375
Zero values: 0
Sparsity: 0.0000 (0.00%)

Top 5 genes (first 5 rows):
  fgfr1op2: min=-0.2939, max=7.9207, mean=-0.0006
  si:dkey-21h14.12: min=-0.0168, max=101.5491, mean=-0.0000
  si:dkey-285e18.2: min=-0.0079, max=182.2779, mean=0.0001
  znf1114: min=-0.1110, max=18.6632, mean=0.0004
  si:dkey-21h14.10: min=-0.0139, max=136.0288, mean=-0.0019

Top 5 cells (first 5 columns):
  AAACCAAAGCTTAACG-1: min=-2.3688, max=124.4653, mean=-0.0308
  AAACCAAAGTACGCCG-1: min=-3.3848, max=44.7331, mean=-0.0216
  AAACCATTCGCTCATT-1: min=-2.0438, max=22.0059, mean=0.0395
  AAACCCGCACCTCACG-1: min=-2.7271, max=35.6094, mean=0.0179
  AAACCCGCATCGGGAC-1: min=-2.5089, max=22.6247, mean=0.0050

adata.raw type: <class 'anndata._core.raw.Raw'>

=== adata.raw.X (Raw Expression Matrix) ===
Shape: (160875, 24597) (cells: 160875, genes: 24597)
Type: sparse
Data type: float32
Min: 0.000000
Max: 8.389225
Mean: 0.113462
Std: 0.443252
Non-zero values: 294862321
Zero values: 3662180054
Sparsity: 0.9255 (92.55%)

Top 5 genes (first 5 rows):
  fgfr1op2: min=0.0000, max=3.0251, mean=0.1080
  si:dkey-21h14.12: min=0.0000, max=1.9556, mean=0.0003
  si:dkey-285e18.2: min=0.0000, max=1.9756, mean=0.0001
  znf1114: min=0.0000, max=2.6757, mean=0.0159
  si:dkey-21h14.10: min=0.0000, max=2.1449, mean=0.0002

Top 5 cells (first 5 columns):
  AAACCAAAGCTTAACG-1: min=0.0000, max=4.7681, mean=0.0959
  AAACCAAAGTACGCCG-1: min=0.0000, max=5.3671, mean=0.0970
  AAACCATTCGCTCATT-1: min=0.0000, max=5.0932, mean=0.1313
  AAACCCGCACCTCACG-1: min=0.0000, max=6.1345, mean=0.1182
  AAACCCGCATCGGGAC-1: min=0.0000, max=5.2011, mean=0.1099

==================================================
COMPARISON: adata.X vs adata.raw.X
==================================================
Same shape: True
Number of different elements: 3957042359
Identical matrices: False

==================================================
CATEGORICAL OBSERVATIONS
==================================================
Found 4 categorical columns:

Column: sample (8 unique values)
Values: ['Zebra', 'TH115', 'TH44', 'TH54', 'TH55', 'TH56', 'TH57', 'TH71']
---
Column: renamed_samples (3 unique values)
Values: ['LD', 'NMDA', 'Control']
---
Column: leiden (92 unique values)
Values: ['3', '57', '13', '0', '56', '51', '12', '25', '35', '11', '87', '8', '30', '18', '28', '10', '33', '52', '61', '6', '5', '89', '16', '4', '42', '24', '44', '32', '72', '48', '70', '41', '9', '1', '45', '47', '93', '17', '86', '26', '29', '59', '94', '23', '15', '31', '78', '14', '34', '20', '27', '22', '92', '91', '54', '7', '85', '67', '21', '60', '37', '64', '2', '55', '19', '46', '75', '68', '43', '38', '36', '65', '63', '39', '53', '88', '90', '84', '83', '50', '69', '40', '82', '71', '77', '62', '49', '73', '74', '79', '81', '95']
---
Column: celltype (15 unique values)
Values: ['PR precursors', 'HC', 'MG', 'MGPC', 'Rod', 'BC', 'Cones', 'AC', 'RGC', 'RPE', 'Microglia_ImmuneCells', 'Oligodenrocyte', 'Perycites', 'Endothelial', 'Melanocyte']
---

ADDITIONAL METADATA
==================================================
adata.obs columns: ['sample', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_genes', 'renamed_samples', 'predicted_doublets', 'leiden', 'celltype']
adata.var columns: ['mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std']

No additional layers found.



## How to run Snakemake 

For dry run to check everythign before actual run:

    snakemake -j1 -p --configfile config.yaml -n 

For Actual run: 
    

##### Update the configfile 

The config files has the sample names with the suffix in samples.tsv 
and the annotations in annotations.txt. Edit both files correspondingly. You can edit the config file to change their names. 




