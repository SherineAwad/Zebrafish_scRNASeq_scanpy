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

### `umap_Harmonyzebrafishes.png`

This UMAP plot shows zebrafish single-cell data **after batch effect correction** using **Harmony**.

By aligning shared biological structure across different batches or samples, Harmony improves the clustering and visualization. Cells are now grouped primarily by **biological similarity** rather than technical batch differences, enabling more accurate downstream analysis.

![umap_Harmonyzebrafishes.png](figures/umap_Harmonyzebrafishes.png)

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

### UMAP after clustering 

<img src="figures/umap_clusters.png" alt="Umap Annotations" width="90%">


## After re-Clustering annotation using marker genes

<img src="figures/Zebrafishes_annotations.png" alt="Umap Annotations" width="90%">
<img src="figures/Zebrafishes_annotationsON.png" alt="Umap Annotations ON" width="90%">


## Remove Celltype Cones_MG_MGPC_PostMitotic as likely doublet 

<img src="figures/Zebrafishes_Nannotations.png" alt="Umap Annotations" width="90%">
<img src="figures/Zebrafishes_NannotationsON.png" alt="Umap Annotations ON" width="90%">
 

## Cell ratio after increasing resolution and re-annotations 

<img src="figures/Restacked_bar_sample_by_celltype.png" width="90%">

## Umap for Control,LD, and NMDA merged 

<img src="figures/umap_merged_annotated_clustered_corrected_doubletRemoved_Zebrafishes.png" width="90%"> 

## Umap per sample 
<img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Control.png" width="33%"><img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_LD.png" width ="33%"><img src="figures/umap_annotated_clustered_corrected_doubletRemoved_Zebrafishes_NMDA.png" width="33%">


### UMAPs of Major Retinal Cell Subtypes

The following UMAP plots illustrate the clustering of subtypes within each major retinal cell type. These visualizations highlight transcriptional diversity captured by single-cell RNA sequencing.

#### Cones Photoreceptor Subtypes before batch correction 
<img src="figures/umapCones.png" width="70%"> 
<img src="figures/umap_merged_Cones.png" width ="70%">

<img src= "figures/umap_Cones_NMDA.png" width ="33%"><img src= "figures/umap_Cones_LD.png" width ="33%"><img src= "figures/umap_Cones_Control.png" width ="33%">


#### Cones Photoreceptor Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/Cones_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%"> 


#### Amacrine Cell (AC) Subtypes
<img src="figures/umapAC.png" width="70%">
<img src="figures/umap_merged_AC.png" width ="70%">

<img src= "figures/umap_AC_NMDA.png" width ="33%"><img src= "figures/umap_AC_LD.png" width ="33%"><img src= "figures/umap_AC_Control.png" width ="33%">


#### AC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/AC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%">


#### Retinal Ganglion Cell (RGC) Subtypes
<img src="figures/umapRGC.png" width="70%">
<img src="figures/umap_merged_RGC.png" width ="70%">

<img src= "figures/umap_RGC_NMDA.png" width ="33%"><img src= "figures/umap_RGC_LD.png" width ="33%"><img src= "figures/umap_RGC_Control.png" width ="33%">


#### RGC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/RGC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%">


#### Horizontal Cell (HC) Subtypes
<img src="figures/umapHC.png" width="70%">
<img src="figures/umap_merged_HC.png" width ="70%">

<img src= "figures/umap_HC_NMDA.png" width ="33%"><img src= "figures/umap_HC_LD.png" width ="33%"><img src= "figures/umap_HC_Control.png" width ="33%">


#### HC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/HC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%">


#### Bipolar Cell (BC) Subtypes
<img src="figures/umapBC.png" width="70%">
<img src="figures/umap_merged_BC.png" width ="70%">

<img src= "figures/umap_BC_NMDA.png" width ="33%"><img src= "figures/umap_BC_LD.png" width ="33%"><img src= "figures/umap_BC_Control.png" width ="33%">

#### BC Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/BC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/BC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/BC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/BC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/BC_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%">

#### MG Subtypes
<img src="figures/umapMG.png" width="70%">
<img src="figures/umap_merged_MG.png" width ="70%">

<img src= "figures/umap_MG_NMDA.png" width ="33%"><img src= "figures/umap_MG_LD.png" width ="33%"><img src= "figures/umap_MG_Control.png" width ="33%">


#### MG Subtypes after batch correction (STRICT PARAMETERS) 

<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony.png" width ="90%">

<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_LD.png" width ="33%"><img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_NMDA.png" width ="33%"><img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_Control.png" width ="33%">

<img src="figures/MG_annotated_clustered_corrected_doubletRemoved_Zebrafishes_Harmony_leiden.png" width ="90%">


## How to run Snakemake 

For dry run to check everythign before actual run:

    snakemake -j1 -p --configfile config.yaml -n 

For Actual run: 
    

##### Update the configfile 

The config files has the sample names with the suffix in samples.tsv 
and the annotations in annotations.txt. Edit both files correspondingly. You can edit the config file to change their names. 




