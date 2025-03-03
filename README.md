# Zebrafish Project 

We analyse Zebrafish LD, NMDA and control scRNASeq samples using scanpy. 

###### How to run Snakemake 

For dry run to check everythign before actual run:

    snakemake -j1 -p --configfile config.yaml -n 

For Actual run: 
    
    snakemake -j1 -p --configfile config.yaml

##### Update the configfile 

The config files has the sample names with the suffix in samples.tsv 
and the annotations in annotations.txt. Edit both files correspondingly. You can edit the config file to change their names. 

