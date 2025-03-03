with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()

SUBSET = ['Cones', 'AC'] 

rule all:
         input:
            expand("{all}.h5ad", all= config['ALL']), 
            expand("renamed_{all}.h5ad", all=config['ALL']),
            expand("doubletRemoved_{all}.h5ad", all=config['ALL']),
            expand("corrected_{all}.h5ad", all=config['ALL']),
            expand("clustered_{all}.h5ad", all=config['ALL']), 
            #expand("annotated_{all}.h5ad", all=config['ALL']), 
            #expand("{subset}"_{all}.h5ad", all=config['ALL'], subset = SUBSET),
 
rule preprocess: 
        input:  
            expand("{sample}_filtered_feature_bc_matrix.h5", sample = samples) 
        output: 
          expand("{all}.h5ad", all= config['ALL']), 
        params: 
          samples = config['SAMPLES'],  
          name = config['ALL']
        shell: 
            """
           python preprocess.py {params.samples}  {params.name}  
           """ 
rule rename: 
       input:
           expand("{all}.h5ad", all= config['ALL']),
       output:
          expand("renamed_{all}.h5ad", all= config['ALL']),
       shell:
           """
           python preprocess.py {input} 
           """

rule remove_doublet: 
      input:
           expand("renamed_{all}.h5ad", all= config['ALL']),
      output:
          expand("doubletRemoved_{all}.h5ad", all= config['ALL']),
      shell:
           """
           python doublet.py {input} 
           """

rule batch: 
     input: 
         expand("doubletRemoved_{all}.h5ad", all=config['ALL'])
     output: 
         expand("corrected_{all}.h5ad", all=config['ALL'])
     shell: 
        """ 
        python batch.py {input} 
        """ 


rule cluster: 
       input:
          expand("corrected_{all}.h5ad", all=config['ALL']) 
       output:
          expand("clustered_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python cluster.py {input}
          """

rule annotate:
       input:
          expand("clustered_{all}.h5ad", all=config['ALL'])
       output:
          expand("annotated_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python annotate.py {input}
          """

rule subset: 
      input: 
         expand("annotated_{all}.h5ad", all=config['ALL'])
      output: 
         expand("{subset}_{all}.h5ad", all=config['ALL'], subset = SUBSET),
      shell: 
          """ 
          python subset.py {input} ConesSubtypes.txt" 
          """
      
