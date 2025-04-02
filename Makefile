REF= /nfs/turbo/umms-thahoang/Tools/Fish11_genome/Drerio_genome11_unfiltered
#REF =/nfs/turbo/umms-thahoang/Tools/Fish11_genome/Drerio_genome11
#REF = /nfs/turbo/umms-thahoang/sherine/zebrafish2/Drerio_genome 
CELLRANGER = /nfs/turbo/umms-thahoang/Tools/cellranger-9.0.0/bin
FASTQS= /nfs/turbo/umms-thahoang/Share/Zebrafish/Lighdamage_old/THoa080921
TH65= /nfs/turbo/umms-thahoang/sherine/scanpy/TH65
TH66= /nfs/turbo/umms-thahoang/sherine/scanpy/TH66
TH67= /nfs/turbo/umms-thahoang/sherine/scanpy/TH67	
TH68= /nfs/turbo/umms-thahoang/sherine/scanpy/TH68
TH69= /nfs/turbo/umms-thahoang/sherine/scanpy/TH69
TH70= /nfs/turbo/umms-thahoang/sherine/scanpy/TH70



pattern='104094654|108234960|94171330|94279840|640099|640109|32582302|27434894|64173947|24777101|84196894|81559376|81645307|86025563'

out.txt:
	zgrep -F ${pattern} test.txt > out.txt

Danio_rerio.GRCz11.105.gtf:
	wget http://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz
	gunzip Danio_rerio.GRCz11.105.gtf.gz


Danio_rerio.GRCz11.dna.primary_assembly.fa:
	wget http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
	gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz


Danio_rerio.GRCz11.105.filtered.gtf:
	/nfs/turbo/umms-thahoang/sherine/tools/cellranger-7.2.0/cellranger mkgtf Danio_rerio.GRCz11.105.gtf Danio_rerio.GRCz11.105.filtered.gtf --attribute=gene_biotype:protein_coding


Drerio_genome:
	/nfs/turbo/umms-thahoang/sherine/tools/cellranger-7.2.0/cellranger mkref --genome=Drerio_genome --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa --genes=Danio_rerio.GRCz11.105.filtered.gtf

TH70/TH70_S72_L004_R2_001.fastq.gz: 
	find ${FASTQS}/TH65* -type f -name "*.fastq.gz" -exec ln -fs {} TH65/ \;
	find ${FASTQS}/TH66* -type f -name "*.fastq.gz" -exec ln -fs {} TH66/ \;
	find ${FASTQS}/TH67* -type f -name "*.fastq.gz" -exec ln -fs {} TH67/ \;
	find ${FASTQS}/TH68* -type f -name "*.fastq.gz" -exec ln -fs {} TH68/ \;
	find ${FASTQS}/TH69* -type f -name "*.fastq.gz" -exec ln -fs {} TH69/ \;
	find ${FASTQS}/TH70* -type f -name "*.fastq.gz" -exec ln -fs {} TH70/ \;

FTH65/outs: 
	${CELLRANGER}/cellranger count --id=FTH65 --transcriptome=${REF} --fastqs=${TH65} --expect-cells=10000 --sample=TH65 --create-bam=false

FTH66/outs: 
	${CELLRANGER}/cellranger count --id=FTH66 --transcriptome=${REF} --fastqs=${TH66} --expect-cells=10000 --sample=TH66 --create-bam=false

FTH67/outs: 
	${CELLRANGER}/cellranger count --id=FTH67 --transcriptome=${REF} --fastqs=${TH67} --expect-cells=10000 --sample=TH67 --create-bam=false

FTH68/outs:
	${CELLRANGER}/cellranger count --id=FTH68 --transcriptome=${REF} --fastqs=${TH68} --expect-cells=10000 --sample=TH68 --create-bam=false

FTH69/outs:
	${CELLRANGER}/cellranger count --id=FTH69 --transcriptome=${REF} --fastqs=${TH69} --expect-cells=10000 --sample=TH69 --create-bam=false

FTH70/outs:
	${CELLRANGER}/cellranger count --id=FTH70 --transcriptome=${REF} --fastqs=${TH70} --expect-cells=10000 --sample=TH70 --create-bam=false
