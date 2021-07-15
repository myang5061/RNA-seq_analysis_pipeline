# RNA-seq_analysis_pipeline
Integrated fastqc multiqc samtools bwa kraken2 bracken subread homer and bedtools
Before running the script, please paste index folder in current path and create a folder named "input".

The first time you use the pipeline, please run the following command to configure the environment:
conda create --yes -n 03162021 star fastqc multiqc samtools bwa kraken2 bracken subread homer bedtools

## running command:
./mypipeline_test.sh singleEndKraken2
input files should be named as *.sorted.unmapped.fastq.gz

./mypipeline_test.sh pairedEndKraken2
input files should be named as *.sorted.unmapped.R1.fastq.gz *.sorted.unmapped.R1.fastq.gz

./mypipeline_test.sh singleEnd
input files should be named as *.fastq.gz

./mypipeline_test.sh pairedEnd
input files should be named as *_1.fastq.gz *_2.fastq.gz

## Before you run, you need to add the following files to the index folder：
hg38.fa
hg38.fa.amb
hg38.fa.ann
hg38.fa.bwt
hg38.fa.pac
hg38.fa.sa
gencode.v37.chr_patch_hapl_scaff.annotation.gtf  	 
minikraken2_v2_8GB_201904_UPDATE

1. for the hg38.fa file please download from http://hgdownload.soe.ucsc.edu/downloads.html#human
You can get hg38.fa.XXX files by running:
bwa mem hg38.fa
2. for the gencode.v37.chr_patch_hapl_scaff.annotation.gtf please download from https://www.gencodegenes.org/human/
3. for the minikraken2_v2_8GB_201904_UPDATE folder please download from https://www.ccb.jhu.edu/software/kraken2/index.shtml?t=downloads

## Running error fix solution:
File "/xxxx/conda_envs/03162021/lib/python3.9/site-packages/networkx/algorithms/dag.py", line 23, in <module>
    from fractions import gcd
change this line with "import math"
and change line 467 "g = gcd(g, levels[u] - levels[v] + 1)" to "g = math.gcd(g, levels[u] - levels[v] + 1)"
