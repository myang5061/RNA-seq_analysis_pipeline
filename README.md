# RNA-seq_analysis_pipeline
Integrated fastqc multiqc samtools bwa kraken2 bracken subread homer and bedtools
Before running the script, please paste index folder in current path and create a folder named "input".

The first time you use the pipeline, please run the following command to configure the environment:
conda create --yes -n 03162021 star fastqc multiqc samtools bwa kraken2 bracken subread homer bedtools

running command:
./mypipeline_test.sh singleEndKraken2
input files should be named as *.sorted.unmapped.fastq.gz

./mypipeline_test.sh pairedEndKraken2
input files should be named as *.sorted.unmapped.R1.fastq.gz *.sorted.unmapped.R1.fastq.gz

./mypipeline_test.sh singleEnd
input files should be named as *.fastq.gz

./mypipeline_test.sh pairedEnd
input files should be named as *_1.fastq.gz *_2.fastq.gz


Running error fix solution:
File "/gpfs/loomis/project/gerstein/my435/conda_envs/03162021/lib/python3.9/site-packages/networkx/algorithms/dag.py", line 23, in <module>
    from fractions import gcd
change this line with "import math"
and change line 467 "g = gcd(g, levels[u] - levels[v] + 1)" to "g = math.gcd(g, levels[u] - levels[v] + 1)"
