#!/bin/bash

script_path=$(pwd)

index=$script_path/index
input=$script_path/input

kraken2_rf=$index/minikraken2_v2_8GB_201904_UPDATE
homer_annotation=$index/gencode.v37.chr_patch_hapl_scaff.annotation.gtf
bwa_rf=$index/hg38.fa

module load miniconda
#conda create --yes -n 03162021 star fastqc multiqc samtools bwa kraken2 bracken subread homer bedtools
conda activate 03162021

######################################paired_end############################################

function fast_qc {
  cd $input
  for f in *_1.fastq.gz; do
      r1=${f};
      r2=${f/_1.fastq.gz/}_2.fastq.gz;
      fastqc -o $fastqc_output $r1 $r2;
  done
  echo "fastqc done!"
  multiqc $fastqc_output -o $multiqc_output
  echo "multiqc done!"
}

function bwa_mapping {
  cd $input
  for f in *_1.fastq.gz; do
      r1=${f};
      r2=${f/_1.fastq.gz/}_2.fastq.gz;
      echo $r1 $r2;
      r3=${f/_1.fastq.gz/}.sam;
      bwa mem -t 6 $bwa_rf $r1 $r2 > $mappings/$r3;
  done
  echo "bwa_mapping done!"
}

function mapping_post_processing {
  cd $mappings
  for file in *.sam; do
     a=$file;
     b=${file/.sam/}.fixmate.bam;
     c=${file/.sam/}.sorted.bam;
     d=${file/.sam/}.sorted.dedup.bam;
     e=${file/.sam/}.sorted.dedup.q20.bam;
     d1=${file/.sam/}.tmps1.bam;
     d2=${file/.sam/}.tmps2.bam;
     d3=${file/.sam/}.tmps3.bam;
     d4=${file/.sam/}.tmps.bam;
     f=${file/.sam/}.sorted.unmapped.bam;
     g=${file/.sam/}.sorted.unmapped.R1.fastq.gz;
     h=${file/.sam/}.sorted.unmapped.R2.fastq.gz;
     i=${file/.sam/};
     samtools sort -n -O sam $mappings/$a | samtools fixmate -m -O bam - $mappings/$b
     #rm $mappings/$a
     ##sort
     samtools sort -O bam -o $mappings/$c $mappings/$b
     rm $mappings/$b
     ##mark duplicates
     samtools markdup -r -S $mappings/$c $mappings/$d
     #rm $mappings/$c
     ##extract q20 mappers
     samtools view -h -b -q 20 $mappings/$d > $mappings/$e;
     ##extract unmapped
     samtools view -u  -f 4 -F 264 $mappings/$d  > $mappings/$d1
     samtools view -u -f 8 -F 260 $mappings/$d  > $mappings/$d2
     samtools view -u -f 12 -F 256 $mappings/$d > $mappings/$d3
     samtools merge -u $mappings/$d4 $mappings/$d1 $mappings/$d2 $mappings/$d3
     samtools merge -u $mappings/$d4 $mappings/$d1 $mappings/$d2 $mappings/$d3
     var1=$(wc -l $mappings/$d1 | cut -f1 -d' ')
     var2=$(wc -l $mappings/$d2 | cut -f1 -d' ')
     var3=$(wc -l $mappings/$d3 | cut -f1 -d' ')
     echo $'Unmapped_summary:' $'\nSample =' $i $'\nUnmapped read whose mate is mapped =' $var1 $'\nMapped read whos mate is unmapped =' $var2 $'\nBoth reads of the pair are unmapped =' $var3
     #rm $mappings/$d $mappings/$d1 $mappings/$d2 $mappings/$d3
     samtools sort -n $mappings/$d4 > $mappings/$f
     rm $mappings/$d4
     ##covert to fastq
     bedtools bamtofastq -i $mappings/$f -fq $mappings/$g -fq2 $mappings/$h
     ##delete not needed files
     #rm $mappings/$f
  done
  echo "mapping_post_processing done!"
}

function kraken_2 {
  cd $mappings
  for file in *.sorted.unmapped.R1.fastq.gz; do
     a=$file;
     b=${file/.sorted.unmapped.R1.fastq.gz/}.sorted.unmapped.R2.fastq.gz;
     c=${file/.sorted.unmapped.R1.fastq.gz/}.output.txt;
     d=${file/.sorted.unmapped.R1.fastq.gz/}.report.txt;
     kraken2 --db $kraken2_rf --threads 20 --minimum-base-quality 20 --output $kraken2_output/$c --report $kraken2_output_report/$d --confidence 0.1 --paired $a $b;
  done
  echo "kraken2 done!"
}

function homer_analyze {
  cd $mappings
  for f in *.sorted.bam; do
      r1=${f};
      r2=${f/.sorted.bam/};
      makeTagDirectory $r2/ $r1
      # raw counts
      analyzeRepeats.pl $homer_annotation hg38 -count exons -d $r2/ -noadj > $r2.homer.counts
      cp $r2.homer.counts $homer_counts 
      # calculate rpkm
      analyzeRepeats.pl $homer_annotation hg38 -count exons -d $r2/ -rpkm > $r2.homer.rpkm
      cp $r2.homer.rpkm $homer_rpkm
      # calculate rpm/cpm
      analyzeRepeats.pl $homer_annotation hg38 -count exons -d $r2/ -norm 1e7 > $r2.homer.rpm
      cp $r2.homer.rpm $homer_rpm
  done
  echo "homer_analyze done!"
}

######################################single_end############################################

function fast_qc_s {
  cd $input
  for f in *.fastq.gz; do
      r1=${f};
      fastqc -o $fastqc_output $r1;
  done
  echo "fastqc done!"
  multiqc $fastqc_output -o $multiqc_output
  echo "multiqc done!"
}

function bwa_mapping_s {
  cd $input
  for f in *.fastq.gz; do
      r1=${f};
      echo $r1 $r2;
      r3=${f/.fastq.gz/}.sam;
      bwa mem -t 6 $bwa_rf $r1 > $mappings/$r3;
  done
  echo "bwa_mapping done!"
}

function mapping_post_processing_s {
  cd $mappings
  for file in *.sam; do
     a=$file;
     b=${file/.sam/}.fixmate.bam;
     c=${file/.sam/}.sorted.bam;
     d=${file/.sam/}.sorted.dedup.bam;
     e=${file/.sam/}.sorted.dedup.q20.bam;
     d1=${file/.sam/}.tmps1.bam;
     d2=${file/.sam/}.tmps2.bam;
     d3=${file/.sam/}.tmps3.bam;
     d4=${file/.sam/}.tmps.bam;
     f=${file/.sam/}.sorted.unmapped.bam;
     g=${file/.sam/}.sorted.unmapped.fastq;

     i=${file/.sam/};
     samtools sort -n -O sam $mappings/$a | samtools fixmate -m -O bam - $mappings/$b
     #rm $mappings/$a
     ##sort
     samtools sort -O bam -o $mappings/$c $mappings/$b
     rm $mappings/$b
     ##mark duplicates
     samtools markdup -r -S $mappings/$c $mappings/$d
     #rm $mappings/$c
     ##extract q20 mappers
     samtools view -h -b -q 20 $mappings/$d > $mappings/$e;
     ##extract unmapped
     samtools view -u  -f 4 -F 264 $mappings/$d  > $mappings/$d1
     samtools view -u -f 8 -F 260 $mappings/$d  > $mappings/$d2
     samtools view -u -f 12 -F 256 $mappings/$d > $mappings/$d3
     samtools merge -u $mappings/$d4 $mappings/$d1 $mappings/$d2 $mappings/$d3
     samtools merge -u $mappings/$d4 $mappings/$d1 $mappings/$d2 $mappings/$d3
     var1=$(wc -l $mappings/$d1 | cut -f1 -d' ')
     var2=$(wc -l $mappings/$d2 | cut -f1 -d' ')
     var3=$(wc -l $mappings/$d3 | cut -f1 -d' ')
     echo $'Unmapped_summary:' $'\nSample =' $i $'\nUnmapped read whose mate is mapped =' $var1 $'\nMapped read whos mate is unmapped =' $var2 $'\nBoth reads of the pair are unmapped =' $var3
     #rm $mappings/$d $mappings/$d1 $mappings/$d2 $mappings/$d3
     samtools sort -n $mappings/$d4 > $mappings/$f
     rm $mappings/$d4
     ##covert to fastq
     bedtools bamtofastq -i $mappings/$f -fq $mappings/$g
     ##delete not needed files
     #rm $mappings/$f
  done
  echo "mapping_post_processing done!"
}

function kraken_2_s {
  cd $mappings
  for file in *.sorted.unmapped.fastq; do
     a=$file;

     c=${file/.sorted.unmapped.fastq/}.output.txt;
     d=${file/.sorted.unmapped.fastq/}.report.txt;
     kraken2 --db $kraken2_rf --threads 20 --minimum-base-quality 20 --output $kraken2_output/$c --report $kraken2_output_report/$d --confidence 0.1 $a;
  done
  echo "kraken2 done!"
}


#######################################main_function###############################################
function singleEnd {
  mkdir $script_path/output $script_path/output/fastqc_result $script_path/output/fastqc_result/fastqc_output $script_path/output/fastqc_result/multiqc_output $script_path/output/mappings $script_path/output/kraken2_result $script_path/output/kraken2_result/kraken2_output $script_path/output/kraken2_result/kraken2_output_report $script_path/output/homer_result $script_path/output/homer_result/homer_counts $script_path/output/homer_result/homer_rpkm $script_path/output/homer_result/homer_rpm

  fastqc_output=$script_path/output/fastqc_result/fastqc_output
  multiqc_output=$script_path/output/fastqc_result/multiqc_output
  mappings=$script_path/output/mappings
  kraken2_output=$script_path/output/kraken2_result/kraken2_output
  kraken2_output_report=$script_path/output/kraken2_result/kraken2_output_report
  homer_counts=$script_path/output/homer_result/homer_counts
  homer_rpkm=$script_path/output/homer_result/homer_rpkm
  homer_rpm=$script_path/output/homer_result/homer_rpm

  fast_qc_s
  bwa_mapping_s
  mapping_post_processing_s
  kraken_2_s
  homer_analyze
}

function pairedEnd {
    mkdir $script_path/output $script_path/output/fastqc_result $script_path/output/fastqc_result/fastqc_output $script_path/output/fastqc_result/multiqc_output $script_path/output/mappings $script_path/output/kraken2_result $script_path/output/kraken2_result/kraken2_output $script_path/output/kraken2_result/kraken2_output_report $script_path/output/homer_result $script_path/output/homer_result/homer_counts $script_path/output/homer_result/homer_rpkm $script_path/output/homer_result/homer_rpm

  fastqc_output=$script_path/output/fastqc_result/fastqc_output
  multiqc_output=$script_path/output/fastqc_result/multiqc_output
  mappings=$script_path/output/mappings
  kraken2_output=$script_path/output/kraken2_result/kraken2_output
  kraken2_output_report=$script_path/output/kraken2_result/kraken2_output_report
  homer_counts=$script_path/output/homer_result/homer_counts
  homer_rpkm=$script_path/output/homer_result/homer_rpkm
  homer_rpm=$script_path/output/homer_result/homer_rpm

  fast_qc
  bwa_mapping
  mapping_post_processing
  kraken_2
  homer_analyze
}
function pairedEndKraken2 {
  mkdir $script_path/output $script_path/output/kraken2_result $script_path/output/kraken2_result/kraken2_output $script_path/output/kraken2_result/kraken2_output_report

  kraken2_output=$script_path/output/kraken2_result/kraken2_output
  kraken2_output_report=$script_path/output/kraken2_result/kraken2_output_report
  
  cd $input
  for file in *.sorted.unmapped.R1.fastq.gz; do
     a=$file;
     b=${file/.sorted.unmapped.R1.fastq.gz/}.sorted.unmapped.R2.fastq.gz;
     c=${file/.sorted.unmapped.R1.fastq.gz/}.output.txt;
     d=${file/.sorted.unmapped.R1.fastq.gz/}.report.txt;
     kraken2 --db $kraken2_rf --threads 20 --minimum-base-quality 20 --output $kraken2_output/$c --report $kraken2_output_report/$d --confidence 0.1 --paired $a $b;
  done
  echo "kraken2 done!"
}

function singleEndKraken2 {
  mkdir $script_path/output $script_path/output/kraken2_result $script_path/output/kraken2_result/kraken2_output $script_path/output/kraken2_result/kraken2_output_report

  kraken2_output=$script_path/output/kraken2_result/kraken2_output
  kraken2_output_report=$script_path/output/kraken2_result/kraken2_output_report

  cd $input
  for file in *.sorted.unmapped.fastq.gz; do
     a=$file;

     c=${file/.sorted.unmapped.fastq.gz/}.output.txt;
     d=${file/.sorted.unmapped.fastq.gz/}.report.txt;
     kraken2 --db $kraken2_rf --threads 20 --minimum-base-quality 20 --output $kraken2_output/$c --report $kraken2_output_report/$d --confidence 0.1 $a;
  done
  echo "kraken2 done!"
}

"$@"

conda deactivate
