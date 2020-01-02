#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --output=FastQC.log
#SBATCH -p intel
#SBATCH --chdir=/bigdata/jialab/rli012/PCa/data/fromSRA/GSE54460/

fastqc=/bigdata/jialab/rli012/software/FastQC/fastqc

N=$SLURM_ARRAY_TASK_ID
CPU=$SLURM_NTASKS

FILE=`ls raw/SRR*\.fastq.gz | grep _1.fastq.gz | head -n $N | tail -n 1`
PREFIX=${FILE%_1.fastq.gz}
PREFIX=${PREFIX#raw/}

fq1=$FILE
fq2=${FILE/_1/_2}

#PREFIX=SRR2973290
#bam=${PREFIX}.bam
#fq1=${PREFIX}_1.fastq.gz
#fq2=${PREFIX}_2.fastq.gz

echo 'Start QC...'
echo $PREFIX

### QC ###

#FILES=`ls raw/SRR*\.fastq.gz`

$fastqc $fq1 --outdir=FastQC/
$fastqc $fq2 --outdir=FastQC/

echo 'Done!'
