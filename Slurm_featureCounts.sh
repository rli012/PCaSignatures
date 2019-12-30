#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=featureCounts.gencode.log
#SBATCH -p intel
#SBATCH --chdir=/bigdata/jialab/rli012/PCa/data/fromSRA/GSE54460/

featureCounts=/rhome/rli012/bigdata/software/subread-2.0.0-Linux-x86_64/bin/featureCounts

annotation=/rhome/rli012/bigdata/PCa/data/Reference/gencode.v32.annotation.gtf # gtf annotation file
genomeFa=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38.primary_assembly.genome.fa # fasta sequence file
genomeDir=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38/ #output directory

N=$SLURM_ARRAY_TASK_ID
CPU=$SLURM_NTASKS

BAMS=`ls BAM/*\.bam`
OUTPUT='GSE54460_featureCounts_Gencode.v32.txt'

### Quantification ###
echo 'Start Counting...'

$featureCounts -T $CPU --primary --ignoreDup -p -t exon -g gene_id -a $annotation -o count.tmp $BAMS
tail -n+2 count.tmp | cut --complement -f2-5 > $OUTPUT

echo 'Done!'
