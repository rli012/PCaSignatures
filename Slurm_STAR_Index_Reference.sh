#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=STAR_index.log
#SBATCH -p batch

#sbatch Slurm_STAR_Index_Reference.sh

STAR=/rhome/rli012/bigdata/software/STAR-2.7.3a/bin/Linux_x86_64/STAR

# annotation
annotation=/rhome/rli012/bigdata/PCa/data/Reference/gencode.v32.annotation.gtf # gtf annotation file
genomeFa=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38.primary_assembly.genome.fa # fasta sequence file
genomeDir=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38/ #output directory

CPU=$SLURM_NTASKS

echo "Indexing..."

$STAR --runThreadN $CPU \
      --runMode genomeGenerate \
      --genomeDir $genomeDir \
      --genomeFastaFiles $genomeFa \
      --sjdbGTFfile $annotation \
      --sjdbOverhang 100

echo 'Done'
