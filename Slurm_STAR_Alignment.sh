#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=Alignment.log
#SBATCH -p intel
#SBATCH --chdir=/bigdata/jialab/rli012/PCa/data/fromSRA/GSE54460/

# sbatch --array 1-106 Slurm_STAR_Alignment.sh

STAR=/rhome/rli012/bigdata/software/STAR-2.7.3a/bin/Linux_x86_64/STAR
samtools=/rhome/rli012/bigdata/software/samtools-1.9/bin/samtools

annotation=/rhome/rli012/bigdata/PCa/data/Reference/gencode.v32.annotation.gtf # gtf annotation file
genomeFa=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38.primary_assembly.genome.fa # fasta sequence file
genomeDir=/rhome/rli012/bigdata/PCa/data/Reference/GRCh38/ #output directory

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

echo 'Start Alignment...'
echo $PREFIX

### Alignment ###

$STAR --runThreadN $CPU \
      --genomeDir $genomeDir \
      --twopassMode Basic \
      --readFilesIn $fq1 $fq2 \
      --readFilesCommand zcat \
      --outSAMtype BAM Unsorted \
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
      --clip5pNbases 5 \
      --clip3pNbases 5 \
      --limitBAMsortRAM 19732153018 \
      --outFileNamePrefix alignment/${PREFIX}

$samtools sort -@ $CPU alignment/${PREFIX}Aligned.out.bam -T alignment/${PREFIX} -o alignment/${PREFIX}.bam
rm alignment/${PREFIX}Aligned.out.bam

echo 'Done!'
