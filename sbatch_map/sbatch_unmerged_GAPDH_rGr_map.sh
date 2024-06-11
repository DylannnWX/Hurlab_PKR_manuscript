#!/bin/bash
#SBATCH -c 6                               # Request one core
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=16G                         # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)


set -e

Group=$1

#module load gcc/6.2.0 bowtie2/2.3.4.3
#Run this first: bowtie2-build Repbase.fasta Repbase
#Run this first: bowtie2-build rRNA.fasta rRNA


module load java/jdk-1.8u112 
module load trimmomatic/0.36

java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -threads 6 -phred33 ${Group}*1.fq.gz ${Group}*2.fq.gz ${Group}_R1_paired.fq ${Group}_R1_unpaired.fq ${Group}_R2_paired.fq ${Group}_R2_unpaired.fq ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:10

module load fastqc/0.11.9

#fastqc ${FASTQ1}_paired.fq
#fastqc ${FASTQ2}_paired.fq

module load gcc/6.2.0 bowtie2/2.3.4.3
bowtie2 -p 6 --end-to-end -x /home/xiw545/bowtie2/rRNA_human -1 ${Group}_R1_paired.fq -2 ${Group}_R2_paired.fq -S ${Group}_rRNA.sam --un-conc ./${Group}_unused_rRNA_reads
bowtie2 -p 6 --end-to-end -x /home/xiw545/bowtie2/Human_Gchr -1 ${Group}_unused_rRNA_reads.1 -2 ${Group}_unused_rRNA_reads.2 -S ${Group}_human.sam --un-conc ./${Group}_unmappable

module load gcc/6.2.0 samtools/1.9
samtools view -bS  ${Group}_rRNA.sam >  ${Group}_rRNA_unsorted.bam
samtools view -bS  ${Group}_human.sam >  ${Group}_human_unsorted.bam

samtools sort ${Group}_rRNA_unsorted.bam -o ${Group}_rRNA_sorted.bam
samtools sort ${Group}_human_unsorted.bam -o ${Group}_human_sorted.bam

samtools view -b -F 260 ${Group}_rRNA_sorted.bam > ${Group}_rRNA.bam
samtools view -b -F 260 ${Group}_human_sorted.bam > ${Group}_human.bam

rm ${Group}_unused_rRNA_reads.1 ${Group}_unused_rRNA_reads.2
rm ${Group}*unsorted.bam
rm ${Group}*sorted.bam
rm ${Group}*unpaired.fq
rm ${Group}*paired.fq
rm ${Group}*.sam

rRNA_count=$(samtools view -c -F 260 ${Group}_rRNA.bam)
Human_count=$(samtools view -c -F 260 ${Group}_human.bam)

NewFileName_rRNA="${Group}_rRNA_${rRNA_count}.bam"
NewFileName_human="${Group}_human_${Human_count}.bam"

mv "${Group}_rRNA.bam" "$NewFileName_rRNA"
mv "${Group}_human.bam" "$NewFileName_human"


#samtools index ${NewFileName_rRNA}
samtools index ${NewFileName_human}

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${NewFileName_human}> ${NewFileName_human}_rev1.bam
#samtools index ${NewFileName_human}_rev1.bam

samtools view -b -f 80 ${NewFileName_human} > ${NewFileName_human}_rev2.bam
#samtools index ${NewFileName_human}_rev2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${NewFileName_human}_minus.bam ${NewFileName_human}_rev1.bam ${NewFileName_human}_rev2.bam
samtools index ${NewFileName_human}_minus.bam

# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${NewFileName_human} > ${NewFileName_human}_fwd1.bam
#samtools index ${NewFileName_human}_fwd1.bam

samtools view -b -f 64 -F 16 ${NewFileName_human} > ${NewFileName_human}_fwd2.bam
#samtools index fwd2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${NewFileName_human}_plus.bam ${NewFileName_human}_fwd1.bam ${NewFileName_human}_fwd2.bam
samtools index ${NewFileName_human}_plus.bam

group1_plusG=$(samtools depth -r Gchr12:6534517-6538371 ${NewFileName_human}_plus.bam |  awk '{sum+=$3} END { print (sum+0.01)/(NR+1)}')
group1_minusG=$(samtools depth -r Gchr12:6534517-6538371 ${NewFileName_human}_minus.bam | awk '{sum+=$3} END { print (sum+0.01)/(NR+1)}')

NewFileName_human_GAPDH="${Group}_human_${Human_count}_plus_${group1_plusG}_minus_${group1_minusG}.bam"

mv "${NewFileName_human}" "$NewFileName_human_GAPDH"

rm ${NewFileName_human}_rev1.bam ${NewFileName_human}_rev2.bam ${NewFileName_human}_fwd1.bam ${NewFileName_human}_fwd2.bam ${NewFileName_human}_plus.bam ${NewFileName_human}_minus.bam ${NewFileName_human}_plus.bam.bai ${NewFileName_human}_minus.bam.bai ${NewFileName_human}.bai