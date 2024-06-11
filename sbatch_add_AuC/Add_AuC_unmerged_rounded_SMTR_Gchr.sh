#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=16G                         # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

subtract_pair=$1 group1=$2

module load gcc/6.2.0 python/2.7.12 samtools/1.9 deeptools/3.0.2 macs2/2.1.1.20160309

samtools merge -f ${group1}_unsorted.bam ${group1}*.bam
samtools sort ${group1}_unsorted.bam -o ${group1}_all.bam
samtools index ${group1}_all.bam
samtools idxstats ${group1}_all.bam | cut -f 1 | grep -v RNA45SN1 | grep -v 5s_rRNA | xargs samtools view -b ${group1}_all.bam > ${group1}_all_normalize.bam
rm ${group1}_unsorted.bam

#starting strand separation, this is unmerged

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${group1}_all.bam > ${group1}_rev1.bam
#samtools index ${group1}_rev1.bam

samtools view -b -f 80 ${group1}_all.bam > ${group1}_rev2.bam
#samtools index ${group1}_rev2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${group1}_all_plus.bam ${group1}_rev1.bam ${group1}_rev2.bam
#samtools index ${group1}_minus.bam

# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${group1}_all.bam > ${group1}_fwd1.bam
#samtools index ${group1}_fwd1.bam

samtools view -b -f 64 -F 16 ${group1}_all.bam > ${group1}_fwd2.bam
#samtools index fwd2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${group1}_all_minus.bam ${group1}_fwd1.bam ${group1}_fwd2.bam
#samtools index ${group1}_all_plus.bam 


rm ${group1}_fwd1.bam
rm ${group1}_fwd2.bam
rm ${group1}_rev1.bam
rm ${group1}_rev2.bam

#ending strand separation

module load gcc/9.2.0 bedtools/2.30.0

bedtools genomecov -ibam ${group1}_all.bam -bg > ${group1}_all.bedgraph
bedtools genomecov -ibam ${group1}_all_plus.bam -bg > ${group1}_all_plus.bedgraph
bedtools genomecov -ibam ${group1}_all_minus.bam -bg > ${group1}_all_minus.bedgraph

module load gcc/6.2.0 python/2.7.12 samtools/1.9 deeptools/3.0.2 macs2/2.1.1.20160309

sed -i 's/^Gchr/chr/' ${group1}_all.bedgraph
sed -i 's/^Gchr/chr/' ${group1}_all_plus.bedgraph
sed -i 's/^Gchr/chr/' ${group1}_all_minus.bedgraph

group1c=$(samtools view -c ${group1}_all_normalize.bam)

rm ${group1}_all_normalize.bam

awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all.bedgraph > ${group1}_all_PM.bedgraph
awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all_plus.bedgraph > ${group1}_all_plus_PM.bedgraph
awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all_minus.bedgraph > ${group1}_all_minus_PM.bedgraph

awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_PM.bedgraph > ${group1}_all_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_plus_PM.bedgraph > ${group1}_all_plus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_minus_PM.bedgraph > ${group1}_all_minus_PM_rounded.bedgraph

module load gcc/9.2.0
module load bedtools/2.30.0

bedtools intersect -a ${subtract_pair}.bed -b ${group1}_all_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${group1}_AuC.bed
bedtools intersect -a ${subtract_pair}_plus.bed -b ${group1}_all_plus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${group1}_plus_AuC.bed
bedtools intersect -a ${subtract_pair}_minus.bed -b ${group1}_all_minus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${group1}_minus_AuC.bed

rm ${group1}*_rounded.bedgraph

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7 }' ${group1}_AuC.bed > ${group1}_AuC_dropped.bed
sort -n -r -k 5 ${group1}_AuC_dropped.bed > ${group1}_AuC_dropped_sorted.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${group1}_AuC" | cat - ${group1}_AuC_dropped_sorted.bed > ${group1}_AuC_dropped_sorted_header.bed
awk '{ gsub(/[\t]/,","); print }' ${group1}_AuC_dropped_sorted_header.bed > ${subtract_pair}_${group1}_AuC.csv
rm ${group1}_AuC.bed
rm ${group1}_AuC_dropped.bed
rm ${group1}_AuC_dropped_sorted.bed
rm ${group1}_AuC_dropped_sorted_header.bed

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7 }' ${group1}_plus_AuC.bed > ${group1}_plus_AuC_dropped.bed
sort -n -r -k 5 ${group1}_plus_AuC_dropped.bed > ${group1}_plus_AuC_dropped_sorted.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${group1}_plus_AuC" | cat - ${group1}_plus_AuC_dropped_sorted.bed > ${group1}_plus_AuC_dropped_sorted_header.bed
awk '{ gsub(/[\t]/,","); print }' ${group1}_plus_AuC_dropped_sorted_header.bed > ${subtract_pair}_${group1}_plus_AuC.csv
rm ${group1}_plus_AuC.bed
rm ${group1}_plus_AuC_dropped.bed
rm ${group1}_plus_AuC_dropped_sorted.bed
rm ${group1}_plus_AuC_dropped_sorted_header.bed

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7 }' ${group1}_minus_AuC.bed > ${group1}_minus_AuC_dropped.bed
sort -n -r -k 5 ${group1}_minus_AuC_dropped.bed > ${group1}_minus_AuC_dropped_sorted.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${group1}_minus_AuC" | cat - ${group1}_minus_AuC_dropped_sorted.bed > ${group1}_minus_AuC_dropped_sorted_header.bed
awk '{ gsub(/[\t]/,","); print }' ${group1}_minus_AuC_dropped_sorted_header.bed > ${subtract_pair}_${group1}_minus_AuC.csv
rm ${group1}_minus_AuC.bed
rm ${group1}_minus_AuC_dropped.bed
rm ${group1}_minus_AuC_dropped_sorted.bed
rm ${group1}_minus_AuC_dropped_sorted_header.bed


