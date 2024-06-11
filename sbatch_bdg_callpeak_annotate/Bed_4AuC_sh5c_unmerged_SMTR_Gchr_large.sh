#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=64G                         # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

subtract_pair=$1 group1=$2 group2=$3 group3=$4 group4=$5 Macs2cutoff=$6 subtract_pair2=$7

module load gcc/6.2.0 python/2.7.12 samtools/1.9 deeptools/3.0.2 macs2/2.1.1.20160309

samtools merge -f ${group1}_unsorted.bam ${group1}*.bam
samtools sort ${group1}_unsorted.bam -o ${group1}_all.bam
samtools index ${group1}_all.bam
samtools idxstats ${group1}_all.bam | cut -f 1 | grep -v RNA45SN1 | grep -v 5s_rRNA | xargs samtools view -b ${group1}_all.bam > ${group1}_all_normalize.bam
rm ${group1}_unsorted.bam

samtools merge -f ${group2}_unsorted.bam ${group2}*.bam
samtools sort ${group2}_unsorted.bam -o ${group2}_all.bam
samtools index ${group2}_all.bam
samtools idxstats ${group2}_all.bam | cut -f 1 | grep -v RNA45SN1 | grep -v 5s_rRNA | xargs samtools view -b ${group2}_all.bam > ${group2}_all_normalize.bam
rm ${group2}_unsorted.bam

samtools merge -f ${group3}_unsorted.bam ${group3}*.bam
samtools sort ${group3}_unsorted.bam -o ${group3}_all.bam
samtools index ${group3}_all.bam
samtools idxstats ${group3}_all.bam | cut -f 1 | grep -v RNA45SN1 | grep -v 5s_rRNA | xargs samtools view -b ${group3}_all.bam > ${group3}_all_normalize.bam
rm ${group3}_unsorted.bam

samtools merge -f ${group4}_unsorted.bam ${group4}*.bam
samtools sort ${group4}_unsorted.bam -o ${group4}_all.bam
samtools index ${group4}_all.bam
samtools idxstats ${group4}_all.bam | cut -f 1 | grep -v RNA45SN1 | grep -v 5s_rRNA | xargs samtools view -b ${group4}_all.bam > ${group4}_all_normalize.bam
rm ${group4}_unsorted.bam

#starting strand separation, this is unmerged SMTR

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
#samtools index ${group1}_plus.bam

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
#samtools index ${group1}_all_minus.bam 

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${group2}_all.bam > ${group2}_rev1.bam
#samtools index ${group2}_rev1.bam

samtools view -b -f 80 ${group2}_all.bam > ${group2}_rev2.bam
#samtools index ${group2}_rev2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${group2}_all_plus.bam ${group2}_rev1.bam ${group2}_rev2.bam
#samtools index ${group2}_all_plus.bam

# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${group2}_all.bam > ${group2}_fwd1.bam
#samtools index ${group2}_fwd1.bam

samtools view -b -f 64 -F 16 ${group2}_all.bam > ${group2}_fwd2.bam
#samtools index fwd2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${group2}_all_minus.bam ${group2}_fwd1.bam ${group2}_fwd2.bam
#samtools index ${group2}_all_minus.bam

rm ${group1}_fwd1.bam
rm ${group1}_fwd2.bam
rm ${group1}_rev1.bam
rm ${group1}_rev2.bam
rm ${group2}_fwd1.bam
rm ${group2}_fwd2.bam
rm ${group2}_rev1.bam
rm ${group2}_rev2.bam

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${group3}_all.bam > ${group3}_rev1.bam
#samtools index ${group3}_rev1.bam

samtools view -b -f 80 ${group3}_all.bam > ${group3}_rev2.bam
#samtools index ${group3}_rev2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${group3}_all_plus.bam ${group3}_rev1.bam ${group3}_rev2.bam
#samtools index ${group3}_plus.bam

# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${group3}_all.bam > ${group3}_fwd1.bam
#samtools index ${group3}_fwd1.bam

samtools view -b -f 64 -F 16 ${group3}_all.bam > ${group3}_fwd2.bam
#samtools index fwd2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${group3}_all_minus.bam ${group3}_fwd1.bam ${group3}_fwd2.bam
#samtools index ${group3}_all_minus.bam 

# Reverse strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${group4}_all.bam > ${group4}_rev1.bam
#samtools index ${group4}_rev1.bam

samtools view -b -f 80 ${group4}_all.bam > ${group4}_rev2.bam
#samtools index ${group4}_rev2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${group4}_all_plus.bam ${group4}_rev1.bam ${group4}_rev2.bam
#samtools index ${group4}_all_plus.bam

# Forward strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${group4}_all.bam > ${group4}_fwd1.bam
#samtools index ${group4}_fwd1.bam

samtools view -b -f 64 -F 16 ${group4}_all.bam > ${group4}_fwd2.bam
#samtools index fwd2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${group4}_all_minus.bam ${group4}_fwd1.bam ${group4}_fwd2.bam
#samtools index ${group4}_all_minus.bam

rm ${group3}_fwd1.bam
rm ${group3}_fwd2.bam
rm ${group3}_rev1.bam
rm ${group3}_rev2.bam
rm ${group4}_fwd1.bam
rm ${group4}_fwd2.bam
rm ${group4}_rev1.bam
rm ${group4}_rev2.bam

#ending strand separation

module load gcc/9.2.0 bedtools/2.30.0

bedtools genomecov -ibam ${group1}_all.bam -bg > ${group1}_all.bedgraph
bedtools genomecov -ibam ${group1}_all_plus.bam -bg > ${group1}_all_plus.bedgraph
bedtools genomecov -ibam ${group1}_all_minus.bam -bg > ${group1}_all_minus.bedgraph

bedtools genomecov -ibam ${group2}_all.bam -bg > ${group2}_all.bedgraph
bedtools genomecov -ibam ${group2}_all_plus.bam -bg > ${group2}_all_plus.bedgraph
bedtools genomecov -ibam ${group2}_all_minus.bam -bg > ${group2}_all_minus.bedgraph

bedtools genomecov -ibam ${group3}_all.bam -bg > ${group3}_all.bedgraph
bedtools genomecov -ibam ${group3}_all_plus.bam -bg > ${group3}_all_plus.bedgraph
bedtools genomecov -ibam ${group3}_all_minus.bam -bg > ${group3}_all_minus.bedgraph

bedtools genomecov -ibam ${group4}_all.bam -bg > ${group4}_all.bedgraph
bedtools genomecov -ibam ${group4}_all_plus.bam -bg > ${group4}_all_plus.bedgraph
bedtools genomecov -ibam ${group4}_all_minus.bam -bg > ${group4}_all_minus.bedgraph

sed -i 's/^Gchr/chr/' ${group1}_all.bedgraph
sed -i 's/^Gchr/chr/' ${group1}_all_plus.bedgraph
sed -i 's/^Gchr/chr/' ${group1}_all_minus.bedgraph
sed -i 's/^Gchr/chr/' ${group2}_all.bedgraph
sed -i 's/^Gchr/chr/' ${group2}_all_plus.bedgraph
sed -i 's/^Gchr/chr/' ${group2}_all_minus.bedgraph
sed -i 's/^Gchr/chr/' ${group3}_all.bedgraph
sed -i 's/^Gchr/chr/' ${group3}_all_plus.bedgraph
sed -i 's/^Gchr/chr/' ${group3}_all_minus.bedgraph
sed -i 's/^Gchr/chr/' ${group4}_all.bedgraph
sed -i 's/^Gchr/chr/' ${group4}_all_plus.bedgraph
sed -i 's/^Gchr/chr/' ${group4}_all_minus.bedgraph

module load gcc/6.2.0 python/2.7.12 samtools/1.9 deeptools/3.0.2 macs2/2.1.1.20160309

group1c=$(samtools view -c ${group1}_all_normalize.bam)
group2c=$(samtools view -c ${group2}_all_normalize.bam)
group3c=$(samtools view -c ${group3}_all_normalize.bam)
group4c=$(samtools view -c ${group4}_all_normalize.bam)
rm ${group1}_all_normalize.bam
rm ${group2}_all_normalize.bam
rm ${group3}_all_normalize.bam
rm ${group4}_all_normalize.bam

awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all.bedgraph > ${group1}_all_PM.bedgraph
awk -v varA="$group2c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group2}_all.bedgraph > ${group2}_all_PM.bedgraph
awk -v varA="$group3c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group3}_all.bedgraph > ${group3}_all_PM.bedgraph
awk -v varA="$group4c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group4}_all.bedgraph > ${group4}_all_PM.bedgraph

awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all_plus.bedgraph > ${group1}_all_plus_PM.bedgraph
awk -v varA="$group2c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group2}_all_plus.bedgraph > ${group2}_all_plus_PM.bedgraph
awk -v varA="$group3c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group3}_all_plus.bedgraph > ${group3}_all_plus_PM.bedgraph
awk -v varA="$group4c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group4}_all_plus.bedgraph > ${group4}_all_plus_PM.bedgraph

awk -v varA="$group1c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group1}_all_minus.bedgraph > ${group1}_all_minus_PM.bedgraph
awk -v varA="$group2c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group2}_all_minus.bedgraph > ${group2}_all_minus_PM.bedgraph
awk -v varA="$group3c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group3}_all_minus.bedgraph > ${group3}_all_minus_PM.bedgraph
awk -v varA="$group4c" -v varB=1000000 -v OFS='\t' '{$4=$4*varB/varA; print}' ${group4}_all_minus.bedgraph > ${group4}_all_minus_PM.bedgraph

macs2 bdgcmp -t ${group1}_all_PM.bedgraph -c ${group2}_all_PM.bedgraph -o ${subtract_pair}.bedgraph -m subtract
macs2 bdgcmp -t ${group1}_all_plus_PM.bedgraph -c ${group2}_all_plus_PM.bedgraph -o ${subtract_pair}_plus.bedgraph -m subtract
macs2 bdgcmp -t ${group1}_all_minus_PM.bedgraph -c ${group2}_all_minus_PM.bedgraph -o ${subtract_pair}_minus.bedgraph -m subtract

macs2 bdgcmp -t ${group3}_all_PM.bedgraph -c ${group4}_all_PM.bedgraph -o ${subtract_pair2}.bedgraph -m subtract
macs2 bdgcmp -t ${group3}_all_plus_PM.bedgraph -c ${group4}_all_plus_PM.bedgraph -o ${subtract_pair2}_plus.bedgraph -m subtract
macs2 bdgcmp -t ${group3}_all_minus_PM.bedgraph -c ${group4}_all_minus_PM.bedgraph -o ${subtract_pair2}_minus.bedgraph -m subtract

macs2 bdgpeakcall -l 30 -i ${subtract_pair}.bedgraph -c ${Macs2cutoff} -o ${subtract_pair}_MACS2_raw.bed
awk 'NR>1' ${subtract_pair}_MACS2_raw.bed > ${subtract_pair}_rawA.bed
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' ${subtract_pair}_rawA.bed > ${subtract_pair}.bed
rm ${subtract_pair}_rawA.bed
rm ${subtract_pair}_MACS2_raw.bed

macs2 bdgpeakcall -l 30 -i ${subtract_pair}_plus.bedgraph -c ${Macs2cutoff} -o ${subtract_pair}_plus_MACS2_raw.bed
awk 'NR>1' ${subtract_pair}_plus_MACS2_raw.bed > ${subtract_pair}_plus_rawA.bed
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' ${subtract_pair}_plus_rawA.bed > ${subtract_pair}_plus.bed
rm ${subtract_pair}_plus_rawA.bed
rm ${subtract_pair}_plus_MACS2_raw.bed

macs2 bdgpeakcall -l 30 -i ${subtract_pair}_minus.bedgraph -c ${Macs2cutoff} -o ${subtract_pair}_minus_MACS2_raw.bed
awk 'NR>1' ${subtract_pair}_minus_MACS2_raw.bed > ${subtract_pair}_minus_rawA.bed
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' ${subtract_pair}_minus_rawA.bed > ${subtract_pair}_minus.bed
rm ${subtract_pair}_minus_rawA.bed
rm ${subtract_pair}_minus_MACS2_raw.bed

awk '$4 > 0' ${subtract_pair}.bedgraph > ${subtract_pair}_PM_subtraction.bedgraph

awk '$4 > 0' ${subtract_pair}_plus.bedgraph > ${subtract_pair}_plus_PM_subtraction.bedgraph

awk '$4 > 0' ${subtract_pair}_minus.bedgraph > ${subtract_pair}_minus_PM_subtraction.bedgraph

awk '$4 > 0' ${subtract_pair2}.bedgraph > ${subtract_pair2}_PM_subtraction.bedgraph

awk '$4 > 0' ${subtract_pair2}_plus.bedgraph > ${subtract_pair2}_plus_PM_subtraction.bedgraph

awk '$4 > 0' ${subtract_pair2}_minus.bedgraph > ${subtract_pair2}_minus_PM_subtraction.bedgraph

awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_PM.bedgraph > ${group1}_all_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group2}_all_PM.bedgraph > ${group2}_all_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group3}_all_PM.bedgraph > ${group3}_all_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group4}_all_PM.bedgraph > ${group4}_all_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair}.bedgraph > ${subtract_pair}_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair2}.bedgraph > ${subtract_pair2}_rounded.bedgraph

awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_plus_PM.bedgraph > ${group1}_all_plus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group2}_all_plus_PM.bedgraph > ${group2}_all_plus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group3}_all_plus_PM.bedgraph > ${group3}_all_plus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group4}_all_plus_PM.bedgraph > ${group4}_all_plus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair}_plus.bedgraph > ${subtract_pair}_plus_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair2}_plus.bedgraph > ${subtract_pair2}_plus_rounded.bedgraph

awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group1}_all_minus_PM.bedgraph > ${group1}_all_minus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group2}_all_minus_PM.bedgraph > ${group2}_all_minus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group3}_all_minus_PM.bedgraph > ${group3}_all_minus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${group4}_all_minus_PM.bedgraph > ${group4}_all_minus_PM_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair}_minus.bedgraph > ${subtract_pair}_minus_rounded.bedgraph
awk -F'\t' '{OFS="\t"; $4=int($4); print}' ${subtract_pair2}_minus.bedgraph > ${subtract_pair2}_minus_rounded.bedgraph


module load gcc/9.2.0
module load bedtools/2.30.0

bedtools intersect -a ${subtract_pair}.bed -b ${subtract_pair}_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${subtract_pair}_AuC.bed
bedtools intersect -a ${subtract_pair}_AuC.bed -b ${group1}_all_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $11} {print $1, $2, $3, $4, $5, $6, $7, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k8,8nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group1}_AuC.bed
bedtools intersect -a ${group1}_AuC.bed -b ${group2}_all_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $12} {print $1, $2, $3, $4, $5, $6, $7, $8, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k9,9nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group2}_AuC.bed
bedtools intersect -a ${group2}_AuC.bed -b ${group3}_all_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $13} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k10,10nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group3}_AuC.bed
bedtools intersect -a ${group3}_AuC.bed -b ${group4}_all_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $14} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k11,11nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group4}_AuC.bed
bedtools intersect -a ${group4}_AuC.bed -b ${subtract_pair2}_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $15} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k12,12nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${subtract_pair2}_AuC.bed
#rm ${subtract_pair2}_all.bedgraph


bedtools intersect -a ${subtract_pair}_plus.bed -b ${subtract_pair}_plus_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${subtract_pair}_plus_AuC.bed
bedtools intersect -a ${subtract_pair}_plus_AuC.bed -b ${group1}_all_plus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $11} {print $1, $2, $3, $4, $5, $6, $7, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k8,8nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group1}_plus_AuC.bed
bedtools intersect -a ${group1}_plus_AuC.bed -b ${group2}_all_plus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $12} {print $1, $2, $3, $4, $5, $6, $7, $8, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k9,9nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group2}_plus_AuC.bed
bedtools intersect -a ${group2}_plus_AuC.bed -b ${group3}_all_plus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $13} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k10,10nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group3}_plus_AuC.bed
bedtools intersect -a ${group3}_plus_AuC.bed -b ${group4}_all_plus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $14} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k11,11nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group4}_plus_AuC.bed
bedtools intersect -a ${group4}_plus_AuC.bed -b ${subtract_pair2}_plus_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $15} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k12,12nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${subtract_pair2}_plus_AuC.bed
#rm ${subtract_pair2}_all_plus.bedgraph


bedtools intersect -a ${subtract_pair}_minus.bed -b ${subtract_pair}_minus_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $10} {print $1, $2, $3, $4, $5, $6, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k7,7nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${subtract_pair}_minus_AuC.bed
bedtools intersect -a ${subtract_pair}_minus_AuC.bed -b ${group1}_all_minus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $11} {print $1, $2, $3, $4, $5, $6, $7, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k8,8nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group1}_minus_AuC.bed
bedtools intersect -a ${group1}_minus_AuC.bed -b ${group2}_all_minus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $12} {print $1, $2, $3, $4, $5, $6, $7, $8, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k9,9nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group2}_minus_AuC.bed
bedtools intersect -a ${group2}_minus_AuC.bed -b ${group3}_all_minus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $13} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k10,10nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group3}_minus_AuC.bed
bedtools intersect -a ${group3}_minus_AuC.bed -b ${group4}_all_minus_PM_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $14} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k11,11nr | awk '!seen[$1,$2,$3]++' | sort -k4,4nr > ${group4}_minus_AuC.bed
bedtools intersect -a ${group4}_minus_AuC.bed -b ${subtract_pair2}_minus_rounded.bedgraph -wa -wb | awk -v OFMT='%f' 'BEGIN {OFS = "\t"} {key = $1 "\t" $2 "\t" $3} {sums[key] += $15} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, sums[key]}' | sort -k1,1 -k2,2n -k3,3n -k12,12nr | awk '!seen[$1,$2,$3]++' | sort -k5,5nr > ${subtract_pair2}_minus_AuC.bed
#rm ${subtract_pair2}_all_minus.bedgraph

rm ${group1}*_rounded.bedgraph
rm ${group2}*_rounded.bedgraph
rm ${group3}*_rounded.bedgraph
rm ${group4}*_rounded.bedgraph
rm ${subtract_pair}*_rounded.bedgraph
rm ${subtract_pair2}*_rounded.bedgraph


bedtools intersect -wao -a ${subtract_pair2}_AuC.bed -b Human_genes.bed > ${subtract_pair}_human_annotated.bed
bedtools intersect -wao -a ${subtract_pair}_human_annotated.bed -b Repbase_IR_annotated.bed > ${subtract_pair}_both_annotated.bed
bedtools intersect -wao -a ${subtract_pair2}_plus_AuC.bed -b Human_genes.bed > ${subtract_pair}_plus_human_annotated.bed
bedtools intersect -wao -a ${subtract_pair}_plus_human_annotated.bed -b Repbase_IR_annotated.bed > ${subtract_pair}_plus_both_annotated.bed
bedtools intersect -wao -a ${subtract_pair2}_minus_AuC.bed -b Human_genes.bed > ${subtract_pair}_minus_human_annotated.bed
bedtools intersect -wao -a ${subtract_pair}_minus_human_annotated.bed -b Repbase_IR_annotated.bed > ${subtract_pair}_minus_both_annotated.bed

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12, $14, $15, $16, $17, $18, $20, $21, $22, $23, $24, $25, $26 }' ${subtract_pair}_both_annotated.bed  > ${subtract_pair}_both_annotated_dropped.bed
sort -n -r -k 5 ${subtract_pair}_both_annotated_dropped.bed > ${subtract_pair}_AuC_annotated.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${subtract_pair}_AuC\t${group1}_AuC\t${group2}_AuC\t${group3}_AuC\t${group4}_AuC\t${subtract_pair2}_AuC\tgene_start\tgene_end\tgene_name\tgene_length\tgene_direction\tgene_overlap_bp\trepeat_start\trepeat_end\trepeat_name\tIR_ID\trepeat_direction\trepeat_overlap_bp" | cat - ${subtract_pair}_AuC_annotated.bed > ${subtract_pair}_AuC_annotated_header.bed
awk '{ gsub(/[\t]/,","); print }' ${subtract_pair}_AuC_annotated_header.bed > ${subtract_pair}_AuC_annotated.csv

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12, $14, $15, $16, $17, $18, $20, $21, $22, $23, $24, $25, $26 }' ${subtract_pair}_plus_both_annotated.bed  > ${subtract_pair}_plus_both_annotated_dropped.bed
sort -n -r -k 5 ${subtract_pair}_plus_both_annotated_dropped.bed > ${subtract_pair}_plus_AuC_annotated.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${subtract_pair}_AuC\t${group1}_AuC\t${group2}_AuC\t${group3}_AuC\t${group4}_AuC\t${subtract_pair2}_AuC\tgene_start\tgene_end\tgene_name\tgene_length\tgene_direction\tgene_overlap_bp\trepeat_start\trepeat_end\trepeat_name\tIR_ID\trepeat_direction\trepeat_overlap_bp" | cat - ${subtract_pair}_plus_AuC_annotated.bed > ${subtract_pair}_plus_AuC_annotated_header.bed
awk '{ gsub(/[\t]/,","); print }' ${subtract_pair}_plus_AuC_annotated_header.bed > ${subtract_pair}_plus_AuC_annotated.csv

awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12, $14, $15, $16, $17, $18, $20, $21, $22, $23, $24, $25, $26 }' ${subtract_pair}_minus_both_annotated.bed  > ${subtract_pair}_minus_both_annotated_dropped.bed
sort -n -r -k 5 ${subtract_pair}_minus_both_annotated_dropped.bed > ${subtract_pair}_minus_AuC_annotated.bed
echo -e "Chr\tpeak_start\tpeak_end\tpeak_name\tpeak_score\t${subtract_pair}_AuC\t${group1}_AuC\t${group2}_AuC\t${group3}_AuC\t${group4}_AuC\t${subtract_pair2}_AuC\tgene_start\tgene_end\tgene_name\tgene_length\tgene_direction\tgene_overlap_bp\trepeat_start\trepeat_end\trepeat_name\tIR_ID\trepeat_direction\trepeat_overlap_bp" | cat - ${subtract_pair}_minus_AuC_annotated.bed > ${subtract_pair}_minus_AuC_annotated_header.bed
awk '{ gsub(/[\t]/,","); print }' ${subtract_pair}_minus_AuC_annotated_header.bed > ${subtract_pair}_minus_AuC_annotated.csv




rm ${subtract_pair}_AuC.bed
rm ${group1}_AuC.bed
rm ${group2}_AuC.bed
rm ${group3}_AuC.bed
rm ${group4}_AuC.bed
rm ${subtract_pair}_plus_AuC.bed
rm ${group1}_plus_AuC.bed
rm ${group2}_plus_AuC.bed
rm ${group3}_plus_AuC.bed
rm ${group4}_plus_AuC.bed
rm ${subtract_pair}_minus_AuC.bed
rm ${group1}_minus_AuC.bed
rm ${group2}_minus_AuC.bed
rm ${group3}_minus_AuC.bed
rm ${group4}_minus_AuC.bed
rm ${subtract_pair}_human_annotated.bed
rm ${subtract_pair}_plus_human_annotated.bed
rm ${subtract_pair}_minus_human_annotated.bed
rm ${subtract_pair}_both_annotated_dropped.bed
rm ${subtract_pair}_AuC_annotated_header.bed
rm ${subtract_pair}_plus_both_annotated_dropped.bed
rm ${subtract_pair}_plus_AuC_annotated_header.bed
rm ${subtract_pair}_minus_both_annotated_dropped.bed
rm ${subtract_pair}_minus_AuC_annotated_header.bed
rm ${subtract_pair}_both_annotated.bed
rm ${subtract_pair}_plus_both_annotated.bed 
rm ${subtract_pair}_minus_both_annotated.bed
rm ${subtract_pair}_AuC_annotated.bed
rm ${subtract_pair}_plus_AuC_annotated.bed
rm ${subtract_pair}_minus_AuC_annotated.bed
rm ${subtract_pair2}_AuC.bed
rm ${subtract_pair2}_plus_AuC.bed
rm ${subtract_pair2}_minus_AuC.bed

mkdir ${subtract_pair}
mv ${subtract_pair}* ${subtract_pair}
mv ${subtract_pair2}* ${subtract_pair}
mv ${group1}*bedgraph ${subtract_pair}
mv ${group2}*bedgraph ${subtract_pair}
mv ${group3}*bedgraph ${subtract_pair}
mv ${group4}*bedgraph ${subtract_pair}




