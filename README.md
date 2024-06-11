# Hurlab_PKR_manuscript
Welcome to Dylan (Xi Wang, Ph.D.)'s Github page on the PKR manuscript! 

Here you will find the custom bioinformatics pipelines &amp; scripts I designed for PKR manuscript. These pipeline include:

1. Slurm sbatch jobs utilizing Harvard Medical School (HMS) O2 cloud computing platform. You can modify the sbatch shell script for your own slurm.
   1a: sbatch_map: trim & pair the fastq reads, map first to human rRNA (RNA45SN1 and RNA5SN1), then map the leftover reads to human genome. A quality control is included by quantifying the plus/minus coverage values of GAPDH region, and the values were saved in the file name for your convenience.
   
   1b: sbatch_bdg_callpeak_annotate: take the bam files as inputs, merge all bams within the same group (rRNA and human genome), and then normalize to 1M total non-rRNA reads. Then, the normalized bedgraph files were used for (pulldown - input) subtraction, and peak calling. The peaks were then further annotated by human genes and repeats. The human genes utilized NCBI Grch38.p14 annotation, while repeats utilized Repbase (https://www.girinst.org/repbase/)
   
   1c: sbatch_add_AuC: take the peaks called in step 1b, and then add the individual normalized AuC values. The results were exported in csv formats.

2. Python script of finding Inverted Repeats (IR): this jupyter notebook scripts takes the input of Repbase repeat names and coordinates, and search nearby repeats by utilizing NCBI's BLAST engine (blastn v2.15.0; reward=2, penalty=-3, gapopen=5, gapextend=2, word_size=11). The blastn results were saved in additional columns, and the result spreadsheet were saved with "IR_Blasted_" as a header from the original datasheet.


