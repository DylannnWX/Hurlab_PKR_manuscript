Here is the sbatch shell script on Harvard Medical School (HMS) O2 cloud computing platform. You can modify the sbatch accordingly for your own slurm.
Here the below steps were executed:
  1. Fastq file trimming & pairing
  2. Trimmed & paired files -> map to rRNA
  3. Leftover paired reads -> map to human genome
For quality control purposes, the human genome bam file were separated into plus and minus stranded bam files. Then, the GAPDH coverage value (chr12: 6534517-6538371) were calculated, and the coverage were saved in the file name. If your RNA-seq samples are strand specific, depending on the kit, the GAPDH should be overwhelmingly plus-stranded (e.g. SSIV) or minus-stranded (e.g. SMARTer).
