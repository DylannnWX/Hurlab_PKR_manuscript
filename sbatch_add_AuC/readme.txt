Here is the sbatch shell script on Harvard Medical School (HMS) O2 cloud computing platform. You can modify the sbatch accordingly for your own slurm.
The steps in this shell script include:
1. Take the input peaks and bam files, merge all bam files in the same group (e.g. rRNA and human), and normalize to 1M total non-rRNA human reads.
2. Calculate the AuC (area under curve) of peak regions after normalization
3. Export the results in a .csv file.
