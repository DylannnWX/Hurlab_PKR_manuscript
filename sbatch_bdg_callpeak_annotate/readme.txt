Here is the sbatch shell script on Harvard Medical School (HMS) O2 cloud computing platform. You can modify the sbatch accordingly for your own slurm.
This shell script executes the below steps:
1. Take the input of four groups: fCLIP PD, fCLIP input, control PD, control input, as well as peak calling parameter (we usually use 20 based on numerous tests)
2. Merge all bam files in the same group (e.g. rRNA and human genome), normalize to 1M total human genomic read mapped, and then subtract fCLIP input over fCLIP PD.
3. Call peak of the subtracted bedgraph track, and then annotate based on human genome and Repbase annotations.
4. Export the result files in both .bed (standard peak file format) and .csv
