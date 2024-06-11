Here we use python and NCBI's BLAST engine to test whether two neighboring repeats are inverted repeats. 
First, all repeat names and sequences were acquired from Repbase (https://www.girinst.org/repbase/).
Then, this program takes the regions of interest (e.g. csv files with chromosome, start and end position) and checks whether the Repbase regions have any inverted repeats nearby.
It then saves all results as additional columns, and exports the new csv file with "IR_BLASTed_" as a header of the original csv file. 
