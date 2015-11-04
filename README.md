# see_SNP_neighborhoods
Visualize read pileup neighborhoods of SNPs to identify false positives


Install:
Compile from source with make, then export the bin dir to use the py script.


Scripts:
(1) see_gSNP_neighborhood.py -- generates a PDF report of SNP pileups for a single sample. SNPs are nucleotides that are different between the sample and the reference genome.

(2) see_sample-pair_SNP_neighborhood.py -- generates a PDF report of SNP pileups for a pair of samples. SNPs are nucleotides that are different between the two samples irrespective of the reference genome's content.