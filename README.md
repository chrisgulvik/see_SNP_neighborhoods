# see SNP neighborhoods
**Purpose:** visualize read pileup neighborhoods of SNPs to identify false positives


**Install:**
compile from source with make and export the bin dir to execute a python script


###### Scripts:
* **see_gSNP_neighborhood.py** -- generates a PDF report of SNP pileups for a single sample. SNPs are nucleotides that are different between the sample and the reference genome.
* **see_sample-pair_SNP_neighborhood.py** -- generates a PDF report of SNP pileups for a pair of samples. SNPs are nucleotides that are different between the two samples irrespective of the reference genome's content.
