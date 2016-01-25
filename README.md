# see SNP neighborhoods
**Purpose:** visualize read pileup neighborhoods of SNPs to identify false positives

![alt tag](https://github.com/chrisgulvik/images/raw/master/see-SNPs.jpg)

**Install:**
compile from source with make and export the bin dir to execute a python script. The samtools binary shouldn't conflict with what you might already have, because I named it 'samtools-0.1.18'. Currently the pyvcf module is outdated in the pip package manager, so you'll likely need to compile from source the bleeding edge version as an egg to avoid errors parsing ChromIDs in some VCF files: `pip install -e git+https://github.com/jamescasbon/PyVCF.git#egg=PyVCF --upgrade`


###### Scripts:
* **see_gSNP_neighborhood.py** -- generates a PDF report of SNP pileups for a single sample. SNPs are nucleotides that are different between the sample and the reference genome.
* **see_sample-pair_SNP_neighborhood.py** -- generates a PDF report of SNP pileups for a pair of samples. SNPs are nucleotides that are different between the two samples irrespective of the reference genome's content.
