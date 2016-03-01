# see SNP neighborhoods
**Purpose:** visualize read pileup neighborhoods of SNPs to identify false positives

![alt tag](https://github.com/chrisgulvik/images/raw/master/see-SNPs.jpg)

###### Scripts:
* **see_gSNP_neighborhood.py** -- generates a PDF report of SNP pileups for a single sample. SNPs are nucleotides that are different between the sample and the reference genome.
* **see_sample-pair_SNP_neighborhood.py** -- generates a PDF report of SNP pileups for a pair of samples. SNPs are nucleotides that are different between the two samples irrespective of the reference genome's content.
* **batch_see_SNP_neighborhood_sample-pairs.py** -- generates a series of PDF reports, each being SNPs between a pair of samples (as see_sample-pair_SNP_neighborhood.py); well-suited for lyve-SET's output.

## Example Install:
Compile from source with make and export the bin dir to execute a python script. The samtools binary shouldn't conflict with what you might already have, because I named it 'samtools-0.1.18'. Currently the pyvcf module is outdated in the pip package manager, so you'll likely need to compile from source the bleeding edge version as an egg to avoid errors parsing ChromIDs in some VCF files.

    pip install -e git+https://github.com/jamescasbon/PyVCF.git#egg=PyVCF --upgrade
    git clone https://github.com/chrisgulvik/see_SNP_neighborhoods.git $HOME/see_SNP_neighborhoods
    cd $HOME/see_SNP_neighborhoods
    make
    echo 'export PATH="$PATH:$HOME/see_SNP_neighborhoods/bin"' >> $HOME/.bash_profile

## Example Execution:
After running lyve-SET, where ~/lSET_out is the parent outdir containing {asm,bam,log,msa,out,reads,reference,vcf}, you'll need to give the same reference file

    batch_see_SNP_neighborhood_sample-pairs.py -S ~/lSET_out -f ~/ref.fna
    gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -q -sOutputFile=SNP_hoods/Summary_SNP_hoods.pdf SNP_hoods/*.pdf
