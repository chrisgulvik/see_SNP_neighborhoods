#!/usr/bin/env python
# usage: python see_gSNP_neighborhood.py -b Input.sorted.bam -f Input.fasta -v Input.vcf[.gz] -o Output.pdf


import argparse
import os
import random
import subprocess
import reportlab.platypus as Platypus
import vcf
from operator import itemgetter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus.paragraph import Paragraph
from reportlab.lib.units import inch
from reportlab.rl_config import defaultPageSize


def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Prints raw read pileups to manually inspect genomic SNPs')
	parser.add_argument('-b', '--bam', required=True, help='indexed BAM input file')
	parser.add_argument('-f', '--fasta', required=True, help='FastA input file')
	parser.add_argument('-v', '--vcf', required=True, help='vcf[.gz] input file')
	parser.add_argument('-o', '--output', required=False, default='SNP_report.pdf', help='output PDF filename')
	parser.add_argument('-n', '--numsites', required=False, type=int, default='12', help='number of SNP sites to investigate')
	#parser.add_argument('-s', '--sites', required=False, default='', help='comma-separated positions of SNP sites to investigate')
	return parser.parse_args()


def get_vcf(infile):
	VCF = vcf.Reader(filename=infile, strict_whitespace=True)
	INDELs = []
	SNPs = []
	for record in VCF:
		if record.is_indel:
			INDELs.append((record.CHROM, record.POS, record.QUAL))
		elif record.is_snp:
			SNPs.append((record.CHROM, record.POS, record.QUAL, str(record.REF) + '->' + str(record.ALT).strip('[]')))
	return SNPs, INDELs


def get_ttview(bam, fasta, chrom, pos, rows):
	x = 100
	viewCmd = ['ttview',
		'-X', str(x), #number of columns to show including the SNP site
		'-g', '%s:%s' % (chrom, pos), #position to inspect
		bam, fasta]
	pileup = subprocess.Popen(viewCmd, stdout=subprocess.PIPE)
	out, err = pileup.communicate()
	return '\n'.join(out.split('\n')[0:rows])


def main():
	args = parseArgs()
	inbam = args.bam
	infasta = args.fasta
	invcf = args.vcf
	numsites = args.numsites
	output = args.output

	wantSites = int(numsites)
	SNPs, INDELs = get_vcf(invcf)
	vcfSNPs = len(SNPs)
	print 'Identified %s putative SNPs' % vcfSNPs
	#print 'Identified %s putative INDEL bases' % len(INDELs)
	if wantSites > vcfSNPs:
		print 'Number of requested SNP sites to investigate exceeds number of identified SNPs'
		print 'Selecting all %s putative SNP sites to print...' % vcfSNPs
		wantSites = vcfSNPs
	randSitesUnsorted = random.sample(SNPs, wantSites) #randomly selected sites to investigate
	randSites = sorted(randSitesUnsorted, key=itemgetter(1), reverse=True)

	#calculate page and frames
	doc = Platypus.BaseDocTemplate(output, topMargin=0, bottomMargin=0, leftMargin=10, rightMargin=0)
	doc.pagesize = landscape(A4) #ISO Code A4
	# A4 Dimensions: 297 mm x 210 mm ; 11.69 in x 8.27 in ; 842 pt x 595 pt

	frameCount = (numsites + (-numsites%6)) // 6
	frameWidth = doc.height / frameCount
	frameHeight = doc.width - .05 * inch
	frames = []

	#construct a frame for each column
	for frame in range(frameCount):
		leftMargin = doc.leftMargin + frame * frameWidth
		column = Platypus.Frame(leftMargin, doc.bottomMargin, frameWidth, frameHeight)
		frames.append(column)
	template = Platypus.PageTemplate(frames=frames)
	doc.addPageTemplates(template)
	
	PAGE_HEIGHT = defaultPageSize[0] #size in points
	PAGE_WIDTH = defaultPageSize[1]
	styles = getSampleStyleSheet()
	style = styles['Normal']
	style.fontName = 'Courier'
	style.fontSize = 6.5

	Title = "Report containing " + str(wantSites) + " of " + str(vcfSNPs) + " putative SNPs for " + os.path.basename(invcf)
	report = [Paragraph(Title, styles['Heading2'])]
	while wantSites > 0:
		pileup = Platypus.Preformatted(get_ttview(inbam, infasta, \
			randSites[wantSites-1][0], randSites[wantSites-1][1] - 50, 10), style)
		if randSites[wantSites-1][2]:
			# print 'mapqual: %s' % mapqual
			mapqual = str(round(randSites[wantSites-1][2], 1))
		else: #avoids error when 'QUAL' in VCF is '.'
			mapqual = '*not available*'
		pos = str(randSites[wantSites-1][1])
		SNP = str(randSites[wantSites-1][3])
		header = Paragraph(SNP + " at position " + pos + " with a base quality score of " + mapqual, styles['Heading6'])
		gap = Platypus.Spacer(0.25, 0.05*inch)
		report.append(Platypus.KeepTogether([gap, header, pileup]))
		wantSites -= 1
	doc.build(report)


if __name__ == '__main__':
	main()
