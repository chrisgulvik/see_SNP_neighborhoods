#!/usr/bin/env python
# usage: python see_gSNP_neighborhood.py -1 Sample1 -2 Sample1 -f reference.fasta -b1 Sample1.Input.sorted.bam -b2 Sample2.Input.sorted.bam -v1 Sample1.Input.vcf[.gz] -v2 Sample2.Input.vcf[.gz] -s out.snpmatrix.tsv
# see_sample-pair_SNP_neighborhood.py -n 4 -f reference/reference.fasta -1 3000113810 -2 3000113846 -b1 bam/3000113810.fastq-reference.sorted.bam -b2 bam/3000113846.fastq-reference.sorted.bam -s msa/out.snpmatrix.tsv -v1 vcf/3000113810.fastq-reference.vcf -v2 vcf/3000113846.fastq-reference.vcf


import argparse
import os
import random
import re
import subprocess
import sys
import reportlab.platypus as Platypus
import vcf
from operator import itemgetter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus.paragraph import Paragraph
from reportlab.lib.units import inch
from reportlab.rl_config import defaultPageSize


def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Prints raw read pileups to manually inspect putative SNPs of a sample pair')
	parser.add_argument('-1', '--Sample1Name', help='First sample name in VCF to compare (pairwise)')
	parser.add_argument('-2', '--Sample2Name', help='Second sample name in VCF to compare (pairwise)')
	parser.add_argument('-b1', '--bam1', required=True, help='indexed BAM input file')
	parser.add_argument('-b2', '--bam2', required=True, help='indexed BAM input file')
	parser.add_argument('-f', '--fasta', required=True, help='FastA input file')
	parser.add_argument('-s', '--snpmatrix', required=True, help='out.snpmatrix.tsv')
	parser.add_argument('-v1', '--vcf1', required=True, help='VCF input file')
	parser.add_argument('-v2', '--vcf2', required=True, help='VCF input file')
	parser.add_argument('-o', '--output', default='SNP_Report.pdf', help='output PDF filename')
	parser.add_argument('-r', '--rows', type=int, default=12, help='maximum number of rows/lines to show per pileup')
	parser.add_argument('-n', '--numsites', type=int, default=10, help='number of SNP sites to investigate')
	parser.add_argument('-p', '--positions', help='CSV-delimited list of SNP positions to investigate')
	return parser.parse_args()


def get_longestSubstr(s1, s2):
	longestSubstr = ''
	str1 = os.path.basename(s1) #avoid matching filepaths
	str2 = os.path.basename(s2)
	len1, len2 = len(str1), len(str2)
	for i in range(len1):
		match = ''
		for j in range(len2):
			if (i+j<len1 and str1[i+j]==str2[j]):
				match += str2[j]
			else:
				if (len(match) > len(longestSubstr)):
					longestSubstr = match
				else:
					match = ''
		if len(match) > len(longestSubstr):
			longestSubstr = match
	return longestSubstr.rstrip('.') #discard trailing period if exists


def get_pairwiseSNPlist(snpMatrix, Sample1Name, Sample2Name):
	pairedSNPs = {}
	with open(snpMatrix) as snpTSV:
		for line in snpTSV:
			#get sample ID names from snpmatrix file header
			if line.startswith('#'):
				cleanedLine = re.sub(r':\w+\t', '\t', line).strip('# ').rstrip()
				header = re.split(r'\t', cleanedLine)[3:]
				sampleIDs = [] #sample ID order is +3 positions with respect to file columns (0-base)
				for i in header:
					sampleIDs.append(re.sub(r'\[[0-9]+\]', '', i))

			#create pairwise SNP list of lists containing chrom ID, pos, Sample1's nucleotide, and Sample2's nucleotide
			if not line.startswith('#'):
				tabs = re.split(r'\t', line.rstrip())
				#locate col num for each sample ID
				Sample1Col = sampleIDs.index(Sample1Name)
				Sample2Col = sampleIDs.index(Sample2Name)

				Chrom = tabs[0]
				SamplePos = tabs[1]
				Sample1Nucl = tabs[Sample1Col + 3]
				Sample2Nucl = tabs[Sample2Col + 3]
				if Sample1Nucl is not Sample2Nucl:
					if Sample1Nucl is '.':
						Sample1Nucl = tabs[2]
					if Sample2Nucl is '.':
						Sample2Nucl = tabs[2]
					pairedSNPs[int(SamplePos)] = [Chrom, SamplePos, Sample1Nucl, Sample2Nucl]
					#pairedSNPs.append([Chrom, SamplePos, Sample1Nucl, Sample2Nucl])
				# else:
				# 	print 'no SNP in position %s: %s & %s are identical' % (SamplePos, Sample1Nucl, Sample2Nucl)
	return pairedSNPs


def get_pairwiseSNP_vcf(pos, infile):
	VCF = vcf.Reader(filename=infile, strict_whitespace=True)
	vcfData = []
	for record in VCF: #class 'vcf.parser.Reader' object
		if not pos:
			#print 'quit parsing on %s (to save time)' % record.POS
			break
		if int(pos[0]) == record.POS:
			if 'DP' in record.INFO: #record.INFO is a dict
				DP = record.INFO['DP']	
			else:
				DP = '*not available*'
			if 'MQ' in record.INFO:
				MQ = record.INFO['MQ']
			else:
				MQ = '*not available*'
			if record.QUAL:
				QUAL = str(round(record.QUAL, 1))
			else:  #avoids error when 'QUAL' in VCF is '.'
				QUAL = '*not available*'
			vcfData.append([int(pos[0]), record.CHROM, QUAL, DP, MQ])
			del pos[0]
		elif int(pos[0]) < record.POS:
			print '%s position not found in %s' % (pos[0], infile)
			vcfData.append([int(pos[0]), record.CHROM, '*POS not in VCF*', '*POS not in VCF*', '*POS not in VCF*'])
			del pos[0]
	return sorted(vcfData, reverse=True)


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
	inbam1 = args.bam1
	inbam2 = args.bam2
	infasta = args.fasta
	invcf1 = args.vcf1
	invcf2 = args.vcf2
	quantityWantSites = int(args.numsites)
	output = args.output
	rows = args.rows
	Sample1Name = args.Sample1Name
	Sample2Name = args.Sample2Name
	snpMatrix = args.snpmatrix

	# attempt to auto-detect sample names (experimental)
	if not Sample1Name:
		Sample1Name = (get_longestSubstr(inbam1, invcf1)).rstrip('.fastq-reference')
	if not Sample2Name:
		Sample2Name = (get_longestSubstr(inbam2, invcf2)).rstrip('.fastq-reference')

	pairwiseSNPs = get_pairwiseSNPlist(snpMatrix, Sample1Name, Sample2Name)
	pairwiseSNPpositions = pairwiseSNPs.keys()

	# enable explicitly given positions (e.g., -p '89,3969181,44,123456') to investigate
	if args.positions:
		posList = [int(n) for n in args.positions.split(',')]
		unsortedListWantSitesVCF = []
		for i in posList:
			if i in pairwiseSNPpositions:
				unsortedListWantSitesVCF.append(i)
			else:
				sys.exit('ERROR: specified position %s not a pairwise SNP' % i)
		quantityWantSites = len(unsortedListWantSitesVCF)
		listWantSitesVCF = sorted(unsortedListWantSitesVCF, reverse=True) #for pileup report building
		SitesVCF = sorted(unsortedListWantSitesVCF)
		Sites1 = SitesVCF[:] #copy because I delete each pos/element as I iterate through the list in the get_pairwiseSNP_vcf function
		Sites2 = SitesVCF[:]

	# randomly select SNP sites if positions (-p) unspecified
	else:
		if quantityWantSites > len(pairwiseSNPpositions):
			print 'Number of requested SNP sites to investigate exceeds number of identified SNPs'
			print 'Selecting all %s putative SNP sites to print...' % len(pairwiseSNPpositions)
			quantityWantSites = len(pairwiseSNPpositions)
		randSitesUnsorted = random.sample(pairwiseSNPpositions, quantityWantSites) #randomly selected sites to investigate
		randSitesUnsortedInts = [int(i) for i in randSitesUnsorted] #convert strings to integers
		listWantSitesVCF = sorted(randSitesUnsortedInts, reverse=True) #for pileup report building
		randSites = sorted(randSitesUnsortedInts)
		Sites1 = randSites[:] #copy because I delete each pos/element as I iterate through the list in the get_pairwiseSNP_vcf function
		Sites2 = randSites[:]

	vcfData1 = get_pairwiseSNP_vcf(Sites1, invcf1)
	vcfData2 = get_pairwiseSNP_vcf(Sites2, invcf2)

	if len(vcfData1) != len(vcfData2):
		sys.exit('ERROR in VCF parsing')
	numSites = len(vcfData1)

	#calculate page and frames
	doc = Platypus.BaseDocTemplate(output, topMargin=0, bottomMargin=0, leftMargin=10, rightMargin=0)
	doc.pagesize = landscape(A4) #ISO Code A4
	# A4 Dimensions: 297 mm x 210 mm ; 11.69 in x 8.27 in ; 842 pt x 595 pt

	frameCount = 2  #(numsites + (-numsites%6)) // 6
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

	Title = ('Report containing ' + str(quantityWantSites) + ' of ' + str(len(pairwiseSNPpositions)) +
			' putative SNPs for ' + os.path.basename(invcf1) + ' and ' +
			os.path.basename(invcf2))
	report = [Paragraph(Title, styles['Heading2'])]

	while numSites > 0:
		#handle SNP sites at beginning of Chrom (<50 bases in)
		posMinusFitty = (listWantSitesVCF[numSites - 1] - 50)
		centeredPos = str(max(posMinusFitty, 0)) 

		pileup1 = Platypus.Preformatted(get_ttview(inbam1, infasta,
		vcfData1[numSites-1][1], centeredPos, rows), style)
		pileup2 = Platypus.Preformatted(get_ttview(inbam2, infasta,
		vcfData2[numSites-1][1], centeredPos, rows), style)

		pos = str(vcfData1[numSites-1][0])
		MQ1 = str(vcfData1[numSites-1][4])
		MQ2 = str(vcfData2[numSites-1][4])
		DP1 = str(vcfData1[numSites-1][3])
		DP2 = str(vcfData2[numSites-1][3])
		qual1 = str(vcfData1[numSites-1][2])
		qual2 = str(vcfData2[numSites-1][2])
		SNP = str(pairwiseSNPs[int(pos)][2] + '-vs-' + pairwiseSNPs[int(pos)][3])

		header1 = Paragraph(SNP + ' at position ' + pos +
			' for ' + Sample1Name +
			' with a QUAL of ' + qual1 + ", MQ of " + MQ1 +
			', and raw read depth (DP) of ' + DP1
			, styles['Heading6'])
		header2 = Paragraph(SNP + ' at position ' + pos +
			' for ' + Sample2Name +
			' with a QUAL of ' + qual2 + ", MQ of " + MQ2 +
			', and raw read depth (DP) of ' + DP2
			, styles['Heading6'])

		gap = Platypus.Spacer(0.25, 0.05*inch)
		report.append(Platypus.KeepTogether([gap, header1, pileup1]))
		report.append(Platypus.KeepTogether([gap, header2, pileup2]))
		numSites -= 1
	doc.build(report)


if __name__ == '__main__':
	main()
