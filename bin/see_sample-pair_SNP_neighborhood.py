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
	parser.add_argument('-1', '--Sample1Name', required=True, help='First sample name in VCF to compare (pairwise)')
	parser.add_argument('-2', '--Sample2Name', required=True, help='Second sample name in VCF to compare (pairwise)')
	parser.add_argument('-b1', '--bam1', required=True, help='indexed BAM input file')
	parser.add_argument('-b2', '--bam2', required=True, help='indexed BAM input file')
	parser.add_argument('-f', '--fasta', required=True, help='FastA input file')
	parser.add_argument('-s', '--snpmatrix', required=True, help='out.snpmatrix.tsv')
	parser.add_argument('-v1', '--vcf1', required=True, help='VCF input file')
	parser.add_argument('-v2', '--vcf2', required=True, help='VCF input file')
	parser.add_argument('-o', '--output', required=False, default='SNP_report2.pdf', help='output PDF filename')
	parser.add_argument('-r', '--rows', required=False, type=int, default='12', help='maximum number of rows/lines to show per pileup')
	parser.add_argument('-n', '--numsites', required=False, type=int, default='12', help='number of SNP sites to investigate')
	return parser.parse_args()


def pairwiseSNPlist(snpMatrix, Sample1Name, Sample2Name):
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
			print 'quit parsing on %s (to save time)' % record.POS
			break
		if int(pos[0]) == record.POS:
			#print pos[0], record
			if record.QUAL:
				qual = str(record.QUAL)
			else:
				qual = '*absent from VCF*'
			if 'ADP' in record.INFO:
				ADP = str(record.INFO['ADP'])  #record.INFO is a dict
			else:
				ADP = '*absent from VCF*'
			sampleDataStr = str(record.samples).strip('CallData()[]')
			uglyList = re.findall(r'\w+=[\d./E]+', sampleDataStr)
			sampleDataDict = dict([pair.split('=', 1) for pair in uglyList])
			if 'SDP' in sampleDataDict:
				SDP = sampleDataDict['SDP']
			else:
				SDP = '*absent from VCF*'
			vcfData.append([int(pos[0]), record.CHROM, qual, ADP, SDP])
			#vcfData[int(pos[0])] = [record.CHROM, qual, ADP, SDP]
			del pos[0]
		elif int(pos[0]) < record.POS:
			print '%s position not found in VCF' % pos[0]
			vcfData.append([int(pos[0]), record.CHROM, '*absent from VCF*', '*absent from VCF*', '*absent from VCF*'])
			#vcfData[int(pos[0])] = [record.CHROM, 'absent from VCF', 'absent from VCF', 'absent from VCF']
			del pos[0]
		#VCF.fetch(record.CHROM, int(i), int(i)))
	return vcfData


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
	wantSites = int(args.numsites)
	output = args.output
	rows = args.rows
	Sample1Name = args.Sample1Name
	Sample2Name = args.Sample2Name
	snpMatrix = args.snpmatrix

	pairwiseSNPs = pairwiseSNPlist(snpMatrix, Sample1Name, Sample2Name)
	pairwiseSNPpositions = pairwiseSNPs.keys()

	if wantSites > len(pairwiseSNPpositions):
		print 'Number of requested SNP sites to investigate exceeds number of identified SNPs'
		print 'Selecting all %s putative SNP sites to print...' % len(pairwiseSNPpositions)
		wantSites = len(pairwiseSNPpositions)
	randSitesUnsorted = random.sample(pairwiseSNPpositions, wantSites) #randomly selected sites to investigate
	randSitesUnsortedInts = [int(i) for i in randSitesUnsorted] #convert strings to integers
	randSites = sorted(randSitesUnsortedInts)
	randSitesRev = sorted(randSitesUnsortedInts, reverse=True)
	randSites1 = randSites[:]
	randSites2 = randSites[:]
	randSitesRev2 = randSitesRev

	vcfData1 = get_pairwiseSNP_vcf(randSites1, invcf1)
	vcfData1Rev = sorted(vcfData1, reverse=True)
	vcfData2 = get_pairwiseSNP_vcf(randSites2, invcf2)
	vcfData2Rev = sorted(vcfData2, reverse=True)
	print 'parsing of VCFs completed...'

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

	Title = ('Report containing ' + str(wantSites) + ' of ' + str(len(pairwiseSNPpositions)) +
			' putative SNPs for ' + os.path.basename(invcf1) + ' and ' +
			os.path.basename(invcf2))
	report = [Paragraph(Title, styles['Heading2'])]

	while numSites > 0:
		if randSitesRev[numSites-1] < 50:
			position = (vcfData1Rev[numSites-1][0])
		else:
			position = (vcfData1Rev[numSites-1][0] - 50)
		pileup1 = Platypus.Preformatted(get_ttview(inbam1, infasta,
		vcfData1Rev[numSites-1][1], position, rows), style)
		pileup2 = Platypus.Preformatted(get_ttview(inbam2, infasta,
		vcfData2Rev[numSites-1][1], position, rows), style)

		pos = str(vcfData1Rev[numSites-1][0])
		SDP1 = vcfData1Rev[numSites-1][4]
		SDP2 = vcfData2Rev[numSites-1][4]
		qual1 = vcfData1Rev[numSites-1][2]
		qual2 = vcfData2Rev[numSites-1][2]
		SNP = str(pairwiseSNPs[int(pos)][2] + '-vs-' + pairwiseSNPs[int(pos)][3])

		header1 = Paragraph(SNP + ' at position ' + pos +
			' for ' + os.path.basename(invcf1) +
			' with a base quality score of ' + qual1 +
			' and raw read depth of ' + SDP1
			, styles['Heading6'])
		header2 = Paragraph(SNP + ' at position ' + pos +
			' for ' + os.path.basename(invcf2) +
			' with a base quality score of ' + qual2 +
			' and raw read depth of ' + SDP2
			, styles['Heading6'])
		gap = Platypus.Spacer(0.25, 0.05*inch)
		report.append(Platypus.KeepTogether([gap, header1, pileup1]))
		report.append(Platypus.KeepTogether([gap, header2, pileup2]))
		numSites -= 1
	doc.build(report)


if __name__ == '__main__':
	main()
