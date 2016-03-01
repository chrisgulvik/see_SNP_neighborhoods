#!/usr/bin/env python

import argparse
import os
import random
import re
import subprocess
import sys
import reportlab.platypus as Platypus
import vcf
import vcf.utils
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
	parser.add_argument('-f', '--fasta', required=True, help='FastA reference file')
	parser.add_argument('-s', '--snpmatrix', required=True, help='out.snpmatrix.tsv')
	parser.add_argument('-v1', '--vcf1', required=True, help='VCF input file')
	parser.add_argument('-v2', '--vcf2', required=True, help='VCF input file')
	parser.add_argument('-o', '--output', default='SNP_Report.pdf', help='output PDF filename')
	parser.add_argument('-r', '--rows', type=int, default=10, help='maximum number of read rows to show per pileup')
	parser.add_argument('-n', '--numsites', type=int, default=10, help='number of SNP sites to investigate')
	parser.add_argument('-p', '--positions', help='CSV-delimited list of SNP positions to investigate')
	return parser.parse_args()

def get_prefix(s1, s2):
	str1 = os.path.basename(s1) #avoid matching filepaths
	str2 = os.path.basename(s2)
	l = [str1, str2]
	return os.path.commonprefix(l).rstrip('.')

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
	return pairedSNPs

def make_data_dict(record, sample_ID):
	''' hack to make a dict from a _Call object '''
	data_dict = {}
	for item in re.split(', [A-Z]', str((record.genotype(sample_ID)).data).lstrip('CallData(').rstrip(')')):
		if '=' in item:  #avoid unusual string containing 'one]' or 'one'
			key, val = item.split('=')
			data_dict[key] = val
	return data_dict

def get_vcf_data(pos, infile1, infile2):
	''' takes in a list of coordinates and uses 2 vcf files
	to return 2 lists of the coordinates with associated quality data '''
	VCF1 = vcf.Reader(filename=infile1, strict_whitespace=True)
	ID1 = VCF1.samples[0]
	VCF2 = vcf.Reader(filename=infile2, strict_whitespace=True)
	ID2 = VCF2.samples[0]
	vcf_iter = vcf.utils.walk_together(VCF1, VCF2, vcf_record_sort_key=lambda r: (r.CHROM, r.POS, r.REF))
	vcfData1 = []  #order is: [pos, chrom, qual, DP, MQ, alt]
	vcfData2 = []
	for (rec1, rec2) in vcf_iter:
		if not pos:
			break
		if rec1 is not None and rec2 is not None:
			# rec1.FORMAT gives 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR:FT:DP4'
			if pos[0] == rec1.POS:
				for (rec, ID, vcfData) in [(rec1, ID1, vcfData1), (rec2, ID2, vcfData2)]:
					single_pos_VCF_dataset = []
					single_pos_VCF_dataset.extend((pos[0], rec.CHROM))
					if rec.QUAL:
						single_pos_VCF_dataset.append(rec.QUAL)
					else:
						single_pos_VCF_dataset.append('-')

					# add metrics of interest
					metrics = ['DP', 'MQ']
					vcf_data_dict = make_data_dict(rec, ID)
					for metric in metrics:
						if metric in vcf_data_dict:
							single_pos_VCF_dataset.append(vcf_data_dict[metric])
						else:
							single_pos_VCF_dataset.append('-')

					# get nucleotide identity [ATCG] even if same as Ref
					if rec.ALT[0] is None:
						nucleotide = rec.REF[0]
					else:
						nucleotide = str(rec.ALT[0])
					single_pos_VCF_dataset.append(nucleotide)
					vcfData.append(single_pos_VCF_dataset)
				del pos[0]

			elif pos[0] < rec1.POS:
				print '{} position not found in {}'.format(pos[0], infile1)
				del pos[0]

	return sorted(vcfData1, reverse=True), sorted(vcfData2, reverse=True)

def get_ttview(bam, fasta, chrom, pos, rows):
	x = 100
	viewCmd = ['ttview',
		'-X', str(x), #number of columns to show including the SNP site
		'-g', '{}:{}'.format(chrom, pos), #position to inspect
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
	rows = (args.rows + 2)
	Sample1Name = args.Sample1Name
	Sample2Name = args.Sample2Name
	snpMatrix = args.snpmatrix

	# Attempt to auto-detect sample names (experimental)
	if not Sample1Name:
		Sample1Name = (get_prefix(inbam1, invcf1)).rstrip('.fastq-reference')
	if not Sample2Name:
		Sample2Name = (get_prefix(inbam2, invcf2)).rstrip('.fastq-reference')

	pairwiseSNPs = get_pairwiseSNPlist(snpMatrix, Sample1Name, Sample2Name)
	pairwiseSNPpositions = pairwiseSNPs.keys()

	# Enable explicitly given positions (e.g., -p '89,3969181,44,123456') to investigate
	if args.positions:
		posList = [int(n) for n in args.positions.split(',')]
		unsortedListWantSitesVCF = []
		for i in posList:
			if i in pairwiseSNPpositions:
				unsortedListWantSitesVCF.append(i)
			else:
				sys.exit('ERROR: specified position {} not a pairwise SNP'.format(i))
		quantityWantSites = len(unsortedListWantSitesVCF)
		listWantSitesVCF = sorted(unsortedListWantSitesVCF, reverse=True) #for pileup report building
		SitesVCF = sorted(unsortedListWantSitesVCF)

	# Randomly select SNP sites if positions (-p) unspecified
	else:
		if quantityWantSites > len(pairwiseSNPpositions):
			print 'Selecting all {} putative SNP sites to print...'.format(len(pairwiseSNPpositions))
			quantityWantSites = len(pairwiseSNPpositions)
		randSitesUnsorted = random.sample(pairwiseSNPpositions, quantityWantSites)  #randomly selected sites to investigate
		randSitesUnsortedInts = [int(i) for i in randSitesUnsorted]  #convert strings to integers
		listWantSitesVCF = sorted(randSitesUnsortedInts, reverse=True)  #for pileup report building
		randSites = sorted(randSitesUnsortedInts)

	vcfData1, vcfData2 = get_vcf_data(randSites, invcf1, invcf2)

	if len(vcfData1) != len(vcfData2):
		print 'ERROR: unequal number of SNPs identified between {} and {}'.format(Sample1Name, Sample2Name)
		sys.exit('ERROR in VCF parsing')
	numSites = len(vcfData1)

	# Calculate page and frames
	doc = Platypus.BaseDocTemplate(output, topMargin=0, bottomMargin=0, leftMargin=10, rightMargin=0)
	doc.pagesize = landscape(A4) #ISO Code A4
	# A4 Dimensions: 297 mm x 210 mm ; 11.69 in x 8.27 in ; 842 pt x 595 pt

	frameCount = 2  #(numsites + (-numsites%6)) // 6
	frameWidth = doc.height / frameCount
	frameHeight = doc.width - .05 * inch
	frames = []

	# Construct a frame for each column
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

	Title = ('Report containing {} of {} putative SNPs for {} and {}'.format(
		str(quantityWantSites), str(len(pairwiseSNPpositions)),
		os.path.basename(invcf1), os.path.basename(invcf2)))
	report = [Paragraph(Title, styles['Heading2'])]

	while numSites > 0:
		# handle SNP sites at beginning of Chrom (<50 bases in)
		posMinusFitty = (listWantSitesVCF[numSites - 1] - 50)
		centeredPos = str(max(posMinusFitty, 0))
		pos = str(vcfData1[numSites-1][0])
		gap = Platypus.Spacer(0.25, 0.05*inch)

		# vcfData lists are ordered: [pos, chrom, qual, DP, MQ, alt]
		for sample in [[vcfData1, inbam1, Sample1Name], [vcfData2, inbam2, Sample2Name]]:
			pileup = Platypus.Preformatted(get_ttview(sample[1], infasta,
							(sample[0][1][1]), centeredPos, rows), style)

			qual = sample[0][numSites-1][2]
			DP   = sample[0][numSites-1][3]
			MQ   = sample[0][numSites-1][4]
			SNP  = str(sample[0][numSites-1][5])
			
			header = Paragraph('{} at position {} for {} with a QUAL of {}, MQ of {}, and raw read depth (DP) of {}'.format(
				SNP, pos, sample[2], qual, MQ, DP), styles['Heading6'])
			report.append(Platypus.KeepTogether([gap, header, pileup]))

		numSites -= 1

	doc.build(report)

if __name__ == '__main__':
	main()
