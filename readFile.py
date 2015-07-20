#!/bin/python2.7

import pysam
import sys

# take a sorted.mapped.bam file
# find matepair that is on desired chromosome and region
# if no region/chromosome is provided
# for each chromosome find out all the matepairs
## count the total number of mate pairs that are on other chromosomes txt/cls?
## count the number of mate pairs that are far away from each other according to threshold

def parseFile(filename, chromosome, region):
	"""
	Return one or several text files that contains the mate pair 
	mapping information
	
	@param
	filename: a sorted.mapped.bam file (CANNOT BE NONE)
	chromosome: chromosome that you want to target (NONE if you are targeting 
		all the chromosomes)
	region: a list with len == 2 that indicates which region on that chromosome i.e [100,200] means from 100bp to 200bp
		(NONE if you are targeting all the chromosomes)

	@return
	chromosome1.chromosome2.matepairs: contains the mate pair that maps to 
		another chromosome 
	chromosome1.distance.matepairs: contains the mate pair that both ends 
		mapped to same chromosome but are far away from each other (threshld 
		indicated by distance)
	"""

	# open a bam file
	mappedBam = pysam.AlignmentFile(filename,"rb")


	# if we want to focus on all the chromosomes
	if chromosome != None and region != None:
		start = region[0]
		end = region[1]
		# fetch the reads within region on chromosome
		for read in samFile.fetch(chromosome, start, end):
			# check if the mate is mapped or not 
			if not read.mate_is_unmapped:
				readsonchrV.append(read)
	else:
		

