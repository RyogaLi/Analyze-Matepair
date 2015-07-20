#!/bin/python2.7

import pysam
import sys
import os

# chr ID = [0,1,2,...,16]
# take a sorted.mapped.bam file
# find matepair that is on desired chromosome and region
# if no region/chromosome is provided
# for each chromosome find out all the matepairs
## count the total number of mate pairs that are on other chromosomes txt/cls?
## count the number of mate pairs that are far away from each other according to threshold

def parseFile(filename, chromosome, region, threshold):
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
	ID_Name = {0:"chrI", 1:"chrII", 2:"chrIII", 3:"chrIV", 4:"chrV", 5:"chrVI", 6:"chrVII", 7:"chrVIII", 8:"chrIX", 9:"chrX", 10:"chrXI", 11:"chrXII", 12:"chrXIII", 13:"chrXIV", 14:"chrXV", 15:"chrXVI", 16:"chrM"}

	# open a bam file
	mappedBam = pysam.AlignmentFile(filename,"rb")


	# if we want to focus on all the chromosomes
	if chromosome != None and region != None:
		start = region[0]
		end = region[1]
		# readPairs = []
		# fName = chromosome+"."+
		# f = open(fName)
		# fetch the reads within region on chromosome
		for read in mappedBam.fetch(chromosome, start, end):
			# check if the mate is mapped or not 
			if not read.mate_is_unmapped:
				# find it's mate pair
				mate = mappedBam.mate(read)
				# if mate pair is on another chromosome
				if mate.reference_id != read.reference_id:
					# make a new file and store the mate pairs 
					fName = chromosome+"."+ID_Name[mate.reference_id]+".matepairs"
					f = open(fName, "a")
					f.write(str(read)+"\n")
					f.write(str(mate)+"\n")
				else:
					fName = chromosome+"."+str(threshold)+".matepairs"
					f = open(fName, "a")
					read = str(read).split()
					mate = str(mate).split()
					if (int(read[3]) - int(mate[3])) >= threshold:

						f.write(str(read)+"\n")
						f.write(str(mate)+"\n")
				# readPairs.append((read,mappedBam.mate(read)))
	# else:

if __name__ == '__main__':
	parseFile("WY1769.mapped.bam", "chrV", [246500,249500], 100)
