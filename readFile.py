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
## count the number of mate pairs that are far away from each other on the same chromosome according to threshold

def parseFile(filename, chromosome, start, end, threshold):
	"""
	Return one or several text file(s) that contains the mate pair 
	mapping information
	
	@param
	filename: a sorted.mapped.bam file (CANNOT BE NONE)
	chromosome: chromosome that you want to target (NONE if you are targeting 
		all the chromosomes)
	region: a list with len == 2 that indicates which region on that chromosome i.e [100,200] means from 100bp to 200bp
		(NONE if you are targeting all the chromosomes)
	threshold: a value which is used to select the distance between a mate pair

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
	# print(chromosome,start, end)
	# if we want to focus on a region on one sepecific chromosome
	if chromosome != "" and start != "" and end != "":
		start = int(start)
		end = int(end)
		# fetch the reads within region on chromosome
		print ("Finding mate pairs .... This step will take a while")
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
				else: # if both mates are on same chromosome
					fName = chromosome+"."+str(threshold)+".matepairs"
					f = open(fName, "a")
					read = str(read).split()
					mate = str(mate).split()
					if (int(read[3]) - int(mate[3])) >= threshold:
						f.write(str(read)+"\n")
						f.write(str(mate)+"\n")
				# readPairs.append((read,mappedBam.mate(read)))
	elif chromosome != "" and start == "" and end == "":
		print ("Finding mate pairs .... This step will take a while")
		# fetch the reads on chromosome
		for read in mappedBam.fetch(chromosome):
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
				else: # if both mates are on same chromosome
					fName = chromosome+"."+str(threshold)+".matepairs"
					f = open(fName, "a")
					read = str(read).split()
					mate = str(mate).split()
					if (int(read[3]) - int(mate[3])) >= threshold:
						f.write(str(read)+"\n")
						f.write(str(mate)+"\n")



if __name__ == '__main__':
	bamFile = sys.argv[1]
	chromosome = sys.argv[2]
	start = sys.argv[3]
	end = sys.argv[4]
	threshold = sys.argv[5]

	parseFile(bamFile, chromosome, start, end, threshold)
	quit()
