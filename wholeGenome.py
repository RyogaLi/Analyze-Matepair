## Analyze the mate pairs in the whole genome
## Generate a histogram for the distribution of distance between two ends of the mate pairs

import pysam
import matplotlib.pyplot as plt


def rFile(filename):
	"""
	Range: 0-200;200-400;400-600;600-800;800-1000;1000-1200;1200-1400;1400-1800;1800-2000;2000-100000;100000+

	"""
	dic = {(0,200):0,(201,400):0,(401,600):0,(601,800):0,
			(801,1000):0,(1001,1200):0,(1201,1400):0,(1401,1600):0,1800:0,2000:0,3000:0,}

	mappedBam = pysam.AlignmentFile(filename,"rb")

	for line in mappedBam:
		line = line.split()
		# check if the mate pair is mapped to some where in the genome
		if not read.mate_is_unmapped:
			# if they are on same chromosome 
			if line[2] == line[6]:
				# calculate mate pair distance
				dis = abs(line[7] - line[3])

			else:
				dis = 2000000


def plotHistogram():
	pass


if __name__ == '__main__':
	main()