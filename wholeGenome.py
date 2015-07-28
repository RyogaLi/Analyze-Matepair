## Analyze the mate pairs in the whole genome
## Generate a histogram for the distribution of distance between two ends of the mate pairs

import pysam
import matplotlib.pyplot as plt


def rFile(filename):
	"""
	Range: 0-200;200-400;400-600;600-800;800-1000;1000-1200;1200-1400;1400-1800;1800-2000;2000-100000;100000+

	"""
	mappedBam = pysam.AlignmentFile(filename,"rb")

	for line in mappedBam:
		line = line.split()
		data = []
		# check if the mate pair is mapped to some where in the genome
		if not read.mate_is_unmapped:
			# if they are on same chromosome 
			if line[2] == line[6]:
				# calculate mate pair distance
				dis = abs(line[7] - line[3])
				data.append(dis)

			else:
				dis = 2000000
				data.append(dis)
	return data


def plotHistogram(data):
	plt.title('Mate pairs analysis')
	plt.hist(data, bins=[0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 100000])


if __name__ == '__main__':
	data = rFile()
	plotHistogram(data)