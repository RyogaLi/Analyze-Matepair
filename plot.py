import matplotlib.pyplot as plt
def coordinates(filename):
	"""
	Take out the mapping position (start position) from a formated text file and save them in (x,y)
	"""
	x = []
	y = []
	with open(filename, "r") as inputfile:
		for line in inputfile:
			line = line.split()
			x.append(line[3])
			y.append(line[7])
			inputfile.next()
		return (x, y)

def generateGraph(tuple):
	"""
	"""


