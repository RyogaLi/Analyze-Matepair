import sys
import matplotlib.pyplot as plt
chrs=["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX",
"chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM"]
def plot(filename):
        """
        Take out the mapping position (start position) from a formated text file and save them in (x,y)
        """
        x = []
        y = []
        with open(filename, "r") as inputfile:
                for line in inputfile:
                        line = line.split()
                        x.append(int(line[7]))
                        y.append(int(line[3]))
                        inputfile.next()
        filename2 = filename.split('/')[1].split('.')
        plt.clf()
        plt.xlabel(filename2[1])
        if filename2[0] not in chrs:
                plt.ylabel(filename2[1])
        else:
                plt.ylabel(filename2[0])
        plt.title('Mate pairs')
        plt.plot(x, y, '.')
        figure = plt.gcf()
        figure.savefig(filename2[0]+'.'+filename2[1] + '.png', dpi=200)
        return
        
if __name__ == '__main__':
        f = sys.argv[1]
        plot(f)

