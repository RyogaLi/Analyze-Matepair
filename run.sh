# ===================TEST CASE ONE====================== #
# No chromosome is provided. Run default parameters
# generate 16 directories and each directory contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

# ===================TEST CASE TWO====================== #
# Only chromosome name is provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

# ===================TEST CASE THREE====================== #
# Chromosome name; start region; end region were provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

# ===================TEST CASE FOUR====================== #
# Chromosome name; start region; end region; threshold were provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs

# ===================TEST CASE FIVE====================== #
# Chromosome name; threshold were provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs

# ===================TEST CASE SIX====================== #
# Only threshold is provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs


## TEST CASE ONE
## ./runme.sh -f WY1769.mapped.bam 

# TEST CASE TWO
./runme.sh -f /users/RyogaLi/Desktop/WY1769.mapped.bam

## TEST CASE THREE
##./runme.sh -f WY1769.mapped.bam -c chrV -rs 20000 -re 25000

## TEST CASE FOUR
## ./runme.sh -f WY1769.mapped.bam -c chrV -t 5000

## TEST CASE FIVE 
## ./runme.sh -f WY1769.mapped.bam -c chrV