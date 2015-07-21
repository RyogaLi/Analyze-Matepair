#!/bin/bash

# default: analyze all the chromosomes
# store all the chromosomes' name in an array
chrs=("chrI" "chrII" "chrIII" "chrIV" "chrV" "chrVI" "chrVII" "chrVIII" "chrIX" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI" "chrM")
DEFAULT_THRESH=1000
# user did not put filename 
if [[ "$#" = 0 ]]
	then
	echo ""
	echo "ERROR: Filename cannot be none."
	echo ""
	echo "usage: runnme.sh -f <full path of bam file>"
	echo "-c [optional]: <chromosome name>"
	echo "-t [optional]: <threshold>"
	echo "-rs [optional]: <start position> "
	echo "-re [optional]: <end position>"
	echo ""
	echo "DEFAULT: Analyze mate pairs on all chromosomes"
	echo ""
	exit
fi

while [[ $# > 0 ]]
do
	key="$1"

	case $key in
		-f|--filepath)
	    FILEPATH="$2"
	    shift # past argumen0t
	    ;;
	    -c|--chromosome)
	    CHROMOSOME="$2"
	    shift # past argument
	    ;;
	    -rs|--regionstart)
	    REGIONSTART="$2"
	    shift # past argument
	    ;;
	    -re|--regionend)
	    REGIONEND="$2"
	    shift # past argument
	    ;;
	    -t|--threshold)
	    THRESHOLD="$2"
	    shift # past argument
	    ;;
	esac
	shift # past argument or value
done
# echo "========================================"
# echo FILE = "${FILEPATH}"
# echo CHROMOSOME  = "${CHROMOSOME}"
# echo REGION START     = "${REGIONSTART}"
# echo REGION END     = "${REGIONEND}"
# echo THRESHOLD    = "${THRESHOLD}"
# echo "========================================"

# case one: only filename is provided. All the chromosmes will be analyzed 
# one directory will be generated for each chromosome
# threshold will be set to 1000

# ===================TEST CASE ONE====================== #
# No chromosome is provided. Run default parameters
# generate 16 directories and each directory contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

if [ -z "$CHROMOSOME" ] && [ -z "$REGIONSTART" ] && [ -z "$REGIONEND" ] && [ -z "$THRESHOLD" ]
 	then
 	echo "========================================"
 	echo "Analyzing all the chromosomes in file ${FILEPATH}............................."
 	echo "Default threshold: 1000"
 	echo "Press control + z to quit the program"
 	echo "========================================"

 	# check if result file exists 
 	if [ -d results_DEFAULT ] # if yes, remove the old result directory
 		then
 		# echo "DIR exists"
 		rm -f -rf results_DEFAULT
 	fi
 	if [ ! -d results_DEFAULT ] # if no, make a new directory
 		then 
 		# echo "DIR not exists"
 		mkdir results_DEFAULT
 	fi

 	cd results_DEFAULT
 	for i in "${chrs[@]}"
 	do
 		# creaate directory for each chromosome
 		mkdir $i
 		cd $i
 		python ../../readFile.py ../../"${FILEPATH}" "$i" "" "" "$DEFAULT_THRESH"
 		cd ..
	done

# # case two: chromosome name is provided. Analyze the chromosome
# # one dir will be generated 

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
elif ! [ -z "$CHROMOSOME" ]
	then 

 	# check if result file exists 
 	if [ -d results_"$CHROMOSOME" ] # if yes, remove the old result directory
 		then
 		# echo "DIR exists"
 		rm -f -rf results_"$CHROMOSOME"
 	fi
 	if [ ! -d results_"$CHROMOSOME" ] # if no, make a new directory
 		then 
 		# echo "DIR not exists"
 		mkdir results_"$CHROMOSOME"
 	fi

	cd results_"$CHROMOSOME"
	# ===================TEST CASE TWO====================== #
	# Chromosome name; start region; end region; threshold were provided 
	# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs

	if ! [ -z "$REGIONSTART"] && ! [ -z "$REGIONEND"] && ! [ -z "$THRESHOLD"]
		then
		echo "========================================"
 		echo "Analyzing ${CHROMOSOME} in file ${FILEPATH}............................."
 		echo "Region start: ${REGIONSTART}"
 		echo "Region end: ${REGIONEND}"
 		echo "Default threshold: ${THRESHOLD}"
 		echo "Press control + z to quit the program"
 		echo "========================================"
		
		python ../readFile.py ../"${FILEPATH}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$THRESHOLD"
	fi
	# elif ! [ -z "$REGIONSTART"] && ! [ -z "$REGIONEND"] && [ -z "$THRESHOLD"]
	# 	then
	# 	python ../readFile.py "${FILENAME}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$DEFAULT_THRESH"
	# else 
	# 	# then 
	# 	python ../readFile.py "${FILENAME}" "" "" "" "$DEFAULT_THRESH"
	# cd ..
	# fi
fi
