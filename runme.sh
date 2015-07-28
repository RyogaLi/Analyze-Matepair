#!/bin/bash

# sort a bam file


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
	echo "usage: runnme.sh -f <FULL PATH of bam file>"
	echo "-c [optional]: <chromosome name>"
	echo "-t [optional]: <threshold>"
	echo "-rs [optional]: <start position> "
	echo "-re [optional]: <end position>"
	echo ""
	echo "DEFAULT(Only filename is provided): Analyze mate pairs on all chromosomes"
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
echo "========================================"
echo "NOTE: FILE NAME CANNOT BE EMPTY"
echo FILE = "${FILEPATH}"
echo CHROMOSOME  = "${CHROMOSOME}"
echo REGION START     = "${REGIONSTART}"
echo REGION END     = "${REGIONEND}"
echo THRESHOLD    = "${THRESHOLD}"
echo "========================================"

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
 	mkdir results_DEFAULT

 	cd results_DEFAULT
 	for i in "${chrs[@]}"
 	do
 		# creaate directory for each chromosome
 		mkdir $i
 		cd $i
 		python ../../readFile.py "${FILEPATH}" "$i" "" "" "$DEFAULT_THRESH"
 		for f in $FILE
 		do
 			# echo "$f"
 			python ../../plot.py "$f"
		done
 		cd ..
	done

# # case two: chromosome name is provided. Analyze the chromosome
# # one dir will be generated 

# ===================TEST CASE TWO====================== #
# Only chromosome name is provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

# ===================TEST CASE FIVE====================== #
# Chromosome name; threshold were provided 
# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs




elif ! [ -z "$CHROMOSOME" ]
	then 

 	# check if result file exists 
 	if [ ! -d results_"$CHROMOSOME" ] # if no, make a new dir
 		then

 		mkdir results_"$CHROMOSOME"
 	fi


	cd results_"$CHROMOSOME"
	# ===================TEST CASE TWO====================== #
	# Chromosome name; start region; end region; threshold were provided 
	# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs

	if ! [ -z "$REGIONSTART" ] && ! [ -z "$REGIONEND" ] && ! [ -z "$THRESHOLD" ]
		then
		# check if result file exists 
	 	if [ -d results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$THRESHOLD" ] # if yes, remove the old result directory
	 		then
	 		# echo "DIR exists"
	 		rm -f -rf results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$THRESHOLD"
	 	fi
		# if no, make a new directory
	 	mkdir results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$THRESHOLD"
	 	cd results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$THRESHOLD"
		echo "========================================"
 		echo "Analyzing ${CHROMOSOME} in file ${FILEPATH}............................."
 		echo "Region start: ${REGIONSTART}"
 		echo "Region end: ${REGIONEND}"
 		echo "Default threshold: 1000"
 		echo "Press control + z to quit the program"
 		echo "========================================"
		
		python ../../readFile.py "${FILEPATH}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$THRESHOLD"
		echo "DONE"
		FILE="./*"
		for f in $FILE
 		do
 			# echo "$f"
 			python ../../plot.py "$f"
		done

	# ===================TEST CASE THREE====================== #
	# Chromosome name; start region; end region were provided 
	# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs

	elif ! [ -z "$REGIONSTART" ] && ! [ -z "$REGIONEND" ] && [ -z "$THRESHOLD" ]
		then
		# check if result file exists 
	 	if [ -d results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$DEFAULT_THRESH" ] # if yes, remove the old result directory
	 		then
	 		# echo "DIR exists"
	 		rm -f -rf results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$DEFAULT_THRESH"
	 	fi
	 	# if no, make a new directory
		mkdir results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$DEFAULT_THRESH"
		cd results_"$CHROMOSOME"_"$REGIONSTART"_"$REGIONEND"_"$DEFAULT_THRESH"
		echo "========================================"
 		echo "Analyzing ${CHROMOSOME} in file ${FILEPATH}............................."
 		echo "Region start: ${REGIONSTART}"
 		echo "Region end: ${REGIONEND}"
 		echo "Default threshold: 1000"
 		echo "Press control + z to quit the program"
 		echo "========================================"
		python ../../readFile.py "${FILEPATH}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$DEFAULT_THRESH" 
		for f in $FILE
 		do
 			# echo "$f"
 			python ../../plot.py "$f"
		done
		

	# ===================TEST CASE FOUR====================== #
	# Chromosome name; threshold were provided 
	# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.threshold.matepairs
	elif [ -z "$REGIONSTART" ] && [ -z "$REGIONEND" ] && ! [ -z "$THRESHOLD" ]
		then
		# check if result file exists 
	 	if [ -d results_"$CHROMOSOME"_"$THRESHOLD" ] # if yes, remove the old result directory
	 		then
	 		rm -f -rf results_"$CHROMOSOME"_"$THRESHOLD"
	 	fi
	 	mkdir results_"$CHROMOSOME"_"$THRESHOLD"
		cd results_"$CHROMOSOME"_"$THRESHOLD"
		echo "========================================"
 		echo "Analyzing ${CHROMOSOME} in file ${FILEPATH}............................."
 		echo "Region is not specified"
 		echo "Threshold: ${THRESHOLD}"
 		echo "Press control + z to quit the program"
 		echo "========================================"
		python ../../readFile.py "${FILEPATH}" "$CHROMOSOME" "" "" "$THRESHOLD"
		for f in $FILE
 		do
 			# echo "$f"
 			python ../../plot.py "$f"
		done

	# ===================TEST CASE FIVE====================== #
	# Only chromosome name is provided 
	# generate one directory which contains chromosome1.chromosome2.matepairs and chromosome.1000.matepairs
	elif [ -z "$REGIONSTART" ] && [ -z "$REGIONEND" ] && [ -z "$THRESHOLD" ]
		then
	 	# then 
	 	# check if result file exists 
	 	if [ -d results_"$CHROMOSOME"_"$DEFAULT_THRESH" ] # if yes, remove the old result directory
	 		then
	 		rm -f -rf results_"$CHROMOSOME"_"$DEFAULT_THRESH"
	 	fi
	 	mkdir results_"$CHROMOSOME"_"$DEFAULT_THRESH"
		cd results_"$CHROMOSOME"_"$DEFAULT_THRESH"
		echo "========================================"
 		echo "Analyzing ${CHROMOSOME} in file ${FILEPATH}............................."
 		echo "Region is not specified"
 		echo "Default threshold: 1000"
 		echo "Press control + z to quit the program"
 		echo "========================================"
	 	python ../../readFile.py "${FILEPATH}" "$CHROMOSOMEf" "" "" "$DEFAULT_THRESH"
		for f in $FILE
 		do
 			# echo "$f"
 			python ../../plot.py "$f"
		done
	
	# ===================TEST CASE SIX====================== #
	# ERROR ONE: only start/end position is entered
	# ERROR TWO: no file name 
	else
		echo "Please provide valid parameters"
		echo ""
		echo "usage: runnme.sh -f <full path of bam file>"
		echo "-c [optional]: <chromosome name>"
		echo "-t [optional]: <threshold>"
		echo "-rs [optional]: <start position> "
		echo "-re [optional]: <end position>"
		echo ""
		echo "DEFAULT(Only filename is provided): Analyze mate pairs on all chromosomes"
		echo ""
		exit

	fi
else
	echo "Please provide valid parameters"
	echo ""
	echo "usage: runnme.sh -f <full path of bam file>"
	echo "-c [optional]: <chromosome name>"
	echo "-t [optional]: <threshold>"
	echo "-rs [optional]: <start position> "
	echo "-re [optional]: <end position>"
	echo ""
	echo "DEFAULT(Only filename is provided): Analyze mate pairs on all chromosomes"
	echo ""
	exit

fi
