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
		-f|--filename)
	    FILENAME="$2"
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
echo FILE = "${FILENAME}"
echo CHROMOSOME  = "${CHROMOSOME}"
echo REGION START     = "${REGIONSTART}"
echo REGION END     = "${REGIONEND}"
echo THRESHOLD    = "${THRESHOLD}"
echo "========================================"

# case one: only filename is provided. All the chromosmes will be analyzed 
# one directory will be generated for each chromosome
# threshold will be set to 1000
if [ -z "$CHROMOSOME" ] && [ -z "$REGIONSTART" ] && [ -z "$REGIONEND" ] && [ -z "$THRESHOLD" ]
 	then
 	echo "Analyzing all the chromosomes............................."
 	echo "Press control + z to quit the program"
 	rm -f -rf results_DEFAULT
 	mkdir results_DEFAULT
 	cd results_DEFAULT
 	for i in "${chrs[@]}"
 	do
 		# creaate directory for each chromosome
 		mkdir $i
 		cd $i
 		python ../../readFile.py "${FILENAME}" "$i" "" "" "$DEFAULT_THRESH"
 		cd ..
	done

# case two: chromosome name is provided. Analyze the chromosome
# one dir will be generated 
elif ! [ -z "$CHROMOSOME" ]
	then 
		rm -f -rf "$CHROMOSOME"
		mkdir "$CHROMOSOME"
 		cd "$CHROMOSOME"
 		if ! [ -z "$REGIONSTART"] && ! [ -z "$REGIONEND"] && ! [ -z "$THRESHOLD"]
 			then
 			python ../readFile.py "${FILENAME}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$THRESHOLD"
 		elif ! [ -z "$REGIONSTART"] && ! [ -z "$REGIONEND"] && [ -z "$THRESHOLD"]
 			then
 			python ../readFile.py "${FILENAME}" "${CHROMOSOME}" "$REGIONSTART" "$REGIONEND" "$DEFAULT_THRESH"
 		else 
 			# then 
 			python ../readFile.py "${FILENAME}" "" "" "" "$DEFAULT_THRESH"
 		cd ..
 		fi
fi
