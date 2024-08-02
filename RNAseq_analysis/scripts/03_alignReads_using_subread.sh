#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify"
    echo "1) path to folder containing fastq files"
    echo "2) path to folder for output files"
    echo "3) path to the subread index folder"
    exit 1
fi

inDir=$1
outDir=$2
subreadIndex=$3

mkdir -p ${outDir}

echo "Doing the read alignment on dm6 assembly using subread"
for FILE in ${inDir}/*"_1.fq.gz";
do
    filename=${FILE##*/}
    filename_1=$filename
    filename_2=${filename/"1.fq.gz"/"2.fq.gz"}
    output=$2${filename/".fq.gz"/".bam"}
    
    echo $filename_1
    echo $filename_2
    echo $output
    
    if [ -f "$output" ];
    then
        echo "$output already exists"
        continue
    fi
    
    if [[ "$FILE" = *".bam" ]];
    then
        echo "$FILE already exists"
        continue
    fi

    subread-align -i subread_index_dmelr634 -t 0 -r ${inDir}/$filename_1 -R ${inDir}/$filename_2 -o $output -i ${subreadIndex} -T 10  --sortReadsByCoordinates

done

