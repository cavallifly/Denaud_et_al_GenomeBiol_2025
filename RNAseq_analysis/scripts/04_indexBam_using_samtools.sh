#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify 1) input folder containing sorted bam files and 2) output folder for _s.bai index files"
    exit 1
fi

inDir=$1
outDir=$2

mkdir -p ${outDir}

echo "index bam"
for FILE in ${inDir}/*".bam";
do
    filename=${FILE##*/}
    filename=${filename/}
    output=${outDir}/$filename".bai"
    
    if [ -f "$output" ];
    then
        echo "$output already exists"
        continue
    fi

    samtools index $FILE $output
done
