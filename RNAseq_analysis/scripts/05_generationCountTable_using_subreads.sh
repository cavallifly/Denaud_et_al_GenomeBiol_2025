#!/bin/bash
if [ "$#" -ne 3 ];
then
        echo "Please specify 1) input folder containing sorted bam files and 2) path of the .gtf file 3) output folder for count files"
	exit 1
fi

inDir=$1
gtfFile=$2
outDir=$3

mkdir -p ${outDir}

echo "Gerating count-table using featureCounts of subread"
for FILE in ${inDir}/*".bam";
do
    filename=${FILE##*/}
    output=${outDir}/${filename/".bam"/".tsv"}
    echo "$FILE"
    echo "$output"
    
    featureCounts="featureCounts -p --countReadPairs -s 2 -t exon -g gene_id -a ${gtfFile} -o $output $FILE "
    echo "Running command $featureCounts"
    $featureCounts
done
