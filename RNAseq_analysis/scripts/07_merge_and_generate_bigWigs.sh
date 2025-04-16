#!/bin/bash

bamFiles=$(ls -1 WT*.bam)
output=WT_allRep.bam

# Merging bam files
samtools merge -o $output $bamFiles
# Indexing merged bam file for WT
samtools index ${output} ${output%.bam}.bai
# Generating .bigWig files
bamCoverage --normalizeUsing BPM -e 0 -bs 10 -b ${output} -o ${output%.bam}bs10.BPM.bigWig


bamFiles=$(ls -1 D*.bam)

output=double_allRep.bam
# Merging bam files
samtools merge -o $output $bamFiles
# Indexing merged bam file for double
samtools index ${output} ${output%.bam}.bai
# Generating .bigWig files
bamCoverage --normalizeUsing BPM -e 0 -bs 10 -b ${output} -o ${output%.bam}bs10.BPM.bigWig
