assembly=dm6
mishaDB=./mishaDB/trackdb/${assembly}/
scriptsDir=${PWD}/scripts/

chunkLength=5000000
echo "Working with chunks of ${chunkLength}bp"

nChromosomes=$(wc -l ${mishaDB}/chrom_sizes.txt | awk '{print $1}')

for track in hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1 hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1 hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1 /
	     hic.hic_larvae_double_Rep1_dm6_BeS.hic_larvae_double_Rep1_dm6_BeS_1 hic.hic_larvae_double_Rep2_dm6_BeS.hic_larvae_double_Rep2_dm6_BeS_1 hic.hic_larvae_double_Rep3_dm6_BeS.hic_larvae_double_Rep3_dm6_BeS_1 /
	     hic.hic_larvae_DPRE1_Rep1_dm6_BeS.hic_larvae_DPRE1_Rep1_dm6_BeS_1 hic.hic_larvae_DPRE1_Rep2_dm6_BeS.hic_larvae_DPRE1_Rep2_dm6_BeS_1 /
	     hic.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS_1 hic.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS_1 ;
do

    wDir=_tmp_wDir_${track}
    if [[ -d ${wDir} ]];
    then
	continue
    fi
    echo $wDir
    mkdir -p ${wDir}
    cd ${wDir}

    name=$(echo $track | sed -e "s,\., ,g" | awk '{print $NF}')
    echo $name
    for nC1 in $(seq 1 1 ${nChromosomes}) ;
    do
	c1=$(awk -v n=${nC1} '{if(NR==n){print $1}}' ${mishaDB}/chrom_sizes.txt)
	chrom1=chr${c1}       
	start1=0
	end1=$(grep -w ${chrom1//chr/} ${mishaDB}/chrom_sizes.txt | awk '{print $2}')
	echo $chrom1 $start1 $end1
	
	for nC2 in $(seq 1 1 ${nChromosomes}) ;
	do
	    if [[ $nC2 -gt $nC1 ]];
	    then
		continue
	    fi
	    c2=$(awk -v n=${nC2} '{if(NR==n){print $1}}' ${mishaDB}/chrom_sizes.txt)
	    chrom2=chr${c2}
	    start2=0
	    end2=$(grep -w ${chrom2//chr/} ${mishaDB}/chrom_sizes.txt | awk '{print $2}')
	    echo $chrom2 $start2 $end2
	    
	    outFile=${name}_${chrom1}_${start1}_${end1}bp_${chrom2}_${start2}_${end2}bp_noDup.ints
	    if [[ -s ${outFile} ]];
	    then
		echo "File ${outFile} exists! Go to the next job!"
		continue
	    fi
	    
	    # Estimate nchunks
	    nchunk=0
	    for s1 in $(seq ${start1} ${chunkLength} ${end1});
	    do
		e1=$(awk -v s=${s1} -v cl=${chunkLength} 'BEGIN{print s+cl}')
		for s2 in $(seq ${start2} ${chunkLength} ${end2});
		do
		    e2=$(awk -v s=${s2} -v cl=${chunkLength} 'BEGIN{print s+cl}')
		    nchunk=$((${nchunk}+1))
		    echo $s1 $e1 $s2 $e2
		done # Close cycle over $s2
	    done # Close cycle over $s1
	    echo "Number of chunks ${nchunk}"

	    # Obtaining the matrices
	    nchunk=0
	    for s1 in $(seq ${start1} ${chunkLength} ${end1});
	    do
		e1=$(awk -v s=${s1} -v cl=${chunkLength} 'BEGIN{print s+cl}')
		for s2 in $(seq ${start2} ${chunkLength} ${end2});
		do
		    e2=$(awk -v s=${s2} -v cl=${chunkLength} 'BEGIN{print s+cl}')

		    if [[ ${chrom1} == ${chrom2} ]];
		    then
			if [[ $s1 -lt $s2 ]];
			then
			    continue
			fi
		    fi
		    
		    nchunk=$((${nchunk}+1))
		    (
			echo ${nchunk} ${chrom1} ${s1} ${e1} ${chrom2} ${s2} ${e2}
			rm -fvr *_chunk_${nchunk}.ints
			
			sed -e "s/XXXchrom1XXX/${chrom1}/g" -e "s/XXXstart1XXX/${s1}/g" -e "s/XXXend1XXX/${e1}/g" -e "s/XXXchrom2XXX/${chrom2}/g" -e "s/XXXstart2XXX/${s2}/g" -e "s/XXXend2XXX/${e2}/g" -e "s/XXXnchunkXXX/${nchunk}/g" -e "s/XXXtrackXXX/${track}/g" ${scriptsDir}/17_retrieve_points_from_mishaDB.R > _tmp_${nchunk}.R
			Rscript _tmp_${nchunk}.R &> /dev/null
			outfile=$(ls -1 *_chunk_${nchunk}.ints)
			echo "Output file ${outfile}"
			
			#awk '{printf("%s %s %s %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7)}' ${outfile} > _tmp_${nchunk} ; mv _tmp_${nchunk} ${outfile}
			awk '{if($1!=$4){print $0}else{if($2>=$5){print $0}}}' ${outfile} > _tmp_${nchunk} ; mv _tmp_${nchunk} ${outfile}
			rm _tmp_${nchunk}.R
		    ) #&
		    check=$(awk -v nc=${nchunk} 'BEGIN{if(nc%10==0){print 1}else{print 0}}')

		    #echo "Check $check"
		    if [[ ${check} -eq 1 ]];
		    then
			echo "Waiting for the first ${nchunk} to finish!"
			wait
		    fi
		done # Close cycle over $s2
	    done # Close cycle over $s1
	    #wait
	    
	    cat *_chunk_*.ints | awk '{if(NR != 1 && $1=="chrom1"){next}; print $0}' | sort -k 2n,2n -k 5n,5n | uniq > ${outFile}
	    rm -fvr Observed_merged_*_chunk_*.ints

	done # Exit cycle over ${chrom2}    
    done # Exit cycle over ${chrom1}
    cd ..
done # Exit cycle over ${track}
#exit

echo "Generating .ints files"
#for sample in hic_larvae_DWT_Rep1_dm6_BeS hic_larvae_DWT_Rep2_dm6_BeS hic_larvae_DWT_Rep3_dm6_BeS hic_larvae_DWT_merge_dm6_BeS ;
#for sample in hic_larvae_double_Rep1_dm6_BeS hic_larvae_double_Rep2_dm6_BeS hic_larvae_double_Rep3_dm6_BeS hic_larvae_double_merge_dm6_BeS ;
for sample in hic_LE_WT_Rep1_dm6_YO hic_LE_WT_Rep2_dm6_YO hic_LE_WT_Rep3_dm6_VL hic_LE_WT_merge_dm6_YOVL ;
#for sample in hic_larvae_DPRE1_Rep1_dm6_BeS hic_larvae_DPRE1_Rep2_dm6_BeS hic_larvae_DPRE1_merge_dm6_BeS ;
#for sample in hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS ;
#for sample in hic_larvae_double_merge_dm6_BeS ;	      
do
    name=$(echo $sample | sed -e "s/_merge/ /g" -e "s/res_dm6/ /g" | awk '{print $1}')
    # NOTE: the (Rep1 and Rep1res) and (Rep2 and Rep2res) are merged together.
    sample=$(echo $sample | sed -e "s/res_dm6/_dm6/g")
    echo $name $sample
    
    outFile=${sample}.ints
    if [[ -e $outFile ]];
    then
	#awk -v s=${sample} '{cnt+=$7}END{print s,"valid-pairs",cnt}' ${outFile}
        continue
    fi
    head -1 ./*${name}_*/*chr2L*chr2L* | grep chrom1 | head -1 | awk '{for(i=1;i<=7;i++){printf("%s\t",$i)}; printf("\n")}' | uniq > ${outFile}
    head ${outFile}
    
    ls -lrtha ./*${name}_*/*ints
    for nC1 in $(seq 1 1 ${nChromosomes}) ;
    do
        c1=$(awk -v n=${nC1} '{if(NR==n){print "chr"$1}}' ${mishaDB}/chrom_sizes.txt)
	for nC2 in $(seq 1 1 ${nChromosomes}) ;
	do
            if [[ $nC2 -gt $nC1 ]];
            then
		continue
            fi
            c2=$(awk -v n=${nC2} '{if(NR==n){print "chr"$1}}' ${mishaDB}/chrom_sizes.txt)
	    echo $c1 $c2
	    if [[ ${nC2} -eq ${nC1} ]];
	    then
		ls -lrtha ./_tmp*${name}_*/*${c1}_*${c2}_*ints
		cat ./_tmp*${name}_*/*${c1}_*${c2}_*ints | grep -v chrom1 | awk '{if(NF<7){next}; pair=$1"_"$2"_"$3"_"$4"_"$5"_"$6; h[pair]+=$7}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | awk '{for(i=1;i<=7;i++){printf("%s\t",$i)}; printf("\n")}' | sort -k 1,1d -k 4,4d -k 2,2n -k 5,5n >> ${outFile}
	    else
		ls -lrtha ./_tmp*${name}_*/*${c1}_*${c2}_*ints ./_tmp*${name}_*/*${c2}_*${c1}_*ints 2> /dev/null		
		cat ./_tmp*${name}_*/*${c1}_*${c2}_*ints ./_tmp*${name}_*/*${c2}_*${c1}_*ints 2> /dev/null | grep -v chrom1 | awk '{if(NF<7){next}; pair=$1"_"$2"_"$3"_"$4"_"$5"_"$6; h[pair]+=$7}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | awk '{for(i=1;i<=7;i++){printf("%s\t",$i)}; printf("\n")}' | sort -k 1,1d -k 4,4d -k 2,2n -k 5,5n >> ${outFile}
	    fi
	done
    done
done
