
for inFile in $(ls -1 *merge*ints);
do

    if [[ -e 01_ints2hic_${inFile%.ints}.out ]];
    then
	continue
    fi

    outName=${inFile%.ints}
    outFile=${outName}.hic
    echo "$inFile -> ${outFile}"

    if [[ ! -e ${outFile} ]];
    then
        touch ${outFile}

	outFileVP=${outName}.allValidPairs
	if [[ ! -e ${outFileVP} ]];
	then
            touch ${outFileVP}
            cat ${inFile} | awk 'BEGIN{n=0}{if($1=="chrom1"){next}; for(i=0;i<=$7;i++){printf("R%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\t0\t1\t42\t42\n",n,$1,$2,"+",$4,$5,"+");n++}}' > ${outFileVP}
	fi

	./scripts/hicpro2juicebox.sh -i ${outFileVP} -g ./scripts/chrom_sizes.txt -j ./scripts/juicer_tools_1.22.01.jar -t . -o .
	mv -v ${outFileVP} ${outFile}
    fi

    echo "#Normalizing ${outFile}"
    java -Xmx300G -jar scripts/juicer_tools_1.22.01.jar addNorm -d -F -k SCALE -j 16 $outFile
    
done
