resolution=3000

for hicFile in $(ls -1 *merge*.hic);
do
    name=$(echo $hicFile | sed -e "s/\.hic//g")
    echo $hicFile $name

    echo "# Dumping the matrix for panel S1d"
    start=16000000
    end=17000000
    java -Xmx300G -jar scripts/juicer_tools_1.22.01.jar dump observed scale ${name}.hic 2L:${start}:${end} 2L:${start}:${end} BP ${resolution} > dac_${name}_chr2L_${start}_${end}bp_scale_${resolution}bp_FigS1D.tab
    
    start=16321000
    end=16512000
    java -Xmx300G -jar scripts/juicer_tools_1.22.01.jar dump observed scale ${name}.hic 2L:${start}:${end} 2L:${start}:${end} BP ${resolution} > dac_${name}_chr2L_${start}_${end}bp_scale_${resolution}bp_FigS1D.tab  

    echo "# Dumping the matrix for panel S2b"
    start=16300000
    end=16600000
    java -Xmx300G -jar scripts/juicer_tools_1.22.01.jar dump observed scale ${name}.hic 2L:${start}:${end} 2L:${start}:${end} BP ${resolution} > dac_${name}_chr2L_${start}_${end}bp_scale_${resolution}bp_FigS2B.tab
    
done # Close cycle over $hicFile

