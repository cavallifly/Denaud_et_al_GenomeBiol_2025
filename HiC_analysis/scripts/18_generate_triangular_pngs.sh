size=800

for inFile in $(ls -1 *toModify*.png);
do
    outFile=${inFile%_toModify.png}.png
    echo "Input file : ${outFile}"
    echo "Cropping..."
    convert -background white -alpha remove ${inFile} -crop ${size}x0 ${outFile}

    if [[ ! -e ${outFile%.png}-0.png ]];
    then
	for file in $(ls -1 ${outFile%.png}.png) ;
	do
	    convert +repage ${file} _tmp ; mv _tmp ${file}
	    identify ${file}
	    echo "Trimming..."
	    convert -trim ${file} ${file%.png}_trimmed.png
	    convert +repage ${file%.png}_trimmed.png _tmp ; mv _tmp ${file%.png}_trimmed.png
	    identify ${file%.png}_trimmed.png
	    echo "Rotating..."
	    convert -rotate 45 -background white -alpha remove ${file%.png}_trimmed.png _tmp_${file}
	    convert +repage _tmp_${file} _tmp; mv _tmp _tmp_${file}
	    identify _tmp_${file}
	    echo "Cropping..."
	    convert -background white -alpha remove _tmp_${file} -crop 1042x521-0-0 ${file%.png}_triang.png
	    convert +repage ${file%.png}_triang.png _tmp ; mv _tmp ${file%.png}_triang.png
	    identify ${file%.png}_triang.png
	    rm -fr ${file} _tmp*png *_trimmed.png
	done
    fi  
    
    if [[ -e ${outFile%.png}-0.png ]];
    then
	for file in $(ls -1 ${outFile%.png}-?.png) ;
	do
	    convert +repage ${file} _tmp ; mv _tmp ${file}
	    identify ${file}
	    echo "Trimming..."
	    convert -trim ${file} ${file%.png}_trimmed.png
	    convert +repage ${file%.png}_trimmed.png _tmp ; mv _tmp ${file%.png}_trimmed.png
	    identify ${file%.png}_trimmed.png
	    echo "Rotating..."
	    convert -rotate 45 -background white -alpha remove ${file%.png}_trimmed.png _tmp_${file}
	    convert +repage _tmp_${file} _tmp; mv _tmp _tmp_${file}
	    identify _tmp_${file}
	    echo "Cropping..."
	    convert -background white -alpha remove _tmp_${file} -crop 1042x521-0-0 ${file%.png}_triang.png
	    convert +repage ${file%.png}_triang.png _tmp ; mv _tmp ${file%.png}_triang.png
	    identify ${file%.png}_triang.png
	    rm -fr ${file} _tmp*png *_trimmed.png
	done
    fi

    rm -fr ${outFile} ${inFile}
done

# Rename triangular matrices
#LD_larvae_DWT_vs_larvae_double_k250_kexp250
mv -v diffMaps_chr2L_16350000_16500000_r3000bp_k250_kexp250_Fig1F_annot_FALSE_triang.png diffMaps_chr2L_16350000_16500000_r3000bp_k250_kexp250_Fig1F_annot_FALSE_DWTvsdouble_triang.png 
mv -v diffMaps_chr2L_16350000_16500000_r3000bp_k250_kexp250_Fig1F_annot_TRUE_triang.png  diffMaps_chr2L_16350000_16500000_r3000bp_k250_kexp250_Fig1F_annot_TRUE_DWTvsdouble_triang.png
#LD_larvae_DPRE1Up_vs_larvae_DPRE1_k250_kexp250
mv -v diffMaps_chr2L_16300000_16600000_r3000bp_k250_kexp250_Fig2G_annot_FALSE_triang.png diffMaps_chr2L_16300000_16600000_r3000bp_k250_kexp250_Fig2G_annot_FALSE_DPRE1UpvsDPRE1_triang.png
mv -v diffMaps_chr2L_16300000_16600000_r3000bp_k250_kexp250_Fig2G_annot_TRUE_triang.png  diffMaps_chr2L_16300000_16600000_r3000bp_k250_kexp250_Fig2G_annot_TRUE_DPRE1UpvsDPRE1_triang.png

mv -v scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_FALSE-0_triang.png scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_FALSE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_FALSE-1_triang.png scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_FALSE-larvae_double_triang.png
mv -v scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_FALSE-0_triang.png scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_FALSE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_FALSE-1_triang.png scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_FALSE-larvae_double_triang.png
mv -v scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_TRUE-0_triang.png  scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_TRUE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_TRUE-1_triang.png  scoreMapsk250kexp500_chr2L_16000000bp_17000000bp_r3000bp_Fig1C_annot_TRUE-larvae_double_triang.png
mv -v scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_TRUE-0_triang.png  scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_TRUE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_TRUE-1_triang.png  scoreMapsk250kexp500_chr2L_16321000bp_16512000bp_r3000bp_Fig1C_annot_TRUE-larvae_double_triang.png

mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-0_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-1_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-larvae_DPRE1_triang.png
mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-2_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_FALSE-larvae_DPRE1Up_triang.png
mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-0_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-larvae_DWT_triang.png
mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-1_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-larvae_DPRE1_triang.png
mv -v scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-2_triang.png scoreMapsk250kexp500_chr2L_16300000bp_16600000bp_r3000bp_Fig2C_annot_TRUE-larvae_DPRE1Up_triang.png



