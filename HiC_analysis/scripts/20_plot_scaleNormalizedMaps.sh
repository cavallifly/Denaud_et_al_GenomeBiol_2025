resolution=3000


for locus in dac ;
do
    
    for matrix in $(ls -1 ${locus}*merge*${resolution}bp*tab);
    do
	
	echo $matrix
	xmin=$(echo $matrix | sed -e "s/_/ /g" -e "s/bp//g" | awk '{print $9 }' | awk '{printf("%.2f\n",$1/1000000)}')
	xmax=$(echo $matrix | sed -e "s/_/ /g" -e "s/bp//g" | awk '{print $10}' | awk '{printf("%.2f\n",$1/1000000)}')
	ymin=$(echo $matrix | sed -e "s/_/ /g" -e "s/bp//g" | awk '{print $9 }' | awk '{printf("%.2f\n",$1/1000000)}')
	ymax=$(echo $matrix | sed -e "s/_/ /g" -e "s/bp//g" | awk '{print $10}' | awk '{printf("%.2f\n",$1/1000000)}')	
	
	
	awk '{if(NF==3){print $1; print $2}}' ${matrix} | grep -vi Warn | uniq | sort | uniq | sort -k 2,2n | uniq | awk '{print NR-1,$1}' > _bins
	#head _bins
	#tail _bins    
	
	awk '{if(NF==3){printf("%s\t%s\t%f\n",$1,$2,$3)}}' ${matrix} > _matrix
	#head _matrix
	#wc -l _matrix
	

	size=$(tail -1 _bins | awk '{print $1}')	
	awk -v s=$size 'BEGIN{for(i=0;i<=s;i++){for(j=0;j<=s;j++){m[i,j]=0.0}}}{if(NF==2){bin[$2]=$1}else{v=$3; m[int(bin[$1]),int(bin[$2])]=v; m[int(bin[$2]),int(bin[$1])]=v}}END{for(i=0;i<=s;i++){for(j=0;j<=s;j++){print i,j,m[i,j]}}}' _bins _matrix | sort -k 1,1n -k 2,2n | uniq > _data    
	head _data
	tail _data
	wc -l _data

	cbmin=$(awk '{if(NR==1){min=$3}; if($3<min){min=$3}}END{print min}' _matrix)
	
	for cbmax in 300;
	do
	    outFile=${matrix%.tab}_cbMax_${cbmax}.pdf
	    if [[ ! -e ${outFile} ]];
	    then
		echo $size $cbmin $cbmax
		sed -e "s/XXXcbminXXX/${cbmin}/g" -e "s/XXXxminXXX/${xmin}/g" -e "s/XXXxmaxXXX/${xmax}/g" -e "s/XXXyminXXX/${ymin}/g" -e "s/XXXymaxXXX/${ymax}/g" -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXsizeXXX/${size}/g" ./scripts/20_plot_scaleNormalizedMaps.gp | gnuplot
		file=data.ps
		ps2pdf ${file}
		./scripts/pdfcrop --ini --margin 0 ${file%.ps}.pdf _tmp_${file%.ps}.pdf #&> /dev/null    
		mv _tmp_${file%.ps}.pdf ${file%.ps}.pdf
		rm $file	
		mv data.pdf ${outFile}
	    fi
	done
	rm -fvr _* data.png
    done # Close cycle over $matrix
done
