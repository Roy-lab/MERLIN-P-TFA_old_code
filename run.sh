#!/bin/sh

if [ $# != 6 ]
then
	echo "Incorrect number of arguments $#"
	echo "Usage: ${0} TFs_file TGTs_file Prior_file Exp_file lambda Out_dir "
	exit
fi

tfs=$1
tgs=$2
net=$3
exp=$4
lam=$5
out=$6

echo "TFs: ${tfs}"
echo "TGTs: ${tgs}"
echo "Prior: ${net}"
echo "Data: ${exp}"
echo "Lambda: ${lam}"
echo "Outdir: ${out}"

if [ -d ${out} ]; 
then 
	echo -e "\e[1mDirectory exist, checking existing files...\e[0m"
	j=-1; 
	for ((i=0;i<10;i++)); 
	do 
		if [ -d ${out}/${i} ]; 
		then 
			j=$i; 
		fi; 
	done 
	echo -e "\e[1mIndex is: $j\e[0m"
	if [ $j -gt -1 ] && [ -f ${out}/${j}/0/fold0/prediction_k300.txt ]
	then
		echo -e "\e[1mWas stopped in middle of MERLIN, resume MERLIN\e[0m"
		time ./MERLINNCA -d ${exp} -r ${tfs} -g ${tgs} -p ${net} -l ${lam} -n10 -t 1 -q $j -a -o ${out} > ${out}/out.txt 
	else
		let j=j-1;
		if [ $j -gt -1 ]
		then 
			echo -e "\e[1mWas stopped in middle of NCA, resume NCA\e[0m"
			time ./MERLINNCA -d ${exp} -r ${tfs} -g ${tgs} -p ${net} -l ${lam} -n10 -t 1 -q $j -o ${out} > ${out}/out.txt 
		else
			echo -e "\e[1mStart over\e[0m"
			mkdir ${out}
			time ./MERLINNCA -d ${exp} -r ${tfs} -g ${tgs} -p ${net} -l ${lam} -n10 -t 1 -o ${out} > ${out}/out.txt
		fi
	fi
else
	echo -e "\e[1mStart from scratch\e[0m"
	mkdir ${out}
	time ./MERLINNCA -d ${exp} -r ${tfs} -g ${tgs} -p ${net} -l ${lam} -n10 -t 1 -o ${out} > ${out}/out.txt
fi

