#!/bin/bash

declare -a file=(desc-fb-24  desc-fr-24  desc-hw-24  desc-kg1-24  desc-kg2-24 desc-lj-24  desc-ok-24  desc-pk-24  desc-rd-24 desc-rm-24  desc-tw-24 desc-wk-24)
#declare -a kernel=(THD_expand_sort WAP_expand_sort CTA_expand_sort THD_bu_expand_sort WAP_bu_expand_sort CTA_bu_expand_sort)
declare -a kernel=(sw_expand_warp sw_expand_thd)
declare -a metric=(gld_transactions gst_transactions)

kens=""
count=0
for ken in ${kernel[@]}
do
	if [ $count -eq 0 ]
	then
		kens=$ken
	else
		kens=$kens"|"$ken
	fi
	count=$((count+1))
done
echo $kens

for met in ${metric[@]}
do
	echo $met
	for fle in ${file[@]}
	do
		grep -E "$kens|$met" $fle.3.6.gld.gst.count | grep -A 1 "Kernel" | sed -e '/^--$/d'|awk -v av=0 -v count=0 '{if(NR%2==0){av+=$NF; count++;}}END{printf "%.2f\n",av/count}' 
	done
done
