#declare -a file=(desc-fb-24  desc-fr-24  desc-hw-24  desc-kg0-24  desc-kg1-24  desc-kg2-24 desc-lj-24  desc-ok-24  desc-pk-24  desc-rd-24 desc-rm-24  desc-tw-24 desc-wk-24)
#declare -a file=(desc-hw-24 desc-lj-24 desc-tw-24 desc-wk-24 desc-pk-24)
declare -a file=(desc-hw-24)

#echo "#Share diff hub"
#for name in ${file[@]}
#do
#	#echo $name
#	i=1
#	for share in 4 16 64 128
#	do
#		fqsz=$(grep "Break-fq-sz" $name-diff-outdegree-"$share".log | awk -F " " -v fq=10000000 '{if(fq>$2)fq=$2}END{print fq}')
#		array[$i]=$(grep -E "Break-fq-sz: $fqsz|Traversal-" $name-diff-outdegree-"$share".log | grep -A 1 Break | sed -e "/--/d;/Break/d" | awk -F " " -v count=0 -v time=0 '{time+=$2;count++}END{print time/count}')
#		i=$((i+1))	
#	done
#	printf "%s," "$name"
#	printf "%s," "${array[@]}"
#	echo
#done



echo "#Share same hub"
for name in ${file[@]}
do
	#echo $name
	i=1
	for ((num_groups=1;num_groups<=1048576;num_groups=num_groups*2)) 
	#for share in 16
	#for num_groups in 1 2 3 4 5 6 7 8 9 10 50 100 200 300 400 500 600 700 800 900 1000
	do
		fqsz=$(grep "Break-fq-sz" "$name"."$num_groups".nogroupby.log | awk -F " " -v fq=10000000 '{if(fq>$2)fq=$2}END{print fq}')
		array[$i]=$(grep -E "Break-fq-sz: $fqsz|Traversal-" "$name"."$num_groups".nogroupby.log | grep -A 1 Break | sed -e "/--/d;/Break/d" | awk -F " " -v count=0 -v time=0 '{time+=$2;count++}END{if(count!=0){print time/count}}')
		i=$((i+1))	
	done
	printf "%s\n" "$name"
	printf "%s\n" "${array[@]}"
	echo
done
