declare -a data=(rm)
for arr in ${data[@]}
do
	echo "Processing $arr ... "
	for ((i=1;i<=32;i=i*2))
	do
		file=$(find ./ -maxdepth 1 -iname "traversal-"$i"gpu-$arr.o*")

		if [ -f $file ]
		then
			tm=$(grep "Traversal-iter" $file | awk -v tm=0 -v count=0 '{tm+=$2;count++}END{print tm}')	
			echo $tm
		fi
	done
	echo "----------------"
done
