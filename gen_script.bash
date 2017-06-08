#!/bin/bash

declare -a name=(fb hw kg1 ok pk rd rm tw wk)
for file in ${name[@]}
do
	for ((i=1;i<=32;i=i*2))
	do
		#cp jobScript"$i"gpu-fb jobScript"$i"gpu-$file
		sed -i "s/TG-CCR120025/TG-CTS100002/g" jobScript"$i"gpu-$file
		#sed -i "s/"$i"gpu/-"$i"gpu-$file/g" jobScript"$i"gpu
	done
done

