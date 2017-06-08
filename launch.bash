#declare -a name=(fb ok kg1 rd hw tw pk wk)
declare -a name=(rm)

for ((gcount=1;gcount<=32;gcount=gcount*2))
do
	for file in ${name[@]}
	do
		file_inst=jobScript"$gcount"gpu-$file
		if [ -f $file_inst ]
		then
			sbatch $file_inst 
		fi
	done
done
