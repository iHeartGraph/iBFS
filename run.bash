#declare -a file=(desc-fb-24  desc-fr-24  desc-hw-24  desc-kg0-24  desc-kg1-24  desc-kg2-24 desc-lj-24  desc-ok-24  desc-pk-24  desc-rd-24 desc-rm-24  desc-tw-24 desc-wk-24)
#declare -a file=(desc-hw-24 desc-lj-24 desc-tw-24 desc-wk-24)
#declare -a file=(desc-hw-24 desc-tw-24)
declare -a file=(desc-hw-24)

for name in ${file[@]}
do
	echo $name
	for ((i=16;i<=1048576;i=i*2)) 
	do
		stdbuf -o 0 ./ibfs-5.0 ../desc-files/$name haha.input $i | tee $name.$i.nogroupby.log
	done
done
