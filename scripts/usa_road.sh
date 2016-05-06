#!/bin/bash
base_path="/home/mzj/mzj/transfer_data"
graphdir=(usa)
graphname=(usa_rn.txt)
mem=(8192)

export GRAPHCHI_ROOT=/home/mzj/block-graphchi
#for i in "${graphname[@]}"
if [ "${#graphdir[@]}" = "${#graphname[@]}" -a "${#graphname[@]}" = "${#mem[@]}" ];then
	for (( i=0; i<${#graphdir[@]}; i++ ))
	do
		echo "running ${graphname[i]}"
		if [ "${graphname[i]}" != "" ];then
			echo edgelist | /home/mzj/block-graphchi/bin/example_apps/partition_hash --file=$base_path/${graphdir[i]}/${graphname[i]}\
				 --membudget_mb=${mem[i]} > $base_path/${graphdir[i]}/${graphname[i]}-hash_output 2>&1 

			echo edgelist | /home/mzj/block-graphchi/bin/example_apps/partition_bfs --file=$base_path/${graphdir[i]}/${graphname[i]}\
				 --membudget_mb=${mem[i]} > $base_path/${graphdir[i]}/${graphname[i]}-bfs_output 2>&1 

			sleep 1m
			#sudo ../clear_cache.sh
		fi
		#echo ${graphdir[i]}	${graphname[i]}
	done
else
	echo "array length is not equal"
	exit 1
fi
