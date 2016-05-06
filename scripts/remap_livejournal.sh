#!/bin/bash
#base_path="/home/mzj/mzj/transfer_data"
base_path="/home/mzj/mzjdata"

graphdir=(livejournal livejournal)
#graphdir=(google)
#graphname=(web-Google_weight.txt web-BerkStan_weight.txt cage15.graph_weight.txt soc-LiveJournal1_weight.txt usa_rn_weight.txt)

graphname=(soc-LiveJournal1.txt soc-LiveJournal1_weight.txt)
#graphname=(web-Google_weight.txt web-BerkStan_weight.txt soc-LiveJournal1_weight.txt pokec_weight.txt ca_rn_weight.txt tx_rn_weight.txt)

mem=(512 512)
#mem=(50 50 1024 256 50 30)
#mem=(50)

export GRAPHCHI_ROOT=/home/mzj/block-graphchi
#for i in "${graphname[@]}"
if [ "${#graphdir[@]}" = "${#graphname[@]}" -a "${#graphname[@]}" = "${#mem[@]}" ];then
	for (( i=0; i<${#graphdir[@]}; i++ ))
	do
		(cd $base_path/${graphdir[i]} && exec ./clean.sh >/dev/null 2>&1)
		echo "running $base_path/${graphdir[i]}/clean.sh"

		if [ "${graphname[i]}" != "" ];then
			echo "hash par on ${graphname[i]}"
			echo edgelist | /home/mzj/block-graphchi/bin/example_apps/partition_hash --file=$base_path/${graphdir[i]}/${graphname[i]}\
				 --membudget_mb=${mem[i]} > $base_path/${graphdir[i]}/${graphname[i]}-hash_output 2>&1 
				# --membudget_mb=${mem[i]} > $base_path/${graphdir[i]}/${graphname[i]}-hash_output 2>&1 
			#$base_path/${graphdir[i]}/clean.sh;	
				
			#echo "clean ${graphdir[i]}"
			(cd $base_path/${graphdir[i]} && exec ./clean.sh>/dev/null 2>&1)
			echo "cleaned hash ${graphdir[i]}/clean.sh"
			
			echo "bfs par on ${graphname[i]}"
			echo edgelist | /home/mzj/block-graphchi/bin/example_apps/partition_bfs --file=$base_path/${graphdir[i]}/${graphname[i]}\
				 --membudget_mb=4096 > $base_path/${graphdir[i]}/${graphname[i]}-bfs_output 2>&1 
				 #--membudget_mb=${mem[i]} > $base_path/${graphdir[i]}/${graphname[i]}-bfs_output 2>&1 
			
			#$base_path/${graphdir[i]}/clean.sh	
			
			(cd $base_path/${graphdir[i]} && exec ./clean.sh > /dev/null 2>&1)
			echo "cleaned bfs ${graphdir[i]}"
			#sudo ../clear_cache.sh
		fi
		#echo ${graphdir[i]}	${graphname[i]}
	done
else
	echo "array length is not equal"
	exit 1
fi
