#!/bin/bash
#base_path="/home/mzj/mzj/transfer_data"
base_path="/home/mzj/mzjdata"

graphdir=(google berkstan livejournal pokec ca_rn tx_rn)
#graphdir=(google)
#graphname=(web-Google.txt web-BerkStan.txt cage15.graph.txt soc-LiveJournal1.txt usa_rn.txt)

graphname=(web-Google.txt web-BerkStan.txt soc-LiveJournal1.txt pokec.txt ca_rn.txt tx_rn.txt)

#graphname=(web-Google.txt)
mem=(50 50 512 256 50 30)
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
