#!/bin/bash
#echo  -e edgelist | bin/example_apps/connectedcomponents --file=/home/mzj/mzjdata/web-Google/web-Google.txt --membudget_mb=50 
if [ "$1" != "" ];then
cat $1 | awk '{for(i=1;i<=NF;i++) if ($i=="niters:" || $i=="membudget_mb:"|| $i=="nshards:" || $i=="app:" || $i=="execute-updates:" || $i=="runtime:" || $i=="file:") print $(i),"\t",$(i+1)}' > $2
else
	echo "input file parameter is empty"
fi
