#!/bin/bash
#@ Yulong Li
# run admixture using different seeds and K.

infile=$1
maxk=$2
reps=$3
T=$4
if [ $# -lt 4 ] || [ $1 = -h ];
then
    echo "run_admix.sh [input ped file] [maximum K] [replicates] [threads]"
    exit 1
fi
i=0
# run admixture
for k in `seq 1 $maxk`;
    do for j in `seq 1 $reps`;   
        do let i++
	pre=`basename -s .bed $infile`
        cur="$pre.run_$i.$k"
        admixture --cv=20 -j$T -s $RANDOM $infile $k | tee $cur.log
	mv $pre.$k.Q $cur.Q
    done
done
echo "Done!"
    
