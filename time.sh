SPACEDIR=~/is2014/space
TIMEDIR=~/is2014/time
QUERYDIR=~/is2014/queries

DATASET=$1
QUERY=$2
TRIALS=$3

TIME=$(date "+%Y-%m-%d_%H:%M:%S")
echo $DATASET $QUERY $TRIALS
for i in $(ls $SPACEDIR/$DATASET/)
do
	for j in $(seq 1 $TRIALS)
	do
		echo "[$TIME] Running $i $QUERY (trial $j)"
		./benchmark $SPACEDIR/$DATASET/$i $QUERYDIR/$DATASET/$QUERY >> $TIMEDIR/$DATASET-$TIME

	done
done
