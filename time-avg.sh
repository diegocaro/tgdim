SPACEDIR=~/is2014/space
TIMEDIR=~/is2014/time-avg
QUERYDIR=~/is2014/queries

DATASET=$1
STRUCT=$2
QUERY=$3
TRIALS=$4

if [ $# -eq 0 ]
  then
    echo "bash time.sh <dataset> <struct> <query> <trials>"
    exit
fi

TIME=$(date "+%Y-%m-%d_%H:%M:%S")
echo $DATASET $STRUCT $QUERY $TRIALS
for i in $(ls $SPACEDIR/$DATASET*$STRUCT)
do
	for j in $(seq 1 $TRIALS)
	do
		echo "[$TIME] Running $i $QUERY (trial $j)"
		./benchmark $i $QUERYDIR/$DATASET/$QUERY >> $TIMEDIR/$DATASET-$TIME

	done
done
