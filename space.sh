TYPEGRAPH=$1
PARAMS=$2
INPUT=$3
OUTPUT=$4

echo Params $PARAMS
echo Typegraph $TYPEGRAPH
echo Reading $INPUT

HEADER=`gzcat $INPUT|head -n1`
LEVELS=`python -c "import math; h = '$HEADER'.split(); print int(math.ceil(math.log(max(int(h[0]),int(h[2]))/math.log(2) )));"`

#for ((lf=4; lf <= $LEVELS; lf+=4 ))
#do
#    OUTFILE=$OUTPUT-$PARAMS.mxf
#	if [ ! -f "$OUTFILE" ]; then
#    echo Creating $OUTFILE
#    gzcat $INPUT | ./create -s MXF -g $TYPEGRAPH -f $PARAMS $OUTFILE
#	fi
#done

for ((l=0; l <= $LEVELS; l+=4 ))
do
    lki=$((l*2))
    
    OUTFILE=$OUTPUT-$PARAMS.mxd
	if [ ! -f "$OUTFILE" ]; then
		echo Creating $OUTFILE
		gzcat $INPUT | ./create -s MXD -g $TYPEGRAPH -f $PARAMS $OUTFILE
	fi
	
    OUTFILE=$OUTPUT-$PARAMS.prb
    if [ ! -f "$OUTFILE" ]; then
		echo Creating $OUTFILE
		gzcat $INPUT | ./create -s PRB -g $TYPEGRAPH -f $PARAMS $OUTFILE
	fi
	
    for((F=2; F <= 16; F+=2 ))
    do
        OUTFILE=$OUTPUT-$PARAMS.prb2
		if [ ! -f "$OUTFILE" ]; then
        	echo Creating $OUTFILE
			gzcat $INPUT | ./create -s PRB2 -g $TYPEGRAPH -f $PARAMS $OUTFILE
		fi
    done
	
#    OUTFILE=$OUTPUT-2,2,0,$lki,0,0-$BITT-$BITB.prw
#	if [ ! -f "$OUTFILE" ]; then
#    echo Creating $OUTFILE
#    gzcat $INPUT | ./create -s PRW -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0,0 $OUTFILE
#	fi
done
