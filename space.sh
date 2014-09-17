TYPEGRAPH=$1
BITT=$2
BITB=$3
INPUT=$4
OUTPUT=$5

echo BitmapT $BITT
echo BitmapB $BITB
echo Typegraph $TYPEGRAPH
echo Reading $INPUT

HEADER=`gzcat $INPUT|head -n1`
LEVELS=`python -c "import math; h = '$HEADER'.split(); print int(math.ceil(math.log(max(int(h[0]),int(h[2]))/math.log(2) )));"`

for ((lf=4; lf <= $LEVELS; lf+=4 ))
do
    OUTFILE=$OUTPUT-2,2,0,0,$lf-$BITT-$BITB.mxf
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s MXF -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,0,$lf $OUTFILE
done

for ((l=0; l <= $LEVELS; l+=4 ))
do
    lki=$((l*2))
    
    OUTFILE=$OUTPUT-2,2,0,$lki,0-$BITT-$BITB.mxd
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s MXD -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0 $OUTFILE
    
    OUTFILE=$OUTPUT-2,2,0,$lki,0-$BITT-$BITB.prb
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s PRB -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0 $OUTFILE
    
    OUTFILE=$OUTPUT-2,2,0,$lki,0-$BITT-$BITB.prw
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s PRW -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0 $OUTFILE
done
