TYPEGRAPH=$1
BITT=$2
BITB=$3
INPUT=$4
OUTPUT=$5

echo BitmapT $BITT
echo BitmapB $BITB
echo Typegraph $TYPEGRAPH
echo Reading $INPUT

OUTFILE=$OUTPUT-2,2,0,0,0-$BITT-$BITB.mxd
echo Creating $OUTFILE
gzcat $INPUT | ./create -s MXD -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,0,0 $OUTFILE

for lf in {1..32}
do
    OUTFILE=$OUTPUT-2,2,0,0,$lf-$BITT-$BITB.mxf
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s MXF -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,0,$lf $OUTFILE
done

for l in {1..32}
do
    lki=$((l*2))
    OUTFILE=$OUTPUT-2,2,0,$lki,0-$BITT-$BITB.prb
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s PRB -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0 $OUTFILE
    
    OUTFILE=$OUTPUT-2,2,0,$lki,0-$BITT-$BITB.prw
    echo Creating $OUTFILE
    gzcat $INPUT | ./create -s PRW -t $BITT -b $BITB -g $TYPEGRAPH -f 2,2,0,$lki,0 $OUTFILE
done
