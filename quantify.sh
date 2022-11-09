#!/bin/bash

# curl -Ls https://github.com/pachterlab/scATAK/blob/main/lib/737K-cratac-v1.txt?raw=true > 737K-cratac-v1.txt
# 10xATAC

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output
    -i, --index
    -x, --technology
    -g, --t2g
    -w, --whitelist
    -c, --cbfastq
    -1, --r1fastq
    -2, --r2fastq
    "
    exit 1
}

while getopts ":o:i:x:g:w:c:1:2:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        i|--INDEX)
            INDEX=$OPTARG
            ;;
        x|--technology)
            TECH=$OPTARG
            ;;
        g|--t2g)
            T2G=$OPTARG
            ;;
        w|--whitelist)
            WL=$OPTARG
            ;;
        c|--cbfastq)
            CBFQ=$OPTARG
            ;;
        1|--r1fastq)
            R1FQ=$OPTARG
            ;;
        2|--r2fastq)
            R2FQ=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid argument"
            usage
            ;;
        :)
            echo "Add arguments"
            usage
            ;;
    esac
done

# check options        
if [ -z "$OUTPUT" -o -z "$INDEX" -o -z "$TECH" -o -z "$T2G" -o -z "$WL" -o -z "$CBFQ" -o -z "$R1FQ" -o -z "$R2FQ" ]
then
    echo "Error"
    usage
fi


kb count \
-i $INDEX \
-g $T2G \
-x $TECH \
-o $OUTPUT \
-w  $WL \
--h5ad \
$CBFQ $R1FQ $R2FQ