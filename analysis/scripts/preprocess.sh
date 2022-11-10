#!/bin/bash

# curl -Ls https://github.com/pachterlab/scATAK/blob/main/lib/737K-cratac-v1.txt?raw=true > 737K-cratac-v1.txt
# 10xATAC

# quantify.sh -o atac/ -i peaks.idx -x 10xATAC -g t2g.txt -w 737K-cratac-v1.txt -c R2.fastq.gz -1 R1.fastq.gz -2 R3.fastq.gz

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output
    -f, --fasta 
    -g, --t2g 
    -x, --technology
    -w, --whitelist
    -c, --cbfastq
    -1, --r1fastq
    -2, --r2fastq
    "
    exit 1
}

while getopts ":o:f:g:x:w:c:1:2:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        f|--fasta)
            FASTA=$OPTARG
            ;;
        g|--genemap)
            T2G=$OPTARG
            ;;
        x|--technology)
            TECH=$OPTARG
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
if [ -z "$OUTPUT" -o -z "$T2G" -o -z "$FASTA" -o -z "$TECH" -o -z "$WL" -o -z "$CBFQ" -o -z "$R1FQ" -o -z "$R2FQ" ]
then
    echo "Error"
    usage
fi

# Build pseudoalignment index
kallisto index -i tmp/peaks.idx $FASTA

# Quantify
kb count \
-i tmp/peaks.idx \
-g $T2G \
-x $TECH \
-o $OUTPUT \
-w  $WL \
--h5ad \
$CBFQ $R1FQ $R2FQ

kb count -i tmp/peaks.idx -g $T2G -x $TECH -o $OUTPUT -w  $WL --h5ad $CBFQ $R1FQ $R2FQ