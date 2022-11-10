#!/bin/bash

# curl -Ls https://github.com/pachterlab/scATAK/blob/main/lib/737K-cratac-v1.txt?raw=true > 737K-cratac-v1.txt
# 10xATAC

# quantify.sh -o atac/ -i peaks.idx -x 10xATAC -g t2g.txt -w 737K-cratac-v1.txt -c R2.fastq.gz -1 R1.fastq.gz -2 R3.fastq.gz

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output
    -p, --peaks 
    -x, --technology
    -w, --whitelist
    -c, --cbfastq
    -1, --r1fastq
    -2, --r2fastq
    "
    exit 1
}

while getopts ":o:p:x:w:c:1:2:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        p|--peaks)
            PEAKS=$OPTARG
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
if [ -z "$OUTPUT" -o -z "$PEAKS" -o -z "$TECH" -o -z "$WL" -o -z "$CBFQ" -o -z "$R1FQ" -o -z "$R2FQ" ]
then
    echo "Error"
    usage
fi

# Merge peaks from multiple samples
cat $PEAKS | bedtools sort | bedtools merge > $OUTPUT/peaks.merged.bed


# Create merged peaks FASTA
bedtools getfasta -fi tmp/genome.fa -bed $OUTPUT/peaks.merged.bed -fo $OUTPUT/peaks.merged.fa
cat $OUTPUT/peaks.merged.fa | awk '{if($1~/>/)print $1"\t"$1"\t"$1}' > $OUTPUT/t2g.txt
sed -i 's/>//g' $OUTPUT/t2g.txt

# Build pseudoalignment index
kallisto index -i $OUTPUT/peaks.idx $OUTPUT/peaks.merged.fa

# Quantify
kb count \
-i $OUTPUT/peaks.idx \
-g $OUTPUT/t2g.txt \
-x $TECH \
-o $OUTPUT \
-w  $WL \
--h5ad \
$CBFQ $R1FQ $R2FQ
