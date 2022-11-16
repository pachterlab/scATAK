#!/bin/bash

# curl -Ls https://github.com/pachterlab/scATAK/blob/main/lib/737K-cratac-v1.txt?raw=true > 737K-cratac-v1.txt
# 10xATAC

# quantify.sh -o atac/ -i peaks.idx -x 10xATAC -g t2g.txt -w 737K-cratac-v1.txt -c R2.fastq.gz -1 R1.fastq.gz -2 R3.fastq.gz

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output
    -p, --peaks
    -g, --genome 
    "
    exit 1
}

while getopts ":o:p:g:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        p|--peaks)
            PEAKS=$OPTARG
            ;;
        g|--genome)
            GENOME=$OPTARG
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
if [ -z "$OUTPUT" -o -z "$PEAKS" -o -z "$GENOME" ]
then
    echo "Error"
    usage
fi

mkdir -p $OUTPUT
mkdir -p tmp

zcat $GENOME | fold -w 80 > tmp/genome.fa

# Merge peaks from multiple samples
cat $PEAKS | bedtools sort | bedtools merge > $OUTPUT/peaks.merged.bed


# Create merged peaks FASTA
bedtools getfasta -fi tmp/genome.fa -bed $OUTPUT/peaks.merged.bed -fo $OUTPUT/peaks.merged.fa
cat $OUTPUT/peaks.merged.fa | awk '{if($1~/>/)print $1"\t"$1"\t"$1}' > $OUTPUT/t2g.txt
sed -i 's/>//g' $OUTPUT/t2g.txt