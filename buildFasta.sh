#!/bin/bash

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output            output FASTA
    -g, --genome            Genome primary assembly (fa.gz)
    -cbc, --cellbc          Fastq for the cell barcode
    -r1, --read1            Fastq for biological paired end 1
    -r2, --read2            Fastq for biological paired end 2
    "
    exit 1
}

while getopts ":o:g:c:1:2:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        g|--genome)
            GENOME=$OPTARG
            ;;
        c|--cellbc)
            CBCFQ=$OPTARG
            ;;
        1|--read1)
            R1FQ=$OPTARG
            ;;
        2|--read2)
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

echo "Options"
echo "-------"
echo "Output: ${OUTPUT}"
echo "Genome: ${GENOME}"
echo "CellFq: ${CBCFQ}"
echo "Read1Fq: ${R1FQ}"
echo "Read2Fq: ${R2FQ}"

echo ""

# check options        
if [ -z "$OUTPUT" -o -z "$GENOME" -o -z "$CBCFQ" -o -z "$R1FQ" -o -z "$R2FQ" ]
then
    echo "Error"
    usage
fi

BLEN=16 # Could be derived from fastq?
mkdir -p tmp

# Build minimap2 index
minimap2 -d tmp/minimap2_index.mmi $GENOME

# Align reads to genome
minimap2 -ax sr -t 8 tmp/minimap2_index.mmi \
   <(paste <(zcat $R1FQ) <(zcat $CBCFQ) | awk -v blen="$BLEN" '{if(NR%4==1) header=$1; if(NR%4==2) print header"_"substr($2, 0, blen)"_\n"$1; if(NR%4==3 || NR%4==0) print $1;}') \
   <(paste <(zcat $R2FQ) <(zcat $CBCFQ) | awk -v blen="$BLEN" '{if(NR%4==1) header=$1; if(NR%4==2) print header"_"substr($2, 0, blen)"_\n"$1; if(NR%4==3 || NR%4==0) print $1;}') \
   > tmp/genome.sam

# Convert SAM to BAM
sambamba view -t 8 -f bam -S -o tmp/genome.bam tmp/genome.sam

# Delete SAM
rm tmp/genome.sam

# Sort BAM
sambamba sort -n -t 8 -m 8GB --tmpdir=./tmp tmp/genome.bam

# Call peaks from BAM
Genrich -t tmp/genome.sorted.bam -o tmp/genome.narrowPeak -f tmp/genome_peaks.log -v

# Sort the peaks (prepare for slicing genome)
cat tmp/genome.narrowPeak | bedtools sort | bedtools merge > tmp/peaks.bed

# Create peak FASTA
bedtools getfasta -fi $GENOME -bed tmp/peaks.bed -fo $OUTPUT/peaks.fa
cat tmp/peaks.fa | awk '{if($1~/>/)print $1"\t"$1"\t"$1}' > tmp/t2g.txt
sed -i 's/>//g' tmp/t2g.txt

# # Get gene TSS from GTF 
# awk '{if($3=="gene" && $7=="+") print $1"\t"$4"\t"$4"\t"$14; else if($3=="gene" && $7=="-") print $1"\t"$5"\t"$5"\t"$14;}' $GTF > tss.bed
# sed -i 's/"//g' tss.bed
# sed -i 's/;//g' tss.bed

# # Sort TSS (used for peak->gene association)
# bedtools sort -i tss.bed > tss_sort.bed
# mv tss_sort.bed  tss.bed