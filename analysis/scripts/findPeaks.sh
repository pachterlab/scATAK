#!/bin/bash
# Input:
# 1. Genome FASTA
# 2. Cell BC FASTQ
# 3. Read 1 FASTQ
# 4. Read 2 FASTQ
# Output:
# 1. peaks.bed


# 2. peaks.fa
# 3. t2g.txt

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output           output folder
    -g, --genome           Genome primary assembly (fa.gz)
    -c, --cellbc           Fastq for the cell barcode
    -1, --read1            Fastq for biological paired end 1
    -2, --read2            Fastq for biological paired end 2
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
mkdir -p $OUTPUT

# Build minimap2 index
minimap2 -d tmp/minimap2_index.mmi $GENOME

# Align reads to genome
minimap2 -ax sr -t 8 tmp/minimap2_index.mmi \
   <(paste <(zcat $R1FQ) <(zcat $CBCFQ) | awk -v blen="$BLEN" '{if(NR%4==1) header=$1; if(NR%4==2) print header"_"substr($2, 0, blen)"_\n"$1; if(NR%4==3 || NR%4==0) print $1;}') \
   <(paste <(zcat $R2FQ) <(zcat $CBCFQ) | awk -v blen="$BLEN" '{if(NR%4==1) header=$1; if(NR%4==2) print header"_"substr($2, 0, blen)"_\n"$1; if(NR%4==3 || NR%4==0) print $1;}') \
   > tmp/genome.sam

# Convert sam to bam and sort
samtools view -@ 8 -o tmp/genome.bam -b tmp/genome.sam
# Delete SAM
rm tmp/genome.sam
# Sort BAM
samtools sort -n -@ 8 -o tmp/genome.sorted.bam -m 8G tmp/genome.bam 

# Call peaks from BAM
Genrich -t tmp/genome.sorted.bam -o tmp/genome.narrowPeak -f tmp/genome_peaks.log -v

# Sort the peaks (prepare for slicing genome)
cat tmp/genome.narrowPeak | bedtools sort | bedtools merge > $OUTPUT/peaks.bed