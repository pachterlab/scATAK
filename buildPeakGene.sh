#!/bin/bash

# For peak to gene association what do we need?
# Input:
# 1. genome.gtf -> tss
# 3. atac.genes.txt (columns of matrix quantified against peaks.fa)
# 4. atac.barcodes.txt (rows of matrix quantified against peaks.fa)
# 5. atac.mtx (matrix generated from quantifying against peaks.fa)
# Output:
# 1. genes.mtx matrix summed on grouped genes 
# 2. genes.genes.txt (columns of matrix)
# 3. genes.barcodes.txt (rows of matrix)

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output           output FASTA
    -g, --gtf              Genome annotation (gtf.gz)
    -p, --atacpeaks        ATAC peaks (atac.genes.txt)
    -b, --atacbarcodes     ATAC barcodes (atac.barcodes.txt)
    -m, --atacmatrix       ATAC matrix (atac.mtx)
    "
    exit 1
}

while getopts ":o:g:p:b:m:" opt; do
    case $opt in
        o|--output)
            OUTPUT=$OPTARG
            ;;
        g|--GTF)
            GTF=$OPTARG
            ;;
        p|--atacpeaks)
            ATACPEAKS=$OPTARG
            ;;
        b|--atacbarcodes)
            ATACBCS=$OPTARG
            ;;
        m|--atacmatrix)
            ATACMTX=$OPTARG
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
echo "GTF: ${GTF}"
echo "ATACPEAKS: ${ATACPEAKS}"
echo "ATACBCS: ${ATACBCS}"
echo "ATACMTX: ${ATACMTX}"
echo "OUTPUT: ${OUTPUT}"

# check options        
if [ -z "$OUTPUT" -o -z "$GTF" -o -z "$ATACPEAKS" -o -z "$ATACBCS" -o -z "$ATACMTX" ]
then
    echo "Error"
    usage
fi

# GTF=Homo_sapiens.GRCh38.104.gtf.gz
# ATACPEAKS=atac.genes.txt
# ATACBCS=atac.barcodes.txt
# ATACMTX=atac.mtx
# OUTPUT=./out

mkdir -p tmp

# Get gene TSS from GTF 
awk '{if($3=="gene" && $7=="+") print $1"\t"$4"\t"$4"\t"$14; else if($3=="gene" && $7=="-") print $1"\t"$5"\t"$5"\t"$14;}' $GTF > tmp/tss.bed
sed -i 's/"//g' tmp/tss.bed
sed -i 's/;//g' tmp/tss.bed

# Sort TSS (used for peak->gene association)
bedtools sort -i tmp/tss.bed > tmp/tss_sort.bed
mv tmp/tss_sort.bed  tmp/tss.bed

# copy the genes as a bed file
cp $ATACPEAKS tmp/regions.bed

# clean up the bed file
sed -i 's/:/\t/g' tmp/regions.bed
sed -i 's/-/\t/g' tmp/regions.bed

# get the center of each peak
awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)"\t"$1":"$2"-"$3}' tmp/regions.bed > tmp/regions.center.bed

# find the closest tss to the center of each peak
bedtools closest -a tmp/regions.center.bed -b tmp/tss.bed -d > tmp/region_tss_association.txt
rm tmp/regions.center.bed

# Give a  weight to the distance of each TSS <-> center association 
awk '{if($9<=2000) print $4"\t"$8"\t"1; if($9>2000 && $9<=5000) print $4"\t"$8"\t"0.7; if($9>5000 && $9<=10000) print $4"\t"$8"\t"0.5; if($9>10000 && $9<=20000) print $4"\t"$8"\t"0.25; if($9>20000 && $9<=50000) print $4"\t"$8"\t"0.03;}' tmp/region_tss_association.txt > tmp/region_tss_score.txt
awk '{print $4}' tmp/tss.bed | sort | uniq > $OUTPUT/gene.genes.txt

# What does this do??
awk 'NR==FNR{j++;Arr[$1]=NR;next} ($1 in Arr){print Arr[$1]"\t"$2"\t"$3}' $ATACPEAKS tmp/region_tss_score.txt > tmp/region_index_tss_score.txt
awk 'NR==FNR{j++;Arr[$1]=NR;next} ($2 in Arr){print $1"\t"Arr[$2]"\t"$3}' $OUTPUT/gene.genes.txt tmp/region_index_tss_score.txt > tmp/region_index_tss_index_score.txt
rm tmp/region_index_tss_score.txt

# Perform weighted sum based on peak <-> gene association 
BC=1
last_line=$(wc -l < $ATACMTX)
awk -v bc=$BC -v last=$last_line 'NR==FNR{gene[$1]=$2;score[$1]=$3;next} ($2 in gene){if(bc!=$1 || NR==last){for (key in sum) {print bc" "key" "sum[key];} delete sum; bc=$1;} sum[gene[$2]]+=score[$2];}' tmp/region_index_tss_index_score.txt $ATACMTX > $OUTPUT/gene.mtx
echo "$(wc -l < $ATACBCS) $(wc -l < $OUTPUT/gene.genes.txt) $(wc -l < $OUTPUT/gene.mtx)" > tmp/gene_mtx_sum
head -3 $ATACMTX | cat - tmp/gene_mtx_sum  $OUTPUT/gene.mtx > $OUTPUT/gene2.mtx
mv $OUTPUT/gene2.mtx $OUTPUT/gene.mtx