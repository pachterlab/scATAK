args=("$@")
input=${args[0]} 
bin=${args[1]}
#echo "Input file is $input"
awk -v BIN=$bin '{print $1":"int($2/BIN)*BIN"-"int($3/BIN+1)*BIN"\t"$4}' ${input} > track.bin
awk 'NR==FNR{Arr[$1]++;next;} ($1 in Arr){sum[$1]=sum[$1]+$2;} END {for (key in Arr) print key"\t"sum[key]}' track.bin track.bin | sed 's/:/\t/g' - | sed 's/-/\t/g' -

