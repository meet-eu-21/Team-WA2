echo -n "" > "normalisation_report.txt";
#echo "normalisation start";

while getopts 'c:i:' opt
do
    case $opt in
        c) chr_n=$OPTARG;;
        i) input_dir=$OPTARG;;
       \?) echo "ERROR: Invalid option: $USAGE"
           exit 1;;
    esac
done

if [ $chr_n == 'all' ]; then
    chromosomes=($(seq 1 22))
    chromosomes+=("X")
else
    chromosomes=($chr_n)
fi


outfiles=()
for i in "${chromosomes[@]}"
do
    #echo "Normalising chromosome" $i "matrix"
    #python3 ./main.py -chr $i >> "normalisation_report.txt";
    outfile=`python3 ./main.py -chr $i -dir $input_dir`
    outfiles+=($outfile)
    #echo "Chromosome" $i "done"
done
echo ${outfiles[*]}
#echo "normalisation done"
