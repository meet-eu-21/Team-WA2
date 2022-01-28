#!/bin/bash

#chr_n=''
normalisation_type="normal"
peak_find_func="find_peaks"
smooth_func="False"
filter_func="filter_peaks"

while getopts 'i:c:n:p:s:f:' opt
do
    case $opt in
        i) input_dir=$OPTARG;;
        c) chr_n=$OPTARG;;
        n) normalisation_type="kr";;
        p) peak_find_func=$OPTARG;;
        s) smooth_func=$OPTARG;;
        f) filter_func=$OPTARG;;
       \?) echo "ERROR: Invalid option: $USAGE"
           exit 1;;
    esac
done

if [ ${chr_n} == "all" ]
then
    chromosomes=($(seq 1 22))
    chromosomes+=("X")
else
    chromosomes=($chr_n)
fi



echo "

-----Summary of input parameters----------
Folder with input files: ${input_dir}
Chromosomes: ${chromosomes[*]}
Normalisation: ${normalisation_type}
Peak finding function: ${peak_find_func}
Smoothing function: ${smooth_func}
Filtering function: ${filter_func}
------------------------------------------

"

### normalisation ###
p=`pwd`
normalised_matrices=()
if [ $normalisation_type == "normal" ]; then
    cd ./normalisation/NORMALnormalisation/
    nm=`./normalise_all_chromosome_matrices.sh -c $chromosomes -i $p/$input_dir`
    normalised_matrices+=(`pwd`/$nm)
else
    cd ./normalisation/KRnormalisation/
    ./KRnormalise_all_chromosome_matrices.sh
fi
cd ../..

### running TopDom ###
echo ${normalised_matrices[*]}
cd scripts/

for i in "${normalised_matrices[@]}"
do
    python3 TopDom.py -i $i -p $peak_find_func -s $smooth_func -f $filter_func
done

pwd=`pwd`
cd ..

echo "Finished! The output files with defined TADs can be found in $pwd"
# ./test.sh -i HiC/GM12878/100kb_resolution_intrachromosomal/ -c 21
