peak_find_funk=("None" "find_min" "detect_local_extrema" "find_peaks" "find_peaks_2")
smooth_func=("False" "savgol_filter" "smooth" "qspline")
filter_func=("False" "statFilter" "filter_peaks")
chromosomes=($(seq 1 2))
chromosomes+=("X")

echo "start"

for i in "${chromosomes[@]}";
do

    echo "Chromosome $i"
    basic_command="python3 TopDom.py -i /home/eu_plus_3/Team-WA2/normalisation/NORMALnormalisation/normalised_matrices/contact_map_chr$i.tsv -r 100k -d ./combination_outputs/"
    for p in "${peak_find_funk[@]}";
    do
        peakfindfunc="-p $p"

        for s in "${smooth_func[@]}";
        do
            smoothfunc="-s $s"

            for f in "${filter_func[@]}";
            do
                filterfunc="-f $f"
                outfile="-o ${p}_${s}_${f}"
                echo "Topdom for chromosome $i contact map."
                echo "$basic_command $peakfindfunc $smoothfunc $filterfunc $outfile"
                eval "$basic_command $peakfindfunc $smoothfunc $filterfunc $outfile"
                echo "Chromosome $i done."
            done
        done
    done

done

echo "Done!"
