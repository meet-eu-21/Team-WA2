peak_find_funk=("find_min" "detect_local_extrema" "find_peaks" "find_peaks_2")
smooth_func=("False" "savgol_filter" "smooth" "qspline")
filter_func=("False" "statFilter" "filter_peaks")
chromosomes=($(seq 1 22))
chromosomes+=("X")

echo "start"

for p in "${peak_find_funk[@]}";
do
    peakfindfunc="-p $p"

    for s in "${smooth_func[@]}";
    do
        smoothfunc="-s $s"

        for f in "${filter_func[@]}";
        do
            filterfunc="-f $f"
            outfile="${p}_${s}_${f}"
            for i in "${chromosomes[@]}";
            do
                basic_command="python3 TopDom.py -i /home/eu_plus_3/Team-WA2/normalisation/NORMALnormalisation/normalised_matrices/contact_map_chr$i.tsv -r 100k -d ./combination_outputs/"
                echo "Topdom for chromosome $i contact map."
                echo "Chromosome $i $outfile"
                eval "$basic_command $peakfindfunc $smoothfunc $filterfunc -o $outfile"
                echo "Chromosome $i done."
            done
            eval "cat ./combination_outputs/for_overlap_scores__${outfile}__chr*.csv >> ./combination_outputs/for_overlap_scores__${outfile}__all.csv"
            eval "cat -n ./combination_outputs/for_overlap_scores__${outfile}__all.csv | sort -uk2 | sort -nk1 | cut -f2- > ./combination_outputs/for_overlap_scores__${outfile}__allchr.csv"
        done
    done
done

echo "Done!"
