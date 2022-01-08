echo "start";


chromosomes=($(seq 1 22))
chromosomes+=("X")
for i in "${chromosomes[@]}"
do
    echo "Topdom for chromosome" $i "contact map."
    python3 TopDom.py -i /home/eu_plus_3/Team-WA2/normalisation/NORMALnormalisation/normalised_matrices/contact_map_chr$i.tsv -r 100k -d ./topdom_outputs/
    echo "Chromosome" $i "done"
done

cat ./topdom_outputs/for_overlap_scores_chr*.csv >> ./topdom_outputs/for_overlap_scores_1.csv
