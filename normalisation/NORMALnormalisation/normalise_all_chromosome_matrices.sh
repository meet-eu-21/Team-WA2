echo -n "" > "normalisation_report.txt";
echo "start";


chromosomes=($(seq 1 22))
chromosomes+=("X")
for i in "${chromosomes[@]}"
do
    echo "Normalising chromosome" $i "matrix"
    python3 main.py -chr $i >> "normalisation_report.txt";
    echo "Chromosome" $i "done"
done

