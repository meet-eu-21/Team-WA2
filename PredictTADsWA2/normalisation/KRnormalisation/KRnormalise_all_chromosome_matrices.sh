echo -n "" > "KRnormalisation_report.txt";
echo "start";


chromosomes=($(seq 3 22))
#chromosomes+=("X")
res=100
for i in "${chromosomes[@]}"
do
    echo "Normalising chromosome" $i "matrix"
    Rscript Normalize.R --method KRnorm --input /home/eu_plus_3/HiC/GM12878/"$res"kb_resolution_intrachromosomal/chr"$i"_"$res"kb.RAWobserved --output /home/eu_plus_3/normalisation/krnormalised_matrices/chr"$i"_"$res"kb_norm.csv;
    echo "Chromosome" $i "done"
done

