echo "start";


chromosomes=($(seq 1 22))
chromosomes+=("X")
for i in "${chromosomes[@]}"
do
    echo "Computing overlap score for chromosome" $i 
    Rscript ComputeTopDomOverlapScores.R -chr $i >> "overlap_scores_chr${i}.txt";
    echo "Chromosome" $i "done"
done
