source ./package_params_rand_row.sh
while read script; do
    echo Submitting Parallel Script $script
    clusterize -l $TIMEUPPERBOUND:00 -m $MEMUPPERBOUND -p $(($PARALLELNUM + 1))  -c "time $script & sleep "$TIMEUPPERBOUND"m"
done <parallel_script_list.txt