source ./package_params.sh
while read script; do
    echo Submitting Parallel Script $script
    clusterize -l 30:00 -m 600 -p $PARALLELNUM  -c "time $script"
done <parallel_script_list.txt