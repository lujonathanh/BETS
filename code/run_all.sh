while read script; do
    echo Submitting Script $script
    clusterize -m 2000 -c "time ./$script"
done <script_list.txt