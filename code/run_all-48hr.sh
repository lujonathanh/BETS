while read script; do
    echo Submitting Script $script
    clusterize -m 2000 -l 48:00:00 -c ./$script
done <script_list.txt