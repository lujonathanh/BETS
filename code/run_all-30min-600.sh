while read script; do
    echo Submitting Script $script
    clusterize -l 30:00 -m 600 -c "time ./$script"
done <script_list.txt