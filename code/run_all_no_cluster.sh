while read script; do
    echo Submitting Script $script
    time ./$script
done <script_list.txt