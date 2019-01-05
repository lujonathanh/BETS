mkdir cleanup
echo "Moving all cleanup files to cleanup/"
while read name; do
    mv $name cleanup
done < cleanup_list.txt