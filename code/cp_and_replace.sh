TESTS=("lasso" "ridge")
REPLACEMENTS=("insilico_size100_1" "insilico_size100_2" "insilico_size100_3" "insilico_size100_4" "insilico_size100_5")

for t in `seq 0 1`;
do
    TEST=${TESTS[t]}
    for r in `seq 0 4`;
    do
	REPLACEMENT=${REPLACEMENTS[r]}
	cp package_params_cpipeline.sh tmp

	sed -i 's/insilico_size100_1/${REPLACEMENT}/g' package_params_cpipeline.sh
	sed -i  's/enet/${TEST}/g' package_params_cpipeline.sh

	echo
	echo FOLDER IS $FOLDER
	echo
	
	source ./package_params_cpipeline.sh
	./package_for_cluster_cpipeline.sh

	scp -r -P 2222 $FOLDER jhlu@127.0.0.1:/tigress/jhlu/DREAM/$REPLACEMENT
		
	cp tmp package_params_cpipeline.sh
    done
done