#!/usr/bin/env bash

LAGS=("1" "2")
OLDNORM=0mean
NORMS=("0mean" "0mean1var")
DATANUMS=("1" "2" "3" "4" "5")
OLDGENE="insilico_size100_1"
GENES=()
for x in "${DATANUMS[@]}"; do
    GENES+=("insilico_size100_$x")
    done
echo "${GENES[@]}"

for d in `seq 0 4`;
do
    for l in `seq 0 1`;
    do
        for n in `seq 0 1`;
        do
            GENE=${GENES[d]}
            LAG=${LAGS[l]}
            NORM=${NORMS[n]}

            cp package_params_cpipeline.sh tmptmp
            sed '/^export/s/export LAG=.*/export LAG='"$LAG"'/g; s/'"$OLDNORM"'/'"$NORM"'/g; s/'"$OLDGENE"'/'"$GENE"'/g ' tmptmp > package_params_cpipeline.sh


            # sed '/^export/s/'"$OLDNORM"'/'"$NORM"'/g'

            echo "*********************"
            echo

            echo $GENE
            echo $LAG
            echo $NORM

            # more package_params_cpipeline.sh

            ./package_for_cluster_cpipeline.sh

            echo "*********************"
            echo

            cp tmptmp package_params_cpipeline.sh

        done
    done
done