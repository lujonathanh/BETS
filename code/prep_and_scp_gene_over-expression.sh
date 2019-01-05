#!/usr/bin/env bash



export OLDGENE=CEBPB
export NEWTEST="enet"
export NEWLAGS=(1 2)
export OLDNORM=0mean1var
export NEWNORM=0mean
export NEWNULL=g


echo "NOTE FOR ALL USERS: THIS SCRIPT IS NOT ROBUST TO ERRORS AND WILL OVERWRITE OLD FILES."
echo "MAKE BACKUPS AND PROCEED WITH CAUTION"

echo "##################################"
echo "ASSUMES THAT CURRENT cpipeline FILE HAS OLD: "
echo "OLDGENE=" $OLDGENE
echo "OLDNORM=" $OLDNORM

while read NEWGENE; do


    export DELLABASEDIR=/tigress/jhlu/cpipeline/over-expression
    cp package_params_cpipeline.sh tmp

    # Set the gene
    sed "/^export/s/$OLDGENE/$NEWGENE/g" tmp > package_params_cpipeline.sh

    # First round: use OLDNORM

    for (( i=0; i < ${#NEWLAGS[@]}; i++));
    do

        echo "-------------------------"
        echo "LAG " ${NEWLAGS[i]}
        echo "-------------------------"
        #### Make each
        cp package_params_cpipeline.sh tmptmp
        sed '/^export/s/export CAUSAL=.*/export CAUSAL='"$NEWTEST"'/g; s/export LAG=.*/export LAG='"${NEWLAGS[i]}"'/g; s/export NULL=.*/export NULL='"$NEWNULL"'/g;' tmptmp > package_params_cpipeline.sh

        echo "First sed complete"
        source ./package_params_cpipeline.sh > tmptmp
        ./package_for_cluster_cpipeline.sh

        if [ $NULL == "l" ]; then
            export NULLFOLDER="local"
        else
             if [ $NULL == "g" ]; then
                export NULLFOLDER="global"
             else
                echo "NULL is not 'g' or 'l'. Error!!!"
                exit 1
             fi
        fi

        # give it the hyperparameter
        if [ "$NORMALIZATION" = "0mean1var" ];
            then
                export NORM="0mean-1var"
            else
                if [ "$NORMALIZATION" = "0mean" ];
                then
                    export NORM="0mean-unnormalized"
                else
                    echo "Error: NORMALIZATION $NORMALIZATION is not 0mean or 0mean1var"
                    exit 1
                fi
            fi

        if [ ! -d "$FOLDER/hyper" ]; then
            echo "MAKING DIR"
            mkdir $FOLDER/hyper
        fi

        export HYPERFILE=~/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/cpipeline/GGR/$NORM"_"$NULL-null_g-fdr/run/$NEWTEST-${NEWLAGS[i]}/hyper/best_hyper.p
        echo "Best hyper from " $HYPERFILE
        cp $HYPERFILE $FOLDER/hyper

        export DELLAOUTDIR=$DELLABASEDIR/$NORMALIZATION/$NULLFOLDER"_"null
        echo "FOLDER IS " $FOLDER
        echo "SCP to " $DELLAOUTDIR
        scp -r -P 2222 $FOLDER jhlu@127.0.0.1:$DELLAOUTDIR

    done

    echo
    echo
    echo "***************************************"
    echo $NEWNORM
    echo "***************************************"





    # 2nd round: replace everything with new NORM
    cp package_params_cpipeline.sh tmptmp
    sed '/^export/s/'"$OLDNORM"'/'"$NEWNORM"'/g' tmptmp > package_params_cpipeline.sh

    for (( i=0; i < ${#NEWLAGS[@]}; i++));
    do
        #### Make each
        cp package_params_cpipeline.sh tmptmp
        sed '/^export/s/export CAUSAL=.*/export CAUSAL='"$NEWTEST"'/g; s/export LAG=.*/export LAG='"${NEWLAGS[i]}"'/g; s/export NULL=.*/export NULL='"$NEWNULL"'/g;' tmptmp > package_params_cpipeline.sh

        source ./package_params_cpipeline.sh > tmptmp
        ./package_for_cluster_cpipeline.sh

        if [ $NULL == "l" ]; then
            export NULLFOLDER="local"
        else
             if [ $NULL == "g" ]; then
                export NULLFOLDER="global"
             else
                echo "NULL is not 'g' or 'l'. Error!!!"
                exit 1
             fi
        fi

        # give it the hyperparameter
        if [ "$NORMALIZATION" = "0mean1var" ];
            then
                export NORM="0mean-1var"
            else
                if [ "$NORMALIZATION" = "0mean" ];
                then
                    export NORM="0mean-unnormalized"
                else
                    echo "Error: NORMALIZATION $NORMALIZATION is not 0mean or 0mean1var"
                    exit 1
                fi
            fi


        export HYPERFILE=~/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/cpipeline/GGR/$NORM"_"$NULL-null_g-fdr/run/$NEWTEST-${NEWLAGS[i]}/hyper/best_hyper.p
        echo "Best hyper from " $HYPERFILE
        if [ ! -d "$FOLDER/hyper" ]; then
            mkdir $FOLDER/hyper
        fi

        cp $HYPERFILE $FOLDER/hyper

        export DELLAOUTDIR=$DELLABASEDIR/$NORMALIZATION/$NULLFOLDER"_"null
        echo "FOLDER IS " $FOLDER
        echo "SCP to " $DELLAOUTDIR
        scp -r -P 2222 $FOLDER jhlu@127.0.0.1:$DELLAOUTDIR

    done

    cp tmp package_params_cpipeline.sh

done < prep_and_scp_gene_list.txt