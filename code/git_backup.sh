#!/bin/bash


export CODEDIR="/Users/jlu96/v-causal-snps/code"
export REPODIR="/Users/jlu96/causal_pipeline"
export GITDIR=$REPODIR/code
export FILELIST="git_files.txt"
export FOLDERLIST="git_folders.txt"

echo "GIT folder is " $CODEDIR

# remove the old folder and update with new, in case there were any dumb
# delets. No problem since last one was good

rm -r $GITDIR
mkdir $GITDIR

cd $CODEDIR

while read file; do
    echo Copy $file to $GITDIR
    cp $file $GITDIR
done < $FILELIST

while read folder; do
    echo Copy Folder $folder to $GITDIR
    cp -r $folder $GITDIR
done < $FOLDERLIST

cd $REPODIR

