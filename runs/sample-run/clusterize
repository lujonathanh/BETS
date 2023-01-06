#!/bin/bash

usage="clusterize [options] -c \"command\"
	Creates a script in the temp directory to run the user specified command
	with parameters defined by [options], submits it to the cluster, then cleans
	the script file.

where:

	-c <command>
		REQUIRED. The command <command> to be executed surrounded by quotes
	-h
		Shows this help text.

	-w <dir>
		Sets the working directory of the executing script to <dir> (default is pwd,
		for your current directory it is: `pwd`).

	-i
		Starts an interactive session. Only supports -l and -m.
		Emulates srun --partition=general --qos=general --pty bash

	-l <length>
		Sets the upper bound on the amount of time required by this job to <length>.
		It's important to have an accurate estimate. Some schedulers might terminate
		your job if you go over this time and all schedulers will give you higher
		priority if your job is requires less time (default is one hour: 01:00:00).
	-g <gpus>
		Enables GPU and set the number of requested GPUs to <gpus>

	-p <ppn>
		Sets the number of processors per node to <ppn> (default 1)

	-x <wt>
		Sets weight time to submission to <wt> (default 1)

	-n <nodes>
		Sets the number of nodes required by this job to <nodes> (default 1)

	-m <mem>
		Sets the amount of memory required per node to <mem> (default 1G)

	-t
		Test this script for your command. This will create the temporary file,
		tell you where it is, and then abort. You should remember to clean up
		your temporary file after done verifying it.

	-o
		Do NOT combine the std out and std err streams (default is to combine them)"


# these parameters you will likely have to change, either defaults or through input to scripts
WORKING_DIR="`pwd`"
LENGTH="01:00:00"
MEMORY="1G"
COMMAND="${*,2}"
NODES="1"
PPN="1"
INTERACTIVE="0"
GPUS="0"
WT="1"
WORKSET="0"
TESTSCRIPT="0"
IFS=' '
read -r PROGRAM string <<< "$1"
COMB="1"
NAME=""

if [ "$COMMAND" == "-h" -o $# -eq 0 ]; then
  echo $usage
  exit
fi

while [[ $# > 1 ]]
do
key="$1"
shift
COMMAND=${*:1}
read -r PROGRAM string <<< "$1"
case $key in
    -w|-d|--workingdir|--dir)
    WORKING_DIR="$1"
    WORKSET="1"
    shift
    ;;
    -l|--length)
    LENGTH="$1"
    shift
    ;;
    -i|--interactive)
    INTERACTIVE="1"
    ;;
    -e|--name)
    NAME="$1"
    shift
    ;;
    -h|--help)
    echo $usage
	exit
	;;
	-p|--ppn)
    PPN="$1"
    shift
    ;;
	-g|-G)
    GPUS="$1"
    shift
    ;;
	-x|--wt)
    WT="$1"
    shift
    ;;
	-n|--nodes)
    NODES="$1"
    shift
    ;;
    -m|--memory)
    MEMORY="$1"
    shift
    ;;
    -t|--test)
    TESTSCRIPT="1"
    ;;
    -o|--nocomb)
    COMB="0"
    ;;
    -h|--help)
    echo $usage
	exit
    ;;
    -c|--command)
    COMMAND="$1"
    read -r PROGRAM string <<< "$1"
    shift
    ;;
    *)  # unknown option
	echo "unknown option " $key
        echo $usage
	return
    ;;
esac
done

BATCH_JOB_NAME=`whoami`"_"$RANDOM
SLURM_JOB_NAME=$BATCH_JOB_NAME

if [[ "$INTERACTIVE" == "1" ]] || [[ $COMMAND == *"-i"* ]]; then
	echo "starting interactive job..."
	srun --qos=general --mem=$MEMORY --time=$LENGTH --pty bash
	#srun --partition=general --qos=general --mem=$MEMORY --time=$LENGTH --pty bash #edit 06/19/2019 since partition shouldn't be specified
fi

echo "Running clusterize..."
echo "Make sure you have correctly set the job submission parameters.. what I'm using for job:"
echo "Working directory = " $WORKING_DIR
echo "Upper bound on time limit (try to be accurate!) = " $LENGTH
echo "Memory per node = " $MEMORY
echo "Number of nodes to use = " $NODES
echo "Number of processors per node to use = " $PPN

TMPFILE=`mktemp` || exit 1

POSTCOMMAND="/usr/bin/time -f '%E elapsed,%U user,%S system,%M memory,%K avgmem, %x status' "$COMMAND

echo "Command to run = " $POSTCOMMAND


echo "#!/bin/sh" >> $TMPFILE
echo "" >> $TMPFILE
echo "# Request runtime" >> $TMPFILE
echo "#SBATCH --time=$LENGTH" >> $TMPFILE
echo "" >> $TMPFILE



if [ "$GPUS" == "0" ]; then
	echo "# Request a number of CPUs per task:" >> $TMPFILE
	echo "#SBATCH --cpus-per-task=$PPN" >> $TMPFILE
	echo "" >> $TMPFILE
elif [ "$GPUS" -gt 1 ]; then
	echo "# Request a number of GPUs per task:" >> $TMPFILE
        echo "#SBATCH --ntasks-per-socket=2" >> $TMPFILE
        echo "#SBATCH --gres=gpu:$GPUS" >> $TMPFILE
        echo "" >> $TMPFILE
else
        echo "# Request a number of GPUs per task:" >> $TMPFILE
        echo "#SBATCH --gres=gpu:$GPUS" >> $TMPFILE
        echo "" >> $TMPFILE
fi
echo "# Request a number of nodes:" >> $TMPFILE
if [[ "${MEMORY: -1}" -eq "G" || "${MEMORY: -1}" -eq "g" ]] && [[ "${MEMORY::-1}" -gt 256 ]]; then
	#echo "#SBATCH --partition=himem" >> $TMPFILE  #edit 06/19/2019 since partition shouldn't be specified
	echo "#SBATCH --qos=himem" >> $TMPFILE
else
	#echo "#SBATCH --partition=general" >> $TMPFILE  #edit 06/19/2019 since partition shouldn't be specified
	echo "#SBATCH --qos=general" >> $TMPFILE
fi
echo "" >> $TMPFILE
echo "# Request an amount of memory per node:" >> $TMPFILE
echo "#SBATCH --mem="$MEMORY >> $TMPFILE
echo "# usage of --mem and --mem-per-cpu is mutually exclusive, do NOT use both" >> $TMPFILE
echo "# by default in this script --mem is used" >> $TMPFILE
echo "" >> $TMPFILE
echo "# Specify a job name:" >> $TMPFILE
echo "#SBATCH -J $BATCH_JOB_NAME" >> $TMPFILE
echo "" >> $TMPFILE

if [ -z "$NAME" ]; then
	echo "#SBATCH -o "$BATCH_JOB_NAME"_%j.out" >> $TMPFILE
else
	echo "#SBATCH -o "$NAME"_%j.out" >> $TMPFILE
	echo "" >> $TMPFILE
fi
echo "" >> $TMPFILE
if [ "$COMB" == "0" ]; then
  echo "# Specify an output file" >> $TMPFILE
  echo "#SBATCH -o "$BATCH_JOB_NAME"_%j_o.out" >> $TMPFILE
  echo "#SBATCH -e "$BATCH_JOB_NAME"_%j_e.out" >> $TMPFILE
fi


if [ "$WORKSET" == "1" ]; then
  echo "cd "$WORKING_DIR >> $TMPFILE
else
  printf "cd %s\n" `pwd` >> $TMPFILE
fi

echo "" >> $TMPFILE
echo "# Set working directory:" >> $TMPFILE
echo "#SBATCH --workdir=$WORKING_DIR" >> $TMPFILE
echo "" >> $TMPFILE
POSTCOMMAND="/usr/bin/time -f '%E elapsed,%U user,%S system,%M memory,%K avgmem, %x status' "$COMMAND


echo "echo \"" $POSTCOMMAND "\" 1>&2"  >> $TMPFILE
echo "" >> $TMPFILE
echo $POSTCOMMAND >> $TMPFILE
echo "" >> $TMPFILE



if [ "$TESTSCRIPT" == "1" ]; then
  echo "file containing your script is in " $TMPFILE
  echo "please clean this file (rm "$TMPFILE") when done"
else
  echo "Executing script $TMPFILE"
  sbatch $TMPFILE
  sleep $WT
fi
