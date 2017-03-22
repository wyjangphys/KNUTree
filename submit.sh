#!/bin/bash

# PURPOSE: This script helps to submit massive amount of jobs
# AUTHOR: Wooyoung Jang
# DATE: 2014. 12. 22

date -u  # print current time

# First check all the arguments are given from the user.
#if [ !$1 -o !$2 ]; then
#  echo "Usage: $ ./submit.sh set[yymm] [queue name]"
#  echo "Exiting the script."
#  exit 1
#fi

USRNAME=`whoami`
#CMDEXE="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/bin/main"
#EXEC="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/exec.sh"
#OUTDIR="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/Output"
#LISTDIR="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/FileList"
PWD=`pwd`
CMDEXE=$PWD/main
EXEC=$PWD/exec.sh
LISTDIR=$PWD/FileList
OUTDIR=/eos/ams/user/w/wyjang/He3He4/Output
LOGDIR=Output/${1}/log
ERRDIR=Output/${1}/err

if [ -e $CMDEXE ]; then
  echo "\$CMDEXE found .. $CMDEXE"
else
  echo "Can not find $CMDEXE"
fi

if [ -e $EXEC ]; then
  echo "\$EXEC found ... $EXEC"
else
  echo "Can not find $EXEC"
fi

if [ -d $OUTDIR ]; then
  echo "Output files will be located at $OUTDIR"
else
  mkdir -p $OUTDIR
  echo "Warning: $OUTDIR is made to store output files."
fi

if [ -d $OUTDIR/$1 ]; then
  echo "Output files will be located at subdirectory $OUTDIR/$1"
else
  mkdir -p $OUTDIR/$1
  echo "Warning: $OUTDIR/$1 subdirectory is made to store output files."
fi

if [ -d $LOGDIR ]; then
  echo "Log files will be stored at $LOGDIR"
else
  mkdir -p $LOGDIR
  echo "Warning: $LOGDIR is made to store log files."
fi

if [ -d $ERRDIR ]; then
  echo "Error files will be stored at $ERRDIR"
else
  mkdir -p $ERRDIR
  echo "Warning: $ERRDIR is made to store error files."
fi

ls -1 $LISTDIR/$1 > templist_$1

while read irun
do
  # example submit command
  # $ /usr/bin/bsub -J /usr/bin/bsub -J 1305853512 -o Output/set_1105/log/1305853512.log -e Output/set_1105/err/1305853512.err -q 8nh -n 1 exec.sh FileList/set_1105/1305853512 Output/set_1105/1305853512.root
  echo "/usr/bin/bsub -J $irun -o $LOGDIR/$irun.log -e $ERRDIR/$irun.err -q 8nh -n 1 $EXEC $LISTDIR/$1/$irun $OUTDIR/$1/$irun.root"
  /usr/bin/bsub -J $irun -o $LOGDIR/${irun}.log -e $ERRDIR/${irun}.err -M 5000 -q 8nh -n 1 $EXEC $LISTDIR/$1/$irun $OUTDIR/$1/${irun}.root
  #sleep 5
done < "templist_$1"

# Delete the temporary file list
rm templist_$1
echo "templist_$1 is deleted."
