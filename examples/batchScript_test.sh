#!/bin/bash

var_path=`pwd`
echo 'current dir : ' 
echo $var_path

read -p 'Enter AraSim directory : ' AraSimDir
read -p 'Enter setup file name : ' SetUpFile
read -p 'Enter run number : ' RunNo
echo $AraSimDir
echo $SetUpFile
echo $RunNo
cd $AraSimDir
    for (( i = 1; i <= $RunNo; i++ ))
    do
#qsub -l walltime=48:00:00 runaraSim_test.sh -v INPUTFILE=$SetUpFile,RUN_NO=$i,RUN_DIR=$AraSimDir
      if [ $i -eq $RunNo ]; then
        qsub -m e runaraSim_test.sh -v INPUTFILE=$SetUpFile,RUN_NO=$i,RUN_DIR=$AraSimDir
      else 
        qsub runaraSim_test.sh -v INPUTFILE=$SetUpFile,RUN_NO=$i,RUN_DIR=$AraSimDir
      fi
    done

