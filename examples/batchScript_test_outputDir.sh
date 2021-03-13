#!/bin/bash

var_path=`pwd`
echo 'current dir : ' 
echo $var_path

read -p 'Enter AraSim directory : ' AraSimDir

if [ "$AraSimDir" = "" ]
then 
    AraSimDir=$var_path
fi

read -p 'Enter setup file name : ' SetUpFile
read -p 'Enter run number : ' RunNo

read -p 'Enter output dir (wo last /) : ' OutputDir

echo $AraSimDir
echo $SetUpFile
echo $RunNo
echo $OutputDir

read -p 'Enter to run!' RUNNOW

cd $AraSimDir

    for (( i = 1; i <= $RunNo; i++ ))
    do
      if [ $i -eq $RunNo ]; then
        qsub -m e runaraSim_test_outputDir.sh -v INPUTFILE=$SetUpFile,RUN_NO=$i,RUN_DIR=$AraSimDir,OUTPUT_DIR=$OutputDir
      else 
        qsub runaraSim_test_outputDir.sh -v INPUTFILE=$SetUpFile,RUN_NO=$i,RUN_DIR=$AraSimDir,OUTPUT_DIR=$OutputDir
      fi
    done

