#/bin/bash
#PBS -l nodes=1:ppn=1,walltime=10:00:00
if [ "$INPUTFILE" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

if [ "$RUN_NO" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

if [ "$RUN_DIR" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

if [ "$OUTPUT_DIR" = "" ]
then
   echo "INPUTFILE environmental variable must be set" 1>&2
   exit 1
fi

cd $RUN_DIR

./AraSim $INPUTFILE $RUN_NO $TMPDIR

pbsdcp $TMPDIR/'*' $OUTPUT_DIR



