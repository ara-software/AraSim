#!/bin/bash

#SBATCH --job-name=getPol
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PAS0654
#SBATCH --mail-type=FAIL
#SBATCH --time=20:00:00
#SBATCH --output=run_Pol.log   # Standard output and error log

eval 'source /users/PAS0654/jflaherty13/.bashrc' #source yours
eval 'source /users/PAS0654/jflaherty13/.bash_profile'
eval 'cvmfs'
export XDG_RUNTIME_DIR=/users/PAS0654/jflaherty13/source/AraSim/temp/
export RUNLEVEL=3
export QT_QPA_PLATFORM='offscreen'
cd /users/PAS0654/jflaherty13/source/AraSim #go to wherever you have the code


# psiMin=0
# psiMax=90
# psiInterval=1
# outDirectory="outputs/updatedIdlPulseModelNoNoise_DebugMode"
# outDirectory="outputs/updatedIdlPulseModelNoNoise"
# outDirectory="outputs/updatedIdlPulseModelHpolGain"

# outDirectory="outputs/updatedIdlPulseModel"
# outDirectory="outputs/alisaIdlPulseModel"

# outDirectory="outputs/updatedIdlPulseModelHpolGainNoNoise"

# outDirectory="outputs/updatedIdlPulseModelHpolGainNoResponse"
outDirectory="outputs/test"
# outDirectory="outputs/testNeutrino"
# outDirectory="outputs/20230710_pulserSims_1000m_noNoise_perChannelGain"
# outDirectory="outputs/20230710_pulserSims_1000m/withNoise"
# outDirectory="outputs/20230710_pulserSims_1000m/noNoise"

# outDirectory="outputs/20230716_pulserSims_1000m_NFOUR_8192/withNoise/"
psi=${SLURM_ARRAY_TASK_ID}

setupfilename="SETUP/setup_variablePsi.txt"
tempSetupFile="SETUP/setup_variablePsi=${psi}.txt"
search="CLOCK_ANGLE=\b[0-9]*"
replace="CLOCK_ANGLE=$psi"
sed  "s/$search/$replace/" $setupfilename > $tempSetupFile
./AraSim $tempSetupFile $psi $outDirectory
rm $tempSetupFile


