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
# outDirectory="outputs/test"
# outDirectory="outputs/testNeutrino"
# outDirectory="outputs/20230710_pulserSims_1000m_noNoise_perChannelGain"
# outDirectory="outputs/20230710_pulserSims_1000m/withNoise"



# rawDirectory="outputs/20230710_pulserSims_1000m/withNoise/"
# recoDirectory="outputs/20230710_pulserSims_1000m/interferometry/"
# outDirectory="outputs/20230710_pulserSims_1000m/withNoise/deconvolution/"

# rawDirectory="outputs/20230710_pulserSims_1000m/noNoise/"
# recoDirectory="outputs/20230710_pulserSims_1000m/interferometry/"
# outDirectory="outputs/20230710_pulserSims_1000m/noNoise/deconvolution/"

# rawDirectory="outputs/20230716_pulserSims_1000m_NFOUR_8192/withNoise/"
# recoDirectory="outputs/20230710_pulserSims_1000m/interferometry/"
# outDirectory="outputs/20230716_pulserSims_1000m_NFOUR_8192/deconvolution/"






#Deconvolution used in ICRC Proceedings
# rawDirectory="outputs/20230710_pulserSims_1000m/withNoise/"
# recoDirectory="outputs/20230710_pulserSims_1000m/interferometry/"
# outDirectory="outputs/20230710_pulserSims_1000m/withNoise/deconvolution/"


#Testing of deconvolution with larger NFOUR to narrow frequency bins to match spice core deconvolution
rawDirectory="outputs/20230710_pulserSims_1000m/withNoise/"
recoDirectory="outputs/20230710_pulserSims_1000m/interferometry/"
outDirectory="outputs/20230710_pulserSims_1000m/withNoise/deconvolution_NFOUR=4096/"

psi=${SLURM_ARRAY_TASK_ID}

./deconvolveWaveform 2 6 $psi ${rawDirectory}AraOut.setup_variablePsi\=$psi.txt.run$psi.root ${recoDirectory}recangle_reco_out_run_$psi.root $outDirectory



