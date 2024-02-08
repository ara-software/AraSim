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


# rawDirectory="outputs/20230710_pulserSims_1000m/noNoise/"
rawDirectory="/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/data/A2/run_012559/split/"
recoDirectory="/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/"
# outDirectory="outputs/20230713_SpiceReco/"
# outDirectory="outputs/20231012_SpiceReco_NFOUR=8192/"
# outDirectory="outputs/20231018_SpiceReco/"

outDirectory="outputs/20231102_SpiceReco/"

# outDirectory="outputs/20231113_SpiceReco/"

subset=${SLURM_ARRAY_TASK_ID}

# Original deconvolution using a fixed butterworth filter
# ./deconvolveWaveform 2 6 12559_$subset ${rawDirectory}event012559__$subset.root ${recoDirectory}recangle_reco_out_run_12559_$subset.root $outDirectory

#New covolution using a user inputted butterworth filter
# rawDirectory="/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/data/A2/run_012559/split/"
# recoDirectory="/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/"
# outDirectory="outputs/20231102_SpiceReco/"

# freqMin=150
# freqMax=200
# ./deconvolveWaveform 2 6 12559_$subset ${rawDirectory}event012559__$subset.root ${recoDirectory}recangle_reco_out_run_12559_$subset.root $outDirectory/${freqMin}to${freqMax}/ $freqMin $freqMax

# freqMin=200
# freqMax=250
# ./deconvolveWaveform 2 6 12559_$subset ${rawDirectory}event012559__$subset.root ${recoDirectory}recangle_reco_out_run_12559_$subset.root $outDirectory/${freqMin}to${freqMax}/ $freqMin $freqMax

# freqMin=250
# freqMax=300
# ./deconvolveWaveform 2 6 12559_$subset ${rawDirectory}event012559__$subset.root ${recoDirectory}recangle_reco_out_run_12559_$subset.root $outDirectory/${freqMin}to${freqMax}/ $freqMin $freqMax

# freqMin=150
# freqMax=300
# ./deconvolveWaveform 2 6 12559_$subset ${rawDirectory}event012559__$subset.root ${recoDirectory}recangle_reco_out_run_12559_$subset.root $outDirectory/${freqMin}to${freqMax}/ $freqMin $freqMax

#Deconvolution of Spice pulser simulations with Alan's Birefrinegence code
rawDirectory="/fs/scratch/PAS0654/SPICE_birefringence_A2/"
recoDirectory="/fs/scratch/PAS0654/SPICE_birefringence_A2_interferometry/"
outDirectory="/fs/scratch/PAS0654/SPICE_birefringence_A2_deconvolution/"


freqMin=150
freqMax=300
./deconvolveWaveform 2 6 ${subset}0 ${rawDirectory}AraOut.setup_variablePsi\=${subset}0.txt.run${subset}0.root ${recoDirectory}recangle_reco_out_run_${subset}0.root $outDirectory/${freqMin}to${freqMax}/ $freqMin $freqMax
