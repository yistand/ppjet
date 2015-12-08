#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=30:00:00
#PBS -r n
#PBS -V
#PBS -q hep
#PBS -j oe


cd $PBS_O_WORKDIR
echo $PWD
echo "Job Start at `date`"

#module load Apps/ROOT/5.34.34
echo source /home/hep/caines/ly247/.bashrc
source /home/hep/caines/ly247/.bashrc

echo source SetEnvironment.sh
source SetEnvironment.sh

## make sure executable exists
#echo make bin/PicoJetUnderlyingActivity || exit
#make bin/PicoJetUnderlyingActivity || exit


echo "./bin/PicoJetUnderlyingActivity "/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/TransCharge0_MatchTrig_ppJP2_${jobid}.root" "ppJP2" "/home/hep/caines/ly247/Scratch/pp12Pico_151018/sum${jobid}.root" "0" "0" &> Pico${jobid}.log"
./bin/PicoJetUnderlyingActivity "/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/TransCharge0_MatchTrig_ppJP2_${jobid}.root" "ppJP2" "/home/hep/caines/ly247/Scratch/pp12Pico_151018/sum${jobid}.root" "0" "0" &> Pico${jobid}.log

echo "Job End at `date`"

