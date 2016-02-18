#!/bin/bash -l

#STAR data
jobid=$1
trg=$2
jcharge=$3
tcharge=$4
match=$5

##pythia
#jcharge=$1
#tcharge=$2

echo $PWD
echo "Job Start at `date`"

echo source /home/fas/caines/ly247/.bash_profile
source /home/fas/caines/ly247/.bash_profile

echo source /home/fas/caines/ly247/code/ppjet/SetEnvironment.sh
source /home/fas/caines/ly247/code/ppjet/SetEnvironment.sh


echo ./bin/PicoJetUnderlyingActivity "/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/_${jcharge}_${tcharge}_${match}pp${trg}_${jobid}.root" "pp${trg}" "/home/fas/caines/ly247/Scratch/run12ppQA/pp200Y12Pico${trg}_sum${jobid}.root" "0" "0" &> /home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/log/_${match}${trg}Pico${jcharge}${tcharge}_${jobid}.log
./bin/PicoJetUnderlyingActivity "/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/_${jcharge}_${tcharge}_${match}pp${trg}_${jobid}.root" "pp${trg}" "/home/fas/caines/ly247/Scratch/run12ppQA/pp200Y12Pico${trg}_sum${jobid}.root" "0" "0" &> /home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/log/_${match}${trg}Pico${jcharge}${tcharge}_${jobid}.log

#echo ./bin/PicoJetUnderlyingActivity "/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/${jcharge}_${tcharge}_${match}pp${trg}_noTOFMatch_Zdcls6000_${jobid}.root" "pp${trg}" "/home/fas/caines/ly247/Scratch/run12ppQA/pp200Y12Pico${trg}_sum${jobid}.root" "0" "0" &> /home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/log/${match}${trg}Pico${jcharge}${tcharge}_${jobid}.log
#./bin/PicoJetUnderlyingActivity "/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/${jcharge}_${tcharge}_${match}pp${trg}_noTOFMatch_Zdcls6000_${jobid}.root" "pp${trg}" "/home/fas/caines/ly247/Scratch/run12ppQA/pp200Y12Pico${trg}_sum${jobid}.root" "0" "0" &> /home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/log/${match}${trg}Pico${jcharge}${tcharge}_${jobid}.log


wait


echo "Job End at `date`"


