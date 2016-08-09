#!/bin/bash -l

##pythia
jcharge=$1
tcharge=$2

echo $PWD
echo "Job Start at `date`"

echo source /home/fas/caines/ly247/.bash_profile
source /home/fas/caines/ly247/.bash_profile

echo source /home/fas/caines/ly247/code/ppjet/SetEnvironment.sh
source /home/fas/caines/ly247/code/ppjet/SetEnvironment.sh


#echo ./bin/STARPythiaJetUnderlyingActivity "/home/fas/caines/ly247/scratch/pythiadata/${jcharge}_${tcharge}_pythiaMB_160206.root" "pythia" "/home/fas/caines/ly247/scratch/pythiadata/pythia_MB_pp200tree_151117.root" "0" "0"
#./bin/STARPythiaJetUnderlyingActivity "/home/fas/caines/ly247/scratch/pythiadata/${jcharge}_${tcharge}_pythiaMB_160206.root" "pythia" "/home/fas/caines/ly247/scratch/pythiadata/pythia_MB_pp200tree_151117.root" "0" "0"

./bin/EmbedPythiaJetUnderlyingActivity "/home/fas/caines/ly247/Scratch/embedPythia/${jcharge}_${tcharge}_embedPythia_160509.root" "All" "/home/fas/caines/ly247/Scratch/embedPythia/embedpp200Y12Pico_160509.root"  "0"  "0"

wait


echo "Job End at `date`"


