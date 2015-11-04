#!/bin/bash -l
#PBS -N PicoJetUnderlyingEvent_JP2
#PBS -l nodes=1:ppn=1,walltime=30:00:00
#PBS -r n
#PBS -V
#PBS -q hep
#PBS -j oe

cd $PBS_O_WORKDIR
echo $PWD
echo "Job Start at `date`"

echo source SetEnvironment.sh
source SetEnvironment.sh


# make sure executable exists
echo make bin/PicoJetUnderlyingActivity || exit
make bin/PicoJetUnderlyingActivity || exit

echo "./bin/PicoJetUnderlyingActivity > Pico.log"
./bin/PicoJetUnderlyingActivity > Pico.log

echo "Job End at `date`"

