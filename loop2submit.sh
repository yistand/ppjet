#!/bin/bash -l


for i in {0..9}
do
	qsub -N PicoJetUnderlyingEvent_JP2_${i} -v jobid=${i} submitjob.sh
done
