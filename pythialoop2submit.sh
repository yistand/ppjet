#!/bin/bash -l

trg=pythia


for jetcharge in FullJet  # ChargeJet 
do
	for undercharge in TransCharged  # TransNeutral
	do 
		bsub -W 24:00 -J "pythia${jetcharge}${undercharge}" "sh  pythiasubmitjob.sh ${jetcharge} ${undercharge}"
		sleep .2s
	done
done



