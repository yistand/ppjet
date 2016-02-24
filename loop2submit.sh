#!/bin/bash -l

trg=MB					# MB or JP2 for STAR data

match=""
if [ $trg == "JP2" ] 
then
	match=MatchTrig_
fi

tag=""			#TranPhi30_  Dijet_ kT_

echo $trg $match 

for jetcharge in FullJet ChargeJet 
do
	for undercharge in TransCharged #TransNeutral 
	do 
		for i in {0..9}
		do
			bsub -q week -W 168:00 -M 120000 -R "rusage[mem=120000]" -J "${match}${jetcharge}${undercharge}_${trg}_${i}" "sh  submitjob.sh ${i} ${trg} ${jetcharge} ${undercharge} ${match}"
			#bsub -q shared -W 20:00 -M 50000 -R "rusage[mem=50000]" -J "${tag}${match}${jetcharge}${undercharge}_${trg}_${i}" "sh  submitjob.sh ${i} ${trg} ${jetcharge} ${undercharge} ${match} ${tag}"
			sleep .2s
		done
	done
done



