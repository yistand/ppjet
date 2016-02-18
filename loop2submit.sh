#!/bin/bash -l

trg=JP2					# MB or JP2 for STAR data

match=""
if [ $trg == "JP2" ] 
then
	match=MatchTrig_
fi


echo $trg $match 

for jetcharge in ChargeJet FullJet 
do
	for undercharge in TransNeutral TransCharged  
	do 
		for i in {0..9}
		do
			bsub -q shared -W 24:00 -J "${match}${jetcharge}${undercharge}_${trg}_${i}" "sh  submitjob.sh ${i} ${trg} ${jetcharge} ${undercharge} ${match}"
			sleep .2s
		done
	done
done



