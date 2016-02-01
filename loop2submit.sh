#!/bin/bash -l

trg=JP2					# MB or JP2
#jetcharge=ChargeJet			# FullJet or ChargeJet
#undercharge=TransNeutral		# TransCharged or TransNeutral


for jetcharge in ChargeJet FullJet
do
	for undercharge in TransCharged TransNeutral
	do 
		for i in {0..9}
		do
			bsub -W 24:00 -J "Match${jetcharge}${undercharge}_${trg}_${i}" "sh  submitjob.sh ${i} ${trg} ${jetcharge} ${undercharge}"
			sleep .2s
		done
	done
done



