#!/bin/bash -l

trg=JP2					# MB or JP2 or pythia

tag="BemcOrTofMatch_"		# TranPhi30_  Dijet_    kT_    NoTofMatch_ or BemcOrTofMatch_	R0.2	OR ""

for xvariable in leadjetpt   #transntrk  multiplicity  
do
	for jetcharge in FullJet  ChargeJet  
	do
		for undercharge in TransCharged TransNeutral	
		do 
			#for i in {0..9}
			#do
				#qsub -N plot${jetcharge}${undercharge}_${trg}_${i} -v jobid=${i},trg=${trg},jcharge=${jetcharge},tcharge=${undercharge},xvariable=${xvariable} submitjob.sh
			#for tag in TranPhi30_  Dijet_ kT_
			#do
				qsub -N plot${jetcharge}${undercharge}_${trg} -v trg=${trg},jcharge=${jetcharge},tcharge=${undercharge},xvariable=${xvariable},tag=${tag} submitjob.sh
				sleep .2s
			#done
		done
	done
done



