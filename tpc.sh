#!/bin/bash -l

for i in  LeadAreaPtSum  LeadAreaNtrk LeadPtAve SubAreaNtrk SubPtAve SubAreaPtSum TranPtAve TranTotPtSum TranTotNtrk #
do
	for j in  5 4
	do
		for pm in 0 1
		do 
			for abs in 0 1
			do
				./bin/MainUnfold2D ${i} ${j} 1 0 1 ${pm} ${abs} 0 		# ToUnfold, iter, jetweight, scaley, tpcsys, tpcsyspm, tpcsysabs, Rc02Mc05
			done
		done
	done
done


wait

