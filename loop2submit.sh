#!/bin/bash -l

trg=JP2
jetcharge=FullJet
undercharge=TransCharged

for i in {0..9}
do
	qsub -N ${jetcharge}${undercharge}_${trg}_${i} -v jobid=${i},trg=${trg},jcharge=${jetcharge},tcharge=${undercharge} submitjob.sh
done



