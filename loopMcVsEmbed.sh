#!/bin/bash -l 

for ip in 2_3  3_4 4_5 5_7 7_9 9_11 11_15 15_20 20_25 25_35 35_-1			# run12  pp
#for ip in 11_15  15_25  25_35  3_4  35_45  4_5  45_55  55_65  5_7  7_9  9_11		# HC run6 pp
do
	#echo bsub -q shared -W 20:00 -M 50000 -R "rusage[mem=50000]" -J "${ip}" " sh submitMc.sh ${ip}"
	#bsub -q shared -W 20:00 -M 50000 -R "rusage[mem=50000]" -J "${ip}" " sh submitMc.sh ${ip}"
	echo bsub -q week -W 100:00 -M 50000 -R "rusage[mem=50000]" -J "${ip}" " sh submitMc.sh ${ip}"
	bsub -q week -W 100:00 -M 50000 -R "rusage[mem=50000]" -J "${ip}" " sh submitMc.sh ${ip}"
done
