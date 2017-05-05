for j in TranTotNtrk TranTotPtSum LeadAreaPtSum SubAreaPtSum LeadAreaNtrk LeadPtAve SubAreaNtrk SubPtAve TranPtAve
do for i in {1..8}
	do for w in 0 1
		do ./bin/MainUnfold2D ${j} ${i} ${w} 0 0 0 0 1		 # ToUnfold, iter, jetweight, scaley, tpcsys, tpcsyspm, tpcsysabs, Rc02Mc05
		done
	done
done


./tpc.sh

./scaley.sh

