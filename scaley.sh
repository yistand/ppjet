#!/bin/bash -l

fun() {

Y=$1

root -l <<EOF
.L src/ScaleByY.cxx+
ScaleByY *sby;
sby = new ScaleByY("${Y}");
sby->Init4Read("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_170418.root","~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_170418.root");
sby->FillYRatios()
sby->WriteY("Scale${Y}_NoTofMatch_FullJet_TransCharged_160811P12id_R06_HadrCorr_VPDcut_170418.root");
.q
EOF

}


for i in  LeadAreaNtrk LeadPtAve LeadAreaPtSum SubAreaNtrk SubPtAve SubAreaPtSum TranTotNtrk TranPtAve TranTotPtSum 
do
	fun ${i}
	for j in 4 5
	do
		./bin/MainUnfold2D ${i} ${j} 0 1 0 0 0 0		 		# ToUnfold, iter, jetweight, scaley, tpcsys, tpcsyspm, tpcsysabs, Rc02Mc05
	done
done


