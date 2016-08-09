#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=168:00:00
#PBS -r n
#PBS -V
#PBS -q hep
#PBS -j oe


cd $PBS_O_WORKDIR
echo $PWD
echo "Job Start at `date`"

#module load Apps/ROOT/5.34.34
echo source /home/hep/caines/ly247/.bashrc
source /home/hep/caines/ly247/.bashrc

#echo root -l -q -b plotTree2Histo.C++"(\"${xvariable}\",\"~/Scratch/pp200Y12_jetunderlying/\",\"${jcharge}_${tcharge}_pp${trg}_R06_HadrCorr_160227\")" # ,0,5)"	# change jet pt min #MB
#root -l -b <<EOF
#.L plotTree2Histo.C+
#plotTree2Histo("${xvariable}","~/Scratch/pp200Y12_jetunderlying/","${jcharge}_${tcharge}_pp${trg}_R06_HadrCorr_160227") 			
#.q
#EOF

jet=10
jetmax=200
echo root -l -q -b plotTree2Histo.C++"(\"${xvariable}\",\"~/Scratch/pp200Y12_jetunderlying/\",\"${tag}${jcharge}_${tcharge}_MatchTrig_pp${trg}_151030P12id_R06_HadrCorr\",${jet},${jetmax},0)" #JP2
root -l -b <<EOF
.L plotTree2Histo.C+
plotTree2Histo("${xvariable}","~/Scratch/pp200Y12_jetunderlying/","${tag}${jcharge}_${tcharge}_MatchTrig_pp${trg}_151030P12id_R06_HadrCorr",${jet},${jetmax},0)
.q
EOF

#jet=10
#jetmax=200
#echo root -l -q -b DetaDphi.C++"(\"~/Scratch/pp200Y12_jetunderlying/\",\"${tag}${jcharge}_${tcharge}_MatchTrig_pp${trg}_151030P12id_R06_HadrCorr_160314\",${jet},${jetmax},0,0,0,200,1)" #JP2
#root -l -b <<EOF
#.L DetaDphi.C+
#DetaDphi("~/Scratch/pp200Y12_jetunderlying/","${tag}${jcharge}_${tcharge}_MatchTrig_pp${trg}_151030P12id_R06_HadrCorr_160314",${jet},${jetmax},0,0,0,200,1)
#.q
#EOF

#echo root -l -q -b plotTree2Histo.C++"(\"${xvariable}\",\"~/Scratch/pp200Y12_jetunderlying/\",\"${jcharge}_${tcharge}_pythiaMB_160206\")"
#root -l -b <<EOF
#.L plotTree2Histo.C+
#plotTree2Histo("${xvariable}","~/Scratch/pp200Y12_jetunderlying/","${jcharge}_${tcharge}_pythiaMB_160206")
#.q
#EOF
#
echo "Job End at `date`"




