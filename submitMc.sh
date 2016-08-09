#!/bin/bash -l  

echo "Job Start at `date`"

ip=$1

#echo ./bin/PicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/pt${ip}_JetMcVsEmbedMatchTrig_Online.root ppJP2 /home/fas/caines/ly247/Scratch/embedPythia/160728/pp12Pico_pt${ip}\*.root
#./bin/PicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/pt${ip}_JetMcVsEmbedMatchTrig_Online.root ppJP2 /home/fas/caines/ly247/Scratch/embedPythia/160728/pp12Pico_pt${ip}\*.root
#echo ./bin/HPicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/HCpt${ip}_JetMcVsEmbedMatchTrig.root ppJP2 /scratch/fas/caines/hlc7/pp2006PythiaBBC/${ip}/picoDst_${ip}\*.root
#./bin/HPicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/HCpt${ip}_JetMcVsEmbedMatchTrig.root ppJP2 /scratch/fas/caines/hlc7/pp2006PythiaBBC/${ip}/picoDst_${ip}\*.root
echo ./bin/PicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/pt${ip}_JetMcVsEmbed_nobbc.root ppJP2 /home/fas/caines/ly247/Scratch/embedPythia/160804/pp12Pico_pt${ip}\*.root
./bin/PicoJetMcVsEmbed /home/fas/caines/ly247/Scratch/embedPythia/pt${ip}_JetMcVsEmbed_nobbc.root ppJP2 /home/fas/caines/ly247/Scratch/embedPythia/160804/pp12Pico_pt${ip}\*.root

wait

echo "Job End at `date`"
