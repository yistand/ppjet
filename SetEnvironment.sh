#!/bin/bash
# call it by "source SetEnvironment.sh"

#ly 2015.07.09 modified from KKauder's SetEnvironment.csh
#ly set environment variables for linking library


# CHANGE the following to suit your environment
declare -x BASEDIR=/home/fas/caines/ly247/Software

#### ROOT
#declare -x ROOTSYS=/home/hep/share/app/root	# already set up in ~/.barshrc

### FastJet
declare -x FASTJETDIR=${BASEDIR}/fastjet-install

### PYTHIA8
declare -x PYTHIA8DIR=${BASEDIR}/pythia8215
declare -x PYTHIA8DATA=${PYTHIA8DIR}/share/Pythia8/xmldoc

### TStarJetPicoDst structure
declare -x STARPICOPATH=${BASEDIR}/PicoCode/eventStructuredAu
#declare -x STARPICOPATH=/home/hep/caines/ly247/Scratch/pp12Pico_150407/code/StRoot/eventStructure

### RooUnfold. package to unfold jet spectra 
declare -x RooUnfold=${BASEDIR}/RooUnfold


##### On rhic21, you can use
#if ( `echo $HOST|grep -c rhic21` ) then
#    declare -x BASEDIR=/Users/putschke
#    declare -x ROOTSYS=/usr/local/root_v5.32_binary_m64/
#    declare -x FASTJETDIR=${BASEDIR}/fastjet3
#    declare -x PYTHIA8DIR=${BASEDIR}/pythia8100.new
#    declare -x PYTHIA8DATA=${PYTHIA8DIR}/xmldoc
#    declare -x STARPICOPATH=/Users/kkauder/eventStructuredAu
#endif
#


############## Done with indivivual settings.

###### Update paths
if [[ -z $LD_LIBRARY_PATH ]] 
 then
	declare -x LD_LIBRARY_PATH
fi
if [[ -z $DYLD_LIBRARY_PATH ]] 
 then
	declare -x DYLD_LIBRARY_PATH
fi

declare -x PATH=./bin:${ROOTSYS}/bin:${PATH}
declare -x CPATH=./bin:${ROOTSYS}/bin:${STARPICOPATH}:${CPATH}	#ly
declare -x LD_LIBRARY_PATH=${ROOTSYS}/lib:${FASTJETDIR}/lib:${STARPICOPATH}:${PYTHIA8DIR}/lib:${RooUnfold}:${LD_LIBRARY_PATH}
declare -x DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${FASTJETDIR}/lib:${STARPICOPATH}:${PYTHIA8DIR}/lib:${RooUnfold}:${DYLD_LIBRARY_PATH}

# ly if [[ $?TERM == 0 || $?prompt == 0 ]] exit 0

echo ''
#echo 'Setup ROOT, ktJet (including FastJet)'
echo 'Setup ktJet, STARJetPico (including FastJet)'
echo '====================================='
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "ROOT: " $ROOTSYS
echo "FastJet: " $FASTJETDIR
echo "STARPICOPATH: " $STARPICOPATH
echo "<I>---------------Info--------------------<I>"
echo ""
