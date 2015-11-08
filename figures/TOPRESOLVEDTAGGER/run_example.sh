#!/bin/bash

#run ./effZFit -999 SingleMu_2012-22Jan2013_allhadrbits_TnP.root muon_selection_data Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Signal_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg1_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg2_TnP.root
#and make sure that Data and MC have more or less the same yields. If not may have a problem estimating parameters



#Signal.root and Bkg1.root Bkg2.root is where I get the templates from (make it a list)

#Signal.root = tt1l w/ isHadronictop=1; Bkg1.root=w/ isHadronictop=0; Bkg2.root all the non tt1l bg. Write a function which starts from runAllHadronicSelection output and reformat the output as descibed above (i.e: simpler version of mergeAllHadronic.root)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`



dir=TOPRESOLVEDTAGGER
maketemplates=1


#TOPMVACUT_min=-1.0
#TOPMVACUT_max=1.0
#divisions=10.0                                                                                                                                                                                                                                
#var2=$(bc -l <<<  "${TOPMVACUT_min}/${divisions}")
#var3=$(bc -l <<<  "${TOPMVACUT_max}/${divisions}")
#step=$(bc -l <<<  "$var3 - $var2")


TOPMVACUTARR="-0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0"




#echo "TOPMVACUT_min=${TOPMVACUT_min}, TOPMVACUT_max=${TOPMVACUT_max}, step=$step, echo division=$divisions"
#for i in `seq 1 10`;

for TOPMVACUT in `echo $TOPMVACUTARR`
do


#TOPMVACUT=$(bc -l <<<  "${TOPMVACUT_min} + $i*$step")



echo "TOPMVACUT=${TOPMVACUT}"

./effZFit SingleMu_2012-22Jan2013_allhadrbits_TnP.root ${dir}_${TOPMVACUT} Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Signal_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg1_TnP.root Summer12_TTJets_SemiLeptMGDecays_allhadrbits_Bkg2_TnP.root ${TOPMVACUT} ${maketemplates} >& ${dir}_${TOPMVACUT}.log

done
echo "done"