#!/bin/zsh
echo started `date` `pwd` on `hostname` with jobid $PBS_JOBID

#source /srv01/agrp/pitt/.zshrc
WorkDir=/storage/agrp/antonc/InOne
echo cd $WorkDir
cd $WorkDir

#Setup enviroment
echo source setup.sh
source setup.sh

#start run
echo cd $DIR
cd $DIR
mkdir -pv run${run_number}
echo cd run${run_number}
cd run${run_number}

echo current exec dir $PWD
echo $WorkDir/build/HepMCEx02 $WorkDir/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat/PiPPi0macro.in ${run_number}
$WorkDir/build/HepMCEx02 $WorkDir/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat/PiPPi0macro.in ${run_number}
echo mv PFlowNtupleFile_QCD.root $DIR/run${run_number}_homdet.root
mv PFlowNtupleFile_QCD.root $DIR/run${run_number}_homdet.root
cd ../
rm -rf run${run_number}

echo finished at `date`

