#!/bin/bash

inputDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/
outdir=/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_`date +%d%d%Y`
#!/bin/bash
cmssw=${1}
input=${2}
maxEvents=${3}
outdir=${4}

name=`basename ${input}`

cd $cmssw/src
eval `scram r -sh`
cd -

cmsRun $cmssw/src/UserCode/HGCElectronicsValidation/test/hgcoccupancyanalysis_cfg.py input=${input}/GSD output=plots.root maxEvents=${maxEvents} 

mkdir -p $outdir
mv -v plots*.root ${outdir}/${name}.root
