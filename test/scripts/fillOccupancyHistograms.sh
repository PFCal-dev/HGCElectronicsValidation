#!/bin/bash

cmssw=${1}
input=${2}
maxEvents=${3}
outdir=${4}
fold=${5}
adcThrMIP=${6}
adcThrMIPbxm1=${7}
geom=Extended2026D49
if [ ! -z "${8}" ]; then
    geom=${8}
fi

name=`basename ${input}`

cd $cmssw/src
eval `scram r -sh`
cd -

cmsRun $cmssw/src/UserCode/HGCElectronicsValidation/test/hgcoccupancyanalysis_cfg.py \
    input=${input}/GSD \
    output=plots.root \
    fold=${fold} \
    adcThrMIP=${adcThrMIP} adcThrMIPbxm1=${adcThrMIPbxm1} \
    maxEvents=${maxEvents} \
    geometry=${geom}

mkdir -p $outdir
thrbxm1=${adcThrMIPbxm1/./p}
thrbxm1=${thrbxm1/-/m}
fl=${fold}
finalname=${name}_thr${adcThrMIP/./p}_thrbxm1${thrbxm1}_fl${fl}
mv -v plots*.root ${outdir}/${finalname}.root

python $cmssw/src/UserCode/HGCElectronicsValidation/test/scripts/prepareOccupancySummary.py ${outdir}/${finalname}.root
