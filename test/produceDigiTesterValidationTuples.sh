#!/bin/bash

setups=(
    orig $CMSSW_BASE/src/23234.103_TTbar_14TeV+2026D94Aging3000_orig,useVanillaCfg=False,byDoseAlgo=0
    pr $CMSSW_BASE/src/23234.103_TTbar_14TeV+2026D94Aging3000,useVanillaCfg=False,byDoseAlgo=0
)

for i in ${setups[@]}; do
    IFS=',' read tag indir useVanilla doseAlgo <<< "${i}"
    echo "*******************************************"
    echo "Starting with ${tag}"
    cmsRun test/hgcdigitester_cfg.py input=${indir}/step3.root output=${tag}.root ${useVanilla} ${doseAlgo} hardProcOnly=False; 
done
