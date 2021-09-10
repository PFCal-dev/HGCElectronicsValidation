#!/bin/bash

scenarios=(
    startup_600V
    3iab_120-600V_200-300-600V
    4iab_120-600V_200-300-800V
)
for s in ${scenarios[@]}; do
    baseOpts="doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt"
    baseOpts="${baseOpts} geometry=Extended2026D86"
    baseOpts="${baseOpts} savePadInfo=True"
    baseOpts="${baseOpts} scenario=${s}"
    baseOpts="${baseOpts} uvmapfile=UserCode/HGCElectronicsValidation/data/geomCMSSW10052021_corrected.txt"

    #cmsRun test/hgcsiopscan_cfg.py ${baseOpts} wafersFromCMSSW=True output=SiOpScan_${scenarios}_cmssw_wafers.root 
    cmsRun test/hgcsiopscan_cfg.py ${baseOpts} wafersFromCMSSW=False output=SiOpScan_${scenarios}_flattxt_wafers.root 
done
