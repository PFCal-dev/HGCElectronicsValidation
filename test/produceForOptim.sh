#!/bin/bash

scenarios=(
    startup_600V
    3iab_120-600V_200-300-600V
    4iab_120-600V_200-300-800V
)

dosemap=doseParams_3000fb_fluka-6.2.0.1.txt
geometry=Extended2026D86
flatfile=geomCMSSW10052021_corrected.txt

dosemap=doseParams_3000fb_fluka-3.7.20.txt
geometry=Extended2026D49
flatfile=geomnew_corrected_withmult_F_rotations_v11.1.txt

for s in ${scenarios[@]}; do
    baseOpts="doseMap=SimCalorimetry/HGCalSimProducers/data/${dosemap}"
    baseOpts="${baseOpts} geometry=${geometry}"
    baseOpts="${baseOpts} savePadInfo=True"
    baseOpts="${baseOpts} scenario=${s}"
    baseOpts="${baseOpts} uvmapfile=UserCode/HGCElectronicsValidation/data/${flatfile}"

    cmsRun test/hgcsiopscan_cfg.py ${baseOpts} wafersFromCMSSW=False output=SiOpScan_${s}_${geometry}.root 
done
