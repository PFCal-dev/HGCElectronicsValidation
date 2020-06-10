#!/bin/bash

scenarios=(
    3iab_120-600V_200-300-600V
)
for s in ${scenarios[@]}; do
    baseOpts="doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt"
    baseOpts="${baseOpts}  savePadInfo=True scenario=${s}"
    cmsRun test/hgcsiopscan_cfg.py ${baseOpts} \
        uvmapfile=UserCode/HGCElectronicsValidation/data/ld2hd_geomnew_corrected_360.txt \
        output=SiOpScan_mnoy_${s}_v5.root 
    cmsRun test/hgcsiopscan_cfg.py ${baseOpts} \
        uvmapfile=UserCode/HGCElectronicsValidation/data/geomnew_corrected_360.txt \
        output=SiOpScan_${s}_v5.root 
done
