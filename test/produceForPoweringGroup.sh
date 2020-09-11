#!/bin/bash

scenarios=(
    3iab_120-600V_200-300-600V
    #3iab_120-600V_200-300-800V
    #3iab_120-800V_200-300-800V
    #4iab_120-600V_200-300-600V
    4iab_120-600V_200-300-800V
    #4iab_120-800V_200-300-800V
)

geometries=(
Extended2026D49
Extended2026D62
)

#uvmap=geomnew_corrected_withmult_F_rotations_v11.1.txt
uvmap=geomnew_corrected.txt
#uvmap=ld2hd_geomnew_corrected_360.txt
#uvmap=geomnew_corrected_360.txt

for g in ${geometries[@]}; do
    for s in ${scenarios[@]}; do
        cmsRun test/hgcsiopscan_cfg.py \
            geometry=${g} \
	    doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt \
	    uvmapfile=UserCode/HGCElectronicsValidation/data/${uvmap} \
	    savePadInfo=True scenario=${s} output=SiOpScan_${s}_${g}_v6.root &
    done
done
