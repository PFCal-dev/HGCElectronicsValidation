#!/bin/bash

baseOpts="doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt"
baseOpts="${baseOpts} savePadInfo=False scenario=3iab_120-600V_200-300-600V"

cmsRun test/hgcsiopscan_cfg.py ${baseOpts} wafersFromCMSSW=True output=SiOpScan_cmssw.root
cmsRun test/hgcsiopscan_cfg.py ${baseOpts} \
    uvmapfile=UserCode/HGCElectronicsValidation/data/ld2hd_geomnew_corrected_360.txt \
    output=SiOpScan_elecgroup.root
