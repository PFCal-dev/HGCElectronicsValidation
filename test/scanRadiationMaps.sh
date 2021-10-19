#!/bin/bash

dosemap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt
geometry=GeometryExtended2026D86Reco
for conditions in TDR_600V CERN21_600V_10m CERN21_600V_30m CERN21_600V_90m CERN21_600V_120m; do
    cmsRun ../../SimCalorimetry/HGCalSimAlgos/test/hgcsiNoiseMapTester_cfg.py \
        doseMap=${dosemap} geometry=${geometry} conditions=${conditions}
    python test/scripts/drawRadiationMapPlots.py dosemap_output_${geometry}_${conditions}.root 
done
