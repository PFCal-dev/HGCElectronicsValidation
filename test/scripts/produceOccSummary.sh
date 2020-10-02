#!/bin/bash

function process {

    sample=$1
    geom=${2}
    scaleByDoseFactor=${3}

    out_dir=/eos/cms/store/cmst3/user/psilva/HGCAL/EOL/2020-09-23
    max_events=5000
    fold=False
    adcThrMIP=0.5
    adcThrMIPbxm1=-1


    executable=${CMSSW_BASE}/src/UserCode/HGCElectronicsValidation/test/scripts/fillOccupancyHistograms.sh
    arguments="${CMSSW_BASE} ${sample} ${max_events} ${out_dir} ${fold} ${adcThrMIP} ${adcThrMIPbxm1} ${geom} ${scaleByDoseFactor}"

    mkdir -p ${out_dir}
    sh ${executable} ${arguments}
}

baseDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/

#process ${baseDir}/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723  Extended2026D49 1.0
process ${baseDir}/ttbar_D49_1120pre1_PU200_eolupdate_1p5ab_six_20200724  Extended2026D49 0.5
process ${baseDir}/ttbar_D49_1120pre1_PU200_eolupdate_startup_qua_20200723  Extended2026D49 0.0
process ${baseDir}/ttbar_D49_1120pre1_PU200_eolupdate_4ab_800V_six_20200724  Extended2026D49 1.3333

#process /eos/cms/store/group/dpg_hgcal/comm_hgcal/franzoni/ttbar_D62_1120pre3IB_PU200_ttb_ter_eolupdate_20200814 Extended2026D62
