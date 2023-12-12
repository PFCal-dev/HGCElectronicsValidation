#!/bin/bash

setups=(
    #SingleK0LGun_eta1p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=1
    #SingleK0LGun_eta1p8_13_1_0_pre4_D99_startup_15ns,useVanillaCfg=False,byDoseAlgo=1
    #SingleK0LGun_eta2p0_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LGun_eta2p5_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LGun_eta2p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonGun_eta1p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=1
    #SinglePhotonGun_eta1p8_13_1_0_pre4_D99_startup_15ns,useVanillaCfg=False,byDoseAlgo=1
    #SinglePhotonGun_eta2p0_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonGun_eta2p5_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonGun_eta2p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePionGun_eta1p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=1
    #SinglePionGun_eta1p8_13_1_0_pre4_D99_startup_15ns,useVanillaCfg=False,byDoseAlgo=1
    #SinglePionGun_eta2p0_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePionGun_eta2p5_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePionGun_eta2p8_13_1_0_pre4_D99_startup,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta1p8_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    SingleK0LGun_eta1p8_13_1_0_pre4_D99_3iab_15ns,useVanillaCfg=False,byDoseAlgo=0
    SinglePionGun_eta1p8_13_1_0_pre4_D99_3iab_15ns,useVanillaCfg=False,byDoseAlgo=0
    SinglePhotonGun_eta1p8_13_1_0_pre4_D99_3iab_15ns,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta1p8_13_0_0_pre4_D99_3iab_15ns,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta2p0_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta2p5_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta2p8_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LEGun_eta1p8_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LEGun_eta1p8_13_0_0_pre4_D99_3iab_15ns,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LEGun_eta2p0_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LEGun_eta2p5_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SingleK0LEGun_eta2p8_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
    #SinglePhotonEGun_eta1p8_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SinglePhotonEGun_eta1p8_13_0_0_pre4_D99_15ns,useVanillaCfg=True,byDoseAlgo=1
    #SinglePhotonEGun_eta2p0_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SinglePhotonEGun_eta2p5_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SinglePhotonEGun_eta2p8_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SingleK0LEGun_eta1p8_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SingleK0LEGun_eta1p8_13_0_0_pre4_D99_15ns,useVanillaCfg=True,byDoseAlgo=1
    #SingleK0LEGun_eta2p0_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SingleK0LEGun_eta2p5_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
    #SingleK0LEGun_eta2p8_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
)

baseDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production
outDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/DigiTester/2023Dec4
mkdir -p ${outDir}
for i in ${setups[@]}; do
    IFS=',' read tag useVanilla doseAlgo <<< "${i}"
    echo "*******************************************"
    echo "Starting with ${tag}"
    cmsRun test/hgcdigitester_cfg.py input=${baseDir}/${tag} output=${tag}.root ${useVanilla} ${doseAlgo} hardProcOnly=True; 
    mv -v ${tag}.root ${outDir}/${tag}.root
done
