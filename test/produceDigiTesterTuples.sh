#!/bin/bash

setups=(
    SingleK0LGun_eta2p0_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
#    SingleK0LGun_eta2p5_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
#    SinglePhotonGun_eta2p0_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
#    SinglePhotonGun_eta2p5_13_0_0_pre4_D99_3iab,useVanillaCfg=False,byDoseAlgo=0
#    SingleK0LGun_eta2p5_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
#    SingleK0LGun_eta2p0_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
#    SinglePhotonGun_eta2p0_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
#    SinglePhotonGun_eta2p5_13_0_0_pre4_D99,useVanillaCfg=True,byDoseAlgo=1
)

baseDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production
outDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/DigiTester

for i in ${setups[@]}; do
    IFS=',' read tag useVanilla doseAlgo <<< "${i}"
    echo "*******************************************"
    echo "Starting with ${tag}"
    cmsRun test/hgcdigitester_cfg.py input=${baseDir}/${tag} output=${tag}.root ${useVanilla} ${doseAlgo} hardProcOnly=True; 
    mv -v ${tag}.root ${outDir}/${tag}.root
done
