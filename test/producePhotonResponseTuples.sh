#!/bin/bash

setups=(
    CloseByParticleGunProducer_16x_vanilla_20201114,useVanillaCfg=True,byDoseAlgo=1
    CloseByParticleGunProducer_realisticStartup_20201123,useVanillaCfg=False,byDoseAlgo=1
    CloseByParticleGunProducer_aged3iab_20201117,useVanillaCfg=False,byDoseAlgo=0
    FlatRandomPtGunProducer_16x_vanilla_20201112,useVanillaCfg=True,byDoseAlgo=1
    FlatRandomPtGunProducer_realisticStartup_20201110,useVanillaCfg=False,byDoseAlgo=1
    FlatRandomPtGunProducer_aged3iab_20201022,useVanillaCfg=False,byDoseAlgo=0
    FlatRandomPtGunProducer_aged3iab_20201023,useVanillaCfg=False,byDoseAlgo=0
)

baseDir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production
outDir=/eos/cms/store/cmst3/user/psilva/HGCAL/emResponse

for i in ${setups[@]}; do
    IFS=',' read tag useVanilla doseAlgo <<< "${i}"
    echo "*******************************************"
    echo "Starting with ${tag}"
    cmsRun test/hgcdigitester_cfg.py input=${baseDir}/${tag}/GSD output=${tag}.root ${useVanilla} ${doseAlgo}; 
    mv -v ${tag}.root ${outDir}/${tag}.root
done
