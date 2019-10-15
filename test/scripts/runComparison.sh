#!/bin/bash

baseDir=/eos/cms/store/cmst3/user/psilva/HGCal/Occupancies
outdir=/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_`date +%d%b%y`

#python prepareOccupancySummary.py --waferPlots 0,3:0,5:0,10 --onlyLayers CEE:5,CEE:15,CEH:2 \
#    noOOT,noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_noNoise_20191010.root:ana \
#    noOOT:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_20191010.root:ana \
#    noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noNoise_20191010.root:ana \
#    aged:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_pu200_20191010.root:ana \
#    -o ${outdir}

python drawOccupancySummary.py -o ${outdir} \
    aged:${outdir}/summary.pck:3 \
    noOOT,noNoise:${outdir}/summary.pck:0 \
    noOOT:${outdir}/summary.pck:1 \
    noNoise:${outdir}/summary.pck:2 \
    
cp /eos/user/p/psilva/www/index.php ${outdir}
for i in `seq 1 4`; do
    cp /eos/user/p/psilva/www/index.php ${outdir}/proc${i};
done


