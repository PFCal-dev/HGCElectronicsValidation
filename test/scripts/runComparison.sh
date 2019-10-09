#!/bin/bash

outdir=/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_8Oct

python test/scripts/prepareOccupancySummary.py \
    --waferPlots -3,0:-5,0:-10,0 \
    pfs:FlatRandomPtGunProducer_NeutrinoGun_v11_aged_20190902_numEvent2500.root:ana \
    franzoni:Franzoni_digiinfo_PU200.root:ana \
    -o ${outdir}

python test/scripts/drawOccupancySummary.py \
    -o ${outdir} \
    pfs:${outdir}/summary.pck:0 franzoni:${outdir}/summary.pck:1 

#python test/scripts/prepareOccupancySummary.py \
#    --waferPlots -3,0:-5,0:-10,0 \
#    aged:FlatRandomPtGunProducer_NeutrinoGun_v11_aged_20190902_numEvent2500.root:ana \
#    noOOT:FlatRandomPtGunProducer_NeutrinoGun_v11_aged_nooot_20190925_numEvent2500.root:ana \
#    noOOT,noNoise:FlatRandomPtGunProducer_NeutrinoGun_v11_aged_nooot_nonoise_20190927_numEvent2500.root:ana \
#    -o /eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct_new/

#python test/scripts/drawOccupancySummary.py \
#    -o /eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct/start \
#    aged:/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct/aged/summary.pck:0 \
#    noOOT:/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct/aged/summary.pck:1 \
#    noOOT,noNoise:/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct_new/aged/summary.pck:2

cp /eos/user/p/psilva/www/index.php ${outdir}
cp /eos/user/p/psilva/www/index.php ${outdir}/proc1
cp /eos/user/p/psilva/www/index.php ${outdir}/proc2

#cp /eos/user/p/psilva/www/index.php /eos/user/p/psilva/www/HGCal/Electronics/Occupancies_3Oct/
