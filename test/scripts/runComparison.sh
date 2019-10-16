#!/bin/bash

baseDir=/eos/cms/store/cmst3/user/psilva/HGCal/Occupancies
outdir=/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_`date +%d%b%y`

function copyphp {

    cp /eos/user/p/psilva/www/index.php ${1}
    for i in `seq 1 4`; do
        cp /eos/user/p/psilva/www/index.php ${1}/proc${i};
    done

}

function dottbar {
#    python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
#        "#nugun":${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_pu200_20191010.root:ana \
#        "t#bar{t}":${baseDir}/ttbar_ttbar_v11_aged_unbiased_20191009.root:ana \
#        -o ${outdir}/ttbar

    python drawOccupancySummary.py -o ${outdir}/ttbar \
        "#nugun":${outdir}/ttbar/summary.pck:0 \
        "t#bar{t}":${outdir}/ttbar/summary.pck:1

    copyphp $outdir/ttbar
}

function donugun {
    #python prepareOccupancySummary.py --waferPlots -1,3:-2,0:-3,0:-4,0:-5,0:-6,0:-7,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
    #    noOOT,noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_noNoise_20191010.root:ana \
    #    -o ${outdir}/nootnonoise
    
    python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
        noOOT,noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_noNoise_20191010.root:ana \
        noOOT:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_20191010.root:ana \
        noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noNoise_20191010.root:ana \
        aged:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_pu200_20191010.root:ana \
        -o ${outdir}
    python drawOccupancySummary.py -o ${outdir} \
        aged:${outdir}/summary.pck:3 \
        noOOT,noNoise:${outdir}/summary.pck:0 \
        noOOT:${outdir}/summary.pck:1 \
        noNoise:${outdir}/summary.pck:2

    copyphp $outdir/nugun
}


dottbar


