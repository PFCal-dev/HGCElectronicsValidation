#!/bin/bash

baseDir=/eos/cms/store/cmst3/user/psilva/HGCal/Occupancies/4Nov/
outdir=/eos/user/p/psilva/www/HGCal/Electronics/Occupancies_`date +%d%b%y`
mkdir -p $outdir

function copyphp {

    cp /eos/user/p/psilva/www/index.php ${1}
    for i in `seq 1 4`; do
        if [ -d ${1}/proc${i} ]; then
            cp /eos/user/p/psilva/www/index.php ${1}/proc${i};
        fi
    done

}

function dottbar {

#    python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
#        "min.bias":${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_Aged_20191101.root:ana \
#        "t#bar{t}":${baseDir}/ttbar_ttbar_v11_aged_unbiased_20191101.root:ana \
#        -o ${outdir}/ttbar

#    python drawOccupancySummary.py -o ${outdir}/ttbar \
#        "min.bias":${outdir}/ttbar/summary.pck:0 \
#        "t#bar{t}":${outdir}/ttbar/summary.pck:1

#    copyphp $outdir/ttbar

    python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
        "t#bar{t}":${baseDir}/ttbar_ttbar_v11_aged_unbiased_20191101.root:ana \
        "t#bar{t}biased":${baseDir}/ttbar_ak8GenJetsNoNu_100_1_ttbar_v11_aged_biased_20191101.root:ana \
        -o ${outdir}/ttbarbias

    copyphp $outdir/ttbarbias
}

function donugun {

    mkdir -p ${outdir}/nugun_noot_nonoise
    python prepareOccupancySummary.py --waferPlots -2,0:-3,0:-4,0:-5,0:-6,0:-7,0:-8,0 --onlyLayers CEE:2,CEE:5,CEE:20,CEH:2 \
        noOOT,noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_noOOT_noNoise_aged_20191101.root:ana \
        -o ${outdir}/nugun_noot_nonoise

    copyphp $outdir/nugun_noot_nonoise

#    python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
#        noOOT,noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_noOOT_noNoise_aged_20191101.root:ana \
#        noOOT:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_noOOT_aged_20191101.root:ana \
#        noNoise:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_noNoise_aged_20191101.root:ana \
#        aged:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_Aged_20191101.root:ana \
#        -o ${outdir}/nugun

#    python drawOccupancySummary.py -o ${outdir}/nugun \
#        aged:${outdir}/nugun/summary.pck:3 \
#        noOOT,noNoise:${outdir}/nugun/summary.pck:0 \
#        noOOT:${outdir}/nugun/summary.pck:1 \
#        noNoise:${outdir}/nugun/summary.pck:2 \
#        --yratioran 0.22,1.08

#    copyphp $outdir/nugun

    #python prepareOccupancySummary.py --waferPlots -3,0:-5,0:-8,0 --onlyLayers CEE:5,CEE:20,CEH:2 \
    #    aged:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_Aged_20191101.root:ana \
    #    agedv2:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_aged_20191104.root:ana \
    #    startup:${baseDir}/FlatRandomPtGunProducer_NeutrinoGun_v11_startup_20191104.root:ana \
    #    -o ${outdir}/nugunv2
    
    #copyphp $outdir/nugunv2
}

#dottbar
donugun


