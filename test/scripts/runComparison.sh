#!/bin/bash

python test/scripts/prepareOccupancySummary.py \
    --waferPlots -3,0:-5,0:-10,0 \
    t#bar{t}:ttbar_ttbar_v11_aged_biased_20191009_numEvent${maxEvents}.root:ana \
    t#bar{t},biased:ttbar_ak8GenJetsNoNu_100_1_ttbar_v11_aged_biased_20191009_numEvent${maxEvents}.root:ana \
    -o ${outdir}

python test/scripts/drawOccupancySummary.py \
    -o ${outdir} \
    t#bar{t}:${outdir}/summary.pck:0 \
    t#bar{t},biased:${outdir}/summary.pck:1

cp /eos/user/p/psilva/www/index.php ${outdir}
cp /eos/user/p/psilva/www/index.php ${outdir}/proc1
cp /eos/user/p/psilva/www/index.php ${outdir}/proc2


