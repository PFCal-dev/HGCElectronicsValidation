# HGCElectronicsValidation

A set of analyzers/analysis scripts to help with the validation of the electronics simulation in CMSSW.
These tools are used to debug the development on the main cmssw repository.
On top of your current CMSSW work area do

```
cmsrel CMSSW_13_1_0_pre4
cd CMSSW_13_1_0_pre4/src/
cmsenv

#OPTIONAL (mostly if you need to do some development in the base code)
git cms-addpkg SimCalorimetry/HGCalSimProducers
git cms-addpkg SimCalorimetry/HGCalSimAlgos
#END OPTIONAL

git clone https://github.com/PFCal-dev/HGCElectronicsValidation.git UserCode/HGCElectronicsValidation
scram b -j 8
```

When working on the packages and committing changes remember to pull and merge locally, before pushing the code.

```
git pull https://github.com/PFCal-dev/cmssw ${work_branch}
git push https://github.com/PFCal-dev/cmssw HEAD:${work_branch}
```

## Digi tester

An extensive comparison of the digis and the sim hits can be performed using the HGCDigiTester.cc class.
The files produced are large (approx. 6M hits per event) so don't be greedy running the tool.
The option 'hardProcOnly' can be used to store only the hits associated to MC truth (i.e. not from pileup) 
which greatly reduces the footprint in particular for particle guns.
An example of how to run it is the following:

```
cmsRun test/hgcdigitester_cfg.py \
       input=/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723/GSD \
       output=test.root \                                                                                               
       useVanillaCfg=False \
       byDoseAlgo=0 \
       hardProcOnly=False \
       maxEvents=10
```

The output ROOT file will contain a tree with information about sim and digitized hits.

## Radiation map analysis

A set of control plots for the noise, charge collection efficiency, fluence is produced based on a radiation map file
and on a geometry. These analyzers can be found in SimCalorimetry/HGCalSimAlgos/test. An example is given below, assuming the current directory location

```
cmsRun ../../SimCalorimetry/HGCalSimAlgos/test/hgcsiNoiseMapTester_cfg.py \
       doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt \
       geometry=GeometryExtended2026D86Reco \
       conditions=TDR_600V
python3 test/scripts/drawRadiationMapPlots.py dosemap_output_GeometryExtended2026D86Reco_TDR_600V.root 
cmsRun ../../SimCalorimetry/HGCalSimAlgos/test/hgchebacksignalscaler_cfg.py \
       doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt \
       sipmMap=SimCalorimetry/HGCalSimProducers/data/sipmParams_geom-10.txt \
       geometry=GeometryExtended2026D86Reco
python3 test/scripts/drawRadiationMapPlots.py sipmontile_dosemap_GeometryExtended2026D86Reco.root sci
```

You can also run `test/scanRadiationMaps.sh` for all the scenarios.

## Wafer map

The class src/HGCGeometryScan.cc dumps a txt file with the positions of wafers at the same radius in each layer.
This file is useful for the design studies and it's used as input for the occupancy analysis.
To run it you can do:

```
cmsRun test/hgcgeometryscan_cfg.py geometry=Extended2026D46
```

The output is a ROOT file with TGraph2D showing the position of the wafers in each layer.
You can plot them with `python test/scripts/drawGeometry.py` and an ascii file called `wafer_pos.dat`.
To use the later in the occupancy analyzer copy it to the data folder.

## Wafer position optimization

Prepare a ntuple with the fluence and S/N in each wafer, layer and sub-detector

```
cmsRun test/hgcsiopscan_cfg.py \
       doseMap=SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt \
       uvmapfile=UserCode/HGCElectronicsValidation/data/geomnew_corrected_360.txt \
       savePadInfo=True \
       scenario=3iab_600V \
       output=SiOpScan_3iab_600V.root
```

Then you can use the SiSensorOptim notebook to interact with its contents.

## Occupancy analysis

A set of control plots is produced globally per layer and per wafer-equivalent position in each layer.
The number of plots is very large and an additional script is used to loop over them
and produce a summary in terms of quantiles of the distributions.
A third script is used to build final profiles (function of radius, layer, etc.)
or to dump an ascii file with the summarized information per wafer equivalent and layer.

The filling of the histograms can be called by running src/HGCOccupancyAnalyzer.cc with

```
cmsRun test/hgcoccupancyanalysis_cfg.py input=input_dir
```

The histograms can be digested into a summary of occupancies etc. by running the following script

```
python test/scripts/prepareOccupancySummary.py occupancies.root
```

After having produced the summary of the occupancies one can compare different results by using the jupyter notebook


## Generator-level jet profiles

May be useful to profile the generator level jets which are being thrown on the HGCal surface area
to get an idea of its characteristics (pT, eta, em fraction spectra etc.).
The first script fills the histograms using the genJets collection, the second one dumps the plots

```
python test/scripts/genJetAnalyzer.py input_dir output_name
python test/scripts/compareJetProfiles.py name1=file1.root ...
```
