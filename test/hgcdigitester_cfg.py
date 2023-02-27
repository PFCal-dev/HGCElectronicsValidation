import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9 as Era_Phase2
process = cms.Process("ANALYSIS", Era_Phase2)


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', 
                 '/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/SingleK0LGun_eta2p5_13_0_0_pre4_D99',
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 "input directory")
options.register('geometry', 
                 'Extended2026D99', 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 'geometry to use')
options.register('sipmMap', 
                 'SimCalorimetry/HGCalSimProducers/data/sipmParams_geom-10.txt', 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 'geometry to use')
options.register('useTDCOnsetAuto',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'use tdcOnset values that were automatically set instead of the externally specified tdcOnset_fC')
options.register('useVanillaCfg', 
                 False, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.bool, 
                 'use vanilla fe parameters from the cfg')
options.register('hardProcOnly', 
                 False, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.bool, 
                 'filter hits for hard process only (matching SimHits')
options.register('onlyROCTree', 
                 False, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.bool, 
                 'save only the ROC summary tree')
options.register('byDoseAlgo', 
                 0, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.int, 
                 'dose algorithm to use')
options.register('pxFiringRate', 
                 -1, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 'SiPM pixel firing rate')

options.parseArguments()


#set geometry/global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.Geometry%sReco_cff'%options.geometry)
process.load('Configuration.Geometry.Geometry%s_cff'%options.geometry)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 500

#source/number events to process
import os
if os.path.isdir(options.input):
    fList = ['file:'+os.path.join(options.input,f) for f in os.listdir(options.input) if '.root' in f]
else:
    fList = ['file:'+x if not x.find('/store')==0 else x for x in options.input.split(',')]
print(fList)


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fList),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                        )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#prepare digitization parameters fo the end of life
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi')

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCal_setEndOfLifeNoise
HGCal_setEndOfLifeNoise(process,byDoseAlgo=options.byDoseAlgo)

process.HGCAL_noise_heback.pxFiringRate  = cms.double(options.pxFiringRate)
process.hgchebackDigitizer.digiCfg.sipmMap = cms.string(options.sipmMap)

#analyzer
process.ana = cms.EDAnalyzer("HGCDigiTester",
                             hgceeDigitizer=process.hgceeDigitizer,                             
                             hgcehDigitizer=process.hgchefrontDigitizer,
                             hgcehsciDigitizer=process.hgchebackDigitizer,
                             hgcee_fCPerMIP=process.HGCalRecHit.HGCEE_fCPerMIP,
                             hgceh_fCPerMIP=process.HGCalRecHit.HGCHEF_fCPerMIP,
                             hgcehsci_keV2DIGI=process.HGCalRecHit.HGCHEB_keV2DIGI,
                             useTDCOnsetAuto=cms.bool(options.useTDCOnsetAuto),
                             useVanillaCfg=cms.bool(options.useVanillaCfg),
                             hardProcOnly=cms.bool(options.hardProcOnly),
                             onlyROCTree=cms.bool(options.onlyROCTree),
                             maxDeltaR=cms.double(0.4),
                             clustJetAlgo=cms.int32(2)
                         )

process.RandomNumberGeneratorService.ana = cms.PSet( initialSeed = cms.untracked.uint32(0),
                                                     engineName = cms.untracked.string('TRandom3')
                                                 )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )



print('*'*50)
print('fcPerMIP or keV2DIGI')
print(process.HGCalRecHit.HGCEE_fCPerMIP)
print(process.HGCalRecHit.HGCHEF_fCPerMIP)
print(process.HGCalRecHit.HGCHEB_keV2DIGI)
#from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import calcWeights,weightsPerLayer_V16
#weights=calcWeights(weightsPerLayer_V16)
#thickness_weights=[0.75, 0.76, 0.75, 0.85, 0.85, 0.84, 0.69]
#print('thick corections')
#print(process.HGCalRecHit.thicknessCorrection)
#print('de/dx weights')
#print(process.HGCalRecHit.layerWeights)
print('*'*50)

process.p = cms.Path(process.ana)
