import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process("ANALYSIS", Phase2C9)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', 
                 '/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer__20201022/GSD/',
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 "input directory")
options.register('geometry', 
                 'Extended2026D49', 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 'geometry to use')
options.register('useVanillaCfg', 
                 False, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.bool, 
                 'use vanilla fe parameters from the cfg')
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
    fList = ['file:'+x for x in options.input.split(',')]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fList),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                        )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#prepare digitization parameters fo the end of life
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi')
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCal_setEndOfLifeNoise
#HGCal_setEndOfLifeNoise(process,byDoseAlgo=6) #no noise, no CCE
#HGCal_setEndOfLifeNoise(process,byDoseAlgo=4) #no noise
#HGCal_setEndOfLifeNoise(process,byDoseAlgo=2) #no CCE
HGCal_setEndOfLifeNoise(process,byDoseAlgo=0) #with noise and CCE


#analyzer
process.ana = cms.EDAnalyzer("HGCDigiTester",
                             hgceeDigitizer=process.hgceeDigitizer,                             
                             hgcehDigitizer=process.hgchefrontDigitizer,
                             hgcehsciDigitizer=process.hgchebackDigitizer,
                             hgcee_fCPerMIP=process.HGCalRecHit.HGCEE_fCPerMIP,
                             hgceh_fCPerMIP=process.HGCalRecHit.HGCHEF_fCPerMIP,
                             hgcehsci_keV2DIGI=process.HGCalRecHit.HGCHEB_keV2DIGI,
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
print('thick corections')
print(process.HGCalRecHit.thicknessCorrection)
print('de/dx weights')
print(process.HGCalRecHit.layerWeights)
print('*'*50)

process.p = cms.Path(process.ana)
