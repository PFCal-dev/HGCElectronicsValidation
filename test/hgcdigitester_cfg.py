import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('geometry', 'Extended2026D49', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'geometry to use')
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

#process.MessageLogger = cms.Service("MessageLogger",
#    debugModules = cms.untracked.vstring('*'),
#    destinations = cms.untracked.vstring('cout'),
#    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#    )
#)


#source/number events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

#prepare digitization parameters fo the end of life
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi')
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCal_setEndOfLifeNoise
HGCal_setEndOfLifeNoise(process)

#analyzer
process.ana = cms.EDAnalyzer("HGCDigiTester",
                             hgceeDigitizer=process.hgceeDigitizer)

process.RandomNumberGeneratorService.ana = cms.PSet( initialSeed = cms.untracked.uint32(0),
                                                     engineName = cms.untracked.string('TRandom3')
                                                 )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )

process.p = cms.Path(process.ana)
