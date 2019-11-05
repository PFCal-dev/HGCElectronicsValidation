import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('geometry', 'Extended2026D46', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'geometry to use')
options.register("doseMap", "",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

#ddFZ Si operation
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_cceParams_toUse as HGCAL_ddfzParams
process.siopddfz = cms.EDAnalyzer("HGCSiOperationScan",
                                  doseMap            = cms.string( options.doseMap ),
                                  doseMapAlgo        = cms.uint32(0),
                                  ileakParam         = HGCAL_ileakParam_toUse,
                                  cceParams          = HGCAL_ddfzParams
                              )

#epi Si operation
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_cceParams_toUse, cceParamFine_epi600, cceParamThin_epi600, cceParamThick_epi600
HGCAL_epiParams = cms.PSet(
    cceParamFine  = cms.vdouble(cceParamFine_epi600),
    cceParamThin  = cms.vdouble(cceParamThin_epi600),
    cceParamThick = cms.vdouble(cceParamThick_epi600)
)
process.siopepi = cms.EDAnalyzer("HGCSiOperationScan",
                                  doseMap            = cms.string( options.doseMap ),
                                  doseMapAlgo        = cms.uint32(0),
                                  ileakParam         = HGCAL_ileakParam_toUse,
                                  cceParams          = HGCAL_epiParams
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )

process.p = cms.Path(process.siopddfz
                     *process.siopepi)
