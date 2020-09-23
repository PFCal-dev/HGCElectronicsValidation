import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', 
                 '/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/ttbar_D49_1120pre1_PU200_eolupdate_ter_20200713/GSD/',
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 "input directory")
options.register('adcThrMIP',
                 0.5, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "threshold (in-time)")
options.register('scaleByDoseFactor',
                 1.0, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "scale fluence by this factor")
options.register('adcThrMIPbxm1',
                 2.5, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "threshold (BX-1)")
options.register('fold',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "fold wafers histos x6/x3 for EE/HE")
options.register('geometry', 
                 'Extended2026D49', 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 'geometry to use')
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


#analyzer
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse,HGCAL_cceParams_toUse
process.ana = cms.EDAnalyzer("HGCOccupancyAnalyzer",
                             adcThrMIP       = cms.double(options.adcThrMIP),
                             adcThrMIPbxm1   = cms.double(options.adcThrMIPbxm1),
                             fold            = cms.bool(options.fold),
                             doseMap         = cms.string('SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt'),
                             scaleByDoseAlgo = cms.uint32(0),
                             scaleByDoseFactor = cms.double(options.scaleByDoseFactor),
                             ileakParam      = HGCAL_ileakParam_toUse,
                             cceParams       = HGCAL_cceParams_toUse
                        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )

process.p = cms.Path(process.ana)
