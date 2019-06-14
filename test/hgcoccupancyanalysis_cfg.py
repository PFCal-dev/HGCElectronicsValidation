import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', '/eos/cms/store/cmst3/user/psilva/CMSSW_10_6_0/TTJets/PU0', VarParsing.multiplicity.singleton, VarParsing.varType.string, "input directory")
options.parseArguments()

#set geometry/global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 500


#source/number events to process
import os
fList = ['file:'+os.path.join(options.input,f) for f in os.listdir(options.input) if '.root' in f]
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fList),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                        )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )


#analyzer
process.ana0p3 = cms.EDAnalyzer("HGCOccupancyAnalyzer",
                                mipEqThr=cms.double(0.3),
                                fudgeFactor=cms.double(1.0))
process.ana0p5 = cms.EDAnalyzer("HGCOccupancyAnalyzer",
                                mipEqThr=cms.double(0.5),
                                fudgeFactor=cms.double(1.0))
process.ana1p0 = cms.EDAnalyzer("HGCOccupancyAnalyzer",
                                mipEqThr=cms.double(1.0),
                                fudgeFactor=cms.double(1.0))


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("occupancies.root")
                               )
process.p = cms.Path(process.ana0p3+process.ana0p5+process.ana1p0)
