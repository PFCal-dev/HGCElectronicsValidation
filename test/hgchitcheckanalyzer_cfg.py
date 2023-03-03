import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9 as Era_Phase2
process = cms.Process("ANALYSIS", Era_Phase2)


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', 
                 '/eos/cms/store/relval/CMSSW_13_0_0_pre4/RelValCloseByPGun_CE_E_Front_120um/GEN-SIM-RECO/130X_mcRun4_realistic_v2_2026D88noPU-v1/00000/0437a2c5-a715-4f03-a8e7-e5d245a23bd0.root',
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 "input directory")
options.register('geometry', 
                 'Extended2026D88', 
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
    fList = ['file:'+x if not x.find('/store')==0 else x for x in options.input.split(',')]
print(fList)


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fList),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                        )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#analyzer
process.ana = cms.EDAnalyzer("HGCHitCheckAnalyzer")
process.p = cms.Path(process.ana)
