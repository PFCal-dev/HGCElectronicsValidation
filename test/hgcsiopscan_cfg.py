import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('geometry', 'Extended2026D46', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'geometry to use')
options.register("doseMap", "",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register("savePadInfo", False,  VarParsing.multiplicity.singleton, VarParsing.varType.bool)
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

#import standard Ileak and CCE parameterizations
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import ileakParam_600V
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import cceParamFine_epi600, cceParamThin_tdr600, cceParamThick_tdr600

#define the Si types to scan
siTypesToScan = cms.VPSet(
    cms.PSet( tag        = cms.string('epi80fine'),
              mipEqfC    = cms.double(65),
              cellVol    = cms.double(0.52*(80.e-4)),
              cellCap    = cms.double(74),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('epi80coarse'),
              mipEqfC    = cms.double(65),
              cellVol    = cms.double(1.18*(80.e-4)),
              cellCap    = cms.double(167),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('epi100fine'),
              mipEqfC    = cms.double(66),
              cellVol    = cms.double(0.52*(100.e-4)),
              cellCap    = cms.double(59),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('epi100coarse'),
              mipEqfC    = cms.double(66),
              cellVol    = cms.double(1.18*(100.e-4)),
              cellCap    = cms.double(134),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('epi120fine'),
              mipEqfC    = cms.double(67),
              cellVol    = cms.double(0.52*(120.e-4)),
              cellCap    = cms.double(49),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('epi120coarse'),
              mipEqfC    = cms.double(67),
              cellVol    = cms.double(1.18*(120.e-4)),
              cellCap    = cms.double(111),
              cceParam   = cms.vdouble(cceParamFine_epi600)),
    cms.PSet( tag        = cms.string('ddfz200fine'),
              mipEqfC    = cms.double(70),
              cellVol    = cms.double(0.52*(200.e-4)),
              cellCap    = cms.double(29),
              cceParam   = cms.vdouble(cceParamThin_tdr600)),
    cms.PSet( tag        = cms.string('ddfz200coarse'),
              mipEqfC    = cms.double(70),
              cellVol    = cms.double(1.18*(200.e-4)),
              cellCap    = cms.double(65),
              cceParam   = cms.vdouble(cceParamThin_tdr600)),
    cms.PSet( tag        = cms.string('ddfz300fine'),
              mipEqfC    = cms.double(73),
              cellVol    = cms.double(0.52*(300.e-4)),
              cellCap    = cms.double(20),
              cceParam   = cms.vdouble(cceParamThick_tdr600)),
    cms.PSet( tag        = cms.string('ddfz300coarse'),
              mipEqfC    = cms.double(73),
              cellVol    = cms.double(1.18*(300.e-4)),
              cellCap    = cms.double(45),
              cceParam   = cms.vdouble(cceParamThick_tdr600)),
)


#analyzer (EOL conditions)
process.siop_eol = cms.EDAnalyzer("HGCSiOperationScan",
                                  savePadInfo = cms.bool( options.savePadInfo ),
                                  doseMap     = cms.string( options.doseMap ),
                                  doseMapAlgo = cms.uint32(0),
                                  ileakParam  = cms.vdouble(ileakParam_600V),
                                  aimMIPtoADC = cms.int32(10),
                                  siTypes     = cms.VPSet(siTypesToScan)
                              )

#analyzer (startup conditions)
process.siop_startup = process.siop_eol.clone( doseMapAlgo        = cms.uint32(1) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )

process.p = cms.Path(process.siop_eol)
#                     *process.siop_startup)

