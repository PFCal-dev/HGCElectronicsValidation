import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('geometry', 'Extended2026D46', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'geometry to use')
options.register("doseMap", "",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register("scenario", "600V",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
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

#tfileservice
process.TFileService = cms.Service("TFileService",fileName = cms.string(options.output))

#import standard Ileak and CCE parameterizations
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import ileakParam_600V,ileakParam_800V
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import cceParamFine_epi600, cceParamFine_tdr600, cceParamThin_tdr600, cceParamThick_tdr600
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import cceParamFine_tdr800, cceParamThin_tdr800, cceParamThick_tdr800

#define the Si types to scan
def defineSiTypesToScan(cceParamEpi, cceParamThin, cceParamThick, qefC=1.60217646e-4):

    """returns a dict of interesting sensors and their main properties"""
    
    return { 'epi80fine': cms.PSet( tag        = cms.string('epi80fine'),
                                    mipEqfC    = cms.double(65*80*qefC),
                                    cellVol    = cms.double(0.52*80e-4),
                                    cellCap    = cms.double(74),
                                    cceParam   = cms.vdouble(cceParamEpi)),
             'epi80coarse': cms.PSet( tag        = cms.string('epi80coarse'),
                                      mipEqfC    = cms.double(65*80*qefC),
                                      cellVol    = cms.double(1.18*80e-4),
                                      cellCap    = cms.double(167),
                                      cceParam   = cms.vdouble(cceParamEpi)),
             'epi100fine': cms.PSet( tag        = cms.string('epi100fine'),
                                     mipEqfC    = cms.double(66*100*qefC),
                                     cellVol    = cms.double(0.52*100e-4),
                                     cellCap    = cms.double(59),
                                     cceParam   = cms.vdouble(cceParamEpi)),
             'epi100coarse': cms.PSet( tag        = cms.string('epi100coarse'),
                                       mipEqfC    = cms.double(66*100*qefC),
                                       cellVol    = cms.double(1.18*100e-4),
                                       cellCap    = cms.double(134),
                                       cceParam   = cms.vdouble(cceParamEpi)),
             'epi120fine': cms.PSet( tag        = cms.string('epi120fine'),
                                     mipEqfC    = cms.double(67*120*qefC),
                                     cellVol    = cms.double(0.52*120e-4),
                                     cellCap    = cms.double(49),
                                     cceParam   = cms.vdouble(cceParamEpi)),
             'epi120coarse': cms.PSet( tag        = cms.string('epi120coarse'),
                                       mipEqfC    = cms.double(67*120*qefC),
                                       cellVol    = cms.double(1.18*120e-4),
                                       cellCap    = cms.double(111),
                                       cceParam   = cms.vdouble(cceParamEpi)),
             'ddfz200fine': cms.PSet( tag        = cms.string('ddfz200fine'),
                                      mipEqfC    = cms.double(70*200*qefC),
                                      cellVol    = cms.double(0.52*200e-4),
                                      cellCap    = cms.double(29),
                                      cceParam   = cms.vdouble(cceParamThin)),
             'ddfz200coarse': cms.PSet( tag        = cms.string('ddfz200coarse'),
                                        mipEqfC    = cms.double(70*200*qefC),
                                        cellVol    = cms.double(1.18*200e-4),
                                        cellCap    = cms.double(65),
                                        cceParam   = cms.vdouble(cceParamThin)),
             'ddfz300fine': cms.PSet( tag        = cms.string('ddfz300fine'),
                                      mipEqfC    = cms.double(73*300*qefC),
                                      cellVol    = cms.double(0.52*300e-4),
                                      cellCap    = cms.double(20),
                                      cceParam   = cms.vdouble(cceParamThick)),
             'ddfz300coarse': cms.PSet( tag        = cms.string('ddfz300coarse'),
                                        mipEqfC    = cms.double(73*300*qefC),
                                        cellVol    = cms.double(1.18*300e-4),
                                        cellCap    = cms.double(45),
                                        cceParam   = cms.vdouble(cceParamThick)),
         }



siat600V=defineSiTypesToScan(cceParamFine_epi600,cceParamThin_tdr600,cceParamThick_tdr600)

from math import sqrt
#analyzer template
process.siop_template = cms.EDAnalyzer("HGCSiOperationScan",
                                       savePadInfo  = cms.bool( options.savePadInfo ),
                                       doseMap      = cms.string( options.doseMap ),
                                       doseMapAlgo  = cms.uint32(0),
                                       ileakParam   = cms.vdouble(),
                                       aimMIPtoADC  = cms.int32(10),                                      
                                       siTypes      = cms.PSet(),                                       
                                       encCommonNoiseSub = cms.double(sqrt(1.0)),
                                       fluenceSF    = cms.double(1.0),
                                   )

#instantiate analyzers for different scenarios
for s in ['epi','ddfz']:
    thicknesses=[80,100,120] if s=='epi' else [200,300]
    for t in thicknesses:
        for d in ['fine','coarse']:
            sitag='{0}{1}{2}'.format(s,t,d)
            print "Init scan for",sitag
            
            setattr(process,
                    "siop_{0}_startup_600V".format(sitag),
                    process.siop_template.clone(doseMapAlgo       = cms.uint32(0),
                                                ileakParam        = cms.vdouble(ileakParam_600V),
                                                siType            = cms.PSet(siat600V[sitag]),
                                                encCommonNoiseSub = cms.double(1.0),
                                                fluenceSF         = cms.double(1.0)))
            
            setattr(process,
                    "siop_{0}_3iab_600V".format(sitag),
                    process.siop_template.clone(doseMapAlgo       = cms.uint32(0),
                                                ileakParam        = cms.vdouble(ileakParam_600V),
                                                siType            = cms.PSet(siat600V[sitag]),
                                                encCommonNoiseSub = cms.double(1.0),
                                                fluenceSF         = cms.double(1.0)))

            setattr(process,
                    "siop_{0}_4iab_600V".format(sitag),
                    process.siop_template.clone(doseMapAlgo       = cms.uint32(0),
                                                ileakParam        = cms.vdouble(ileakParam_600V),
                                                siType            = cms.PSet(siat600V[sitag]),
                                                encCommonNoiseSub = cms.double(1.0),
                                                fluenceSF         = cms.double(4.0/3.0)))
            
#define path according to the scenario to simulate
if options.scenario=='startup_600V':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_600V
        *process.siop_epi80coarse_startup_600V
        *process.siop_epi100fine_startup_600V
        *process.siop_epi100coarse_startup_600V
        *process.siop_epi120fine_startup_600V
        *process.siop_epi120coarse_startup_600V
        *process.siop_ddfz200fine_startup_600V    
        *process.siop_ddfz200coarse_startup_600V
        *process.siop_ddfz300fine_startup_600V
        *process.siop_ddfz300coarse_startup_600V)

if options.scenario=='i3ab_600V':
    process.p=cms.Path( 
        process.siop_epi80fine_i3ab_600V
        *process.siop_epi80coarse_i3ab_600V
        *process.siop_epi100fine_i3ab_600V
        *process.siop_epi100coarse_i3ab_600V
        *process.siop_epi120fine_i3ab_600V
        *process.siop_epi120coarse_i3ab_600V
        *process.siop_ddfz200fine_i3ab_600V    
        *process.siop_ddfz200coarse_i3ab_600V
        *process.siop_ddfz300fine_i3ab_600V
        *process.siop_ddfz300coarse_i3ab_600V)

if options.scenario=='i4ab_600V':
    process.p=cms.Path( 
        process.siop_epi80fine_i4ab_600V
        *process.siop_epi80coarse_i4ab_600V
        *process.siop_epi100fine_i4ab_600V
        *process.siop_epi100coarse_i4ab_600V
        *process.siop_epi120fine_i4ab_600V
        *process.siop_epi120coarse_i4ab_600V
        *process.siop_ddfz200fine_i4ab_600V    
        *process.siop_ddfz200coarse_i4ab_600V
        *process.siop_ddfz300fine_i4ab_600V
        *process.siop_ddfz300coarse_i4ab_600V)

