import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('geometry', 'Extended2026D86', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'geometry to use')
options.register("doseMap", "SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register("scenario", "3iab_120-600V_200-300-600V",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register("uvmapfile", "UserCode/HGCElectronicsValidation/data/geomCMSSW10052021_corrected.txt",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register("savePadInfo", False,  VarParsing.multiplicity.singleton, VarParsing.varType.bool)
options.register("wafersFromCMSSW", False,  VarParsing.multiplicity.singleton, VarParsing.varType.bool)
options.parseArguments()

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('demo',Phase2C11I13M9)

#set geometry/global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.Geometry%sReco_cff'%options.geometry)
#process.load('Configuration.Geometry.Geometry%s_cff'%options.geometry)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#source/number events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

#tfileservice
process.TFileService = cms.Service("TFileService",fileName = cms.string(options.output))

#import standard Ileak and CCE parameterizations
from SimCalorimetry.HGCalSimAlgos.hgcSensorOpParams_cfi import hgcSiSensorIleak,hgcSiSensorCCE 

#define the Si types to scan
def defineSiTypesToScan(conditions='TDR_600V', qefC=1.60217646e-4):

    """returns a dict of interesting sensors and their main properties"""
    
    return { 'epi80fine': cms.PSet( tag        = cms.string('epi80fine'),
                                    mipEqfC    = cms.double(65*80*qefC),
                                    cellVol    = cms.double(0.56*80e-4),
                                    cellCap    = cms.double(74),
                                    cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'epi80coarse': cms.PSet( tag        = cms.string('epi80coarse'),
                                      mipEqfC    = cms.double(65*80*qefC),
                                      cellVol    = cms.double(1.26*80e-4),
                                      cellCap    = cms.double(167),
                                      cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'epi100fine': cms.PSet( tag        = cms.string('epi100fine'),
                                     mipEqfC    = cms.double(66*100*qefC),
                                     cellVol    = cms.double(0.56*100e-4),
                                     cellCap    = cms.double(59),
                                     cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'epi100coarse': cms.PSet( tag        = cms.string('epi100coarse'),
                                       mipEqfC    = cms.double(66*100*qefC),
                                       cellVol    = cms.double(1.26*100e-4),
                                       cellCap    = cms.double(134),
                                       cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'epi120fine': cms.PSet( tag        = cms.string('epi120fine'),
                                     mipEqfC    = cms.double(67*120*qefC),
                                     cellVol    = cms.double(0.56*120e-4),
                                     cellCap    = cms.double(49),
                                     cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'epi120coarse': cms.PSet( tag        = cms.string('epi120coarse'),
                                       mipEqfC    = cms.double(67*120*qefC),
                                       cellVol    = cms.double(1.26*120e-4),
                                       cellCap    = cms.double(111),
                                       cceParam   = cms.vdouble(hgcSiSensorCCE(120,conditions))),
             'ddfz200fine': cms.PSet( tag        = cms.string('ddfz200fine'),
                                      mipEqfC    = cms.double(70*200*qefC),
                                      cellVol    = cms.double(0.56*200e-4),
                                      cellCap    = cms.double(29),
                                      cceParam   = cms.vdouble(hgcSiSensorCCE(200,conditions))),
             'ddfz200coarse': cms.PSet( tag        = cms.string('ddfz200coarse'),
                                        mipEqfC    = cms.double(70*200*qefC),
                                        cellVol    = cms.double(1.26*200e-4),
                                        cellCap    = cms.double(65),
                                        cceParam   = cms.vdouble(hgcSiSensorCCE(200,conditions))),
             'ddfz300fine': cms.PSet( tag        = cms.string('ddfz300fine'),
                                      mipEqfC    = cms.double(73*300*qefC),
                                      cellVol    = cms.double(0.56*300e-4),
                                      cellCap    = cms.double(20),
                                      cceParam   = cms.vdouble(hgcSiSensorCCE(300,conditions))),
             'ddfz300coarse': cms.PSet( tag        = cms.string('ddfz300coarse'),
                                        mipEqfC    = cms.double(73*300*qefC),
                                        cellVol    = cms.double(1.26*300e-4),
                                        cellCap    = cms.double(45),
                                        cceParam   = cms.vdouble(hgcSiSensorCCE(300,conditions))),
        }


conditions_list=['TDR_600V','TDR_800V',
                 'CERN21_600V_10m','CERN21_600V_30m','CERN21_600V_90m','CERN21_600V_120m',
                 'CERN21_800V_10m','CERN21_800V_30m','CERN21_800V_90m','CERN21_800V_120m']
siconds={}
siconds=dict( [(conditions,defineSiTypesToScan(conditions)) for conditions in conditions_list] )

from math import sqrt
#analyzer template
process.siop_template = cms.EDAnalyzer("HGCSiOperationScan",
                                       uvmapfile    = cms.string(options.uvmapfile),
                                       setPreassignedWafersFromCMSSW = cms.bool(options.wafersFromCMSSW),
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
import itertools
for s in ['epi','ddfz']:
    thicknesses=[80,100,120] if s=='epi' else [200,300]

    for t,d,conds in itertools.product(thicknesses,['fine','coarse'],conditions_list):

        sitag='{0}{1}{2}'.format(s,t,d)
        print('Init scan for',sitag,'@',conds)
             
        setattr(process,
                "siop_{0}_startup_{1}".format(sitag,conds),
                process.siop_template.clone(doseMapAlgo       = cms.uint32(3),
                                            ileakParam        = cms.vdouble(hgcSiSensorIleak(conds)),
                                            siType            = cms.PSet(siconds[conds][sitag]),
                                            encCommonNoiseSub = cms.double(1.0),
                                            fluenceSF         = cms.double(1.0)))
             
        setattr(process,
                "siop_{0}_3iab_{1}".format(sitag,conds),
                process.siop_template.clone(doseMapAlgo       = cms.uint32(0),
                                            ileakParam        = cms.vdouble(hgcSiSensorIleak(conds)),
                                            siType            = cms.PSet(siconds[conds][sitag]),
                                            encCommonNoiseSub = cms.double(1.0),
                                            fluenceSF         = cms.double(1.0)))

        setattr(process,
                "siop_{0}_4iab_{1}".format(sitag,conds),
                process.siop_template.clone(doseMapAlgo       = cms.uint32(0),
                                            ileakParam        = cms.vdouble(hgcSiSensorIleak(conds)),
                                            siType            = cms.PSet(siconds[conds][sitag]),
                                            encCommonNoiseSub = cms.double(1.0),
                                            fluenceSF         = cms.double(4.0/3.0)))            



            
#define path according to the scenario to simulate
#
#TDR
#
if options.scenario=='startup_600V_TDR':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_TDR_600V
        *process.siop_epi100fine_startup_TDR_600V
        *process.siop_epi120fine_startup_TDR_600V
        *process.siop_epi120coarse_startup_TDR_600V
        *process.siop_ddfz200fine_startup_TDR_600V    
        *process.siop_ddfz200coarse_startup_TDR_600V
        *process.siop_ddfz300fine_startup_TDR_600V
        *process.siop_ddfz300coarse_startup_TDR_600V)

if options.scenario=='3iab_120-600V_200-300-600V_TDR':
    process.p=cms.Path( 
        process.siop_epi80fine_3iab_TDR_600V
        *process.siop_epi100fine_3iab_TDR_600V
        *process.siop_epi120fine_3iab_TDR_600V
        *process.siop_epi120coarse_3iab_TDR_600V
        *process.siop_ddfz200fine_3iab_TDR_600V    
        *process.siop_ddfz200coarse_3iab_TDR_600V
        *process.siop_ddfz300fine_3iab_TDR_600V
        *process.siop_ddfz300coarse_3iab_TDR_600V
    )

if options.scenario=='4iab_120-600V_200-300-800V_TDR':
    process.p=cms.Path( 
        process.siop_epi80fine_4iab_TDR_600V
        *process.siop_epi100fine_4iab_TDR_600V
        *process.siop_epi120fine_4iab_TDR_600V
        *process.siop_epi120coarse_4iab_TDR_600V
        *process.siop_ddfz200fine_4iab_TDR_800V    
        *process.siop_ddfz200coarse_4iab_TDR_800V
        *process.siop_ddfz300fine_4iab_TDR_800V
        *process.siop_ddfz300coarse_4iab_TDR_800V)

#
# CERN21 with 10min annealing
#
if options.scenario=='startup_600V_CERN21_10m':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_CERN21_600V_10m
        *process.siop_epi100fine_startup_CERN21_600V_10m
        *process.siop_epi120fine_startup_CERN21_600V_10m
        *process.siop_epi120coarse_startup_CERN21_600V_10m
        *process.siop_ddfz200fine_startup_CERN21_600V_10m    
        *process.siop_ddfz200coarse_startup_CERN21_600V_10m
        *process.siop_ddfz300fine_startup_CERN21_600V_10m
        *process.siop_ddfz300coarse_startup_CERN21_600V_10m)

if options.scenario=='3iab_120-600V_200-300-600V_CERN21_10m':
    process.p=cms.Path( 
        process.siop_epi80fine_3iab_CERN21_600V_10m
        *process.siop_epi100fine_3iab_CERN21_600V_10m
        *process.siop_epi120fine_3iab_CERN21_600V_10m
        *process.siop_epi120coarse_3iab_CERN21_600V_10m
        *process.siop_ddfz200fine_3iab_CERN21_600V_10m    
        *process.siop_ddfz200coarse_3iab_CERN21_600V_10m
        *process.siop_ddfz300fine_3iab_CERN21_600V_10m
        *process.siop_ddfz300coarse_3iab_CERN21_600V_10m
    )

if options.scenario=='4iab_120-600V_200-300-800V_CERN21_10m':
    process.p=cms.Path( 
        process.siop_epi80fine_4iab_CERN21_600V_10m
        *process.siop_epi100fine_4iab_CERN21_600V_10m
        *process.siop_epi120fine_4iab_CERN21_600V_10m
        *process.siop_epi120coarse_4iab_CERN21_600V_10m
        *process.siop_ddfz200fine_4iab_CERN21_800V_10m
        *process.siop_ddfz200coarse_4iab_CERN21_800V_10m
        *process.siop_ddfz300fine_4iab_CERN21_800V_10m
        *process.siop_ddfz300coarse_4iab_CERN21_800V_10m)

#
# CERN21 with 30min annealing
#
if options.scenario=='startup_600V_CERN21_30m':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_CERN21_600V_30m
        *process.siop_epi100fine_startup_CERN21_600V_30m
        *process.siop_epi120fine_startup_CERN21_600V_30m
        *process.siop_epi120coarse_startup_CERN21_600V_30m
        *process.siop_ddfz200fine_startup_CERN21_600V_30m    
        *process.siop_ddfz200coarse_startup_CERN21_600V_30m
        *process.siop_ddfz300fine_startup_CERN21_600V_30m
        *process.siop_ddfz300coarse_startup_CERN21_600V_30m)

if options.scenario=='3iab_120-600V_200-300-600V_CERN21_30m':
    process.p=cms.Path( 
        process.siop_epi80fine_3iab_CERN21_600V_30m
        *process.siop_epi100fine_3iab_CERN21_600V_30m
        *process.siop_epi120fine_3iab_CERN21_600V_30m
        *process.siop_epi120coarse_3iab_CERN21_600V_30m
        *process.siop_ddfz200fine_3iab_CERN21_600V_30m    
        *process.siop_ddfz200coarse_3iab_CERN21_600V_30m
        *process.siop_ddfz300fine_3iab_CERN21_600V_30m
        *process.siop_ddfz300coarse_3iab_CERN21_600V_30m
    )

if options.scenario=='4iab_120-600V_200-300-800V_CERN21_30m':
    process.p=cms.Path( 
        process.siop_epi80fine_4iab_CERN21_600V_30m
        *process.siop_epi100fine_4iab_CERN21_600V_30m
        *process.siop_epi120fine_4iab_CERN21_600V_30m
        *process.siop_epi120coarse_4iab_CERN21_600V_30m
        *process.siop_ddfz200fine_4iab_CERN21_800V_30m
        *process.siop_ddfz200coarse_4iab_CERN21_800V_30m
        *process.siop_ddfz300fine_4iab_CERN21_800V_30m
        *process.siop_ddfz300coarse_4iab_CERN21_800V_30m)

#
# CERN21 with 90min annealing
#
if options.scenario=='startup_600V_CERN21_90m':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_CERN21_600V_90m
        *process.siop_epi100fine_startup_CERN21_600V_90m
        *process.siop_epi120fine_startup_CERN21_600V_90m
        *process.siop_epi120coarse_startup_CERN21_600V_90m
        *process.siop_ddfz200fine_startup_CERN21_600V_90m    
        *process.siop_ddfz200coarse_startup_CERN21_600V_90m
        *process.siop_ddfz300fine_startup_CERN21_600V_90m
        *process.siop_ddfz300coarse_startup_CERN21_600V_90m)

if options.scenario=='3iab_120-600V_200-300-600V_CERN21_90m':
    process.p=cms.Path( 
        process.siop_epi80fine_3iab_CERN21_600V_90m
        *process.siop_epi100fine_3iab_CERN21_600V_90m
        *process.siop_epi120fine_3iab_CERN21_600V_90m
        *process.siop_epi120coarse_3iab_CERN21_600V_90m
        *process.siop_ddfz200fine_3iab_CERN21_600V_90m    
        *process.siop_ddfz200coarse_3iab_CERN21_600V_90m
        *process.siop_ddfz300fine_3iab_CERN21_600V_90m
        *process.siop_ddfz300coarse_3iab_CERN21_600V_90m
    )

if options.scenario=='4iab_120-600V_200-300-800V_CERN21_90m':
    process.p=cms.Path( 
        process.siop_epi80fine_4iab_CERN21_600V_90m
        *process.siop_epi100fine_4iab_CERN21_600V_90m
        *process.siop_epi120fine_4iab_CERN21_600V_90m
        *process.siop_epi120coarse_4iab_CERN21_600V_90m
        *process.siop_ddfz200fine_4iab_CERN21_800V_90m
        *process.siop_ddfz200coarse_4iab_CERN21_800V_90m
        *process.siop_ddfz300fine_4iab_CERN21_800V_90m
        *process.siop_ddfz300coarse_4iab_CERN21_800V_90m)

#
# CERN21 with 120min annealing
#
if options.scenario=='startup_600V_CERN21_120m':
    process.p=cms.Path( 
        process.siop_epi80fine_startup_CERN21_600V_120m
        *process.siop_epi100fine_startup_CERN21_600V_120m
        *process.siop_epi120fine_startup_CERN21_600V_120m
        *process.siop_epi120coarse_startup_CERN21_600V_120m
        *process.siop_ddfz200fine_startup_CERN21_600V_120m    
        *process.siop_ddfz200coarse_startup_CERN21_600V_120m
        *process.siop_ddfz300fine_startup_CERN21_600V_120m
        *process.siop_ddfz300coarse_startup_CERN21_600V_120m)

if options.scenario=='3iab_120-600V_200-300-600V_CERN21_120m':
    process.p=cms.Path( 
        process.siop_epi80fine_3iab_CERN21_600V_120m
        *process.siop_epi100fine_3iab_CERN21_600V_120m
        *process.siop_epi120fine_3iab_CERN21_600V_120m
        *process.siop_epi120coarse_3iab_CERN21_600V_120m
        *process.siop_ddfz200fine_3iab_CERN21_600V_120m    
        *process.siop_ddfz200coarse_3iab_CERN21_600V_120m
        *process.siop_ddfz300fine_3iab_CERN21_600V_120m
        *process.siop_ddfz300coarse_3iab_CERN21_600V_120m
    )

if options.scenario=='4iab_120-600V_200-300-800V_CERN21_120m':
    process.p=cms.Path( 
        process.siop_epi80fine_4iab_CERN21_600V_120m
        *process.siop_epi100fine_4iab_CERN21_600V_120m
        *process.siop_epi120fine_4iab_CERN21_600V_120m
        *process.siop_epi120coarse_4iab_CERN21_600V_120m
        *process.siop_ddfz200fine_4iab_CERN21_800V_120m
        *process.siop_ddfz200coarse_4iab_CERN21_800V_120m
        *process.siop_ddfz300fine_4iab_CERN21_800V_120m
        *process.siop_ddfz300coarse_4iab_CERN21_800V_120m)
