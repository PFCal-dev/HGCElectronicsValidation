import FWCore.ParameterSet.Config as cms
import os

file_list={'pu':[],'sig':[]}
in_dir='/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/NANOGEN/'
for f in os.listdir(in_dir):
    key='sig'
    if 'MinBias' in f: key='pu'
    file_list[key].append(os.path.join(in_dir,f))
print(file_list)

mixandcluster = cms.PSet( 
    pu = cms.vstring(file_list['pu']),
    sig = cms.vstring(file_list['sig']),
    toaThr = cms.int32(50),
    avgpu = cms.int32(140),
    maxevts = cms.int32(9000),
    jetAlgo = cms.int32(2), #anti-kT=2
    jetR = cms.double(0.4),
    effurl = cms.string('UserCode/HGCElectronicsValidation/bin/timeeff.root')
)
