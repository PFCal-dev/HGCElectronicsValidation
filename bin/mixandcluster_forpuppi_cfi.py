import FWCore.ParameterSet.Config as cms
import os

file_list=[]
in_dir='/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/NANOGEN/'
for f in os.listdir(in_dir):
    file_list.append(os.path.join(in_dir,f))
print(file_list)

mixandcluster = cms.PSet( 
    pu = cms.vstring(file_list),
    avgpu = cms.int32(140),
    maxevts = cms.int32(9000)
)
