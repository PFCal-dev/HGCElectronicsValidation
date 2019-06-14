#! /usr/bin/env python

import ROOT
from DataFormats.FWLite import Events, Handle
import os,sys
        
indir=sys.argv[1]
outname=sys.argv[2]
files = [os.path.join(indir,f) for f in os.listdir(indir) if '.root' in f]
print len(files),'files to analyze'

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

# use Varparsing object
events = Events(files)

# create handle outside of loop
handles = { 'ak4' : ('ak4GenJetsNoNu', Handle("std::vector<reco::GenJet>")),
            'ak8' : ('ak8GenJetsNoNu', Handle("std::vector<reco::GenJet>")) }

jetProfiles={}
for k in handles:
    jetProfiles[k]={ 
        'pt' : ROOT.TH1F('pt_'+k,';Transverse momentum [GeV]; Jets',50,20,250),
        'ee' : ROOT.TProfile('ee_'+k,';Transverse momentum [GeV]; #sigma_{#eta#eta};',50,0,250),
        'pp' : ROOT.TProfile('pp_'+k,';Transverse momentum [GeV]; #sigma_{#phi#phi};',50,0,250),
        'emf': ROOT.TProfile('emf_'+k,';Transverse momentum [GeV];Fraction',50,20,250),
        'hmf': ROOT.TProfile('hmf_'+k,';Transverse momentum [GeV];Fraction',50,20,250),
        'n'  : ROOT.TProfile('n_'+k,';Transverse momentum [GeV];<Constituents>',50,20,250),
        'eta': ROOT.TH1F('eta_'+k,';Pseudo-rapidity; Jets',50,1.5,3.0),
        'eeeta' : ROOT.TProfile('eeeta_'+k,';Pseudo-rapidity; #sigma_{#eta#eta};',50,1.5,3),
        'ppeta' : ROOT.TProfile('ppeta_'+k,';Pseudo-rapidity; #sigma_{#phi#phi};',50,1.5,3),
        'emfeta': ROOT.TProfile('emfeta_'+k,';Pseudo-rapidity;Fraction',50,1.5,3),
        'hmfeta': ROOT.TProfile('hmfeta_'+k,';Pseudo-rapidity;Fraction',50,1.5,3),
        'neta'  : ROOT.TProfile('neta_'+k,';Pseudo-rapidity;<Constituents>',50,1.5,3)
    }
    for hname in jetProfiles[k]:
        jetProfiles[k][hname].SetDirectory(0)
        jetProfiles[k][hname].Sumw2(0)
                     


#loop over event
nini=0
npass={'ak4':0,'ak8':0}
for event in events:
    
    nini=nini+1

    for k in handles:
        label,handle=handles[k]
        event.getByLabel(label, handle)
        jets = handle.product()

        ptThr=30 if k=='ak4' else 100
        minJets=2 if k=='ak4' else 1

        jetInfo=[]
        for j in jets:
            abseta=abs(j.eta())
            if abseta<1.5 or abseta>3.0 :
                continue
            pt=j.pt()
            if pt<20 : 
                continue

            en=j.energy()
            emf=j.emEnergy()/en
            hmf=j.hadEnergy()/en
            ee=j.etaetaMoment()
            pp=j.phiphiMoment()
            nconst=j.nConstituents()
            jetInfo.append((pt,abseta,emf,hmf,ee,pp,nconst))
 
        #filter event
        if len(jetInfo)<minJets : continue

        #fill histos
        npass[k]=npass[k]+1
        for pt,abseta,emf,hmf,ee,pp,nconst in jetInfo:
            jetProfiles[k]['pt'].Fill(pt)
            jetProfiles[k]['eta'].Fill(abseta)
            jetProfiles[k]['emf'].Fill(pt,emf)
            jetProfiles[k]['emfeta'].Fill(abseta,emf)
            jetProfiles[k]['hmf'].Fill(pt,hmf)
            jetProfiles[k]['hmfeta'].Fill(abseta,hmf)
            jetProfiles[k]['n'].Fill(pt,nconst)
            jetProfiles[k]['neta'].Fill(abseta,nconst)
            jetProfiles[k]['ee'].Fill(pt,ee)
            jetProfiles[k]['eeeta'].Fill(abseta,ee)
            jetProfiles[k]['pp'].Fill(pt,ee)
            jetProfiles[k]['ppeta'].Fill(abseta,ee)




print 'Analyzed',nini,'events'
print 'Endcap jet efficiency'
for k in npass:
    print '\t',k,float(npass[k])/float(nini)
fOut=ROOT.TFile.Open(outname,'RECREATE')
for k in jetProfiles:
  for hname in jetProfiles[k]:
        jetProfiles[k][hname].SetDirectory(fOut)
        jetProfiles[k][hname].Write()                     
fOut.Close()
