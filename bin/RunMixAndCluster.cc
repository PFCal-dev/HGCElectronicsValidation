#include <iostream>
#include <algorithm>

#include "Rivet/Tools/ParticleIdUtils.hh"

#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"

#include "fastjet/ClusterSequence.hh"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "Math/Vector4D.h"
#include "TH2F.h"
#include "TGraph2D.h"

using namespace std;
using namespace fastjet;

std::map<TString, TH1F *> histos;

/**
   @short reads the required entries from the chain and returns a vector of particles for clustering
 */
std::vector<PseudoJet> getParticlesFrom(TChain *t,int ievt,int nevts,bool ispu, TString tag="GenPart") {

  //atach variables to read from
  UInt_t maxPseudoParticles(10000);
  UInt_t nPseudoParticle;
  Int_t PseudoParticle_pdgId[maxPseudoParticles], PseudoParticle_status[maxPseudoParticles];
  Float_t PseudoParticle_eta[maxPseudoParticles], PseudoParticle_mass[maxPseudoParticles], PseudoParticle_phi[maxPseudoParticles], PseudoParticle_pt[maxPseudoParticles];
  Bool_t PseudoParticle_crossedBoundary[maxPseudoParticles];
  t->SetBranchAddress("n"+tag,&nPseudoParticle);
  t->SetBranchAddress(tag+"_pdgId",PseudoParticle_pdgId);
  bool requireStatus(tag.Contains("GenPart"));
  if(requireStatus)
    t->SetBranchAddress(tag+"_status",PseudoParticle_status);
  t->SetBranchAddress(tag+"_pt",PseudoParticle_pt);
  t->SetBranchAddress(tag+"_eta",PseudoParticle_eta);
  t->SetBranchAddress(tag+"_phi",PseudoParticle_phi);
  t->SetBranchAddress(tag+"_mass",PseudoParticle_mass);
  bool requireCrossing(tag.Contains("SimClus"));
  if(requireCrossing)
    t->SetBranchAddress(tag+"_crossedBoundary",PseudoParticle_crossedBoundary);

  std::vector<PseudoJet> pseudoParticles;

  //loop over the required events
  UInt_t ientries(t->GetEntries());
  for(int i=ievt; i<ievt+nevts; i++){

    int ientry(i);
    if(ispu) ientry=gRandom->Integer(ientries-1);
    t->GetEntry(ientry);
    
    //convert kinematics to pseudojets
    for(UInt_t n=0; n<nPseudoParticle; n++) {

      if(requireStatus && PseudoParticle_status[n]!=1) continue;
      if(Rivet::PID::isNeutrino(PseudoParticle_pdgId[n])) continue;
      if(requireCrossing && PseudoParticle_crossedBoundary[n]==false) continue;

      ROOT::Math::PtEtaPhiMVector p4(PseudoParticle_pt[n],PseudoParticle_eta[n],PseudoParticle_phi[n],PseudoParticle_mass[n]);
      auto ip = PseudoJet(p4.px(),p4.py(),p4.pz(),p4.energy());

      //check nature of this particle
      bool isHadron( Rivet::PID::isHadron(PseudoParticle_pdgId[n]) );
      bool isCharged( Rivet::PID::isCharged(PseudoParticle_pdgId[n]) );
      bool isMuon( Rivet::PID::isMuon(PseudoParticle_pdgId[n]) );
      int pfid(0);
      if(isMuon) pfid=13;
      else {
        if(isHadron) {
          if(isCharged) pfid=211;
          else pfid=130;
        } else {
          pfid=22;
        }
      }
      
      //identify as pileup or signal and type of particle (e/g, charged hadron, neutral hadron or muon)
      ip.set_user_index(ispu ? -pfid : pfid);
      pseudoParticles.push_back( ip );
    }
  }
  

  return pseudoParticles;
}


/**
   @fill jet constituents based on the energy fraction
*/
struct JetConstituents_t {
  JetConstituents_t() : 
    nhf(0), chf(0), emf(0), muf(0), 
    punhf(0), puchf(0), puemf(0), pumuf(0)
  {
  };
  float nhf, chf, emf, muf;
  float punhf, puchf, puemf, pumuf;
};
JetConstituents_t fillJetConstituents(const PseudoJet &j,TString tag="",bool fillHistos=false) {

  float en(j.e());
  JetConstituents_t jc;
  for(auto c : j.constituents()) {

    float cen(c.e());
    float cenf(cen/en);
    int pfid=c.user_index();

    //compute the time tagging weights

    if(pfid==22)   {
      jc.emf += cenf;
      if(fillHistos) histos[tag+"emf_en"]->Fill(cen);
    }
    else if(pfid==-22)  {
      jc.puemf += cenf;
      if(fillHistos) histos[tag+"puemf_en"]->Fill(cen);
    }
    else if(pfid==130)  {
      jc.nhf+= cenf;
      if(fillHistos) histos[tag+"nhf_en"]->Fill(cen);
    }
    else if(pfid==-130) {
      jc.punhf += cenf;
      if(fillHistos) histos[tag+"punhf_en"]->Fill(cen);
    }
    else if(pfid==13)   jc.muf+= cenf;
    else if(pfid==-13)  jc.pumuf += cenf;
    else if(pfid==211)  {
      jc.chf+= cenf;
      if(fillHistos) histos[tag+"chf_en"]->Fill(cen);
    }
    else if(pfid==-211) {
      jc.puchf += cenf;
      if(fillHistos) histos[tag+"puchf_en"]->Fill(cen);
    }
    else 
      std::cout << "PFid = " << pfid << " ?? " << std::endl;
  }
  
  return jc;
};


int main(int argc, char** argv) {

  //if passed at command line use new cfi
  std::string url("UserCode/HGCElectronicsValidation/bin/mixandcluster_cfi.py");
  if (argc > 1) url = argv[1];
  url = edm::FileInPath(url).fullPath();

  //get configuration
  const std::shared_ptr<edm::ParameterSet> &pset = edm::readPSetsFrom(url);
  const edm::ParameterSet &cfg = pset->getParameter<edm::ParameterSet>("mixandcluster");
  std::vector<std::string> pufiles = cfg.getParameter<std::vector<std::string> >("pu");
  std::vector<std::string> sigfiles = cfg.getParameter<std::vector<std::string> >("sig");
  bool neutonly  = cfg.exists("neutonly") ? cfg.getParameter<bool>("neutonly") : false;
  Int_t toaThr = cfg.getParameter<int>("toaThr");
  int avgpu = cfg.getParameter<int>("avgpu");
  int maxevts = cfg.getParameter<int>("maxevts");
  std::string effurl=cfg.getParameter<std::string>("effurl");
  effurl = edm::FileInPath(effurl).fullPath();

  //fast jet definition
  int jetAlgo = cfg.getParameter<int>("jetAlgo");
  double jetR = cfg.getParameter<double>("jetR");

  //
  TString foutName=Form("jets_ak%d_vbfhgg_%dfC.root",int(jetR*10),toaThr);
  if(neutonly) foutName=foutName.ReplaceAll(".root","_neutonly.root");
  TFile *fout=TFile::Open(foutName,"RECREATE");
  TTree *tree = new TTree("data","data");
  UInt_t nPU,puMode;
  Int_t toaWgtCat,toaThrApplied;
  Float_t GenJet_pt,GenJet_en, GenJet_eta,GenJet_phi;
  Float_t PrunedGenJet_pt,PrunedGenJet_en, PrunedGenJet_eta,PrunedGenJet_phi;
  Float_t Jet_pt,Jet_eta,Jet_en,Jet_phi;
  Float_t PuJet_en,PuJet_eta,PuJet_phi;
  Float_t PuJet_chf, PuJet_puchf, PuJet_nhf, PuJet_punhf, PuJet_emf, PuJet_puemf;
         
  tree->Branch("nPU",&nPU);
  tree->Branch("GenJet_pt",&GenJet_pt);
  tree->Branch("GenJet_en",&GenJet_en);
  tree->Branch("GenJet_eta",&GenJet_eta);
  tree->Branch("GenJet_phi",&GenJet_phi);
  tree->Branch("PrunedGenJet_pt",&PrunedGenJet_pt);
  tree->Branch("PrunedGenJet_en",&PrunedGenJet_en);
  tree->Branch("PrunedGenJet_eta",&PrunedGenJet_eta);
  tree->Branch("PrunedGenJet_phi",&PrunedGenJet_phi);
  tree->Branch("Jet_pt",&Jet_pt);
  tree->Branch("Jet_eta",&Jet_eta);
  tree->Branch("Jet_phi",&Jet_phi);
  tree->Branch("Jet_en",&Jet_en);
  tree->Branch("PuJet_en",&PuJet_en);
  tree->Branch("PuJet_eta",&PuJet_eta);
  tree->Branch("PuJet_phi",&PuJet_phi);
  tree->Branch("toaWgtCat",&toaWgtCat);
  tree->Branch("toaThr",&toaThrApplied);
  tree->Branch("puMode",&puMode);
  tree->Branch("PuJet_chf",&PuJet_chf);
  tree->Branch("PuJet_puchf",&PuJet_puchf);
  tree->Branch("PuJet_nhf",&PuJet_nhf);
  tree->Branch("PuJet_punhf",&PuJet_punhf);
  tree->Branch("PuJet_emf",&PuJet_emf);
  tree->Branch("PuJet_puemf",&PuJet_puemf);
  
  TString parts[]={"chf","puchf","nhf","punhf","emf","puemf"};
  for(size_t i=0; i<sizeof(parts)/sizeof(TString); i++) {
    for(Int_t ithr=-1; ithr<2; ithr++) {
      TString tag(Form("pu%d_",ithr));
      float maxx(50);
      if(parts[i].Contains("emf")) maxx=25;
      histos[tag+parts[i]+"_en"] = new TH1F(tag+parts[i]+"_en",";Energy [GeV]; Particles",100,0,maxx);
    }
  }

  //read time tag and resolution efficiencies
  TFile *fIn=TFile::Open(effurl.c_str(),"READ");
  std::vector<TGraph2D *> eff_map(2),resol_map(2);
  for(size_t i=0; i<2; i++){
    TString key(i==0 ? "had" : "em");
    eff_map[i]=(TGraph2D *)fIn->Get(Form("%s_%d_timetagwgt",key.Data(),toaThr))->Clone();
    eff_map[i]->SetDirectory(fout);
    resol_map[i]=(TGraph2D *)fIn->Get(Form("%s_%d_timeresolwgt",key.Data(),toaThr))->Clone();
    resol_map[i]->SetDirectory(fout);
  }
  fIn->Close();
  fout->cd();
  
  //fastjet definition
  JetDefinition jet_def( (JetAlgorithm)(jetAlgo), jetR);

  //start chains
  TChain *pu=new TChain("Events");
  for(size_t i=0; i<pufiles.size(); i++) pu->AddFile(pufiles[i].c_str());
  int npuEvts(pu->GetEntries());

  TChain *sig=new TChain("Events");
  for(size_t i=0; i<sigfiles.size(); i++) sig->AddFile(sigfiles[i].c_str());
  int nsigEvts(sig->GetEntries());

  std::cout << "Generating an average pileup of " << avgpu 
            << " from " << pufiles.size() << " pileup files with " << npuEvts << " events,"
            << " and signal injected from " << sigfiles.size() << " files with " << nsigEvts << " events." << std::endl;

  //random number generator
  CLHEP::HepJamesRandom *hre = new CLHEP::HepJamesRandom();
  hre->setSeed(0);
  
  maxevts = maxevts<0 || maxevts>nsigEvts ? nsigEvts : maxevts;
  for(int i=0; i<maxevts; i++) {

    nPU = UInt_t( CLHEP::RandPoisson::shoot(hre,avgpu) );

    //pure signal from gen particles
    auto sigParticles = getParticlesFrom(sig,i,1,false);
    ClusterSequence sigCS(sigParticles, jet_def);
    const auto sigGenJets = sorted_by_pt(sigCS.inclusive_jets());

    //purse signal weighted by timetag probability
    std::vector<PseudoJet> sel_sigParticles;
    for(size_t ipart=0; ipart<sigParticles.size(); ipart++) {
      PseudoJet c(sigParticles[ipart]);
      int pfid=c.user_index();
      float en(c.e()),abseta(fabs(c.eta()));

      if(neutonly && abs(pfid)!=130 && abs(pfid)!=22) continue;
      if(en<1000 && abseta>=1.5 && abseta<=3.0){
        std::pair<int,int> map_key(abs(pfid)==22,0);
        c*=eff_map[abs(pfid)==22]->Interpolate(en,abseta);
      }
      sel_sigParticles.push_back(c);
    }
    ClusterSequence sigPrunedCS(sel_sigParticles, jet_def);
    const auto sigPrunedGenJets = sorted_by_pt(sigPrunedCS.inclusive_jets());
    
    //pure signal from sim clusters crossing HCAL boundary
    auto sigSimClusters = getParticlesFrom(sig,i,1,false,"SimCluster");
    ClusterSequence sigsimCS(sigSimClusters, jet_def);
    const auto sigJets = sorted_by_pt(sigsimCS.inclusive_jets());
    
    //pileup
    auto puParticles = getParticlesFrom(pu,0,nPU,true,"SimCluster");

    //pileup+signal under different time-tagging hypothesis
    //puMode=0 all particles
    //puMode=1 only neutrals
    //puMode=2 only neutral hadrons
    auto sigpuParticles(puParticles);
    sigpuParticles.insert(sigpuParticles.end(), sigParticles.begin(), sigParticles.end());
    for(puMode=0; puMode<3; puMode++) {
      
      for(Int_t ithr=-1; ithr<2; ithr++) {

        toaWgtCat=(ithr==-1 ? -1: (ithr==0 ? 0 : 1));
        toaThrApplied=(ithr==-1? -1 : toaThr);
        
        std::vector<PseudoJet> sel_sigpuParticles;
        for(size_t ipart=0; ipart<sigpuParticles.size(); ipart++) {
          
          PseudoJet c(sigpuParticles[ipart]);
          int pfid=c.user_index();
          float en(c.e()),abseta(fabs(c.eta()));
          float timeTagWgt(1.0),resolWgt(1.0);
          if(ithr>=0 && en<1000 && abseta>=1.5 && abseta<=3.0){
            timeTagWgt=eff_map[abs(pfid)==22]->Interpolate(en,abseta);
            resolWgt=resol_map[abs(pfid)==22]->Interpolate(en,abseta);
          }
          float wgt(1.0);
          if(ithr==0) wgt=timeTagWgt;
          else if(ithr>0) wgt=resolWgt*timeTagWgt;
          if(neutonly && abs(pfid)!=130 && abs(pfid)!=22) wgt=1.0;

          //we assume that with some magic (e.g. leading hadron, shower time distribution)
          //the reference signal toa is known          
          if(ithr>=0) {
            if(wgt<0.5) wgt=0; //minimal combined efficiency required
            else if(pfid<0) wgt=1-wgt; //emulate subtraction of pileup
          }
          
          if(puMode==0) c *= wgt; //weight all 
          if(puMode==1) {
            if(abs(pfid)==22 || abs(pfid)==130) c*=wgt; //weight neutrals
            else if(pfid<0) c *= 0.; //charged pileup is magically subtracted
          }
          if(puMode==2) {
            if(abs(pfid)==130) c*=wgt; //weight neutral hadrons
            else if(pfid<0) c*=0;      //other pileup is magically subtracted
          }
          sel_sigpuParticles.push_back(c);
        }
        
        //run clustering
        ClusterSequence sigpuCS(sel_sigpuParticles, jet_def);
        const auto sigpuJets = sorted_by_pt(sigpuCS.inclusive_jets());
        
        //select up to two signal jets and pair with the closest pileup jet
        int nSelJets(0);
        for(size_t ij=0; ij<sigGenJets.size(); ij++) {
          
          //select
          if(sigGenJets[ij].pt()<30) continue;
          if(fabs(sigGenJets[ij].eta())>3) continue;
          if(fabs(sigGenJets[ij].eta())<1.5) continue;
          if(sigGenJets[ij].constituents().size()<3) continue;
          nSelJets+=1;
          if(nSelJets>1) break;
          
          //find closest pruned jet
          auto ijpruned = min_element(begin(sigPrunedGenJets), end(sigPrunedGenJets), [=] (PseudoJet x, PseudoJet y)
          {
            return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
          });
          auto prunedidx = std::distance(begin(sigPrunedGenJets),ijpruned);
                                         
          //find sim cluster jet closest in eta-phi space
          auto ijsc = min_element(begin(sigJets), end(sigJets), [=] (PseudoJet x, PseudoJet y)
          {
            return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
          });
          auto scidx = std::distance(begin(sigJets),ijsc);

          //find pileup closest in eta-phi space
          auto ijpu = min_element(begin(sigpuJets), end(sigpuJets), [=] (PseudoJet x, PseudoJet y)
          {
            return sigJets[scidx].plain_distance(x) < sigJets[scidx].plain_distance(y);
          });
          auto puidx = std::distance(begin(sigpuJets),ijpu);
          JetConstituents_t sigpujc = fillJetConstituents(sigpuJets[puidx],Form("pu%d_",ithr),puMode==0);

          //fill tree variables
          GenJet_pt=sigGenJets[ij].pt();
          GenJet_phi=sigGenJets[ij].phi();
          GenJet_en=sigGenJets[ij].e();
          GenJet_eta=sigGenJets[ij].eta();
          PrunedGenJet_pt=sigPrunedGenJets[prunedidx].pt();
          PrunedGenJet_en=sigPrunedGenJets[prunedidx].e();
          PrunedGenJet_eta=sigPrunedGenJets[prunedidx].eta();
          PrunedGenJet_phi=sigPrunedGenJets[prunedidx].phi();
          Jet_pt=sigJets[scidx].pt();
          Jet_eta=sigJets[scidx].eta();
          Jet_phi=sigJets[scidx].phi();
          Jet_en=sigJets[scidx].e();          
          PuJet_eta=sigpuJets[puidx].eta();
          PuJet_phi=sigpuJets[puidx].phi();
          PuJet_en=sigpuJets[puidx].e();
          PuJet_chf=sigpujc.chf;
          PuJet_puchf=sigpujc.puchf;
          PuJet_nhf=sigpujc.nhf;
          PuJet_punhf=sigpujc.punhf;
          PuJet_emf=sigpujc.emf;
          PuJet_puemf=sigpujc.puemf;
          tree->Fill();          
        }

        //debug
        if(i%10==0 && puMode==0 && ithr==0) 
          {
            std::cout << "Event: " << i << " pu=" << nPU 
                      << " signal/pu particles: " << sigParticles.size() << "/" << puParticles.size()
                      << " signal/signal+pu jets: " << sigJets.size() << "/" << sigpuJets.size() << std::endl;
          }
        
      }
    }

  }

  //save TTree
  tree->Write();
  for(auto it : histos) it.second->Write();
  fout->Close();

  return 0;
}
