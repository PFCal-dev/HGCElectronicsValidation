#include <iostream>
#include <algorithm>

#include "UserCode/HGCElectronicsValidation/bin/RunMixAndClusterCommon.h"

#include "fastjet/contrib/SoftKiller.hh"

#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TFile.h"
#include "TH2F.h"
#include "TGraph2D.h"

using namespace std;
using namespace fastjet;

std::map<TString, TH1F *> histos;

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
JetConstituents_t fillJetConstituents(const PseudoJet &j,std::vector<PseudoJetProperties> &properties,TString tag="",bool fillHistos=false) {

  float en(j.e());
  JetConstituents_t jc;
  for(auto c : j.constituents()) {

    float cen(c.e());
    float cenf(cen/en);
    int idx=c.user_index();
    int pfid=properties[idx].pfid;
    bool ispu=properties[idx].ispu;
    
    //compute the time tagging weights

    if(pfid==22) {
      if(!ispu) {
        jc.emf += cenf;
        if(fillHistos) histos[tag+"emf_en"]->Fill(cen);
      }
      else {
        jc.puemf += cenf;
        if(fillHistos) histos[tag+"puemf_en"]->Fill(cen);
      }
    }
    else if(pfid==130)  {
      if(!ispu) {
        jc.nhf+= cenf;
        if(fillHistos) histos[tag+"nhf_en"]->Fill(cen);
      }
      else {
        jc.punhf += cenf;
        if(fillHistos) histos[tag+"punhf_en"]->Fill(cen);
      }
    }
    else if(pfid==13) {
      if(!ispu) jc.muf+= cenf;
      else      jc.pumuf += cenf;
    }
    else if(pfid==211)  {
      if(!ispu) {
        jc.chf+= cenf;
        if(fillHistos) histos[tag+"chf_en"]->Fill(cen);
      }
      else {
        jc.puchf += cenf;
        if(fillHistos) histos[tag+"puchf_en"]->Fill(cen);
      }
    }
    else 
      std::cout << "PFid = " << pfid << " ?? " << std::endl;
  }
  
  return jc;
};


//
int main(int argc, char** argv) {

  //if passed at command line use new cfi
  std::string url("UserCode/HGCElectronicsValidation/bin/mixandcluster_12fC_cfi.py");
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
  int minevts = cfg.getParameter<int>("minevts");
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
  Double_t sk_thr;
  Float_t GenJet_pt,GenJet_eta,GenJet_phi,GenJet_en;
  Float_t NoPuJet_pt,NoPuJet_eta,NoPuJet_phi,NoPuJet_en;
  Float_t Jet_pt,Jet_eta,Jet_phi,Jet_en;
  Float_t PuJet_pt,PuJet_eta,PuJet_phi,PuJet_en;
  //  Float_t PuJet_chf, PuJet_puchf, PuJet_nhf, PuJet_punhf, PuJet_emf, PuJet_puemf;
         
  tree->Branch("nPU",&nPU);
  tree->Branch("sk_thr",&sk_thr);
  tree->Branch("GenJet_pt",&GenJet_pt);
  tree->Branch("GenJet_en",&GenJet_en);
  tree->Branch("GenJet_eta",&GenJet_eta);
  tree->Branch("GenJet_phi",&GenJet_phi);
  tree->Branch("NoPuJet_pt",&NoPuJet_pt);
  tree->Branch("NoPuJet_eta",&NoPuJet_eta);
  tree->Branch("NoPuJet_phi",&NoPuJet_phi);
  tree->Branch("NoPuJet_en",&NoPuJet_en);
  tree->Branch("PuJet_pt",&PuJet_pt);
  tree->Branch("PuJet_eta",&PuJet_eta);
  tree->Branch("PuJet_phi",&PuJet_phi);
  tree->Branch("PuJet_en",&PuJet_en);
  tree->Branch("Jet_pt",&Jet_pt);
  tree->Branch("Jet_eta",&Jet_eta);
  tree->Branch("Jet_phi",&Jet_phi);
  tree->Branch("Jet_en",&Jet_en);
  
  /*
  tree->Branch("toaWgtCat",&toaWgtCat);
  tree->Branch("toaThr",&toaThrApplied);
  tree->Branch("puMode",&puMode);
  tree->Branch("PuJet_chf",&PuJet_chf);
  tree->Branch("PuJet_puchf",&PuJet_puchf);
  tree->Branch("PuJet_nhf",&PuJet_nhf);
  tree->Branch("PuJet_punhf",&PuJet_punhf);
  tree->Branch("PuJet_emf",&PuJet_emf);
  tree->Branch("PuJet_puemf",&PuJet_puemf);
  */
  
  TString parts[]={"chf","puchf","nhf","punhf","emf","puemf"};
  for(size_t i=0; i<sizeof(parts)/sizeof(TString); i++) {
    histos[parts[i]+"_alphaF"] = new TH1F(parts[i]+"_alphaF",";#alpha^{F};Particles",250,0,20);
    histos[parts[i]+"_alphapF"] = new TH1F(parts[i]+"_alphapF",";#alpha^{F}';Particles",250,0,20);
    histos[parts[i]+"_wF"] = new TH1F(parts[i]+"_wF",";PUPPI Weight^{F};Particles",250,0,1);
    histos[parts[i]+"_wpF"] = new TH1F(parts[i]+"_wpF",";PUPPI Weight^{F}';Particles",250,0,1);
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

  minevts = minevts<0 ? 0 : minevts;
  maxevts = maxevts<0 || maxevts>nsigEvts ? nsigEvts : maxevts;
  for(int i=minevts; i<maxevts; i++) {
    
    nPU = UInt_t( CLHEP::RandPoisson::shoot(hre,avgpu) );

    //
    // MC TRUTH GEN JETS
    //
    std::vector<PseudoJetProperties> gen_properties;
    auto sigParticles = getParticlesFrom(sig,i,1,false,gen_properties);
    ClusterSequence sigCS(sigParticles, jet_def);
    Selector sel_vbfjets = SelectorPtMin(30.) * SelectorAbsRapMax(3.0) * SelectorAbsRapMin(1.5);
    const auto sigGenJets = sel_vbfjets(sigCS.inclusive_jets());

    
    //
    // JETS FROM SIMCLUSTERS CROSSING HGCAL BOUNDARY
    //

    //no pileup case
    std::vector<PseudoJetProperties> properties;
    auto sigSimClusters = getParticlesFrom(sig,i,1,false,properties,"SimCluster");
    ClusterSequence sigsimCS(sigSimClusters, jet_def);
    const auto sigJets = sorted_by_pt(sigsimCS.inclusive_jets());
    
    //pileup+signal case
    auto puSimClusters = getParticlesFrom(pu,0,nPU,true,properties,"SimCluster");
    auto sigpuSimClusters(puSimClusters);
    sigpuSimClusters.insert(sigpuSimClusters.end(), sigSimClusters.begin(), sigSimClusters.end());
    ClusterSequence sigpuCS(sigpuSimClusters, jet_def);
    const auto sigpuJets = sorted_by_pt(sigpuCS.inclusive_jets());
        
    //Soft Killer algorithm
    contrib::SoftKiller soft_killer(3.0,0.4);
    vector<PseudoJet> sk_sigpuSimClusters;
    soft_killer.apply(sigpuSimClusters, sk_sigpuSimClusters, sk_thr);
    ClusterSequence sk_sigpuCS(sk_sigpuSimClusters,jet_def);
    const auto sk_sigpuJets = sk_sigpuCS.inclusive_jets();

    //check numbers add up at this point
    assert(sigpuSimClusters.size()==properties.size());
    assert(sigSimClusters.size()+puSimClusters.size()==sigpuSimClusters.size());
    

    //fill tree variables for jets matching the gen jets
    int nSelJets(0);
    for(size_t ij=0; ij<sigGenJets.size(); ij++) {
      
      //select further based on number of constituents (kinematics were applied before with fastjet::selectors)
      if(sigGenJets[ij].constituents().size()<3) continue;
      nSelJets+=1;
      if(nSelJets>1) break;
 
      GenJet_pt=sigGenJets[ij].pt();
      GenJet_phi=sigGenJets[ij].phi();
      GenJet_en=sigGenJets[ij].e();
      GenJet_eta=sigGenJets[ij].eta();

      //matched sim cluster level jet (PU=0)
      auto ijsigjets = min_element(begin(sigJets), end(sigJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sigjetidx = std::distance(begin(sigJets),ijsigjets);
      NoPuJet_pt=sigJets[sigjetidx].pt();
      NoPuJet_eta=sigJets[sigjetidx].eta();
      NoPuJet_phi=sigJets[sigjetidx].phi();
      NoPuJet_en=sigJets[sigjetidx].e();          

      //matched sim cluster level jet (PU!=0)
      auto ijsigpujets = min_element(begin(sigpuJets), end(sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sigpujetidx = std::distance(begin(sigpuJets),ijsigpujets);
      PuJet_pt=sigpuJets[sigpujetidx].pt();
      PuJet_eta=sigpuJets[sigpujetidx].eta();
      PuJet_phi=sigpuJets[sigpujetidx].phi();
      PuJet_en=sigpuJets[sigpujetidx].e();
      //PuJet_chf=sigpujc.chf;
      //    PuJet_puchf=sigpujc.puchf;
      //    PuJet_nhf=sigpujc.nhf;
      //   PuJet_punhf=sigpujc.punhf;
      //    PuJet_emf=sigpujc.emf;
      //    PuJet_puemf=sigpujc.puemf;

      //matched sim cluster after soft killer (PU!=0)
      auto ijsk_sigpujets = min_element(begin(sk_sigpuJets), end(sk_sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sk_sigpujetidx = std::distance(begin(sk_sigpuJets),ijsk_sigpujets);
      Jet_pt=sk_sigpuJets[sk_sigpujetidx].pt();
      Jet_eta=sk_sigpuJets[sk_sigpujetidx].eta();
      Jet_phi=sk_sigpuJets[sk_sigpujetidx].phi();
      Jet_en=sk_sigpuJets[sk_sigpujetidx].e();
      
      std::cout << nSelJets << ") "
                << GenJet_pt << " " << GenJet_eta << " | "
                << NoPuJet_pt << " " << NoPuJet_eta << " | "
                << PuJet_pt << " " << PuJet_eta << " | "
                << Jet_pt << " " << Jet_eta << std::endl;
        
      tree->Fill();          
    }
  
    //debug
    if(i%10==0)  {
      std::cout << "Event: " << i << " pu=" << nPU 
                << " signal/pu/soft kill: " << sigJets.size() << "/" << sigpuJets.size() << "/" << sk_sigpuJets.size() 
                << " soft kill threshold was: " << sk_thr << "GeV" << std::endl;
    }
    
    /*

    

    //purse signal weighted by timetag probability
    std::vector<PseudoJet> sel_sigParticles;
    for(size_t ipart=0; ipart<sigParticles.size(); ipart++) {
      PseudoJet c(sigParticles[ipart]);
      int ipart_idx=c.user_index();
      int pfid=gen_properties[ipart_idx].pfid;
      float en(c.e()),abseta(fabs(c.eta()));

      if(neutonly && pfid!=130 && pfid!=22) continue;
      if(en<1000 && abseta>=1.5 && abseta<=3.0){
        std::pair<int,int> map_key(pfid==22,0);
        c*=eff_map[pfid==22]->Interpolate(en,abseta);
      }
      sel_sigParticles.push_back(c);
    }
    ClusterSequence sigPrunedCS(sel_sigParticles, jet_def);
    const auto sigPrunedGenJets = sorted_by_pt(sigPrunedCS.inclusive_jets());
    

    //pileup+signal under different time-tagging hypothesis
    //puMode=0 all particles
    //puMode=1 only neutrals
    //puMode=2 only neutral hadrons

    //pileup mode loop
    for(puMode=0; puMode<3; puMode++) {
      
      for(Int_t ithr=-1; ithr<2; ithr++) {

        toaWgtCat=(ithr==-1 ? -1: (ithr==0 ? 0 : 1));
        toaThrApplied=(ithr==-1? -1 : toaThr);
        
        std::vector<PseudoJet> sel_sigpuSimClusters;
        for(size_t ipart=0; ipart<sigpuSimClusters.size(); ipart++) {
          
          PseudoJet c(sigpuSimClusters[ipart]);
          int ipart_idx=c.user_index();
          int pfid=properties[ipart_idx].pfid;
          bool ispu=properties[ipart_idx].ispu;
          float en(c.e()),abseta(fabs(c.eta()));
          float timeTagWgt(1.0),resolWgt(1.0);
          if(ithr>=0 && en<1000 && abseta>=1.5 && abseta<=3.0){
            timeTagWgt=eff_map[pfid==22]->Interpolate(en,abseta);
            resolWgt=resol_map[pfid==22]->Interpolate(en,abseta);
          }
          float wgt(1.0);
          if(ithr==0) wgt=timeTagWgt;
          else if(ithr>0) wgt=resolWgt*timeTagWgt;
          if(neutonly && pfid!=130 && pfid!=22) wgt=1.0;

          //we assume that with some magic (e.g. leading hadron, shower time distribution)
          //the reference signal toa is known          
          if(ithr>=0) {
            if(wgt<0.5) wgt=0; //minimal combined efficiency required
            else if(ispu) wgt=1-wgt; //emulate subtraction of pileup
          }
          
          if(puMode==0) c *= wgt; //weight all 
          if(puMode==1) {
            if(pfid==22 || pfid==130) c*=wgt; //weight neutrals
            else if(ispu) c *= 0.; //charged pileup is magically subtracted
          }
          if(puMode==2) {
            if(pfid==130) c*=wgt; //weight neutral hadrons
            else if(ispu) c*=0;      //other pileup is magically subtracted
          }
          sel_sigpuSimClusters.push_back(c);
        }
        
        //run clustering
        ClusterSequence sigpuCS(sel_sigpuSimClusters, jet_def);
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
          JetConstituents_t sigpujc = fillJetConstituents(sigpuJets[puidx],properties,Form("pu%d_",ithr),puMode==0);

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
      }
    }
    */

  }

  //save TTree
  tree->Write();
  for(auto it : histos) it.second->Write();
  fout->Close();

  return 0;
}


/*
attic
//before clustering the jets we comput PUPPI weights
    std::vector<float> puppiWgts(sigpuSimClusters.size(),0),puppiWgtsp(sigpuSimClusters.size(),0);
    for(size_t ij=0; ij<sigpuSimClusters.size(); ij++) {

      //local shape : inclusive and same ID versions
      float alpha=getLocalPuppiShape(ij,sigpuSimClusters,properties,false);
      float alphap=getLocalPuppiShape(ij,sigpuSimClusters,properties,true);

      //compute puppi weights
      size_t ij_idx=sigpuSimClusters[ij].user_index();
      if(ij_idx>properties.size()) {
        std::cout << ij_idx << " " << ij << " " << properties.size() << std::endl;
        continue;
      }
      int pfid=properties[ij_idx].pfid;
      bool ispu=properties[ij_idx].ispu;
      float abseta(sigpuSimClusters[ij].eta());
      puppiWgts[ij]=getPuppiWgt(alpha,pfid,abseta,false);
      puppiWgtsp[ij]=getPuppiWgt(alphap,pfid,abseta,true);

      //fill some control histograms for the weights
      if(abs(pfid)==13) continue;

      TString part("chf");
      if(pfid==130) part="nhf";
      if(pfid==22) part="emf";
      if(ispu) part="pu"+part;
      histos[part+"_alphaF"]->Fill(alpha);
      histos[part+"_wF"]->Fill(puppiWgts[ij]);
      histos[part+"_alphapF"]->Fill(alphap);
      histos[part+"_wpF"]->Fill(puppiWgtsp[ij]);
    }

*/
