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
JetConstituents_t fillJetConstituents(const PseudoJet &j,
                                      std::vector<PseudoJetProperties> &properties,
                                      TString tag="",
                                      bool fillHistos=true) {

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
  TTree *tree = new TTree("Events","Events");  
  UInt_t nPU,nGenJet;
  Double_t sk_thr;
  Float_t GenJet_pt[2],GenJet_eta[2],GenJet_phi[2],GenJet_en[2];
  Float_t NoPuJet_pt[2],NoPuJet_eta[2],NoPuJet_phi[2],NoPuJet_en[2];
  Float_t PuJet_pt[2],PuJet_eta[2],PuJet_phi[2],PuJet_en[2];
  Float_t Jet_pt[2],Jet_eta[2],Jet_phi[2],Jet_en[2];
  Float_t NeutralTimeTagJet_pt[2],NeutralTimeTagJet_eta[2],NeutralTimeTagJet_phi[2],NeutralTimeTagJet_en[2];
  Float_t FullTimeTagJet_pt[2],FullTimeTagJet_eta[2],FullTimeTagJet_phi[2],FullTimeTagJet_en[2];
  
  tree->Branch("nPU",&nPU,"nPU/I");
  tree->Branch("toaThr",&toaThr,"toaThr/I");
  tree->Branch("sk_thr",&sk_thr,"sk_thr/D");
  tree->Branch("nGenJet",&nGenJet,"nGenJet/I");
  tree->Branch("GenJet_pt",GenJet_pt,"GenJet_pt[nGenJet]/F");
  tree->Branch("GenJet_eta",GenJet_eta,"GenJet_eta[nGenJet]/F");
  tree->Branch("GenJet_phi",GenJet_phi,"GenJet_phi[nGenJet]/F");
  tree->Branch("GenJet_en",GenJet_en,"GenJet_en[nGenJet]/F");
  tree->Branch("NoPuJet_pt",NoPuJet_pt,"NoPuJet_pt[nGenJet]/F");
  tree->Branch("NoPuJet_eta",NoPuJet_eta,"NoPuJet_eta[nGenJet]/F");
  tree->Branch("NoPuJet_phi",NoPuJet_phi,"NoPuJet_phi[nGenJet]/F");
  tree->Branch("NoPuJet_en",NoPuJet_en,"NoPuJet_en[nGenJet]/F");
  tree->Branch("PuJet_pt",PuJet_pt,"PuJet_pt[nGenJet]/F");
  tree->Branch("PuJet_eta",PuJet_eta,"PuJet_eta[nGenJet]/F");
  tree->Branch("PuJet_phi",PuJet_phi,"PuJet_phi[nGenJet]/F");
  tree->Branch("PuJet_en",PuJet_en,"PuJet_en[nGenJet]/F");
  tree->Branch("Jet_pt",Jet_pt,"Jet_pt[nGenJet]/F");
  tree->Branch("Jet_eta",Jet_eta,"Jet_eta[nGenJet]/F");
  tree->Branch("Jet_phi",Jet_phi,"Jet_phi[nGenJet]/F");
  tree->Branch("Jet_en",Jet_en,"Jet_en[nGenJet]/F");
  tree->Branch("NeutralTimeTagJet_pt",NeutralTimeTagJet_pt,"NeutralTimeTagJet_pt[nGenJet]/F");
  tree->Branch("NeutralTimeTagJet_eta",NeutralTimeTagJet_eta,"NeutralTimeTagJet_eta[nGenJet]/F");
  tree->Branch("NeutralTimeTagJet_phi",NeutralTimeTagJet_phi,"NeutralTimeTagJet_phi[nGenJet]/F");
  tree->Branch("NeutralTimeTagJet_en",NeutralTimeTagJet_en,"NeutralTimeTagJet_en[nGenJet]/F");
  tree->Branch("FullTimeTagJet_pt",FullTimeTagJet_pt,"FullTimeTagJet_pt[nGenJet]/F");
  tree->Branch("FullTimeTagJet_eta",FullTimeTagJet_eta,"FullTimeTagJet_eta[nGenJet]/F");
  tree->Branch("FullTimeTagJet_phi",FullTimeTagJet_phi,"FullTimeTagJet_phi[nGenJet]/F");
  tree->Branch("FullTimeTagJet_en",FullTimeTagJet_en,"FullTimeTagJet_en[nGenJet]/F");

  //control histograms for the energy of the constituents
  TString parts[]={"chf","puchf","nhf","punhf","emf","puemf"};
  TString algos[]={"nopu","pu","sk","ntt_sk","fulltt_sk"};
  for(size_t i=0; i<sizeof(parts)/sizeof(TString); i++) {
    float xmax(50);
    if(parts[i].Contains("emf")) xmax=25;
    for(size_t j=0; j<sizeof(algos)/sizeof(TString); j++) {
      TString name(algos[j]+parts[i]+"_en");
      histos[name] = new TH1F(name,";Energy [GeV]; Particles",100,0,xmax);
    }
  }

  //read time tag and resolution efficiencies
  TFile *fIn=TFile::Open(effurl.c_str(),"READ");
  std::vector<TGraph2D *> eff_map(2),resol_map(2);
  for(size_t i=0; i<2; i++){
    TString key(i==0 ? "pdg130" : "pdg22");
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
    float t0=properties[sigSimClusters.at(0).user_index()].gvt;
    
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

    std::vector<PseudoJet> ntt_sk_sigpuSimClusters;
    std::vector<PseudoJet> fulltt_sk_sigpuSimClusters;
    for(size_t ipart=0; ipart<sk_sigpuSimClusters.size(); ipart++) {
          
      PseudoJet c(sk_sigpuSimClusters[ipart]);
      int ipart_idx=c.user_index();
      int pfid=properties[ipart_idx].pfid;
      bool ispu=properties[ipart_idx].ispu;
      float p_t0=properties[ipart_idx].gvt;
      bool vtxTimeWithinTimeResol( fabs(p_t0-t0)<0.09 );
        
      float en(c.e()),abseta(fabs(c.eta())),pt(c.pt());
      float timeTagWgt(1.0),resolWgt(1.0);
      if(en<1000 && abseta>=1.5 && abseta<=3.0){
        timeTagWgt=eff_map[pfid==22]->Interpolate(en,abseta);
        resolWgt=resol_map[pfid==22]->Interpolate(en,abseta);
      }

      //std::cout << pfid << " " << ispu << " "
      //          << en << " " 
      //          << vtxTimeWithinTimeResol
      //         << " " << timeTagWgt << " " << resolWgt << std::endl;

      //neutral time tagging
      if(pfid==130 || pfid==22) {
        float ntt_wgt(resolWgt*timeTagWgt);
        if(!vtxTimeWithinTimeResol) ntt_wgt*=0.;
        c*=ntt_wgt;
        ntt_sk_sigpuSimClusters.push_back(c);
        fulltt_sk_sigpuSimClusters.push_back(c);
      }
      //charged time tagging
      else {
        ntt_sk_sigpuSimClusters.push_back(c);
        float chtt_wgt( (pt<0.7 || (pt>0.7 && vtxTimeWithinTimeResol)) ? 1. : 0. );
        c*=chtt_wgt;
        fulltt_sk_sigpuSimClusters.push_back(c);
      }
    }
    ClusterSequence ntt_sk_sigpuCS(ntt_sk_sigpuSimClusters,jet_def);
    const auto ntt_sk_sigpuJets = ntt_sk_sigpuCS.inclusive_jets();
    ClusterSequence fulltt_sk_sigpuCS(fulltt_sk_sigpuSimClusters,jet_def);
    const auto fulltt_sk_sigpuJets = fulltt_sk_sigpuCS.inclusive_jets();       

    //fill tree variables for jets matching the gen jets
    nGenJet=0;
    for(size_t ij=0; ij<sigGenJets.size(); ij++) {
      
      //select further based on number of constituents (kinematics were applied before with fastjet::selectors)
      if(sigGenJets[ij].constituents().size()<3) continue;
      if(nGenJet>1) break;
 
      GenJet_pt[nGenJet]=sigGenJets[ij].pt();
      GenJet_phi[nGenJet]=sigGenJets[ij].phi();
      GenJet_en[nGenJet]=sigGenJets[ij].e();
      GenJet_eta[nGenJet]=sigGenJets[ij].eta();

      //matched sim cluster level jet (PU=0)
      auto ijsigjets = min_element(begin(sigJets), end(sigJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sigjetidx = std::distance(begin(sigJets),ijsigjets);
      NoPuJet_pt[nGenJet]=sigJets[sigjetidx].pt();
      NoPuJet_eta[nGenJet]=sigJets[sigjetidx].eta();
      NoPuJet_phi[nGenJet]=sigJets[sigjetidx].phi();
      NoPuJet_en[nGenJet]=sigJets[sigjetidx].e();          
      fillJetConstituents(sigJets[sigjetidx],properties,"nopu");
          
      //matched sim cluster level jet (PU!=0)
      auto ijsigpujets = min_element(begin(sigpuJets), end(sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sigpujetidx = std::distance(begin(sigpuJets),ijsigpujets);
      PuJet_pt[nGenJet]=sigpuJets[sigpujetidx].pt();
      PuJet_eta[nGenJet]=sigpuJets[sigpujetidx].eta();
      PuJet_phi[nGenJet]=sigpuJets[sigpujetidx].phi();
      PuJet_en[nGenJet]=sigpuJets[sigpujetidx].e();
      fillJetConstituents(sigpuJets[sigpujetidx],properties,"pu");

      //matched sim cluster after soft killer (PU!=0)
      auto ijsk_sigpujets = min_element(begin(sk_sigpuJets), end(sk_sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto sk_sigpujetidx = std::distance(begin(sk_sigpuJets),ijsk_sigpujets);
      Jet_pt[nGenJet]=sk_sigpuJets[sk_sigpujetidx].pt();
      Jet_eta[nGenJet]=sk_sigpuJets[sk_sigpujetidx].eta();
      Jet_phi[nGenJet]=sk_sigpuJets[sk_sigpujetidx].phi();
      Jet_en[nGenJet]=sk_sigpuJets[sk_sigpujetidx].e();
      fillJetConstituents(sk_sigpuJets[sk_sigpujetidx],properties,"sk");

      //matched sim cluster after soft killer + neutrals time tag (PU!=0)
      auto ijntt_sk_sigpujets = min_element(begin(ntt_sk_sigpuJets), end(ntt_sk_sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto ntt_sk_sigpujetidx = std::distance(begin(ntt_sk_sigpuJets),ijntt_sk_sigpujets);
      NeutralTimeTagJet_pt[nGenJet]=ntt_sk_sigpuJets[ntt_sk_sigpujetidx].pt();
      NeutralTimeTagJet_eta[nGenJet]=ntt_sk_sigpuJets[ntt_sk_sigpujetidx].eta();
      NeutralTimeTagJet_phi[nGenJet]=ntt_sk_sigpuJets[ntt_sk_sigpujetidx].phi();
      NeutralTimeTagJet_en[nGenJet]=ntt_sk_sigpuJets[ntt_sk_sigpujetidx].e();
      fillJetConstituents(ntt_sk_sigpuJets[ntt_sk_sigpujetidx],properties,"ntt_sk");

      //matched sim cluster after soft killer + full time tag (PU!=0)
      auto ijfulltt_sk_sigpujets = min_element(begin(fulltt_sk_sigpuJets), end(fulltt_sk_sigpuJets), [=] (PseudoJet x, PseudoJet y)
      {
        return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
      });
      auto fulltt_sk_sigpujetidx = std::distance(begin(fulltt_sk_sigpuJets),ijfulltt_sk_sigpujets);
      FullTimeTagJet_pt[nGenJet]=fulltt_sk_sigpuJets[fulltt_sk_sigpujetidx].pt();
      FullTimeTagJet_eta[nGenJet]=fulltt_sk_sigpuJets[fulltt_sk_sigpujetidx].eta();
      FullTimeTagJet_phi[nGenJet]=fulltt_sk_sigpuJets[fulltt_sk_sigpujetidx].phi();
      FullTimeTagJet_en[nGenJet]=fulltt_sk_sigpuJets[fulltt_sk_sigpujetidx].e();
      fillJetConstituents(fulltt_sk_sigpuJets[fulltt_sk_sigpujetidx],properties,"fulltt_sk");

      std::cout << nGenJet+1 << ") "
                << GenJet_pt[nGenJet] << " " << GenJet_eta[nGenJet] << " | "
                << NoPuJet_pt[nGenJet] << " " << NoPuJet_eta[nGenJet] << " | "
                << PuJet_pt[nGenJet] << " " << PuJet_eta[nGenJet] << " | "
                << Jet_pt[nGenJet] << " " << Jet_eta[nGenJet] << " | "
                << NeutralTimeTagJet_pt[nGenJet] << " " << NeutralTimeTagJet_eta[nGenJet] << " | "
                << FullTimeTagJet_pt[nGenJet] << " " << FullTimeTagJet_eta[nGenJet] << std::endl;

      nGenJet+=1;
    }
    
    tree->Fill();
  
    //debug
    if(i%10==0)  {
      std::cout << "Event: " << i << " pu=" << nPU 
                << " signal/pu/soft kill: " << sigJets.size() << "/" << sigpuJets.size() << "/" << sk_sigpuJets.size() 
                << " soft kill threshold was: " << sk_thr << "GeV" << std::endl;
    }
    
  }

  //save TTree
  tree->Write();
  for(auto it : histos) it.second->Write();
  fout->Close();

  return 0;
}


/*
attic

    // histos[parts[i]+"_alphaF"] = new TH1F(parts[i]+"_alphaF",";#alpha^{F};Particles",250,0,20);
    // histos[parts[i]+"_alphapF"] = new TH1F(parts[i]+"_alphapF",";#alpha^{F}';Particles",250,0,20);
    // histos[parts[i]+"_wF"] = new TH1F(parts[i]+"_wF",";PUPPI Weight^{F};Particles",250,0,1);
    // histos[parts[i]+"_wpF"] = new TH1F(parts[i]+"_wpF",";PUPPI Weight^{F}';Particles",250,0,1);


//before clustering the jets we compute PUPPI weights
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
