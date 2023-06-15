#include <iostream>
#include <algorithm>

#include "UserCode/HGCElectronicsValidation/bin/RunMixAndClusterCommon.h"

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

std::map<TString, TH2F *> histos;

int main(int argc, char** argv) {

  //if passed at command line use new cfi
  std::string url("UserCode/HGCElectronicsValidation/bin/mixandcluster_forpuppi_cfi.py");
  if (argc > 1) url = argv[1];
  url = edm::FileInPath(url).fullPath();

  //get configuration
  const std::shared_ptr<edm::ParameterSet> &pset = edm::readPSetsFrom(url);
  const edm::ParameterSet &cfg = pset->getParameter<edm::ParameterSet>("mixandcluster");
  std::vector<std::string> pufiles = cfg.getParameter<std::vector<std::string> >("pu");
  int maxevts = cfg.getParameter<int>("maxevts");
  int avgpu = cfg.getParameter<int>("avgpu");

  //
  TString foutName=Form("puppi_localshapes_%dPU.root",avgpu);
  TFile *fout=TFile::Open(foutName,"RECREATE");
  std::map<TString,int> parts = { {"ch",211}, {"em",22}, {"mu",13}, {"nh",130} };
  std::map<int,TString> parts_inv;
  for(auto p : parts) {

    //create the inverse map
    TString pname(p.first);
    int pfid(p.second);
    parts_inv[pfid]=pname;

    //create the histograms
    histos[pname+"_alphaF"] = new TH2F(pname+"_alphaF",";#alpha^{F};Pseudo-rapidity;Particles",250,0,20,10,1.5,3.0);
    histos[pname+"_alphapF"] = new TH2F(pname+"_alphapF",";#alpha^{F}';Pseudo-rapidity;Particles",250,0,20,10,1.5,3.0);
    histos[pname+"_wF"] = new TH2F(pname+"_wF",";w^{F};Pseudo-rapidity;Particles",250,0,1,10,1.5,3.0);
    histos[pname+"_wpF"] = new TH2F(pname+"_wpF",";w^{F}';Pseudo-rapidity;Particles",250,0,1,10,1.5,3.0);
  }
  fout->cd();
  

  //start chain
  TChain *pu=new TChain("Events");
  for(size_t i=0; i<pufiles.size(); i++) pu->AddFile(pufiles[i].c_str());
  int npuEvts(pu->GetEntries());

  std::cout << "Generating an average pileup of " << avgpu 
            << " from " << pufiles.size() << " pileup files with " << npuEvts << " events" << std::endl;

  //random number generator
  CLHEP::HepJamesRandom *hre = new CLHEP::HepJamesRandom();
  hre->setSeed(0);
  
  maxevts = maxevts<0 ? npuEvts : maxevts;
  for(int i=0; i<maxevts; i++) {

    //stores the properties of all the particles we'll look at in the event
    std::vector<PseudoJetProperties> properties;
    
    //get particles for a given pileup
    UInt_t nPU = UInt_t( CLHEP::RandPoisson::shoot(hre,avgpu) );
    auto puParticles = getParticlesFrom(pu,0,nPU,true,properties,"SimCluster");
    
    if(i%20==0) std::cout << "At event #" << i << " with PU=" << nPU << " and " << puParticles.size() << " particles (" << properties.size() << " properties)" << std::endl;

    //build R=jetR cones around each particle and compute the alpha
    for(size_t ij=0; ij<puParticles.size(); ij++) {

      int ij_idx=puParticles[ij].user_index();
      int pfid=properties[ij_idx].pfid;      
      TString t(parts_inv[pfid]);
      float abseta(puParticles[ij].eta());
      
      //puppi metric (local shape)
      float alpha=getLocalPuppiShape(ij,puParticles,properties,false);
      float alphap=getLocalPuppiShape(ij,puParticles,properties,true);

      //fill the histogram (if alpha is reasonable)
      if(alpha>0) {
        histos[t+"_alphaF"]->Fill(log(alpha),abseta);
      }
      if(alphap>0) {
        histos[t+"_alphapF"]->Fill(log(alphap),abseta);
      }
      
      //puppi weight
      float wgt=getPuppiWgt(alpha,pfid,abseta,false);
      histos[t+"_wF"]->Fill(wgt,abseta);
      float wgtp=getPuppiWgt(alphap,pfid,abseta,true);
      histos[t+"_wpF"]->Fill(wgtp,abseta);
      
    }
  }
  
  for(auto it : histos) it.second->Write();
  fout->Close();
  
  return 0;
}
