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
  std::map<int,TString> parts = { {211,"ch"}, {22,"em"}, {13,"mu"}, {130,"nh"} };
  for(auto p : parts) {
    histos[p.second+"_alphaF"] = new TH2F(p.second+"_alphaF",";#alpha^{F};Pseudo-rapidity;Particles",250,0,20,10,1.5,3.0);
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

    //get particles for a given pileup
    UInt_t nPU = UInt_t( CLHEP::RandPoisson::shoot(hre,avgpu) );
    auto puParticles = getParticlesFrom(pu,0,nPU,true,"SimCluster");
    
    if(i%50==0) std::cout << "At event #" << i << " with PU=" << nPU << " and " << puParticles.size() << " particles" << std::endl;

    //build R=jetR cones around each particle and compute the alpha
    for(size_t ij=0; ij<puParticles.size(); ij++) {

      //puppi metric
      float alpha(0.);
      for(size_t jj=0; jj<puParticles.size(); jj++) {
        if(jj==ij) continue;
        float dRij(puParticles[ij].plain_distance(puParticles[jj]));
        if(dRij>0.3 || dRij<0.02) continue;
        alpha += puParticles[jj].e()/dRij;
      }

      //fill the histogram (if alpha is reasonable)
      if(alpha>0) {
        int pfid=abs(puParticles[ij].user_index());
        float abseta(puParticles[ij].eta());
        TString t(parts[pfid]);
        histos[t+"_alphaF"]->Fill(log(alpha),abseta);
      }
    }
  }
  
  for(auto it : histos) it.second->Write();
  fout->Close();
  
  return 0;
}
