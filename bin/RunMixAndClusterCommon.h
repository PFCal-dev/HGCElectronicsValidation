#ifndef _runmixandclustercommon_h_
#define _runmixandclustercommon_h_

#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "Math/Vector4D.h"

#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/ClusterSequence.hh"

/**
   @short reads the required entries from the chain and returns a vector of particles for clustering
 */
std::vector<fastjet::PseudoJet> getParticlesFrom(TChain *t,int ievt,int nevts,bool ispu, TString tag="GenPart") {

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

  std::vector<fastjet::PseudoJet> pseudoParticles;

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
      auto ip = fastjet::PseudoJet(p4.px(),p4.py(),p4.pz(),p4.energy());

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


#endif
