#ifndef _runmixandclustercommon_h_
#define _runmixandclustercommon_h_

#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "Math/Vector4D.h"
#include "Math/ProbFuncMathCore.h"

#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/ClusterSequence.hh"

struct PseudoJetProperties {
  int pfid;
  bool ispu;
  float gvz,gvt;
};

/**
   @short reads the required entries from the chain and returns a vector of particles for clustering
 */
std::vector<fastjet::PseudoJet> getParticlesFrom(TChain *t,int ievt,int nevts,bool ispu, std::vector<PseudoJetProperties> &properties, TString tag="GenPart") {

  //atach variables to read from
  Float_t GenVtx_t0,GenVtx_z;
  UInt_t maxPseudoParticles(10000);
  UInt_t nPseudoParticle;
  Int_t PseudoParticle_pdgId[maxPseudoParticles], PseudoParticle_status[maxPseudoParticles];
  Float_t PseudoParticle_eta[maxPseudoParticles], PseudoParticle_mass[maxPseudoParticles], PseudoParticle_phi[maxPseudoParticles], PseudoParticle_pt[maxPseudoParticles];
  Bool_t PseudoParticle_crossedBoundary[maxPseudoParticles];
  t->SetBranchAddress("GenVtx_t0",&GenVtx_t0);
  t->SetBranchAddress("GenVtx_z",&GenVtx_z);
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

      //save the properties of this particle
      size_t idx=properties.size();
      PseudoJetProperties p;
      p.pfid=pfid;
      p.ispu=ispu;
      p.gvt=GenVtx_t0;
      p.gvz=GenVtx_z;
      properties.push_back(p);
      
      //add particle with properties index
      ip.set_user_index(idx);
      pseudoParticles.push_back( ip );
    }
  }
  

  return pseudoParticles;
}

/**
   @short computes the local puppi shape
*/
float getLocalPuppiShape(size_t ij,std::vector<fastjet::PseudoJet> &particles, std::vector<PseudoJetProperties> &properties,bool requireSameID=false,float maxDR=0.3, float minDR=0.02) {

  float alpha(0.);

  size_t ij_idx=particles[ij].user_index();
  int ii_pfid=properties[ij_idx].pfid;

  for(size_t jj=0; jj<particles.size(); jj++) {

    if(jj==ij) continue;

    //require within cone of interest
    float dRij(particles[ij].plain_distance(particles[jj]));
    if(dRij>maxDR || dRij<minDR) continue;
    
    //require particles of the same type
    if(requireSameID) {
      size_t jj_idx=particles[jj].user_index();
      int jj_pfid=properties[jj_idx].pfid;
      if(jj_pfid!=ii_pfid) continue;
    }
        
    alpha += pow(particles[jj].pt()/dRij,2);
  }

  return alpha>0 ? TMath::Max((float)log(alpha),(float)0.) : 0.;
}

/**
   @short computes the puppi weight for a given particle at a given pseudo-rapidity, known the local shape variable
 */
float getPuppiWgt(float alpha,int pid,float eta,bool useSameID=false) {

  float a_avg,b_avg,c_avg,a_sigma,b_sigma,c_sigma;

  pid=abs(pid);
  if(useSameID){
    
    if(pid==211) {
       a_avg=-0.1117062893367535; b_avg=0.3943703174537112; c_avg=6.464279330262225;
       a_sigma=0.021366203180838022; b_sigma=-0.10521730154550395; c_sigma=0.35814223346744983;
    }
    else if(pid==130) {
      a_avg=-0.11311657421435714; b_avg=0.35944361198517905; c_avg=4.86627303307938; 
      a_sigma=0.09482800182025537; b_sigma=-0.2720427168624265; c_sigma=0.5507493177662346;
    }
    else if(pid==22) {
      a_avg=-0.1554091322378609; b_avg=0.4412464639527004; c_avg=5.461334143487884;
      a_sigma=0.02604973317036911; b_sigma=-0.08991248893733884; c_sigma=0.3162693838653237;
    }
    else if(pid==13) {
      a_avg=-0.08827534544223824; b_avg=-0.3455310214701993; c_avg=4.111067540984826;
      a_sigma=-0.04679968263635605; b_sigma=0.6775003896822301; c_sigma=-0.3515760863095199;
    }
    else {
      return 0.;
    }
    
  } else {
    
    if(pid==211) {
      a_avg=-0.11087621589022852; b_avg=0.3467466375259241; c_avg=6.998355408597446;
      a_sigma=0.023131569741663336; b_sigma=-0.11079591809323978; c_sigma=0.34923503941503337;
    }
    else if(pid==130) {
      a_avg=-0.30236821330615354; b_avg=1.5055583734473736; c_avg=5.326297392679808;
      a_sigma=0.20766277659475557; b_sigma=-1.1499494971481288; c_sigma=1.8305735481866885; 
    }
    else if(pid==22) {
      a_avg=-0.04591190786432857; b_avg=-0.014611940357667871; c_avg=7.458206699764853;
      a_sigma=-0.004409162133472597; b_sigma=0.02509765469664976; c_sigma=0.1771960704737882;
    }
    else if(pid==13) {
      a_avg=-0.14728351881302346; b_avg=0.5159745068218985; c_avg=6.825325373738728; 
      a_sigma=0.001685939556387035; b_sigma=0.005428267885577016; c_sigma=0.19072405101898207;
    }
    else {
      return 0.;
    }

  }
  
  //compute the median and width expected for the distribution of log(alpha)
  eta=fabs(eta);
  float avg=a_avg*eta*eta+b_avg*eta+c_avg;
  float sigma=a_sigma*eta*eta+b_sigma*eta+c_sigma;

  //compute the chi2
  float chi2 = (alpha>a_avg)*pow((alpha-avg)/sigma,2);
  double pvalue = ROOT::Math::chisquared_cdf_c(chi2,1);
  return pvalue;
} 
  


#endif
