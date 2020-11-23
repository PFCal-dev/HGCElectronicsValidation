#ifndef _HGCDigiTester_h_
#define _HGCDigiTester_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerBase.h"  
#include "SimCalorimetry/HGCalSimProducers/interface/HGCEEDigitizer.h"
#include "SimCalorimetry/HGCalSimProducers/interface/HGCHEfrontDigitizer.h"
#include "SimCalorimetry/HGCalSimProducers/interface/HGCHEbackDigitizer.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "TTree.h"

/**
   @class HGCDigiTester
   @short to test energy -> digi -> energy conversion
*/

class HGCDigiTester : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCDigiTester( const edm::ParameterSet& );
  ~HGCDigiTester() {}
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitsCEE_,simHitsCEH_,simHitsCEHSci_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_,digisCEHSci_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;

  std::vector<HGCalSiNoiseMap *>scal_;
  HGCalSciNoiseMap *scalSci_;

  uint32_t mipTarget_[3];
  double tdcLSB_[3],vanilla_adcLSB_fC_[3];
  std::vector<double> avg_mipfC_[3];
  double sci_keV2MIP_;
  double tdcOnset_fC_[3];
  bool useVanillaCfg_,scaleByTileArea_,scaleBySipmArea_;

  Int_t event_,layer_,thick_,isSci_,isToT_;
  Float_t gpt_,geta_,gphi_,genergy_,gvradius_,gvz_;
  Float_t qsim_,qrec_,mipsim_,avgmipsim_,miprec_,avgmiprec_,cce_,eta_,radius_,z_;
  TTree *tree_;
};
 

#endif
