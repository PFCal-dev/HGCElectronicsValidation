#ifndef _HGCDigiTester_h_
#define _HGCDigiTester_h_

#include <tuple>

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

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"  
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSciNoiseMap.h"  
#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerBase.h"  
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

  typedef std::tuple<int, bool, int, int, int> rocKey_t; //not the Balboa, though...
  typedef std::pair<int,float> rocSummary_t;        
  typedef std::map<const rocKey_t,rocSummary_t> rocDeposits_t;

 private:
  
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitsColls_[3];
  edm::EDGetTokenT<HGCalDigiCollection> digisColls_[3];
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;

  std::vector< HGCalSiNoiseMap<HGCSiliconDetId> *> scal_;
  HGCalSciNoiseMap *scalSci_;

  float refSpeed_;
  double bxTime_[3];
  double tofDelay_[3];

  uint32_t mipTarget_[3];
  double tdcLSB_[3],vanilla_adcLSB_fC_[3];
  std::vector<double> avg_mipfC_[2];
  double sci_keV2MIP_;
  double tdcOnset_fC_[3];
  bool useTDCOnsetAuto_;
  bool useVanillaCfg_,scaleByTileArea_,scaleBySipmArea_;
  double pxFiringRate_;

  uint32_t detid_;
  Int_t event_,layer_,u_,v_,roc_,thick_,isSci_,isToT_,isSat_;
  Float_t gpt_,geta_,gphi_,genergy_,gvradius_,gvz_;
  uint32_t adc_, gain_;
  Float_t qsim_,qrec_,mipsim_,avgmipsim_,miprec_,avgmiprec_,cce_,eta_,radius_,z_;
  Float_t qsimInBX_,mipsimInBX_; // in-time BX
  Float_t qsimPreBX_,mipsimPreBX_; // before BX
  Float_t qsimPostBX_,mipsimPostBX_; // after BX
  Int_t nhits_; //only for the rocTree
  Bool_t side_; //only for the rocTree
  TTree *tree_,*rocTree_;
    
  bool hardProcOnly_;
  bool onlyROCTree_;
};
 

#endif
