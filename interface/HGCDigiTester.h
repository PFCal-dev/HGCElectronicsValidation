#ifndef _HGCDigiTester_h_
#define _HGCDigiTester_h_

#include <tuple>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"  
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSciNoiseMap.h"  
#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerBase.h"  
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"


#include "TTree.h"

/**
   @class HGCDigiTester
   @short to test energy -> digi -> energy conversion
*/

class ModuleToBE{
public:
  ModuleToBE(int32_t _l,int32_t _u,int32_t _v,int32_t _c,int32_t _s) : layer(_l), u(_u), v(_v), crate(_c), slot(_s) { }
  ModuleToBE(const ModuleToBE& m) : layer(m.layer), u(m.u), v(m.v), crate(m.crate), slot(m.slot) { }
  bool operator==(const ModuleToBE& m) const
  {
    if(m.layer != this->layer) return false;
    if(m.u != this->u) return false;
    if(m.v != this->v) return false;
    return true;    
  }
  int32_t layer,u,v,crate,slot;
};



class HGCDigiTester : public edm::one::EDAnalyzer<edm::one::SharedResources>
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
  
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitsCEE_,simHitsCEH_,simHitsCEHSci_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_,digisCEHSci_;
  edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  edm::EDGetTokenT<float> genT0_;

  std::vector< HGCalSiNoiseMap<HGCSiliconDetId> *> scal_;
  HGCalSciNoiseMap *scalSci_;

  uint32_t mipTarget_[3];
  double tdcLSB_[3],vanilla_adcLSB_fC_[3];
  std::vector<double> avg_mipfC_[2];
  double sci_keV2MIP_,maxDeltaR_;
  double tdcOnset_fC_[3],toaLSB_ns_[3],bxTime_[3],tofDelay_[3];
  bool useTDCOnsetAuto_;
  bool useVanillaCfg_;
  double pxFiringRate_;

  uint32_t detid_,crate_,slot_;
  Int_t event_,layer_,u_,v_,roc_,thick_,isSci_,isToT_,isSat_,crossCalo_,clustJetAlgo_,inShower_;
  Float_t gpt_,geta_,gphi_,genergy_,gvradius_,gvz_,gvt_,gbeta_;
  uint32_t adc_, gain_, toa_;
  Float_t qsim_,qrec_,mipsim_,avgmipsim_,miprec_,avgmiprec_,cce_,eta_,radius_,z_, toarec_,toasim_;
  Int_t nhits_; //only for the rocTree
  Bool_t side_; //only for the rocTree
  TTree *tree_,*rocTree_;
    
  bool hardProcOnly_;
  bool onlyROCTree_;

  std::vector<ModuleToBE> module2be_map_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
};
 

#endif
