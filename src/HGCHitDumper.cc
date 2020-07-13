#include "UserCode/HGCElectronicsValidation/interface/HGCHitDumper.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <fstream>
#include <iostream>

#include "TLorentzVector.h"

using namespace std;


//
// PLUGIN IMPLEMENTATION
//


//
HGCHitDumper::HGCHitDumper( const edm::ParameterSet &iConfig ) :   
  geo_("HGCalEESensitive"),
  digis_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
  recHits_( consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits")) ),
  genJets_( consumes<std::vector<reco::GenJet> >(edm::InputTag("ak4GenJetsNoNu")) ),
  event_(0)
{  
  edm::Service<TFileService> fs;
  t_ = fs->make<TTree>("hits","hits");
  t_->Branch("ng",     &ng_,     "ng/I");
  t_->Branch("en",     &en_,     "en/F");
  t_->Branch("eta",    &eta_,    "eta/F");
  t_->Branch("layer",  &layer_,  "layer/I");
  t_->Branch("waferL", &waferL_, "waferL/I");
  t_->Branch("gain",   &gain_,   "gain/I");
  t_->Branch("isTDC",  &isTDC_,  "isTDC/O");
  t_->Branch("threshold", &threshold_, "threshold/O");
  t_->Branch("hasTOA", &hasTOA_, "hasTOA/O");
  t_->Branch("digiT",  &digiT_,  "digiT/I");
  t_->Branch("digiE",  &digiE_,  "digiE/I");
  t_->Branch("recT",   &recT_,   "recT/F");
  t_->Branch("recE",   &recE_,   "recE/F");
}

//_
HGCHitDumper::~HGCHitDumper()
{
}

//
void HGCHitDumper::endJob()
{
}

  
//
void HGCHitDumper::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  event_++;

  edm::Handle<std::vector<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJets_,genJetsHandle);
  std::vector<TLorentzVector> photons(2);
  for(auto &j : *genJetsHandle){
    if(fabs(j.eta())<1.5 || fabs(j.eta())>3) continue;   
    size_t phoIdx( j.eta()<0 ? 0 : 1);
    photons[phoIdx].SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),0.);
  }


  //read geometry from event setup
  edm::ESHandle<HGCalGeometry> geoHandle;
  iSetup.get<IdealGeometryRecord>().get(geo_,geoHandle);
  const HGCalDDDConstants &ddd=geoHandle->topology().dddConstants();

  //map digis to detIds
  edm::Handle<HGCalDigiCollection> digisHandle;
  iEvent.getByToken(digis_,digisHandle);
  const auto &digisColl=digisHandle.product();
  std::map<HGCSiliconDetId, size_t> detId2Digi;
  for(size_t idigi=0; idigi<digisColl->size(); idigi++)
    {
      const auto &hit = (*digisColl)[idigi];
      if(hit.size()==0) continue;
      HGCSiliconDetId detId(hit.id());
      detId2Digi[detId]=idigi;
    }
  
  //analyze rec hits matching to digis
  edm::Handle<HGCRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHits_,recHitsHandle);
  for(auto &recHit : *recHitsHandle) {
    
    HGCSiliconDetId detId(recHit.id());
    int phoOffset( detId.zside()<0 ? 0 : 1);
    size_t digisIdx=detId2Digi[detId];
    const auto &digi = (*digisColl)[digisIdx];

    //save info
    ng_        = event_*2+phoOffset;
    en_        = photons[phoOffset].E();
    eta_       = photons[phoOffset].Eta();
    layer_     = detId.layer();
    waferL_    = ddd.waferType(recHit.id());
    threshold_ = digi.sample(2).threshold();
    gain_      = digi.sample(2).gain();
    isTDC_     = digi.sample(2).mode() ;
    hasTOA_    = digi.sample(2).getToAValid();
    digiT_     = digi.sample(2).toa();
    digiE_     = digi.sample(2).data();
    recT_      = recHit.time();
    recE_      = recHit.energy();
    t_->Fill();
  }

}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HGCHitDumper);
