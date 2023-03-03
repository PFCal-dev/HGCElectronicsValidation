//cmsRun test/hgchitcheckanalyzer_cfg.py #check cfg for options

#include <iostream>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"


/**
   @short this is a simple structure to store combined information of geometry, simhits, digis and rechits
*/
struct CustomHGCALDetIdAccumulator {
  CustomHGCALDetIdAccumulator() :
        hasValidDetId(false), hasValidSimHit(false), hasValidDigi(false), hasValidRecHit(false) { }
  CustomHGCALDetIdAccumulator(bool vdid, bool vsim, bool vdigi, bool vrec) :
    hasValidDetId(vdid), hasValidSimHit(vsim), hasValidDigi(vdigi), hasValidRecHit(vrec) { }
  CustomHGCALDetIdAccumulator(const CustomHGCALDetIdAccumulator& o) :
    hasValidDetId(o.hasValidDetId), hasValidSimHit(o.hasValidSimHit), hasValidDigi(o.hasValidDigi), hasValidRecHit(o.hasValidRecHit) { }
  bool hasValidDetId, hasValidSimHit, hasValidDigi, hasValidRecHit;
};


/**
   @short this is an example analyzer to get started with the matching of SIM-DIGI-REC hits
 */
class HGCHitCheckAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  
public:
  
  explicit HGCHitCheckAnalyzer( const edm::ParameterSet& );
  ~HGCHitCheckAnalyzer() {}
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();
  
private:
  
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitsCEE_,simHitsCEH_,simHitsCEHSci_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_,digisCEHSci_;
  edm::EDGetTokenT<HGCRecHitCollection> hitsCEE_,hitsCEH_,hitsCEHSci_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
};
 



//
HGCHitCheckAnalyzer::HGCHitCheckAnalyzer( const edm::ParameterSet &iConfig ) 
  : simHitsCEE_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE")) ),
    simHitsCEH_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront")) ),
    simHitsCEHSci_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback")) ),
    digisCEE_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
    digisCEH_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEfront")) ),
    digisCEHSci_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEback")) ),
    hitsCEE_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits")) ),
    hitsCEH_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits")) ),
    hitsCEHSci_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits")) ),
    caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>())
{   
}

//
void HGCHitCheckAnalyzer::endJob()
{
}

//
void HGCHitCheckAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  //key is the detId, content is the accumulator of information for a given detId
  std::map<uint32_t,CustomHGCALDetIdAccumulator> hgcalDetIdMapper;
  
  //read geometry components for HGCAL
  // 0 - CE-E
  // 1 - CE-H Si
  // 2 - CE-H SiPM-on-tile
  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  assert(geom.isValid());
  
  for(int i=0; i<3; i++) {

    //read the valid detIds from the geometry of the detector region
    auto detcode = (i==0? DetId::HGCalEE : (i==1? DetId::HGCalHSi : DetId::HGCalHSc));
    const HGCalGeometry *geo = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(detcode, ForwardSubdetector::ForwardEmpty));
    const std::vector<DetId>& validIds = geo->getValidDetIds();
    
    //get simHits, digis and recHits from the event
    edm::Handle<edm::PCaloHitContainer> simHits;
    edm::Handle<HGCalDigiCollection> digis;
    edm::Handle<HGCRecHitCollection> recHits;
    if(i==0) {
      iEvent.getByToken(simHitsCEE_, simHits);
      iEvent.getByToken(digisCEE_,   digis);
      iEvent.getByToken(hitsCEE_,    recHits);
    }
    if(i==1) {
      iEvent.getByToken(simHitsCEH_, simHits);
      iEvent.getByToken(digisCEH_,   digis);
      iEvent.getByToken(hitsCEH_,    recHits);
    }
    if(i==2) {
      iEvent.getByToken(simHitsCEHSci_, simHits);
      iEvent.getByToken(digisCEHSci_,  digis);
      iEvent.getByToken(hitsCEHSci_,   recHits);
    }
    assert(simHits.isValid());
    assert(digis.isValid());
    assert(recHits.isValid());
    
    std::cout << "@Sub-detector " << detcode
              << " #det-ids: " << validIds.size()
              << " #simHits: " << simHits->size()
              << " #digis: " << digis->size()
              << " #recHits: " << recHits->size()
              << std::endl;

    //flag valid detIds
    for(auto id : validIds) hgcalDetIdMapper[id.rawId()] = CustomHGCALDetIdAccumulator(true,false,false,false);

    //flag detIds with valid simhits
    for(auto sh : *simHits) {
      uint32_t rawid(sh.id());
      auto it=hgcalDetIdMapper.find(rawid);
      if(it==hgcalDetIdMapper.end()) {
        hgcalDetIdMapper[rawid] = CustomHGCALDetIdAccumulator(false,true,false,false);
      } else {
        it->second.hasValidSimHit=true;
      }
   }

    //flag detIds with valid digis
    for(auto d : *digis) {
      uint32_t rawid(d.id().rawId());
      auto it=hgcalDetIdMapper.find(rawid);
      if(it==hgcalDetIdMapper.end()) {        
        hgcalDetIdMapper[rawid] = CustomHGCALDetIdAccumulator(false,false,true,false);
      } else {
        it->second.hasValidDigi=true;
      }
    }

    //flag detIds with valid rechits
    for(auto r: *recHits) {
      uint32_t rawid(r.id().rawId());
      auto it=hgcalDetIdMapper.find(rawid);
      if(it==hgcalDetIdMapper.end()) {
        hgcalDetIdMapper[rawid] = CustomHGCALDetIdAccumulator(false,false,false,true);
      } else {
        it->second.hasValidRecHit=true;
      }
    }
  }


  //analyze the accumulators, maybe fill some plots etc. :)
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCHitCheckAnalyzer);
