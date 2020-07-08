#include "UserCode/HGCElectronicsValidation/interface/HGCOccupancyAnalyzer.h"

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

using namespace std;


//
// PLUGIN IMPLEMENTATION
//


//
HGCOccupancyAnalyzer::HGCOccupancyAnalyzer( const edm::ParameterSet &iConfig ) :   
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"))),
  genJets_( consumes<std::vector<reco::GenJet> >(edm::InputTag("ak8GenJetsNoNu")) ),
  geoCEE_("HGCalEESensitive"),
  geoCEH_("HGCalHESiliconSensitive"),
  digisCEE_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
  digisCEH_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEfront")) ),
  nevts_(0),
  adcThrMIP_( iConfig.getParameter<double>("adcThrMIP") ),
  adcThrMIPbxm1_(iConfig.getParameter<double>("adcThrMIPbxm1") ),
  doseMap_("SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt")
{  

  edm::Service<TFileService> fs;

  //some global histograms
  cellCount_=fs->make<TH1F>("cellcount",";Layer;Number of cells/layer",50,0.5,50.5);
  tdcCountProf_=fs->make<TProfile>("tdccount",";Layer;Number of hits/layer",50,0.5,50.5);
  toaCountProf_=fs->make<TProfile>("toacount",";Layer;Number of hits/layer",50,0.5,50.5);
  adcCountProf_=fs->make<TProfile>("adccount",";Layer;Number of hits/layer",50,0.5,50.5);
  adcHitsVsPU_=fs->make<TH2F>("adchitsvspu",";Number of PU interactions; Number of ADC hits",100,100,300,1000,0.5e5,5e5);
  tdcHits_.resize(50,0);
  toaHits_.resize(50,0);
  adcHits_.resize(50,0);

  //start the noise map tools
  for(size_t subdet=0; subdet<2; subdet++){
    noiseMaps_[subdet]=new HGCalSiNoiseMap;
    noiseMaps_[subdet]->setDoseMap(iConfig.getParameter<std::string>("doseMap"),
                                   iConfig.getParameter<uint32_t>("scaleByDoseAlgo"));
    noiseMaps_[subdet]->setIleakParam(iConfig.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
    noiseMaps_[subdet]->setCceParam(iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"),
                                    iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"),
                                    iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));

  }
}

//
HGCOccupancyAnalyzer::~HGCOccupancyAnalyzer()
{
}

//
void HGCOccupancyAnalyzer::endJob()
{
  for(auto &it : waferHistos_) it.second->endJob();
}

  
//
void HGCOccupancyAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  nevts_+=1;

  edm::Service<TFileService> fs;

  //reset global counters
  std::fill(tdcHits_.begin(), tdcHits_.end(), 0);
  std::fill(toaHits_.begin(), toaHits_.end(), 0);
  std::fill(adcHits_.begin(), adcHits_.end(), 0);

  //read geometry from event setup
  std::string subdets[]={"CEE","CEH"};

  for(size_t subdet=0; subdet<2; subdet++) {

    //get the geometry
    std::string sd=subdets[subdet];
    edm::ESHandle<HGCalGeometry> ceeGeoHandle;
    iSetup.get<IdealGeometryRecord>().get(subdet==0 ? geoCEE_ : geoCEH_,ceeGeoHandle);
    const HGCalGeometry *geo=ceeGeoHandle.product();
    const HGCalDDDConstants &ddd=geo->topology().dddConstants();

    //configure noise map
    noiseMaps_[subdet]->setGeometry(geo, HGCalSiNoiseMap::AUTO, 10);

    //instantiate occupancy histograms first time around
    if(nevts_==1) {

      std::cout << "[ HGCOccupancyAnalyzer::analyze ] starting occupancy histos for " << subdets[subdet] << std::endl;

      const auto &validDetIds = geo->getValidDetIds();
      int newWafers(0);
      for(const auto &didIt : validDetIds) {
            
        HGCSiliconDetId detId(didIt.rawId());
            
        //use only positive side
        if(detId.zside()<0) continue;
        int layer=detId.layer();
        int layidx(layer);
        std::pair<int,int> waferUV=detId.waferUV();
        
        WaferOccupancyHisto::WaferKey_t key( subdet, layer, waferUV.first, waferUV.second);
        if(waferHistos_.find(key)==waferHistos_.end()) {
          newWafers++;
          waferHistos_[key]=new WaferOccupancyHisto(key);
        }
        waferHistos_[key]->addPad( ddd.waferType(detId) );

        //global pad counter
        if(subdet==1) layidx+=28;
        cellCount_->Fill(layidx);
      }
      
      std::cout << "\t" << validDetIds.size() << " valid detIds => " << newWafers << " wafer histos" << endl;
            
      //init histos (some of them are already done but the call will be ignored)
      for(auto &it : waferHistos_) it.second->bookHistos(&fs);
    }
  

    //analyze digi collections
    edm::Handle<HGCalDigiCollection> digisHandle;
    iEvent.getByToken(subdet==0 ? digisCEE_ : digisCEH_, digisHandle);
    analyzeDigis(subdet,digisHandle,geo);
   
    //get pileup information for profile
    edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(puToken_,PupInfo);
    int npu(0);
    if(PupInfo.isValid()){
      std::vector<PileupSummaryInfo>::const_iterator ipu;
      for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) 
        {
          if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
          npu=ipu->getPU_NumInteractions();
        }
    }

    int totalADC(0);
    for(size_t idx=0; idx<tdcHits_.size(); idx++){
      tdcCountProf_->Fill(idx+1,tdcHits_[idx]);
      toaCountProf_->Fill(idx+1,toaHits_[idx]);
      adcCountProf_->Fill(idx+1,adcHits_[idx]);
      totalADC+=adcHits_[idx];
    }
    adcHitsVsPU_->Fill(npu,totalADC);
  }  
   
  //all done, reset counters
  for(auto &it : waferHistos_) {
    it.second->analyze();
    it.second->resetCounters();    
  }

}

//
void HGCOccupancyAnalyzer::analyzeDigis(int subdet,edm::Handle<HGCalDigiCollection> &digiColl, const HGCalGeometry *geom)
{
  //check inputs
  if(!digiColl.isValid() || geom==NULL) return;
  //analyze hits
  const int itSample(2); //in-time sample
  for(auto &hit : *digiColl)
    {
      if(hit.size()==0) continue;

      //detid inf
      HGCSiliconDetId detId(hit.id());

      if(detId.zside()<0) continue;

      //wafer id
      int layer=detId.layer();
      std::pair<int,int> waferUV=detId.waferUV();
      WaferOccupancyHisto::WaferKey_t key(std::make_tuple(subdet,layer,waferUV.first,waferUV.second));

      //re-compute the thresholds
      HGCalSiNoiseMap::SiCellOpCharacteristics siop=noiseMaps_[subdet]->getSiCellOpCharacteristics(detId);
      int mipADC=siop.mipADC;
      
      //in-time BX info
      uint32_t rawData(hit.sample(itSample).data() );
      bool isTOA( hit.sample(itSample).getToAValid() );
      bool isTDC( hit.sample(itSample).mode() );
      bool isBusy( isTDC && rawData==0 );
      uint32_t thr( std::floor(mipADC*adcThrMIP_) );

      //BX-1 info
      uint32_t rawDatabxm1(hit.sample(itSample-1).data() );
      bool isTOAbxm1( hit.sample(itSample-1).getToAValid() );
      bool isTDCbxm1( hit.sample(itSample-1).mode() );
      bool isBusybxm1( isTDCbxm1 && rawDatabxm1==0 );
      uint32_t thrbxm1( std::floor(mipADC*adcThrMIPbxm1_) );

      //ZX algo
      bool passZS(true);
      uint32_t bshift(4);
      if(siop.core.gain==HGCalSiNoiseMap::q80fC)   bshift=3;
      if(siop.core.gain==HGCalSiNoiseMap::q160fC)  bshift=4;
      if(siop.core.gain==HGCalSiNoiseMap::q320fC)  bshift=3;
      uint32_t zsCorr( (rawDatabxm1>>bshift) );
      passZS=( rawData > zsCorr+thr );
    
      waferHistos_[key]->count(rawData,     rawData-zsCorr, passZS, isTOA,     isTDC,     isBusy,     thr);
      waferHistos_[key]->count(rawDatabxm1, rawDatabxm1,    true,   isTOAbxm1, isTDCbxm1, isBusybxm1, thrbxm1, "_bxm1");
      
      //global counts for the in-time bunch
      size_t layidx(layer-1);
      if(subdet==1) layidx+=28;
      if(!isBusy) {
        if(isTDC)  tdcHits_[layidx]+=1;
        if(isTOA)  toaHits_[layidx]+=1;
        if(!isTDC) adcHits_[layidx]+=1;
      }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCOccupancyAnalyzer);
