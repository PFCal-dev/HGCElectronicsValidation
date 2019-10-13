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
  genJets_( consumes<std::vector<reco::GenJet> >(edm::InputTag("ak8GenJetsNoNu")) ),
  geoCEE_("HGCalEESensitive"),
  geoCEH_("HGCalHESiliconSensitive"),
  digisCEE_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
  digisCEH_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEfront")) )
{  

  edm::Service<TFileService> fs;

  //parse u-v equivalence map file and start histograms
  edm::FileInPath uvmapF("UserCode/HGCElectronicsValidation/data/wafer_pos.dat");
  std::ifstream inF(uvmapF.fullPath());  
  while(inF) {

    std::string buf;
    getline(inF,buf);

    std::stringstream ss(buf);
    std::vector<std::string> tokens;
    while (ss >> buf)
      if(buf.size()>0)
        tokens.push_back(buf);
    if(tokens.size()<14) continue;

    std::string subdet(tokens[0]);
    int ilay(atoi(tokens[1].c_str()));
    int waferU(atoi(tokens[9].c_str()));
    int waferV(atoi(tokens[10].c_str()));
    int ncells(atoi(tokens[11].c_str()));
    // float phi(atof(tokens[6].c_str()));
    // float eta(atof(tokens[8].c_str()));

    WaferEquivalentId_t key(subdet,ilay,waferU,waferV);
    waferHistos_[key]=new WaferOccupancyHisto(subdet,ilay,waferU,waferV,ncells,&fs);
    for(size_t i=13; i<tokens.size(); i+=2) {
      int u=atoi(tokens[i].c_str());
      int v=atoi(tokens[i+1].c_str());
      waferHistos_[key]->addWaferEquivalent(u,v);

      WaferEquivalentId_t ikey(subdet,ilay,u,v);
      uvEqMap_[ikey]=std::pair<int,int>(waferU,waferV);
    }


    std::vector<TH1F *> hotoccHistos;
    for(int iwaf=0; iwaf<7; iwaf++) {
      TString name(Form("%s_lay%d_hottestwafer%d",subdet.c_str(),ilay,iwaf));
      hotoccHistos.push_back( fs->make<TH1F>(name,";Occupancy;",500,0,1) );
      hotoccHistos[iwaf]->Sumw2();
    }
    std::pair<std::string,int> hotkey(subdet,ilay);
    hottestWaferH_[hotkey]=hotoccHistos;
  }

  inF.close();
}

//
HGCOccupancyAnalyzer::~HGCOccupancyAnalyzer()
{
}

//
void HGCOccupancyAnalyzer::endJob()
{
  /*
    for(auto &it : waferHistos_) it.second->endJob();
  */
}

  
//
void HGCOccupancyAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //read geometry from event setup
  edm::ESHandle<HGCalGeometry> ceeGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEE_,ceeGeoHandle);
  hgcGeometries_["CEE"]=ceeGeoHandle.product();
  edm::ESHandle<HGCalGeometry> cehGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEH_,cehGeoHandle);
  hgcGeometries_["CEH"]=cehGeoHandle.product();

  /*
  //best matching wafer to gen (if dR>0.2 no match was found)
  std::map<std::pair<std::string,int>, std::set<WaferOccupancyHisto::UVKey_t> > wafersOfInterest;
  edm::Handle<std::vector<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJets_,genJetsHandle);
  for(auto &j : *genJetsHandle){
    if(j.pt()<100) continue;
    if(fabs(j.eta())<1.5 || fabs(j.eta())>3) continue;   
       
    float minDR(9999.);
    WaferEquivalentId_t bestMatchedWafer;
    for(auto &w : waferEtaPhi_) {
      float dR(deltaR(j.eta(),j.phi(),w.second.first,w.second.second));      
      if(dR>minDR) continue;
      minDR=dR;
      bestMatchedWafer=w.first;
    }

    //add new best matched wafer
    if(minDR>0.2) continue;
    std::string sd=std::get<0>(bestMatchedWafer);
    int lay=std::get<1>(bestMatchedWafer);
    std::pair<std::string,int> key(sd,lay);
    if(wafersOfInterest.find(key)==wafersOfInterest.end())
      wafersOfInterest[key]=std::set<WaferOccupancyHisto::UVKey_t>();

    int u=std::get<2>(bestMatchedWafer);
    int v=std::get<3>(bestMatchedWafer);
    WaferOccupancyHisto::UVKey_t uv(u,v);
    wafersOfInterest[key].insert(uv);
  }
  */
  

  //analyze digi collections
  edm::Handle<HGCalDigiCollection> ceeDigisHandle;
  iEvent.getByToken(digisCEE_,ceeDigisHandle);
  analyzeDigis("CEE",ceeDigisHandle);
  
  edm::Handle<HGCalDigiCollection> cehDigisHandle;
  iEvent.getByToken(digisCEH_,cehDigisHandle);
  analyzeDigis("CEH",cehDigisHandle);

  //fill wafer histos and save the max. found in each layer
  std::map<std::pair<std::string,int>,std::pair<WaferOccupancyHisto::UVKey_t, float> > hotWaferOccPerLayer;
  for(auto &it : waferHistos_) {

    std::string sd=std::get<0>(it.first);
    int lay=std::get<1>(it.first);
    std::pair<std::string,int> key(sd,lay);

    /*
      //filter wafers matched to gen objects by UV coordinates
      std::set<WaferOccupancyHisto::UVKey_t> genMatchedUVs;
      if(wafersOfInterest.find(key)!=wafersOfInterest.end()) {
      for(auto &w:wafersOfInterest[key]) {
      if(!it.second->isUVEquivalent(w)) continue;
      genMatchedUVs.insert(w);
      }
      }
    */

    //analyze counts
    std::set<WaferOccupancyHisto::UVKey_t> genMatchedUVs;
    it.second->analyze(genMatchedUVs);

    //get hottest wafer
    WaferOccupancyHisto::UVKey_t hotWaferUV=it.second->getHotWaferUV();
    float hotWaferOcc=float(it.second->getHotWaferCounts())/float(it.second->getCells());    

    if(hotWaferOccPerLayer.find(key)==hotWaferOccPerLayer.end() || hotWaferOccPerLayer[key].second<hotWaferOcc)       
      hotWaferOccPerLayer[key]=std::pair<WaferOccupancyHisto::UVKey_t,float>(hotWaferUV,hotWaferOcc);
  }
    
  //fill the max. counts per layer and in the neighboring cells
  for(auto &it : hotWaferOccPerLayer) {

    std::string sd=it.first.first;
    int lay=it.first.second;
    WaferOccupancyHisto::UVKey_t hotWaferUV=it.second.first;
    if(hotWaferUV.first==0 && hotWaferUV.second==0) continue; //empty layer
    float hotOcc=it.second.second;

    std::pair<std::string,int> key(sd,lay);
    hottestWaferH_[key][0]->Fill(hotOcc);

    //find neighbors
    std::vector<float> neighborOccs;
    for(auto &jt : waferHistos_) {
      std::string isd=std::get<0>(jt.first);
      if(isd!=sd) continue;
      int ilay=std::get<1>(jt.first);
      if(ilay!=lay) continue;
      int neighborCts( jt.second->getNeighborCounts(hotWaferUV) );
      if(neighborCts<0) continue;
      float occ=(float)neighborCts;
      float ncells=(float)jt.second->getCells();
      occ/=ncells;
      neighborOccs.push_back(occ);
    }

    //sort and fill neighbor histos
    std::sort(neighborOccs.begin(),neighborOccs.end(),std::greater<int>());
    for(size_t i=0; i<neighborOccs.size(); i++) {
      hottestWaferH_[key][i+1]->Fill( neighborOccs[i] );
    }
  }
   
  //all done, reset counters
  for(auto &it : waferHistos_) it.second->resetCounters();  
}


//
void HGCOccupancyAnalyzer::analyzeDigis(std::string sd,edm::Handle<HGCalDigiCollection> &digiColl)
{
  //check inputs
  if(!digiColl.isValid()) return;

  const HGCalGeometry *geom=hgcGeometries_[sd];
  
  //get ddd constants to get cell properties
  const HGCalDDDConstants &ddd=geom->topology().dddConstants();

  //analyze hits
  for(auto &hit : *digiColl)
    {
      if(hit.size()==0) continue;
      int idx=2;

      //detid info
      HGCSiliconDetId detId(hit.id());
      int layer=detId.layer();

      if(detId.zside()<0) continue;

      int waferTypeL = ddd.waferType(detId);

      std::pair<int,int> waferUV=detId.waferUV();
      WaferEquivalentId_t ikey(std::make_tuple(sd,layer,waferUV.first,waferUV.second));
      std::pair<int,int> uvEq=uvEqMap_[ikey];
      WaferEquivalentId_t key(std::make_tuple(sd,layer,uvEq.first,uvEq.second));

      uint32_t rawData(hit.sample(idx).data() );
      bool isTOA( hit.sample(idx).getToAValid() );
      bool isTDC( hit.sample(idx).mode() );
      bool isBusy( isTDC && rawData==0 );
      waferHistos_[key]->count(waferUV.first,waferUV.second,rawData,isTOA,isTDC,isBusy);
    }
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCOccupancyAnalyzer);
