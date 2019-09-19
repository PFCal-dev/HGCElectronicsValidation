#include "UserCode/HGCElectronicsValidation/interface/HGCGeometryScan.h"

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

#include "TGraph2D.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <fstream>
#include <iostream>

using namespace std;


//
// PLUGIN IMPLEMENTATION
//


//
HGCGeometryScan::HGCGeometryScan( const edm::ParameterSet &iConfig ) :   
  geoCEE_("HGCalEESensitive"),
  geoCEH_("HGCalHESiliconSensitive")
{  
  //parse u-v equivalence map file
  edm::FileInPath uvmapF("UserCode/HGCElectronicsValidation/data/uvequiv.dat");
  std::ifstream inF(uvmapF.fullPath());  
  int u,v,sec,ueq,veq;
  while(inF >> u >> v >> sec >> ueq >> veq ) {
    std::pair<int,int> key(u,v),keyeq(ueq,veq);
    uvEqSet_.insert(keyeq);
    uvEqMap_[key] = keyeq;
    uvSectorMap_[key]=sec;
  }
}

//
HGCGeometryScan::~HGCGeometryScan()
{
}

//
void HGCGeometryScan::endJob()
{
}

//
void HGCGeometryScan::prepareAnalysis()
{

  ofstream cell_count("cell_count.dat");

  Int_t colors1[5]={kRed,kRed-7,kRed-10,kRed-8,kRed+4};
  Int_t colors2[5]={kGreen,kGreen+3, kGreen-5, kGreen-10,kYellow-10};
  Int_t colors3[5]={kAzure,kAzure-5, kAzure+2, kAzure+8, kCyan-10};

  //init histograms
  edm::Service<TFileService> fs;
  std::map<TString, TGraph2D *> histos;
  for(auto &it : hgcGeometries_ )
    {

      std::map<std::pair<int,int>, std::vector<int> > cellCounter;

      int subdet(it.first=="CEE" ? 0 : 1);
      const std::vector<DetId> &validDetIds = it.second->getValidDetIds();
      const HGCalDDDConstants &ddd=it.second->topology().dddConstants();      

      for(auto &didIt : validDetIds) {

        HGCSiliconDetId detId(didIt.rawId());

        if(detId.zside()<0) continue;

        int layer=detId.layer();
        int waferTypeL = ddd.waferType(detId);
        std::pair<int,int> waferUV=detId.waferUV();        

        TString key(Form("sd%d_lay%d_xy",subdet,layer));
        if(histos.find(key)==histos.end()) {
          histos[key]=fs->make<TGraph2D>();
          histos[key]->SetName(key);
        }
        
        Int_t cIdx((abs(waferUV.first)*abs(waferUV.second))%5);
        Int_t ci=colors1[cIdx];
        if(waferTypeL==2) ci=colors2[cIdx];
        if(waferTypeL==3) ci=colors3[cIdx];

        GlobalPoint pt=it.second->getPosition(detId);
        int npts=histos[key]->GetN();
        histos[key]->SetPoint(npts,pt.x(),pt.y(),ci);

        //save information just for a sector
        for(auto &iuv : uvEqSet_){
          //std::pair<int,int> iuv=uv.first;
          if(iuv != waferUV) continue;

          if( cellCounter.find(iuv) == cellCounter.end() ) {
            unsigned int nlay(ddd.layers(true));
            cellCounter[iuv]=std::vector<int>(nlay,0);
          }
          cellCounter[iuv][layer-1] += 1;         
        }
      }

      for(auto &it : cellCounter) {

        std::pair<int,int> uv(it.first);
        std::vector<int> counts(it.second);
        for(size_t ilay=0; ilay<counts.size(); ilay++) {
          cell_count << subdet << " " << ilay+1 << " " << uv.first << " " << uv.second << " " << counts[ilay] << endl; 
        }

      }
        
    }
  
  cell_count.close();
}

  
//
void HGCGeometryScan::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //read geometry from event setup
  edm::ESHandle<HGCalGeometry> ceeGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEE_,ceeGeoHandle);
  hgcGeometries_["CEE"]=ceeGeoHandle.product();
  edm::ESHandle<HGCalGeometry> cehGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEH_,cehGeoHandle);
  hgcGeometries_["CEH"]=cehGeoHandle.product();

  prepareAnalysis();

}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCGeometryScan);
