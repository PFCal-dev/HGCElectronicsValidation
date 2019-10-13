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
#include "TVector2.h"

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

  ofstream wafer_pos("wafer_pos.dat");

  Int_t colors1[5]={kRed,kRed-7,kRed-10,kRed-8,kRed+4};
  Int_t colors2[5]={kGreen,kGreen+3, kGreen-5, kGreen-10,kYellow-10};
  Int_t colors3[5]={kAzure,kAzure-5, kAzure+2, kAzure+8, kCyan-10};

  //init histograms
  edm::Service<TFileService> fs;
  std::map<TString, TGraph2D *> histos;
  for(auto &it : hgcGeometries_ )
    {

      std::map<int, std::map< std::pair<int,int>, WaferEquivalentInfo_t > >uvEqMap;
      std::map< std::pair<int,int>, WaferEquivalentInfo_t > uvTemplate;

      int subdet(it.first=="CEE" ? 0 : 1);
      const std::vector<DetId> &validDetIds = it.second->getValidDetIds();
      const HGCalDDDConstants &ddd=it.second->topology().dddConstants();      

      for(auto &didIt : validDetIds) {

        HGCSiliconDetId detId(didIt.rawId());

        if(detId.zside()<0) continue;

        int layer=detId.layer();
        int waferTypeL = ddd.waferType(detId);
        std::pair<int,int> waferUV=detId.waferUV();        
        int ncells= ddd.numberCellsHexagon(layer,waferUV.first,waferUV.second,true);

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

        std::pair<double, double> xy=ddd.waferPosition(waferUV.first,waferUV.second,true);
        double radius=hypot(xy.first,xy.second);
        double phi=atan2(xy.second,xy.first);
        double z(ddd.waferZ(layer,true));
        double eta(TMath::ASinH(z/radius));

        if(uvEqMap.find(layer)==uvEqMap.end()) uvEqMap[layer]=uvTemplate;
        std::pair<int,int> uvkey(waferTypeL,trunc(100*radius));
        if(uvEqMap[layer].find(uvkey)==uvEqMap[layer].end()) {
          WaferEquivalentInfo_t newWafer;
          newWafer.layer=layer;
          newWafer.phi=phi;
          newWafer.radius=radius;
          newWafer.z=z;
          newWafer.eta=eta;
          newWafer.x=xy.first;
          newWafer.y=xy.second;
          newWafer.ncells=ncells;
          newWafer.u=waferUV.first;
          newWafer.v=waferUV.second;
          uvEqMap[layer][uvkey]=newWafer;
        }
        uvEqMap[layer][uvkey].uvList.insert(waferUV);
        if(phi>0 && phi<uvEqMap[layer][uvkey].phi ) {
          uvEqMap[layer][uvkey].phi=phi;
          uvEqMap[layer][uvkey].x=xy.first;
          uvEqMap[layer][uvkey].y=xy.second;
          uvEqMap[layer][uvkey].u=waferUV.first;
          uvEqMap[layer][uvkey].v=waferUV.second;
        }
      }

      //dump info to output
      for(auto &jt : uvEqMap) {
        for(auto &kt : jt.second) {
          wafer_pos.width(5);
          wafer_pos << it.first << " " << jt.first << " ";
          wafer_pos.width(5);
          wafer_pos << kt.first.first << " ";
          wafer_pos.width(12);
          wafer_pos << kt.second.radius << " " << kt.second.x << " " << kt.second.y << " " << kt.second.phi << " "
                    << kt.second.z << " " << kt.second.eta << " ";
          wafer_pos.width(2);
          wafer_pos << kt.second.u << " " << kt.second.v << " "
                    << kt.second.ncells << " " << kt.second.uvList.size() << " ";
          for(auto &lt : kt.second.uvList){
            wafer_pos.width(2);
            wafer_pos << lt.first << " " << lt.second << " ";
          }
          wafer_pos << endl;
        }
      }

    }
  
  wafer_pos.close();
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
