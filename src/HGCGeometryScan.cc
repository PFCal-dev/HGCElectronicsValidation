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
#include <iomanip>

using namespace std;


//
// PLUGIN IMPLEMENTATION
//


//
HGCGeometryScan::HGCGeometryScan( const edm::ParameterSet &iConfig ) :   
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>())
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

      const std::vector<DetId> &validDetIds = it.second->getValidDetIds();
      const HGCalDDDConstants &ddd=it.second->topology().dddConstants();      

      for(auto &didIt : validDetIds) {

        HGCSiliconDetId detId(didIt.rawId());

        if(detId.zside()<0) continue;

        int layer=detId.layer();
        std::pair<int,int> waferUV=detId.waferUV();        
        int waferTypeL = ddd.waferType(layer,waferUV.first,waferUV.second,false);
        
        TString key(Form("%s_lay%d_xy",it.first.c_str(),layer));
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
 
        //find rotation to equivalence sector
        std::pair<double, double> waferPos=ddd.waferPosition(detId,false);
        double waferPhi=atan2(waferPos.second,waferPos.first);
        double waferRadius=hypot(waferPos.first,waferPos.second);
        double z(ddd.waferZ(layer,true));
        double eta(TMath::ASinH(z/waferRadius));

        if(uvEqMap.find(layer)==uvEqMap.end()) uvEqMap[layer]=uvTemplate;
        if(uvEqMap[layer].find(waferUV)==uvEqMap[layer].end()) {
          WaferEquivalentInfo_t newWafer;
          newWafer.layer=layer;
          newWafer.x=waferPos.first;
          newWafer.y=waferPos.second;
          newWafer.phi=waferPhi;
          newWafer.radius=waferRadius;
          newWafer.z=z;
          newWafer.eta=eta;
          newWafer.ncells=0;
          newWafer.u=waferUV.first;
          newWafer.v=waferUV.second;
          uvEqMap[layer][waferUV]=newWafer;
          newWafer.waferType=waferTypeL;
        }
        uvEqMap[layer][waferUV].ncells+=1;
        uvEqMap[layer][waferUV].uvList.insert(waferUV);
      }

      //dump info to output
      for(auto &jt : uvEqMap) {
        for(auto &kt : jt.second) {

          kt.second.ncells /= kt.second.uvList.size();

          wafer_pos.width(5);
          wafer_pos << it.first << " " << jt.first << " ";
          wafer_pos.width(12);
          wafer_pos << std::setprecision(9);
          wafer_pos << kt.second.radius << " " << kt.second.x << " " << kt.second.y << " " << kt.second.phi << " "
                    << kt.second.z << " " << kt.second.eta << " ";
          wafer_pos << std::fixed;
          wafer_pos.width(2);
          wafer_pos << kt.second.u << " " << kt.second.v << " " << kt.second.ncells << " " << kt.second.waferType;
          wafer_pos << endl;
        }
      }

    }
  
  wafer_pos.close();
}

  
//
void HGCGeometryScan::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  //read geometry components for HGCAL
  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  hgcGeometries_["CEE"] = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty));
  hgcGeometries_["CEH"] = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty));
  prepareAnalysis();

}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCGeometryScan);
