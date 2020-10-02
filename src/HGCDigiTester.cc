#include "UserCode/HGCElectronicsValidation/interface/HGCDigiTester.h"

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
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <fstream>
#include <iostream>

#include <cassert>

using namespace std;


//
// PLUGIN IMPLEMENTATION
//


//
HGCDigiTester::HGCDigiTester( const edm::ParameterSet &iConfig ) 
{ 

  //configure digitizer
  edm::ParameterSet cfg(iConfig.getParameter<edm::ParameterSet>("hgceeDigitizer"));
  digitizer_ = std::make_unique<HGCEEDigitizer>(cfg);
  digitizationType_ = cfg.getParameter<uint32_t>("digitizationType");

  //configure noise map
  edm::ParameterSet digiCfg(cfg.getParameter<edm::ParameterSet>("digiCfg"));
  scal_.setDoseMap(digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<std::string>("doseMap"),
                   digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<uint32_t>("scaleByDoseAlgo"));
  scal_.setFluenceScaleFactor(digiCfg.getParameter<edm::ParameterSet>("noise_fC").getParameter<double>("scaleByDoseFactor"));
  scal_.setIleakParam(digiCfg.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
  scal_.setCceParam(digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"),
                    digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"),
                    digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));
  mipTarget_=digiCfg.getParameter<edm::ParameterSet>("feCfg").getParameter<uint32_t>("targetMIPvalue_ADC");
}

//
void HGCDigiTester::endJob()
{
}

//
void HGCDigiTester::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //read geometry
  edm::ESHandle<HGCalGeometry> geoHandle;
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandle);
  const HGCalGeometry *geo=geoHandle.product();


  //get a valid DetId from the geometry
  DetId rawId(geo->getValidDetIds().begin()->rawId());
  std::unordered_set<DetId> validIds;  
  validIds.insert(rawId);
  HGCSiliconDetId cellId(rawId);
  const auto &xy = (geo->topology()).dddConstants().locateCell(cellId.layer(), cellId.waferU(), cellId.waferV(), cellId.cellU(), cellId.cellV(), true, true);
  std::cout << cellId.rawId() << " @ " << xy.first << " " << xy.second << " lay=" << cellId.layer() << std::endl;

  //re-config noise map and retrieve si-operation mode for this detId
  scal_.setGeometry(geo, HGCalSiNoiseMap::AUTO, mipTarget_);
  HGCalSiNoiseMap::SiCellOpCharacteristics siop=scal_.getSiCellOpCharacteristics(cellId);
  HGCalSiNoiseMap::GainRange_t gain((HGCalSiNoiseMap::GainRange_t)siop.core.gain);
  std::cout <<"\t lsb ADC=" << scal_.getLSBPerGain()[gain] << " max ADC=" << scal_.getMaxADCPerGain()[gain] << std::endl;

  //prepare a sim hit data accumulator to be filled with charge-only information
  hgc::HGCSimHitDataAccumulator simData;
  auto simIt = simData.emplace(rawId,hgc::HGCCellInfo()).first;

  //prepare the inputs for the digitization
  auto digiResult = std::make_unique<HGCalDigiCollection>();
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine *engine = &rng->getEngine(iEvent.streamID());

  //loop to digitize different values
  for(float q=0; q<10000; q+=100) { 
    (simIt->second).hit_info[0][0]=q;
    digitizer_->runSimple(digiResult,simData,geo,validIds,engine);
    
    //cout << q << " " << ((*digiResult)[0])[0].data() << endl;
  }



}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
