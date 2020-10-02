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
  //configure noise map
  edm::ParameterSet cfg(iConfig.getParameter<edm::ParameterSet>("hgceeDigitizer"));
  edm::ParameterSet digiCfg(cfg.getParameter<edm::ParameterSet>("digiCfg"));
  edm::ParameterSet feCfg(digiCfg.getParameter<edm::ParameterSet>("feCfg"));
  scal_.setDoseMap(digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<std::string>("doseMap"),
                   digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<uint32_t>("scaleByDoseAlgo"));
  scal_.setFluenceScaleFactor(digiCfg.getParameter<edm::ParameterSet>("noise_fC").getParameter<double>("scaleByDoseFactor"));
  scal_.setIleakParam(digiCfg.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
  scal_.setCceParam(digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"),
                    digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"),
                    digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));
  mipTarget_=feCfg.getParameter<uint32_t>("targetMIPvalue_ADC");

  //configure digitizer
  digitizer_ = std::make_unique<HGCEEDigitizer>(cfg);
  digitizationType_ = cfg.getParameter<uint32_t>("digitizationType");
  tdcLSB_=feCfg.getParameter<double>("tdcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("tdcNbits") );
  tdcOnset_fC_=feCfg.getParameter<double>("tdcOnset_fC");
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

  //re-config noise map and retrieve si-operation mode for this detId
  scal_.setGeometry(geo, HGCalSiNoiseMap::AUTO, mipTarget_);
  HGCSiliconDetId cellId(rawId);
  HGCalSiNoiseMap::SiCellOpCharacteristics siop=scal_.getSiCellOpCharacteristics(cellId);
  HGCalSiNoiseMap::GainRange_t gain((HGCalSiNoiseMap::GainRange_t)siop.core.gain);
  double adcLSB=scal_.getLSBPerGain()[gain];
  double cce=siop.core.cce;
  std::cout << "ADC lsb=" << adcLSB 
            << " TDC lsb=" << tdcLSB_ 
            << " noise=" << siop.core.noise 
            << " mip=" << siop.mipfC << std::endl;

  //prepare a sim hit data accumulator to be filled with charge-only information
  hgc::HGCSimHitDataAccumulator simData;
  simData.reserve(1);
  auto simIt = simData.emplace(rawId,hgc::HGCCellInfo()).first;

  //prepare the inputs for the digitization
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine *engine = &rng->getEngine(iEvent.streamID());

  std::ofstream ofile;
  ofile.open ("digitest.dat");
  ofile << "qsim mode ADC qrec qreccce" << std::endl;
  
  //loop to digitize different values
  float deltaq(0.1);
  for(float q=0; q<10000; q+=deltaq) { 
    
    if(q>10)   deltaq=1;
    if(q>100)  deltaq=10;
    if(q>1000) deltaq=100;
    
    for(int i=0; i<100; i++) {
      auto digiResult = std::make_unique<HGCalDigiCollection>();
      (simIt->second).hit_info[0][9]=q; //fill in-time index only
      digitizer_->run(digiResult,simData,geo,validIds,digitizationType_,engine);
      
      //if a digi was not produced move to next value
      if(digiResult->size()==0) continue;
      
      //read digi
      uint32_t mode=((*digiResult)[0])[2].mode();
      uint32_t adc=((*digiResult)[0])[2].data();
      
      //convert back to charge
      double qrec( adc*adcLSB );
      if(mode){
        qrec=(std::floor(tdcOnset_fC_/adcLSB)+1)*adcLSB +(adc+0.5)*tdcLSB_ ;
        //qrec=tdcOnset_fC_*0.5*tdcLSB_+-0.5*adcLSB+adc*tdcLSB_; //(adc+0.5)*tdcLSB_ ;
      }
      double qrec_cce(qrec/cce);
      
      ofile << q << " " << mode << " " << adc << " " << qrec << " " << qrec_cce << std::endl;
    }
  }
  
  ofile.close();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
