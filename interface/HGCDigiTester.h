#ifndef _HGCDigiTester_h_
#define _HGCDigiTester_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerBase.h"  
#include "SimCalorimetry/HGCalSimProducers/interface/HGCEEDigitizer.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

/**
   @class HGCDigiTester
   @short to test energy -> digi -> energy conversion
*/

class HGCDigiTester : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCDigiTester( const edm::ParameterSet& );
  ~HGCDigiTester() {}
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  std::unique_ptr<HGCEEDigitizer> digitizer_;
  uint32_t digitizationType_;
  HGCalSiNoiseMap scal_;
  uint32_t mipTarget_;
};
 

#endif
