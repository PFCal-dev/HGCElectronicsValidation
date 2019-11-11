#ifndef _HGCSiOperationScan_h_
#define _HGCSiOperationScan_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"

#include <map>
#include <vector>
#include <string>

#include "TString.h"
#include "TH1F.h"
#include "TNtuple.h"


/**
   @class HGCSiOperationScan
   @short an EDAnalyzer to create a summary of the operation mode
*/

class HGCSiOperationScan : public edm::EDAnalyzer 
{

 public:

  explicit HGCSiOperationScan( const edm::ParameterSet& );
  ~HGCSiOperationScan();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  
  typedef std::pair<std::string,int> layKey_t;
  typedef std::pair<int,int> waferKey_t;
  typedef std::vector<HGCalSiNoiseMap::SiCellOpCharacteristics> cellOp_t;
  typedef std::map<waferKey_t,cellOp_t> waferOp_t;
  typedef std::map<waferKey_t, std::pair<float,float> > waferPos_t;
  typedef std::map<layKey_t, waferOp_t> layerOp_t;
  std::map<std::string, layerOp_t> layerOpColl_;
  std::map<layKey_t,waferPos_t> waferPos_;
  
  std::string geoCEE_,geoCEH_;
  std::map<std::string,const HGCalGeometry *> hgcGeometries_;
  std::map<std::string, std::unique_ptr<HGCalSiNoiseMap>> noiseMaps_;
  
  TNtuple *summaryTuple_;
};
 

#endif
