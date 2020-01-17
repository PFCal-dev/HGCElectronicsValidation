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
  typedef std::map<waferKey_t, waferKey_t> waferGeom_t;
  typedef std::map<waferKey_t, std::pair<double,double> > waferPos_t;
  typedef std::map<waferKey_t, int > waferChoice_t;
  typedef std::map<layKey_t, waferOp_t> layerOp_t;
  std::map<std::string, layerOp_t> layerOpColl_;
  std::map<layKey_t,waferPos_t> waferPos_;
  std::map<layKey_t,waferChoice_t> waferPreChoice_;
  std::map<layKey_t,waferGeom_t> waferGeom_;

  std::map<layKey_t,std::map<waferKey_t,std::vector<waferKey_t> > > layerCellUVColl_;
  std::map<layKey_t,std::map<waferKey_t,std::vector<std::pair<double,double> > > > layerCellXYColl_;

  //geometry
  std::string geoCEE_,geoCEH_;
  std::map<std::string,const HGCalGeometry *> hgcGeometries_;

  //radiation map
  std::map<std::string, std::unique_ptr<HGCalSiNoiseMap>> noiseMaps_;
  
  //si characteristics to scan
  std::vector<edm::ParameterSet> siTypes_;
  int aimMIPtoADC_;

  //summary tree
  TTree *data_;
  Int_t t_section,t_layer,t_waferU,t_waferV,t_waferPreChoice,t_waferShape,t_waferRot,t_npads;
  Float_t t_waferX,t_waferY,t_minf,t_medf,t_maxf;
  Bool_t t_isHDWafer;
  Int_t t_padU[500],t_padV[500];
  Float_t t_padX[500],t_padY[500],t_padf[500];
  std::vector<Float_t> t_minsn,t_q10sn,t_medsn,t_meds,t_medn;
  Float_t t_pads[25][500],t_padn[25][500],t_padsn[25][500];
  bool savePadInfo_;
};
 

#endif
