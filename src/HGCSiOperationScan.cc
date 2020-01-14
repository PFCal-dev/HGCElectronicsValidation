#include "UserCode/HGCElectronicsValidation/interface/HGCSiOperationScan.h"

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
#include <algorithm>

#include "TMath.h"

using namespace std;

//
// PLUGIN IMPLEMENTATION
//

//
HGCSiOperationScan::HGCSiOperationScan( const edm::ParameterSet &iConfig ) :   
  geoCEE_("HGCalEESensitive"),
  geoCEH_("HGCalHESiliconSensitive"),
  siTypes_(iConfig.getParameter<std::vector<edm::ParameterSet> >("siTypes")),
  aimMIPtoADC_(iConfig.getParameter<int>("aimMIPtoADC")),
  savePadInfo_(iConfig.getParameter<bool>("savePadInfo"))
{  
  edm::Service<TFileService> fs;

  //start the summary ntuple with information for each si sensor type hypothesis
  data_ = fs->make<TTree>("data","data");
  data_->Branch("section",        &t_section,        "section/I");
  data_->Branch("layer",          &t_layer,          "layer/I");
  data_->Branch("waferUV",        &t_waferUV,        "waferUV/I");
  data_->Branch("waferX",         &t_waferX,         "waferX/F");
  data_->Branch("waferY",         &t_waferY,         "waferY/F");
  data_->Branch("waferPreChoice", &t_waferPreChoice, "waferPreChoice/I");
  data_->Branch("minf",           &t_minf,           "minf/F");
  data_->Branch("medf",           &t_medf,           "medf/F");
  data_->Branch("maxf",           &t_maxf,           "maxf/F");
  data_->Branch("npads",          &t_npads,          "npads/I");
  if(savePadInfo_) {
    data_->Branch("padUV",    t_padUV,     "padUV[npads]/I");
    data_->Branch("padX",     t_padX,      "padX[npads]/F");
    data_->Branch("padY",     t_padY,      "padY[npads]/F");
    data_->Branch("padf",     t_padf,      "padf[npads]/F");
  }

  size_t nSiTypes(siTypes_.size());
  t_minsn.resize(nSiTypes);
  t_q10sn.resize(nSiTypes);
  t_medsn.resize(nSiTypes);
  t_meds.resize(nSiTypes);
  t_medn.resize(nSiTypes);

  for(size_t i=0; i<nSiTypes; i++) {

    edm::ParameterSet si( siTypes_[i] );
    TString tag( (si.getParameter<std::string>("tag")).c_str() );

    data_->Branch("minsn_"+tag,  &t_minsn[i], "minsn_"+tag+"/F");
    data_->Branch("q10sn_"+tag,  &t_q10sn[i], "q10sn_"+tag+"/F");
    data_->Branch("medsn_"+tag,  &t_medsn[i], "medsn_"+tag+"/F");
    data_->Branch("meds_"+tag,   &t_meds[i],  "meds_"+tag+"/F");
    data_->Branch("medn_"+tag,   &t_medn[i],  "medn_"+tag+"/F");

    if(savePadInfo_) {

      t_pads.push_back( new Float_t(500) );
      t_padn.push_back( new Float_t(500) );
      t_padsn.push_back( new Float_t(500) );

      data_->Branch("pads_"+tag,  t_pads[i],  "pads_"+tag+"[npads]/F");
      data_->Branch("padn_"+tag,  t_padn[i],  "padn_"+tag+"[npads]/F");
      data_->Branch("padsn_"+tag, t_padsn[i], "padsn_"+tag+"[npads]/F");
    }
  }
  
  //start noise maps
  std::string doseMapURL(iConfig.getParameter<std::string>("doseMap"));  
  unsigned int doseMapAlgo(iConfig.getParameter<unsigned int>("doseMapAlgo"));
  std::vector<double> ileakParam(iConfig.getParameter<std::vector<double>>("ileakParam"));

  noiseMaps_["CEE"] = std::unique_ptr<HGCalSiNoiseMap>(new HGCalSiNoiseMap);
  noiseMaps_["CEE"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEE"]->setIleakParam(ileakParam);

  noiseMaps_["CEH"] = std::unique_ptr<HGCalSiNoiseMap>(new HGCalSiNoiseMap);
  noiseMaps_["CEH"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEH"]->setIleakParam(ileakParam);

  //parse u-v equivalence map file and start the layer operation map
  edm::FileInPath uvmapF("UserCode/HGCElectronicsValidation/data/geomnew_corrected_360.txt");
  std::ifstream inF(uvmapF.fullPath());  
  cellOp_t baseCellOp;
  waferOp_t baseWaferOp;
  waferPos_t baseWaferPos;
  waferChoice_t baseWaferChoice;
  layerOp_t baseLayerOp;

  while(inF) {

    std::string buf;
    getline(inF,buf);

    std::stringstream ss(buf);
    std::vector<std::string> tokens;
    while (ss >> buf)
      if(buf.size()>0)
        tokens.push_back(buf);
    if(tokens.size()!=8) continue;
  
    int ilay(atoi(tokens[0].c_str()));
    std::string subdet(ilay<=28 ? "CEE" : "CEH");
    if(ilay>28) ilay-=28;
    int sens(atoi(tokens[2].c_str()));
    int wafType(0);
    if(sens==120) wafType=0; // 120, 0.5cm^2
    if(sens==200) wafType=3; // 200, 1 cm^2
    if(sens==300) wafType=5; // 300, 1 cm^2
    float waferX(atof(tokens[3].c_str()));
    float waferY(atof(tokens[4].c_str()));
    int waferU(atoi(tokens[6].c_str()));
    int waferV(atoi(tokens[7].c_str()));

    layKey_t layKey(subdet,ilay);
    if(baseLayerOp.find(layKey)==baseLayerOp.end()){
      baseLayerOp[layKey]     = baseWaferOp;
      waferPos_[layKey]       = baseWaferPos;
      waferPreChoice_[layKey] = baseWaferChoice;
    }
    waferKey_t waferKey(waferU,waferV);
    baseLayerOp[layKey][waferKey]=baseCellOp;
    waferPos_[layKey][waferKey]=std::pair<double,double>(waferX,waferY);
    waferPreChoice_[layKey][waferKey]=wafType;

    //start the helper map to save uv of the pads belonging to each wafer
    if(layerCellUVColl_.find(layKey)==layerCellUVColl_.end()) {
      std::map<waferKey_t,std::vector<waferKey_t> > layerBase;
      layerCellUVColl_[layKey]=layerBase;
    }
    if(layerCellUVColl_[layKey].find(waferKey)==layerCellUVColl_[layKey].end()) {
      std::vector<waferKey_t> waferBase;
      layerCellUVColl_[layKey][waferKey]=waferBase;
    }
  }

  inF.close();

  cout << "Base layer operation has " << baseLayerOp.size() << " layer keys" << endl;

  //replicate layer operation map for the different modes
  for(size_t i=0; i<nSiTypes; i++) {

    edm::ParameterSet si( siTypes_[i] );
    std::string tag( si.getParameter<std::string>("tag") );
    layerOpColl_[tag]=baseLayerOp;
  }
}

//
HGCSiOperationScan::~HGCSiOperationScan()
{
}

//
void HGCSiOperationScan::endJob()
{
  //helper arrays
  double fvals[500], svals[500], nvals[500], snvals[500];
  
  for(auto lit : layerCellUVColl_) {

    layKey_t layKey(lit.first);
    std::string subdet(layKey.first);
    int sdcode=0;
    if(subdet=="CEH") sdcode=1;

    for(auto wit : lit.second ) {

      waferKey_t waferKey(wit.first);
      std::vector< waferKey_t > cellUVs(wit.second);
      std::vector< std::pair<double,double> > cellXYs(layerCellXYColl_[layKey][waferKey]);

      t_section        = sdcode;
      t_layer          = layKey.second;
      t_waferUV        = waferKey.first*1000+waferKey.second;
      t_waferPreChoice = waferPreChoice_[layKey][waferKey];
      t_npads          = cellUVs.size();
      t_waferX         = waferPos_[layKey][waferKey].first;
      t_waferY         = waferPos_[layKey][waferKey].second;      

      for(size_t i=0; i<siTypes_.size(); i++) {

        std::string tag( siTypes_[i].getParameter<std::string>("tag") );
        const cellOp_t &cellOp=layerOpColl_[tag][layKey][waferKey];

        for(int icell=0; icell<t_npads; icell++) {

          fvals[icell]     = cellOp[icell].fluence;

          if(i==0 && savePadInfo_) {
            waferKey_t wafUV=cellUVs[icell];
            std::pair<double,double> wafXY=cellXYs[icell];
            t_padUV[icell] = wafUV.first*1000+wafUV.second;           
            t_padX[icell]  = wafXY.first;
            t_padY[icell]  = wafXY.second;
            t_padf[icell]  = fvals[icell];
          }

          svals[icell]      = cellOp[icell].cce * cellOp[icell].mipfC;
          t_pads[i][icell]  = svals[icell];
          nvals[icell]      = cellOp[icell].noise;
          t_padn[i][icell]  = nvals[icell];
          snvals[icell]     = svals[icell] / nvals[icell];
          t_padsn[i][icell] = snvals[icell];
        }
      
        //only needs to be computed for the first Si type
        if(i==0 && savePadInfo_) {
          t_minf = TMath::MinElement<double>(t_npads,fvals);
          t_medf = TMath::Median<double>    (t_npads,fvals);
          t_maxf = TMath::MaxElement<double>(t_npads,fvals);
        }

        size_t iq=TMath::Floor(0.1*t_npads);
        std::sort(snvals,snvals+t_npads);
        t_minsn[i] = TMath::MinElement<double>(t_npads,snvals);      
        t_q10sn[i] = snvals[iq];
        t_medsn[i] = TMath::Median<double>(t_npads,snvals);
        t_meds[i]  = TMath::Median<double>(t_npads,svals);
        t_medn[i]  = TMath::Median<double>(t_npads,nvals);
      
      }//end Si types

      //fill the ntuple
      data_->Fill();
    
    }//end wafer 
    
  }//end layer

}

  
//
void HGCSiOperationScan::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::ESHandle<HGCalGeometry> ceeGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEE_,ceeGeoHandle);
  hgcGeometries_["CEE"]=ceeGeoHandle.product();
  edm::ESHandle<HGCalGeometry> cehGeoHandle;
  iSetup.get<IdealGeometryRecord>().get(geoCEH_,cehGeoHandle);
  hgcGeometries_["CEH"]=cehGeoHandle.product();

  for(auto &it : hgcGeometries_ )
    {
      noiseMaps_[it.first]->setGeometry( it.second );
      const std::vector<DetId> &validDetIds = it.second->getValidDetIds();
      for(auto &didIt : validDetIds) {
        unsigned int rawId(didIt.rawId());
        HGCSiliconDetId detId(rawId);
        if(detId.zside()<0) continue;

        //decode detid info
        int layer( detId.layer() );
        int subdet( detId.subdet() );
        std::pair<std::string,int> layKey(it.first,layer);
        std::pair<int,int> waferUV=detId.waferUV();
        std::pair<int,int> cellUV=detId.cellUV();
        GlobalPoint pt=it.second->getPosition(detId);
        double radius2 = std::pow(pt.x(), 2) + std::pow(pt.y(), 2);  //in cm
        HGCalSiNoiseMap::GainRange_t gain(HGCalSiNoiseMap::AUTO);
        
        layerCellUVColl_[layKey][waferUV].push_back( cellUV );
        layerCellXYColl_[layKey][waferUV].push_back( std::pair<double,double>(pt.x(),pt.y()) );
                                                     
        //iterate over si types to test
        for(size_t i=0; i<siTypes_.size(); i++) {
          
          edm::ParameterSet si( siTypes_[i] );
          std::string tag( si.getParameter<std::string>("tag") );
          double mipEqfC( si.getParameter<double>("mipEqfC") );
          double cellVol( si.getParameter<double>("cellVol") );
          double cellCap( si.getParameter<double>("cellCap") );
          std::vector<double> cceParam( si.getParameter<std::vector<double> >("cceParam") );
          
          layerOpColl_[tag][layKey][waferUV].push_back(
                                                       noiseMaps_[it.first]->getSiCellOpCharacteristics(cellCap,cellVol,mipEqfC,cceParam,
                                                                                                        subdet,layer,radius2,
                                                                                                        gain,aimMIPtoADC_) );
        }
      }
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSiOperationScan);
