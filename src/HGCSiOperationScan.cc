#include "UserCode/HGCElectronicsValidation/interface/HGCSiOperationScan.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"

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
  setPreassignedWafersFromCMSSW_(iConfig.getParameter<bool>("setPreassignedWafersFromCMSSW")),
  siType_(iConfig.getParameter<edm::ParameterSet>("siType")),
  aimMIPtoADC_(iConfig.getParameter<int>("aimMIPtoADC")),  
  encCommonNoiseSub_(iConfig.getParameter<double>("encCommonNoiseSub")),
  fluenceSF_(iConfig.getParameter<double>("fluenceSF")),
  savePadInfo_(iConfig.getParameter<bool>("savePadInfo"))
{  
  edm::Service<TFileService> fs;

  //start the summary ntuple with information for each si sensor type hypothesis
  data_ = fs->make<TTree>("data","data");
  data_->Branch("section",        &t_section,        "section/I");
  data_->Branch("layer",          &t_layer,          "layer/I");
  data_->Branch("waferU",         &t_waferU,         "waferU/I");
  data_->Branch("waferV",         &t_waferV,         "waferV/I");
  data_->Branch("waferX",         &t_waferX,         "waferX/F");
  data_->Branch("waferY",         &t_waferY,         "waferY/F");
  data_->Branch("waferPreChoice", &t_waferPreChoice, "waferPreChoice/I");
  data_->Branch("waferShape",     &t_waferShape,     "waferShape/I");
  data_->Branch("waferRot",       &t_waferRot,       "waferRot/I");
  data_->Branch("isHDWafer",      &t_isHDWafer,      "isHDWafer/O");
  data_->Branch("minf",           &t_minf,           "minf/F");
  data_->Branch("medf",           &t_medf,           "medf/F");
  data_->Branch("maxf",           &t_maxf,           "maxf/F");
  data_->Branch("npads",          &t_npads,          "npads/I");

  TString tag( (siType_.getParameter<std::string>("tag")).c_str() );
  data_->Branch("minsn_"+tag,  &t_minsn,  "minsn_"+tag+"/F");
  data_->Branch("q10sn_"+tag,  &t_q10sn,  "q10sn_"+tag+"/F");
  data_->Branch("medsn_"+tag,  &t_medsn,  "medsn_"+tag+"/F");
  data_->Branch("meds_"+tag,   &t_meds,   "meds_"+tag+"/F");
  data_->Branch("medn_"+tag,   &t_medn,   "medn_"+tag+"/F");
  data_->Branch("medencs_"+tag,&t_medencs,"medencs_"+tag+"/F");
  data_->Branch("medileak_"+tag,&t_medileak,"medileak_"+tag+"/F");

  if(savePadInfo_) {
    data_->Branch("padU",         t_padU,     "padU[npads]/I");
    data_->Branch("padV",         t_padV,     "padV[npads]/I");
    data_->Branch("padROC",       t_padROC,   "padROC[npads]/I");
    data_->Branch("padX",         t_padX,     "padX[npads]/F");
    data_->Branch("padY",         t_padY,     "padY[npads]/F");
    data_->Branch("padf",         t_padf,     "padf[npads]/F");
    data_->Branch("pads_"+tag,    t_pads,     "pads_"+tag+"[npads]/F");
    data_->Branch("padn_"+tag,    t_padn,     "padn_"+tag+"[npads]/F");
    data_->Branch("padsn_"+tag,   t_padsn,    "padsn_"+tag+"[npads]/F");
    data_->Branch("padencs_"+tag, t_padencs,  "padencs_"+tag+"[npads]/F");  
    data_->Branch("padileak_"+tag, t_padileak,  "padileak_"+tag+"[npads]/D");  
  }
  
  //start noise maps
  std::string doseMapURL(iConfig.getParameter<std::string>("doseMap"));  
  unsigned int doseMapAlgo(iConfig.getParameter<unsigned int>("doseMapAlgo"));
  std::vector<double> ileakParam(iConfig.getParameter<std::vector<double>>("ileakParam"));

  noiseMaps_["CEE"] = std::unique_ptr<HGCalSiNoiseMap<HGCSiliconDetId>>(new HGCalSiNoiseMap<HGCSiliconDetId>);
  noiseMaps_["CEE"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEE"]->setIleakParam(ileakParam);
  noiseMaps_["CEE"]->setENCCommonNoiseSubScale(encCommonNoiseSub_);
  noiseMaps_["CEE"]->setFluenceScaleFactor(fluenceSF_);

  noiseMaps_["CEH"] = std::unique_ptr<HGCalSiNoiseMap<HGCSiliconDetId>>(new HGCalSiNoiseMap<HGCSiliconDetId>);
  noiseMaps_["CEH"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEH"]->setIleakParam(ileakParam);
  noiseMaps_["CEH"]->setENCCommonNoiseSubScale(encCommonNoiseSub_);
  noiseMaps_["CEH"]->setFluenceScaleFactor(fluenceSF_);

  //parse u-v equivalence map file and start the layer operation map
  edm::FileInPath uvmapF(iConfig.getParameter<std::string>("uvmapfile")); 
  std::ifstream inF(uvmapF.fullPath());  
  cellOp_t baseCellOp;
  waferOp_t baseWaferOp;
  waferPos_t baseWaferPos;
  waferChoice_t baseWaferChoice;
  waferGeom_t baseWaferGeom;
  layerOp_t baseLayerOp;

  //converts the shape string to an integer
  std::map<std::string,int> shapeToInt;
  shapeToInt["F"]=0;
  shapeToInt["a"]=1;
  shapeToInt["b"]=2;
  shapeToInt["c"]=3;
  shapeToInt["d"]=4;
  shapeToInt["dm"]=5;
  shapeToInt["gm"]=6;

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
    if(sens==200) wafType=1; // 200, 1 cm^2
    if(sens==300) wafType=2; // 300, 1 cm^2
    int waferShape( shapeToInt[tokens[1].c_str()] );
    int waferRot(atoi(tokens[5].c_str()));
    float waferX(atof(tokens[3].c_str()));
    float waferY(atof(tokens[4].c_str()));
    int waferU(atoi(tokens[6].c_str()));
    int waferV(atoi(tokens[7].c_str()));
                 
    layKey_t layKey(subdet,ilay);
    if(baseLayerOp.find(layKey)==baseLayerOp.end()){
      baseLayerOp[layKey]     = baseWaferOp;
      waferPos_[layKey]       = baseWaferPos;
      waferPreChoice_[layKey] = baseWaferChoice;
      waferGeom_[layKey]      = baseWaferGeom;
    }
    waferKey_t waferKey(waferU,waferV);
    baseLayerOp[layKey][waferKey]=baseCellOp;
    waferPos_[layKey][waferKey]=std::pair<double,double>(waferX,waferY);
    waferPreChoice_[layKey][waferKey]=wafType;
    waferGeom_[layKey][waferKey]=waferKey_t(waferShape,waferRot);

    //start the helper map to save uv of the pads belonging to each wafer
    if(layerCellUVColl_.find(layKey)==layerCellUVColl_.end()) {
      std::map<waferKey_t,std::vector<waferKey_t> > layerBase;
      layerCellUVColl_[layKey]=layerBase;
      std::map<waferKey_t,std::vector<std::pair<double,double> > > xyBase;
      layerCellXYColl_[layKey]=xyBase;
      std::map<waferKey_t,std::vector<int> > rocBase;
      layerCellROCColl_[layKey]=rocBase;
    }
    if(layerCellUVColl_[layKey].find(waferKey)==layerCellUVColl_[layKey].end()) {
      std::vector<waferKey_t> waferBase;
      layerCellUVColl_[layKey][waferKey]=waferBase;
      std::vector<std::pair<double,double> > xyBase;
      layerCellXYColl_[layKey][waferKey]=xyBase;
      std::vector<int> rocBase;
      layerCellROCColl_[layKey][waferKey]=rocBase;
    }
  }

  inF.close();
}

//
HGCSiOperationScan::~HGCSiOperationScan()
{
}

//
void HGCSiOperationScan::endJob()
{
  std::string tag( siType_.getParameter<std::string>("tag") );

  for(auto lit : layerCellUVColl_) {

    layKey_t layKey(lit.first);
    std::string subdet(layKey.first);
    int sdcode=0;
    if(subdet=="CEH") sdcode=1;
    
    for(auto wit : lit.second ) {

      waferKey_t waferKey(wit.first);
      std::vector< waferKey_t > cellUVs(wit.second);
      std::vector< std::pair<double,double> > cellXYs(layerCellXYColl_[layKey][waferKey]);
      std::vector< int> cellROCs(layerCellROCColl_[layKey][waferKey]);

      if(cellUVs.size()==0) {
        std::cout << "Missing (u,v)=(" << waferKey.first << "," << waferKey.second << ") " 
                  << " " << subdet << " " << layKey.second << std::endl;
        continue;
      }
      
      //fill wafer header
      t_section        = sdcode;
      t_layer          = layKey.second;
      t_waferU         = waferKey.first;
      t_waferV         = waferKey.second;
      t_waferPreChoice = waferPreChoice_[layKey][waferKey];
      t_waferShape     = waferGeom_[layKey][waferKey].first;
      t_waferRot       = waferGeom_[layKey][waferKey].second;
      t_isHDWafer      = (t_waferPreChoice==0);
      t_npads          = cellUVs.size();
      t_waferX         = waferPos_[layKey][waferKey].first;
      t_waferY         = waferPos_[layKey][waferKey].second;      
      
      //reset value arrays
      memset(t_padU,     0, t_npads);
      memset(t_padV,     0, t_npads);
      memset(t_padROC,   0, t_npads);
      memset(t_padX,     0, t_npads);
      memset(t_padY,     0, t_npads);
      memset(t_padf,     0, t_npads);
      memset(t_pads,     0, t_npads);
      memset(t_padn,     0, t_npads);        
      memset(t_padsn,    0, t_npads); 
      memset(t_padencs,  0, t_npads); 
      memset(t_padileak, 0, t_npads); 
      
      const cellOp_t &cellOp=layerOp_[layKey][waferKey];
      for(int icell=0; icell<t_npads; icell++) {
        waferKey_t padUV=cellUVs[icell];
        std::pair<double,double> padXY=cellXYs[icell];
        t_padU[icell]    = padUV.first;
        t_padV[icell]    = padUV.second;
        t_padROC[icell]  = cellROCs[icell];
        t_padX[icell]    = padXY.first;
        t_padY[icell]    = padXY.second;                     
        t_padf[icell]    = cellOp[icell].fluence;          
        t_pads[icell]    = cellOp[icell].mipfC;
        t_padn[icell]    = cellOp[icell].core.noise;          
        t_padencs[icell] = cellOp[icell].enc_s; 
        t_padileak[icell] = cellOp[icell].ileak;
        t_padsn[icell]   = t_padn[icell] > 0 ? t_pads[icell]/t_padn[icell] : -1;
      }
      
      t_minf     = TMath::MinElement<float>(t_npads,t_padf);
      t_medf     = TMath::Median<float>    (t_npads,t_padf);
      t_maxf     = TMath::MaxElement<float>(t_npads,t_padf);

      t_minsn    = TMath::MinElement<float>(t_npads,t_padsn);
     
      //for the 10% quantile we need to order the pads
      //use another vector in order not to mess up with the one which is going to be stored
      std::vector<float> ordered_t_padsn(t_padsn,t_padsn+t_npads);
      std::sort(ordered_t_padsn.begin(),ordered_t_padsn.begin()+t_npads);      
      size_t iq    = TMath::Floor(0.1*t_npads);
      t_q10sn    = ordered_t_padsn[iq];

      t_medsn    = TMath::Median<float>(t_npads,t_padsn);
      t_meds     = TMath::Median<float>(t_npads,t_pads);
      t_medn     = TMath::Median<float>(t_npads,t_padn);        
      t_medencs  = TMath::Median<float>(t_npads,t_padencs);        
      t_medileak = TMath::Median<double>(t_npads,t_padileak);        

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
  HGCSiliconDetIdToROC d2roc;

  std::string tag( siType_.getParameter<std::string>("tag") );
  double mipEqfC( siType_.getParameter<double>("mipEqfC") );
  double cellVol( siType_.getParameter<double>("cellVol") );
  double cellCap( siType_.getParameter<double>("cellCap") );
  std::vector<double> cceParam( siType_.getParameter<std::vector<double> >("cceParam") );

  for(auto &it : hgcGeometries_ )
    {

      //get ddd constants to get cell properties     
      const HGCalDDDConstants &ddd=it.second->topology().dddConstants();

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
        double radius = sqrt(std::pow(pt.x(), 2) + std::pow(pt.y(), 2));  //in cm
        HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t gain(HGCalSiNoiseMap<HGCSiliconDetId>::AUTO);
        
        layerCellUVColl_[layKey][waferUV].push_back( cellUV );
        layerCellXYColl_[layKey][waferUV].push_back( std::pair<double,double>(pt.x(),pt.y()) );
        layerCellROCColl_[layKey][waferUV].push_back( d2roc.getROCNumber(detId) );
               
        //override default assignment with the one from CMSSW
        if(setPreassignedWafersFromCMSSW_){
          int waferTypeL = std::get<0>(ddd.waferType(detId));
          waferPreChoice_[layKey][waferUV] = waferTypeL;
        }

        layerOp_[layKey][waferUV].push_back(
                                            noiseMaps_[it.first]->getSiCellOpCharacteristics(cellCap,cellVol,mipEqfC,cceParam,
                                                                                             subdet,layer,radius,
                                                                                             gain,aimMIPtoADC_) );
      }
    }    
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSiOperationScan);
