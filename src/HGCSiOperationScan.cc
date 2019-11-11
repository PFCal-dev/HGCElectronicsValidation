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
  geoCEH_("HGCalHESiliconSensitive")
{  
  edm::Service<TFileService> fs;
  TString vars("section:layer:u:v:x:y:ncells");
  vars+=":minf:medf:maxf";
  vars+=":minsn_fine:q10sn_fine:medsn_fine:meds_fine:medn_fine";
  vars+=":minsn_thin:q10sn_thin:medsn_thin:meds_thin:medn_thin";
  vars+=":minsn_thick:q10sn_thick:medsn_thick:meds_thick:medn_thick";
  summaryTuple_=fs->make<TNtuple>("data","data",vars);

  std::string doseMapURL(iConfig.getParameter<std::string>("doseMap"));
  unsigned int doseMapAlgo(iConfig.getParameter<unsigned int>("doseMapAlgo"));
  std::vector<double> ileakParam(iConfig.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
  std::vector<double> cceParamFine(iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"));
  std::vector<double> cceParamThin(iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"));
  std::vector<double> cceParamThick(iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));

  noiseMaps_["CEE"] = std::unique_ptr<HGCalSiNoiseMap>(new HGCalSiNoiseMap);
  noiseMaps_["CEE"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEE"]->setIleakParam(ileakParam);
  noiseMaps_["CEE"]->setCceParam(cceParamFine, cceParamThin, cceParamThick);

  noiseMaps_["CEH"] = std::unique_ptr<HGCalSiNoiseMap>(new HGCalSiNoiseMap);
  noiseMaps_["CEH"]->setDoseMap(doseMapURL,doseMapAlgo);
  noiseMaps_["CEH"]->setIleakParam(ileakParam);
  noiseMaps_["CEH"]->setCceParam(cceParamFine, cceParamThin, cceParamThick);

  //parse u-v equivalence map file and start the layer operation map
  edm::FileInPath uvmapF("UserCode/HGCElectronicsValidation/data/wafer_pos.dat");
  std::ifstream inF(uvmapF.fullPath());  
  cellOp_t baseCellOp;
  waferOp_t baseWaferOp;
  waferPos_t baseWaferPos;
  layerOp_t baseLayerOp;
  while(inF) {

    std::string buf;
    getline(inF,buf);

    std::stringstream ss(buf);
    std::vector<std::string> tokens;
    while (ss >> buf)
      if(buf.size()>0)
        tokens.push_back(buf);
    if(tokens.size()<11) continue;
  
    std::string subdet(tokens[0]);
    int ilay(atoi(tokens[1].c_str()));
    float waferX(atof(tokens[3].c_str()));
    float waferY(atof(tokens[4].c_str()));
    int waferU(atoi(tokens[8].c_str()));
    int waferV(atoi(tokens[9].c_str()));

    layKey_t layKey(subdet,ilay);
    if(baseLayerOp.find(layKey)==baseLayerOp.end()){
      baseLayerOp[layKey]=baseWaferOp;
      waferPos_[layKey]=baseWaferPos;
    }
    waferKey_t waferKey(waferU,waferV);
    baseLayerOp[layKey][waferKey]=baseCellOp;
    waferPos_[layKey][waferKey]=std::pair<float,float>(waferX,waferY);
  }
  inF.close();

  cout << "Base layer operation has" << baseLayerOp.size() << " layer keys" << endl;

  //replicate layer operation map for the different modes
  layerOpColl_["fine"]=baseLayerOp;
  layerOpColl_["thin"]=baseLayerOp;
  layerOpColl_["thick"]=baseLayerOp;
}

//
HGCSiOperationScan::~HGCSiOperationScan()
{
}

//
void HGCSiOperationScan::endJob()
{
  Float_t *summaryVals=new Float_t(25);

  for(layerOp_t::iterator lit=layerOpColl_["fine"].begin();
      lit!=layerOpColl_["fine"].end();
      lit++) {

    layKey_t layKey(lit->first);
    std::string subdet(layKey.first);
    int sdcode=0;
    if(subdet=="CEH") sdcode=1;
    waferOp_t waferOp(lit->second);
    for(waferOp_t::iterator wit=waferOp.begin();
        wit!=waferOp.end();
        wit++) {
      waferKey_t waferKey(wit->first);
      cellOp_t cellOp(wit->second);
      cellOp_t cellOpThin=layerOpColl_["thin"][layKey][waferKey];
      cellOp_t cellOpThick=layerOpColl_["thick"][layKey][waferKey];

      size_t ncells(cellOp.size());
      double fluence[ncells], 
        snfine[ncells], sfine[ncells], nfine[ncells],  
        snthin[ncells], sthin[ncells], nthin[ncells], 
        snthick[ncells], sthick[ncells], nthick[ncells];
      for(size_t icell=0; icell<ncells; icell++) {
        fluence[icell] = cellOp[icell].fluence;
        sfine[icell]   = cellOp[icell].cce      * cellOp[icell].mipfC;
        nfine[icell]   = cellOp[icell].noise;
        snfine[icell]  = sfine[icell] / nfine[icell];
        sthin[icell]   = cellOp[icell].cce      * cellOp[icell].mipfC;
        nthin[icell]   = cellOp[icell].noise;
        snthin[icell]  = sthin[icell] / nthin[icell];
        sthick[icell]  = cellOp[icell].cce      * cellOp[icell].mipfC;
        nthick[icell]  = cellOp[icell].noise;
        snthick[icell] = sthick[icell] / nthick[icell];
      }     
      
      summaryVals[0]  = sdcode;
      summaryVals[1]  = layKey.second;
      summaryVals[2]  = waferKey.first;
      summaryVals[3]  = waferKey.second;
      summaryVals[4]  = waferPos_[layKey][waferKey].first;
      summaryVals[5]  = waferPos_[layKey][waferKey].second;
      summaryVals[6]  = ncells;
      summaryVals[7]  = TMath::MinElement<double>(ncells,fluence);
      summaryVals[8]  = TMath::Median<double>    (ncells,fluence);
      summaryVals[9]  = TMath::MaxElement<double>(ncells,fluence);
      
      size_t iq=TMath::Floor(0.1*ncells);

      std::sort(snfine,snfine+ncells);
      summaryVals[10]  = TMath::MinElement<double>(ncells,snfine);      
      summaryVals[11]  = snfine[iq];
      summaryVals[12]  = TMath::Median<double>    (ncells,snfine);
      summaryVals[13] = TMath::Median<double>    (ncells,sfine);
      summaryVals[14] = TMath::Median<double>    (ncells,nfine);

      std::sort(snthin,snthin+ncells);
      summaryVals[15] = TMath::MinElement<double>(ncells,snthin);
      summaryVals[16] = snthin[iq]; 
      summaryVals[17] = TMath::Median<double>    (ncells,snthin);
      summaryVals[18] = TMath::Median<double>    (ncells,sthin);
      summaryVals[19] = TMath::Median<double>    (ncells,nthin);

      std::sort(snthick,snthick+ncells);
      summaryVals[20] = TMath::MinElement<double>(ncells,snthick);
      summaryVals[21] = snthick[iq]; 
      summaryVals[22] = TMath::Median<double>    (ncells,snthick);
      summaryVals[23] = TMath::Median<double>    (ncells,sthick);
      summaryVals[24] = TMath::Median<double>    (ncells,nthick);

      summaryTuple_->Fill( summaryVals );
    }
  }
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

        //detid info
        int layer=detId.layer();
        std::pair<std::string,int> layKey(it.first,layer);
        std::pair<int,int> waferUV=detId.waferUV();
       
        HGCalSiNoiseMap::GainRange_t gainToSet(HGCalSiNoiseMap::AUTO);

        //fine hypothesis
        unsigned int masked_word( (HGCSiliconDetId::HGCalFine & HGCSiliconDetId::kHGCalTypeMask) << HGCSiliconDetId::kHGCalTypeOffset);
        rawId &= ~(HGCSiliconDetId::kHGCalTypeMask << HGCSiliconDetId::kHGCalTypeOffset);
        rawId |= (masked_word); 
        layerOpColl_["fine"][layKey][waferUV].push_back(
                                                        noiseMaps_[it.first]->getSiCellOpCharacteristics(HGCSiliconDetId(rawId), gainToSet, 10)
                                                        );

        //thin hypothesis
        masked_word=( (HGCSiliconDetId::HGCalCoarseThin & HGCSiliconDetId::kHGCalTypeMask) << HGCSiliconDetId::kHGCalTypeOffset);
        rawId &= ~(HGCSiliconDetId::kHGCalTypeMask << HGCSiliconDetId::kHGCalTypeOffset);
        rawId |= (masked_word); 
        layerOpColl_["thin"][layKey][waferUV].push_back(
                                                        noiseMaps_[it.first]->getSiCellOpCharacteristics(HGCSiliconDetId(rawId), gainToSet, 10)
                                                        );

        //thick hypothesis
        masked_word=( (HGCSiliconDetId::HGCalCoarseThick & HGCSiliconDetId::kHGCalTypeMask) << HGCSiliconDetId::kHGCalTypeOffset);
        rawId &= ~(HGCSiliconDetId::kHGCalTypeMask << HGCSiliconDetId::kHGCalTypeOffset);
        rawId |= (masked_word); 
        layerOpColl_["thick"][layKey][waferUV].push_back(
                                                         noiseMaps_[it.first]->getSiCellOpCharacteristics(HGCSiliconDetId(rawId), gainToSet, 10)
                                                         );
      }
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSiOperationScan);
