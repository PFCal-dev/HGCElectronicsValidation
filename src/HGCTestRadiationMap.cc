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


/**
   @class HGCTestRadiationMap
   @short an EDAnalyzer to test features in the radiation map class
*/

class HGCTestRadiationMap : public edm::EDAnalyzer 
{

 public:

  explicit HGCTestRadiationMap( const edm::ParameterSet& ) {}
  ~HGCTestRadiationMap(){}
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob() {}

 private:
};
 
//
void HGCTestRadiationMap::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup) {

  HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t gain(HGCalSiNoiseMap<HGCSiliconDetId>::AUTO);
  int aimMIPtoADC(10);
  std::string doseMapURL("SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt");
  std::vector<double> ileakParam={0.993,-42.668};
  double encCommonNoiseSub(1.0);
  double fluenceSF(1.0);
  std::vector<double> cceParsFine={1.5e+15, -3.00394e-17, 0.318083};
  std::vector<double> cceParsThin={1.5e+15, -3.09878e-16, 0.211207};
  std::vector<double> cceParsThick={6e+14,   -7.96539e-16, 0.251751};

  std::string geoNames[2]={"HGCalEESensitive","HGCalHESiliconSensitive"};
  for(size_t i=0; i<sizeof(geoNames)/sizeof(std::string); i++) {

    std::cout << "***************************************************" << std::endl
              << "Checking sector equivalence for " << geoNames[i] << std::endl;

    //get geometry and valid detIds
    edm::ESHandle<HGCalGeometry> geoHandle;
    iSetup.get<IdealGeometryRecord>().get(geoNames[i],geoHandle);
    const HGCalGeometry *geo=geoHandle.product();
    const std::vector<DetId> &validDetIds = geo->getValidDetIds();

    //start alternative noise maps
    HGCalSiNoiseMap<HGCSiliconDetId> sectorEqMap;
    uint32_t doseMapAlgo(0);
    sectorEqMap.setDoseMap(doseMapURL,doseMapAlgo);
    sectorEqMap.setIleakParam(ileakParam);
    sectorEqMap.setENCCommonNoiseSubScale(encCommonNoiseSub);
    sectorEqMap.setFluenceScaleFactor(fluenceSF);
    sectorEqMap.setGeometry(geo);
    sectorEqMap.setCceParam(cceParsFine,cceParsThin,cceParsThick);

    HGCalSiNoiseMap<HGCSiliconDetId> noSectorEqMap;
    doseMapAlgo |= (1 << HGCalSiNoiseMap<HGCSiliconDetId>::CACHEDOP);
    noSectorEqMap.setDoseMap(doseMapURL,doseMapAlgo);
    noSectorEqMap.setIleakParam(ileakParam);
    noSectorEqMap.setENCCommonNoiseSubScale(encCommonNoiseSub);
    noSectorEqMap.setFluenceScaleFactor(fluenceSF);
    noSectorEqMap.setGeometry(geo);
    noSectorEqMap.setCceParam(cceParsFine,cceParsThin,cceParsThick);

    //loop over detIds and compare
    int nfaulty(0);
    for(auto &didIt : validDetIds) {
      HGCalSiNoiseMap<HGCSiliconDetId>::SiCellOpCharacteristics a(sectorEqMap.getSiCellOpCharacteristics(didIt,gain,aimMIPtoADC));
      HGCalSiNoiseMap<HGCSiliconDetId>::SiCellOpCharacteristics b(noSectorEqMap.getSiCellOpCharacteristics(didIt,gain,aimMIPtoADC));
      if(a.fluence!=b.fluence || a.ileak!=b.ileak || a.core.cce!=b.core.cce || a.core.gain!=b.core.gain) {
        nfaulty++;
        std::cout << std::hex << "[WARN] difference found for 0x" << didIt.rawId() << std::dec << std::endl;       
        break;
      }
    }
    
    std::cout << nfaulty << "/" << validDetIds.size() << " wrong detId equivalence assignments were found" << std::endl;
  }
}


//define this as a plug-in                                                                                              
DEFINE_FWK_MODULE(HGCTestRadiationMap);
