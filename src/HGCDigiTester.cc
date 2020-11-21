#include "UserCode/HGCElectronicsValidation/interface/HGCDigiTester.h"
#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <fstream>
#include <iostream>

#include <cassert>

using namespace std;

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// PLUGIN IMPLEMENTATION
//


//
HGCDigiTester::HGCDigiTester( const edm::ParameterSet &iConfig ) 
  : simHitsCEE_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE")) ),
    simHitsCEH_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront")) ),
    simHitsCEHSci_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback")) ),
    digisCEE_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
    digisCEH_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEfront")) ),
    digisCEHSci_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEback")) ),
    genParticles_( consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles")) )
{   
  //configure noise map
  std::string digitizers[]={"hgceeDigitizer","hgcehDigitizer","hgcehsciDigitizer"};
  for(size_t i=0; i<3; i++) {
    edm::ParameterSet cfg(iConfig.getParameter<edm::ParameterSet>(digitizers[i]));
    edm::ParameterSet digiCfg(cfg.getParameter<edm::ParameterSet>("digiCfg"));
    edm::ParameterSet feCfg(digiCfg.getParameter<edm::ParameterSet>("feCfg"));
    if(i<2) {
      HGCalSiNoiseMap *scal=new HGCalSiNoiseMap;
      scal->setDoseMap(digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<std::string>("doseMap"),
                             digiCfg.getParameter<edm::ParameterSet>("noise_fC").template getParameter<uint32_t>("scaleByDoseAlgo"));
      scal->setFluenceScaleFactor(digiCfg.getParameter<edm::ParameterSet>("noise_fC").getParameter<double>("scaleByDoseFactor"));
      scal->setIleakParam(digiCfg.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
      scal->setCceParam(digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"),
                        digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"),
                        digiCfg.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));
      scal_.push_back( scal );
     }

    //other digitizer configs    
    vanilla_mipfC_[i]=feCfg.getParameter<double>("mipfC");
    mipTarget_[i]=feCfg.getParameter<uint32_t>("targetMIPvalue_ADC");
    tdcOnset_fC_[i]=feCfg.getParameter<double>("tdcOnset_fC");
    tdcLSB_[i]=feCfg.getParameter<double>("tdcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("tdcNbits") );
    vanilla_adcLSB_fC_[i] = feCfg.getParameter<double>("adcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("adcNbits") );
  }
  useVanillaCfg_ = iConfig.getParameter<bool>("useVanillaCfg");
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("hits","hits");
  event_=0;
  tree_->Branch("genergy",&genergy_,"genergy/F");
  tree_->Branch("gpt",&gpt_,"gpt/F");
  tree_->Branch("geta",&geta_,"geta/F");
  tree_->Branch("gphi",&gphi_,"gphi/F");
  tree_->Branch("gvradius",&gvradius_,"gvradius/F");
  tree_->Branch("gvz",&gvz_,"gvz/F");
  tree_->Branch("event",&event_,"event/I");
  tree_->Branch("qsim",&qsim_,"qsim/F");
  tree_->Branch("qrec",&qrec_,"qrec/F");
  tree_->Branch("mipsim",&mipsim_,"mipsim/F");
  tree_->Branch("miprec",&miprec_,"miprec/F");
  tree_->Branch("avgmiprec",&avgmiprec_,"avgmiprec/F");
  tree_->Branch("cce",&cce_,"cce/F");
  tree_->Branch("eta",&eta_,"eta/F");
  tree_->Branch("radius",&radius_,"radius/F");
  tree_->Branch("z",&z_,"z/F");
  tree_->Branch("isTOT",&isToT_,"isTOT/I");
  tree_->Branch("layer",&layer_,"layer/I");
  tree_->Branch("thick",&thick_,"thick/I");
  tree_->Branch("isSci",&isSci_,"isSci/I");
}

//
void HGCDigiTester::endJob()
{
}

//
void HGCDigiTester::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  event_++;
 
  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticles_, genParticlesHandle);
  std::map<int,LorentzVector> photons;
  std::map<int,Point> photonVertex;
  for(size_t i = 0; i < genParticlesHandle->size(); ++i )  {    
    const reco::GenParticle &p = (*genParticlesHandle)[i];
    int idx(p.eta()>0 ? 1 : -1);
    photons[idx]=p.p4();
    photonVertex[idx]=p.vertex();
  }

  //read sim hits
  //acumulate total energy deposited in each DetId
  edm::Handle<edm::PCaloHitContainer> simHitsCEE,simHitsCEH,simHitsCEHSci;
  iEvent.getByToken(simHitsCEE_,    simHitsCEE);
  iEvent.getByToken(simHitsCEH_,    simHitsCEH);
  iEvent.getByToken(simHitsCEHSci_, simHitsCEHSci);

  std::map<uint32_t,double> simE;
  for(size_t i=0; i<3; i++) {
    const std::vector<PCaloHit> &simHits((i==0 ? *simHitsCEE : (i==1 ? *simHitsCEH : *simHitsCEHSci)));
    for(auto sh : simHits) {
      //assign a reco-id
      uint32_t key(sh.id());
      if(simE.find(key)==simE.end()) simE[key]=0.;
      simE[key]=simE[key]+sh.energy();
    }
  }    

  //read digis
  edm::Handle<HGCalDigiCollection> digisCEE,digisCEH,digisCEHSci;
  iEvent.getByToken(digisCEE_,digisCEE);
  iEvent.getByToken(digisCEH_,digisCEH);
  iEvent.getByToken(digisCEHSci_,digisCEHSci);

  //read geometry components for HGCAL
  edm::ESHandle<HGCalGeometry> geoHandleCEE,geoHandleCEH,geoHandleCEHSci;

  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandleCEE);
  const HGCalGeometry *geoCEE=geoHandleCEE.product();
  const HGCalTopology &topoCEE=geoCEE->topology();
  const HGCalDDDConstants &dddConstCEE=topoCEE.dddConstants();

  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",geoHandleCEH);
  const HGCalGeometry *geoCEH=geoHandleCEH.product();
  const HGCalTopology &topoCEH=geoCEH->topology();
  const HGCalDDDConstants &dddConstCEH=topoCEH.dddConstants();

  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",geoHandleCEHSci);
  const HGCalGeometry *geoCEHSci=geoHandleCEHSci.product();
  const HGCalTopology &topoCEHSci=geoCEHSci->topology();
  const HGCalDDDConstants &dddConstCEHSci=topoCEHSci.dddConstants();

  //set the geometry
  scal_[0]->setGeometry(geoCEE, HGCalSiNoiseMap::AUTO, mipTarget_[0]);
  scal_[1]->setGeometry(geoCEH, HGCalSiNoiseMap::AUTO, mipTarget_[1]);
  
  //loop over digis
  size_t itSample(2);
  for(size_t i=0; i<2; i++) {

    const HGCalDigiCollection &digis(i==0 ? *digisCEE : *digisCEH);
    const HGCalDDDConstants &dddConst(i==0 ? dddConstCEE : dddConstCEH);

    for(auto d : digis) {

      HGCSiliconDetId cellId(d.id());
      //read digi (in-time sample only)
      uint32_t adc(d.sample(itSample).data() );
      isToT_ = d.sample(itSample).mode();
      HGCalSiNoiseMap::GainRange_t gain = (HGCalSiNoiseMap::GainRange_t)d.sample(itSample).gain();

      //get the conditions for this det id
      HGCalSiNoiseMap::SiCellOpCharacteristics siop=scal_[i]->getSiCellOpCharacteristics(cellId);

      //convert back to charge
      double adcLSB=useVanillaCfg_ ? vanilla_adcLSB_fC_[i] : scal_[i]->getLSBPerGain()[gain];
      qrec_=(adc+0.5)*adcLSB ;
      if(isToT_) 
        qrec_=tdcOnset_fC_[i]+(adc+0.5)*tdcLSB_[i];          

      //get the simulated charge for this hit
      uint32_t key( cellId.rawId() );
      if(simE.find(key)==simE.end()) continue; //ignore for the moment

      qsim_ = simE[key] * 1.0e6 * 0.044259; // GeV -> fC  (1000 eV / 3.62 (eV per e) / 6.24150934e3 (e per fC))
      
      //convert fC to #MIPs
      thick_  = cellId.type();
      double mipEqfC( scal_[i]->getMipEqfC()[thick_] );
      miprec_=qrec_/mipEqfC;
      avgmiprec_=qrec_/vanilla_mipfC_[i];
      mipsim_=qsim_/mipEqfC;
      
      cce_=useVanillaCfg_ ? 1 : siop.core.cce;

      //additional info
      layer_ = cellId.layer();
      const auto &xy(dddConst.locateCell(cellId.layer(), cellId.waferU(), cellId.waferV(), cellId.cellU(), cellId.cellV(), true, true));
      radius_ = sqrt(std::pow(xy.first, 2) + std::pow(xy.second, 2));  //in cm
      int zside(cellId.zside());
      z_      = zside*dddConst.waferZ(layer_,true);
      eta_    = TMath::ATanH(z_/sqrt(radius_*radius_+z_*z_));

      isSci_  = false;
      
      genergy_  = photons[zside].energy();
      gpt_      = photons[zside].pt();
      geta_     = photons[zside].eta();
      gphi_     = photons[zside].phi();
      gvradius_ = sqrt(pow(photons[zside].x(),2)+pow(photons[zside].y(),2));
      gvz_      = photons[zside].z();
      tree_->Fill();
    }
  }

}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
