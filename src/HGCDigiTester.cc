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

//
// PLUGIN IMPLEMENTATION
//


//
HGCDigiTester::HGCDigiTester( const edm::ParameterSet &iConfig ) 
  : simHitsCEE_( consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE")) ),
    digisCEE_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","EE")) ),
    genParticles_( consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles")) )
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
  tdcOnset_fC_=feCfg.getParameter<double>("tdcOnset_fC");
  tdcLSB_=feCfg.getParameter<double>("tdcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("tdcNbits") );
  vanilla_adcLSB_fC_ = feCfg.getParameter<double>("adcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("adcNbits") );
  useVanillaCfg_ = iConfig.getParameter<bool>("useVanillaCfg");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("hits","hits");
  event_=0;
  tree_->Branch("genergy",&genergy_,"genergy/F");
  tree_->Branch("gpt",&gpt_,"gpt/F");
  tree_->Branch("geta",&geta_,"geta/F");
  tree_->Branch("gphi",&gphi_,"gphi/F");
  tree_->Branch("event",&event_,"event/I");
  tree_->Branch("qsim",&qsim_,"qsim/F");
  tree_->Branch("qrec",&qrec_,"qrec/F");
  tree_->Branch("mipsim",&mipsim_,"mipsim/F");
  tree_->Branch("miprec",&miprec_,"miprec/F");
  tree_->Branch("cce",&cce_,"cce/F");
  tree_->Branch("eta",&eta_,"eta/F");
  tree_->Branch("radius",&radius_,"radius/F");
  tree_->Branch("z",&z_,"z/F");
  tree_->Branch("isTOT",&isToT_,"isTOT/I");
  tree_->Branch("layer",&layer_,"layer/I");
  tree_->Branch("thick",&thick_,"thick/I");
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
  for(size_t i = 0; i < genParticlesHandle->size(); ++i )  {    
    const reco::GenParticle &p = (*genParticlesHandle)[i];
    int idx(p.eta()>0 ? 1 : -1);
    photons[idx]=p.p4();
  }

  //read sim hits
  edm::Handle<edm::PCaloHitContainer> simHits;
  iEvent.getByToken(simHitsCEE_,simHits);

  //read digis
  edm::Handle<HGCalDigiCollection> digis;
  iEvent.getByToken(digisCEE_,digis);

  //read geometry
  edm::ESHandle<HGCalGeometry> geoHandle;
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandle);
  const HGCalGeometry *geo=geoHandle.product();
  const HGCalTopology &topo=geo->topology();
  const HGCalDDDConstants &dddConst=topo.dddConstants();

  //set the geometry
  scal_.setGeometry(geo, HGCalSiNoiseMap::AUTO, mipTarget_);
  
  //acumulate total energy deposited in each DetId
  std::map<uint32_t,double> simE;
  for(auto sh : *simHits) {

    //assign a reco-id
    uint32_t key(sh.id());
    if(simE.find(key)==simE.end()) simE[key]=0.;
    simE[key]=simE[key]+sh.energy();
  }

  //loop over digis
  size_t itSample(2);
  for(auto d : *digis) {

    HGCSiliconDetId cellId(d.id());

    //read digi (in-time sample only)
    uint32_t adc(d.sample(itSample).data() );
    isToT_ = d.sample(itSample).mode();
    HGCalSiNoiseMap::GainRange_t gain = (HGCalSiNoiseMap::GainRange_t)d.sample(itSample).gain();

    //get the conditions for this det id
    HGCalSiNoiseMap::SiCellOpCharacteristics siop=scal_.getSiCellOpCharacteristics(cellId);                                          


    //convert back to charge
    double adcLSB=useVanillaCfg_ ? vanilla_adcLSB_fC_ : scal_.getLSBPerGain()[gain];
    qrec_=(adc+0.5)*adcLSB ;
    if(isToT_) 
      qrec_=tdcOnset_fC_+(adc+0.5)*tdcLSB_;          

    //get the simulated charge for this hit
    uint32_t key( cellId.rawId() );
    if(simE.find(key)==simE.end()) continue; //ignore for the moment
    qsim_ = simE[key] * 1.0e6 * 0.044259; // GeV -> fC  (1000 eV / 3.62 (eV per e) / 6.24150934e3 (e per fC))

    //convert fC to #MIPs
    thick_  = cellId.type(); //dddConst.waferType(layer_,cellId.waferU(),cellId.waferV());
    double mipEqfC( scal_.getMipEqfC()[thick_] );
    miprec_=qrec_/mipEqfC;
    mipsim_=qsim_/mipEqfC;

    cce_=useVanillaCfg_ ? 1 : siop.core.cce;

    //additional info
    layer_ = cellId.layer();
    const auto &xy(dddConst.locateCell(cellId.layer(), cellId.waferU(), cellId.waferV(), cellId.cellU(), cellId.cellV(), true, true));
    radius_ = sqrt(std::pow(xy.first, 2) + std::pow(xy.second, 2));  //in cm
    int zside(cellId.zside());
    z_      = zside*dddConst.waferZ(layer_,true);
    eta_    = TMath::ATanH(z_/sqrt(radius_*radius_+z_*z_));

    genergy_=photons[zside].energy();
    gpt_=photons[zside].pt();
    geta_=photons[zside].eta();
    gphi_=photons[zside].phi();

    tree_->Fill();
  }


    /*
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
  double noise=siop.core.noise;
  std::cout << "ADC lsb=" << adcLSB 
            << " TDC lsb=" << tdcLSB_ 
            << " noise=" <<  noise
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
  //use a uniform distribution of charges in the two ranges
  std::vector<float> qinj;
  for(int i=0; i<pow(2,10); i++){
    for(float x=0; x<=1; x+=0.1) {
      float qval=i*x*adcLSB;
      if(qval>tdcOnset_fC_) continue;
      qinj.push_back(qval);
    }
  }
  for(int i=0; i<pow(2,12); i++){
    for(float x=0; x<=1; x+=0.1) {
      float qval=i*x*tdcLSB_;
      if(qval<tdcOnset_fC_) continue;
      qinj.push_back(qval);
    }
  }   

  for(auto &q:qinj) {

    auto digiResult = std::make_unique<HGCalDigiCollection>();
    (simIt->second).hit_info[0][9]=q; //fill in-time index only
    digitizer_->run(digiResult,simData,geo,validIds,digitizationType_,engine);
      
    //if a digi was not produced move to next value
    if(digiResult->size()==0) continue;
      
    //read digi
    uint32_t mode=((*digiResult)[0])[2].mode();
    uint32_t adc=((*digiResult)[0])[2].data();
      
    //convert back to charge
    //double k=1./sqrt(8*3.1415);
    //double nbias=k*noise/adcLSB;
    double nbias(0.);
    double qrec( (adc+0.5+nbias)*adcLSB );
    if(mode){        
      //nbias=k*noise/tdcLSB_;
      nbias=0.;
      qrec=tdcOnset_fC_+(adc+0.5+nbias)*tdcLSB_;      
    }
    qrec=max(0., qrec);
    double qrec_cce(qrec/cce);

    ofile << q << " " << mode << " " << adc << " " << qrec << " " << qrec_cce << std::endl;    
  }
  
  ofile.close();
    */
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
