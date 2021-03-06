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

    //Si specific
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

      std::vector<double> avg_mip=iConfig.getParameter<std::vector<double> >(i==0? "hgcee_fCPerMIP" : "hgceh_fCPerMIP");
      for(auto v:avg_mip)avg_mipfC_[i].push_back(v);
    }

    //Sci-specific
    if(i==2) {

      scaleByTileArea_ = digiCfg.getParameter<bool>("scaleByTileArea");
      scaleBySipmArea_ = digiCfg.getParameter<bool>("scaleBySipmArea");

      scalSci_ = new HGCalSciNoiseMap;
      scalSci_->setDoseMap(digiCfg.getParameter<edm::ParameterSet>("noise").template getParameter<std::string>("doseMap"),
                           digiCfg.getParameter<edm::ParameterSet>("noise").template getParameter<uint32_t>("scaleByDoseAlgo"));
      scalSci_->setFluenceScaleFactor(digiCfg.getParameter<edm::ParameterSet>("noise").getParameter<double>("scaleByDoseFactor"));
      scalSci_->setSipmMap(digiCfg.getParameter<std::string>("sipmMap"));
      sci_keV2MIP_ = iConfig.getParameter<double>("hgcehsci_keV2DIGI");
    }

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
  //tree_->Branch("qsim",&qsim_,"qsim/F");
  // tree_->Branch("qrec",&qrec_,"qrec/F");
  tree_->Branch("mipsim",&mipsim_,"mipsim/F");
  tree_->Branch("miprec",&miprec_,"miprec/F");
  tree_->Branch("avgmiprec",&avgmiprec_,"avgmiprec/F");
  tree_->Branch("avgmipsim",&avgmipsim_,"avgmipsim/F");
  tree_->Branch("cce",&cce_,"cce/F");
  tree_->Branch("eta",&eta_,"eta/F");
  tree_->Branch("radius",&radius_,"radius/F");
  tree_->Branch("z",&z_,"z/F");
  tree_->Branch("isSat",&isSat_,"isSat/I");
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
    int idx(p.pz()>0 ? 1 : -1);
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
  // const HGCalTopology &topoCEHSci=geoCEHSci->topology();
  // const HGCalDDDConstants &dddConstCEHSci=topoCEHSci.dddConstants();

  //set the geometry
  scal_[0]->setGeometry(geoCEE, HGCalSiNoiseMap::AUTO, mipTarget_[0]);
  scal_[1]->setGeometry(geoCEH, HGCalSiNoiseMap::AUTO, mipTarget_[1]);
  scalSci_->setGeometry(geoCEHSci);

  //loop over digis
  size_t itSample(2);
  for(size_t i=0; i<3; i++) {

    const HGCalDigiCollection &digis(i==0 ? *digisCEE : (i==1 ? *digisCEH : *digisCEHSci) );

    for(auto d : digis) {

      //check if it's matched to a simId
      uint32_t key( d.id().rawId() );
      if(simE.find(key)==simE.end()) continue;

      //read digi (in-time sample only)
      uint32_t adc(d.sample(itSample).data() );
      isToT_=d.sample(itSample).mode();

      isSat_=false;
      //if(!isToT_ && adc>1020) isSat_=true; // never happens
      if(isToT_ && adc==4095) isSat_=true;
      
      //reset common variables
      double adcLSB(vanilla_adcLSB_fC_[i]),mipEqfC(1.0),avgMipEqfC(1.0);
      qsim_=0;
      qrec_=0.;
      cce_=1.;
      thick_=-1;
      layer_=0;
      radius_=0;
      z_=0;
      int zside=0;

      HGCalSiNoiseMap::GainRange_t gain = (HGCalSiNoiseMap::GainRange_t)d.sample(itSample).gain();

      //scintillator does not yet attribute different gains
      if(!useVanillaCfg_ && i!=2) {
        adcLSB=scal_[0]->getLSBPerGain()[gain];
      }

      //Si-specific
      if(i<2) {

        //simulated charge
        qsim_ = simE[key] * 1.0e6 * 0.044259; // GeV -> fC  (1000 eV / 3.62 (eV per e) / 6.24150934e3 (e per fC))

        //get the conditions for this det id
        HGCSiliconDetId cellId(d.id());
        HGCalSiNoiseMap::SiCellOpCharacteristics siop=scal_[i]->getSiCellOpCharacteristics(cellId);
        cce_=siop.core.cce;
        
        thick_     = cellId.type();
        mipEqfC     = scal_[i]->getMipEqfC()[thick_];
        avgMipEqfC = avg_mipfC_[i][thick_];
        isSci_     = false;

        //additional info
        layer_ = cellId.layer();
        const HGCalDDDConstants &dddConst(i==0 ? dddConstCEE : dddConstCEH);
        const auto &xy(dddConst.locateCell(cellId.layer(), cellId.waferU(), cellId.waferV(), cellId.cellU(), cellId.cellV(), true, true));
        radius_ = sqrt(std::pow(xy.first, 2) + std::pow(xy.second, 2));  //in cm
        zside=cellId.zside();
        z_ = zside*dddConst.waferZ(layer_,true);
      }

      //Sci-specific
      if(i==2) {

        //simulated "charge" (in reality this is in MIP units)
        qsim_ = simE[key] *1.0e+6 * sci_keV2MIP_; // keV to mip

        //get the conditions for this det id
        //signal scaled by tile and sipm area + dose
        HGCScintillatorDetId cellId(d.id());
        GlobalPoint global = geoCEHSci->getPosition(cellId);
        radius_ = scalSci_->computeRadius(cellId);
        if(!useVanillaCfg_ && scalSci_->algo()==0) {
          double signal_scale( scalSci_->scaleByDose(cellId,radius_).first );
          if(scaleBySipmArea_) signal_scale *= scalSci_->scaleBySipmArea(cellId,radius_);
          if(scaleByTileArea_) signal_scale *= scalSci_->scaleByTileArea(cellId,radius_);
          cce_ = signal_scale;
        }

        mipEqfC    = 1.0;     //the digis are already in MIP units
        avgMipEqfC = 1.0;

        //additional info
        layer_ = cellId.layer();
        z_     = global.z();
        zside  = cellId.zside();
        //zside  = (z_<0 ? -1 : 1);
        thick_ = -1;
        isSci_ = true;
      }
       

      //convert back to charge
      qrec_=(adc+0.5)*adcLSB ;
      if(isToT_) 
        qrec_=tdcOnset_fC_[i]+(adc+0.5)*tdcLSB_[i];          

      //convert charge to #MIPs
      miprec_     = qrec_/mipEqfC;
      mipsim_     = qsim_/mipEqfC;
      avgmiprec_  = qrec_/avgMipEqfC;
      avgmipsim_  = qsim_/avgMipEqfC;

      //additional info
      eta_    = TMath::ATanH(z_/sqrt(radius_*radius_+z_*z_));

      //for CEH shift by CEE layers
      if(i>0) layer_+=28;

      //MC truth
      genergy_  = photons[zside].energy();
      gpt_      = photons[zside].pt();
      geta_     = photons[zside].eta();
      gphi_     = photons[zside].phi();
      gvradius_ = sqrt(pow(photonVertex[zside].x(),2)+pow(photonVertex[zside].y(),2));
      gvz_      = photonVertex[zside].z();

      //store hit
      tree_->Fill();
    }
  }

}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
