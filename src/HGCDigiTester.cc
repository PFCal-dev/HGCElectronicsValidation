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
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
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
#include <math.h> 

using namespace std;

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// PLUGIN IMPLEMENTATION
//


//
HGCDigiTester::HGCDigiTester( const edm::ParameterSet &iConfig )
{   
  refSpeed_ = 0.1 * CLHEP::c_light;

  hardProcOnly_=iConfig.getParameter<bool>("hardProcOnly");
  onlyROCTree_=iConfig.getParameter<bool>("onlyROCTree");
  genParticles_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticleSrc"));

  //configure noise map
  std::string digitizers[]={"hgceeDigitizer","hgcehDigitizer","hgcehsciDigitizer"};
  for(size_t i=0; i<3; i++) {
    edm::ParameterSet cfg(iConfig.getParameter<edm::ParameterSet>(digitizers[i]));
    edm::ParameterSet digiCfg(cfg.getParameter<edm::ParameterSet>("digiCfg"));
    edm::ParameterSet feCfg(digiCfg.getParameter<edm::ParameterSet>("feCfg"));

    //Si specific
    if(i<2) {
      HGCalSiNoiseMap<HGCSiliconDetId> *scal=new HGCalSiNoiseMap<HGCSiliconDetId>;
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
      //pxFiringRate_    = digiCfg.getParameter<edm::ParameterSet>("noise").template getParameter<double>("pxFiringRate");

      scalSci_ = new HGCalSciNoiseMap;
      scalSci_->setDoseMap(digiCfg.getParameter<edm::ParameterSet>("noise").template getParameter<std::string>("doseMap"),
                           digiCfg.getParameter<edm::ParameterSet>("noise").template getParameter<uint32_t>("scaleByDoseAlgo"));
      scalSci_->setFluenceScaleFactor(digiCfg.getParameter<edm::ParameterSet>("noise").getParameter<double>("scaleByDoseFactor"));
      scalSci_->setSipmMap(digiCfg.getParameter<std::string>("sipmMap"));
      sci_keV2MIP_ = iConfig.getParameter<double>("hgcehsci_keV2DIGI");
    }

    simHitsColls_[i] = consumes<std::vector<PCaloHit>>(cfg.getParameter<edm::InputTag>("hitCollection"));
    digisColls_[i] = consumes<HGCalDigiCollection>(cfg.getParameter<edm::InputTag>("digiCollection"));

    bxTime_[i]=cfg.getParameter<double>("bxTime");
    tofDelay_[i]=cfg.getParameter<double>("tofDelay");

    mipTarget_[i]=feCfg.getParameter<uint32_t>("targetMIPvalue_ADC");
    tdcOnset_fC_[i]=feCfg.getParameter<double>("tdcOnset_fC");
    tdcLSB_[i]=feCfg.getParameter<double>("tdcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("tdcNbits") );
    vanilla_adcLSB_fC_[i] = feCfg.getParameter<double>("adcSaturation_fC") / pow(2., feCfg.getParameter<uint32_t>("adcNbits") );
  }
  useTDCOnsetAuto_ = iConfig.getParameter<bool>("useTDCOnsetAuto");
  useVanillaCfg_ = iConfig.getParameter<bool>("useVanillaCfg");

  event_=0;
  edm::Service<TFileService> fs;
  if(!onlyROCTree_) {
    tree_ = fs->make<TTree>("hits","hits");
    tree_->Branch("genergy",&genergy_,"genergy/F");
    tree_->Branch("gpt",&gpt_,"gpt/F");
    tree_->Branch("geta",&geta_,"geta/F");
    tree_->Branch("gphi",&gphi_,"gphi/F");
    tree_->Branch("gvradius",&gvradius_,"gvradius/F");
    tree_->Branch("gvz",&gvz_,"gvz/F");
    tree_->Branch("event",&event_,"event/I");
    tree_->Branch("detid",&detid_,"detid/i");
    tree_->Branch("layer",&layer_,"layer/I");
    tree_->Branch("u",&u_,"u/I");  //u or iphi
    tree_->Branch("v",&v_,"v/I");  //v or ieta
    tree_->Branch("roc",&roc_,"roc/I");
    tree_->Branch("adc",&adc_,"adc/I");
    tree_->Branch("gain",&gain_,"gain/I");
    //tree_->Branch("qsim",&qsim_,"qsim/F");
    // tree_->Branch("qrec",&qrec_,"qrec/F");
    tree_->Branch("mipsim",&mipsim_,"mipsim/F");
    tree_->Branch("mipsimInBX",&mipsimInBX_,"mipsimInBX/F");
    tree_->Branch("mipsimPreBX",&mipsimPreBX_,"mipsimPreBX/F");
    tree_->Branch("mipsimPostBX",&mipsimPostBX_,"mipsimPostBX/F");
    tree_->Branch("mipsimPerBX", &mipsimPerBX_);
    tree_->Branch("miprec",&miprec_,"miprec/F");
    tree_->Branch("avgmiprec",&avgmiprec_,"avgmiprec/F");
    tree_->Branch("avgmipsim",&avgmipsim_,"avgmipsim/F");
    tree_->Branch("cce",&cce_,"cce/F");
    tree_->Branch("eta",&eta_,"eta/F");
    tree_->Branch("radius",&radius_,"radius/F");
    tree_->Branch("z",&z_,"z/F");
    tree_->Branch("isSat",&isSat_,"isSat/I");
    tree_->Branch("isTOT",&isToT_,"isTOT/I");
    tree_->Branch("thick",&thick_,"thick/I");
    tree_->Branch("isSci",&isSci_,"isSci/I");
  }

  rocTree_ = fs->make<TTree>("rocs","rocs");
  rocTree_->Branch("event",&event_,"event/I");
  rocTree_->Branch("layer",&layer_,"layer/I");
  rocTree_->Branch("side",&side_,"side/O");
  rocTree_->Branch("u",&u_,"u/I");
  rocTree_->Branch("v",&v_,"v/I");
  rocTree_->Branch("roc",&roc_,"roc/I");
  rocTree_->Branch("nhits",&nhits_,"nhits/I");
  rocTree_->Branch("miprec",&miprec_,"miprec/F");
}

//
void HGCDigiTester::endJob()
{
}

//
void HGCDigiTester::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  event_++;
  rocDeposits_t rocs;

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

  //read geometry components for HGCAL
  edm::ESHandle<HGCalGeometry> geoHandleCEE,geoHandleCEH,geoHandleCEHSci;

  const HGCalGeometry *geoList[3];

  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandleCEE);
  //const HGCalGeometry *geoCEE=geoHandleCEE.product();
  geoList[0]=geoHandleCEE.product();
  const HGCalTopology &topoCEE=geoList[0]->topology();
  const HGCalDDDConstants &dddConstCEE=topoCEE.dddConstants();

  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",geoHandleCEH);
  //const HGCalGeometry *geoCEH=geoHandleCEH.product();
  geoList[1]=geoHandleCEH.product();
  const HGCalTopology &topoCEH=geoList[1]->topology();
  const HGCalDDDConstants &dddConstCEH=topoCEH.dddConstants();
  HGCSiliconDetIdToROC sid2roc;

  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",geoHandleCEHSci);
  //const HGCalGeometry *geoCEHSci=geoHandleCEHSci.product();
  geoList[2]=geoHandleCEHSci.product();
  // const HGCalTopology &topoCEHSci=geoCEHSci->topology();
  // const HGCalDDDConstants &dddConstCEHSci=topoCEHSci.dddConstants();

  //set the geometry
  scal_[0]->setGeometry(geoList[0], HGCalSiNoiseMap<HGCSiliconDetId>::AUTO, mipTarget_[0]);
  scal_[1]->setGeometry(geoList[1], HGCalSiNoiseMap<HGCSiliconDetId>::AUTO, mipTarget_[1]);
  scalSci_->setGeometry(geoList[2]);

  //read sim hits
  //acumulate total energy deposited in each DetId
  edm::Handle<edm::PCaloHitContainer> simHitsCEE,simHitsCEH,simHitsCEHSci;
  iEvent.getByToken(simHitsColls_[0], simHitsCEE);
  iEvent.getByToken(simHitsColls_[1], simHitsCEH);
  iEvent.getByToken(simHitsColls_[2], simHitsCEHSci);

  std::map<uint32_t,double> simE;
  std::map<uint32_t,double> simEinBX;
  std::map<uint32_t,double> simEpreBX;
  std::map<uint32_t,double> simEpostBX;
  std::map<uint32_t,std::vector<double>> simEperBX;
  const int indexInBX = 9;
  const int nBXs = 15;
  for(size_t i=0; i<3; i++) {
    const std::vector<PCaloHit> &simHits((i==0 ? *simHitsCEE : (i==1 ? *simHitsCEH : *simHitsCEHSci)));
    for(auto sh : simHits) {
      //assign a reco-id
      uint32_t key(sh.id());
      if(simE.find(key)==simE.end()) simE[key]=0.;
      double hitEnergy = sh.energy();
      simE[key]=simE[key]+hitEnergy;

      // Following the setup in https://github.com/cms-sw/cmssw/blob/8fdeed6c37e3a53c53f7ff2d8f9303867f37a2c1/SimCalorimetry/HGCalSimProducers/plugins/HGCDigitizer.cc#L601-L602
      float toa = (float)sh.time();
      float dist2center = geoList[i]->getPosition(key).mag();
      float tof = toa - dist2center / refSpeed_ + tofDelay_[i];
      int itime = std::floor(tof / bxTime_[i]);

      if (itime == 0) { // in-time BX
        if(simEinBX.find(key)==simEinBX.end()) simEinBX[key]=0.;
        simEinBX[key]=simEinBX[key]+hitEnergy;
      }
      else if (itime >= -1*indexInBX && itime < 0) { // preceding BX
        if(simEpreBX.find(key)==simEpreBX.end()) simEpreBX[key]=0.;
        simEpreBX[key]=simEpreBX[key]+hitEnergy;
      }
      else if (itime > 0 && itime <= 5) { // following BX
        if(simEpostBX.find(key)==simEpostBX.end()) simEpostBX[key]=0.;
        simEpostBX[key]=simEpostBX[key]+hitEnergy;
      }

      int iBX = itime + indexInBX;
      if (iBX < nBXs)
      {
        if (simEperBX.find(key) == simEperBX.end()) {
          simEperBX[key].clear();
          simEperBX[key].resize(nBXs, 0);
        }
        simEperBX[key][iBX] += hitEnergy;
      }
    }
  }

  //read digis
  edm::Handle<HGCalDigiCollection> digisCEE,digisCEH,digisCEHSci;
  iEvent.getByToken(digisColls_[0],digisCEE);
  iEvent.getByToken(digisColls_[1],digisCEH);
  iEvent.getByToken(digisColls_[2],digisCEHSci);

  //loop over digis
  size_t itSample(2);
  for(size_t i=0; i<3; i++) {

    const HGCalDigiCollection &digis(i==0 ? *digisCEE : (i==1 ? *digisCEH : *digisCEHSci) );

    for(auto d : digis) {

      //check if it's matched to a simId
      uint32_t key( d.id().rawId() );
      detid_=key;
      bool simEexists(true); 
      if(simE.find(key)==simE.end()){
        if(hardProcOnly_) continue;
        simEexists=false;
      }
      //read digi (in-time sample only)
      uint32_t adc(d.sample(itSample).data() );
      adc_=adc;
      isToT_=d.sample(itSample).mode();

      isSat_=false;
      //if(!isToT_ && adc>1020) isSat_=true; // never happens
      if(isToT_ && adc==4095) isSat_=true;
      
      //reset common variables
      double adcLSB(vanilla_adcLSB_fC_[i]),mipEqfC(1.0),avgMipEqfC(1.0);
      qsim_=0;
      qsimInBX_=0;
      qsimPreBX_=0;
      qsimPostBX_=0;
      qsimPerBX_.clear();
      qrec_=0.;
      cce_=1.;
      thick_=-1;
      layer_=0;
      u_=0;
      v_=0;
      roc_=0;
      radius_=0;
      z_=0;
      int zside=0;

      HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t gain = (HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t)d.sample(itSample).gain();
      gain_ = gain;

      //scintillator does not yet attribute different gains
      if(!useVanillaCfg_ && i!=2) {
        adcLSB=scal_[0]->getLSBPerGain()[gain];
      }

      //Si-specific
      if(i<2) {

        //simulated charge
        qsim_ = simEexists ? simE[key] * 1.0e6 * 0.044259 : 0.; // GeV -> fC  (1000 eV / 3.62 (eV per e) / 6.24150934e3 (e per fC))
        qsimInBX_ = simEexists ? simEinBX[key] * 1.0e6 * 0.044259 : 0.;
        qsimPreBX_ = simEexists ? simEpreBX[key] * 1.0e6 * 0.044259 : 0.;
        qsimPostBX_ = simEexists ? simEpostBX[key] * 1.0e6 * 0.044259 : 0.;
        qsimPerBX_.clear();
        for (std::vector<double>::iterator itBX = simEperBX[key].begin(); itBX != simEperBX[key].end(); ++itBX) {
          float qtmp = simEexists ? (*itBX) * 1.0e6 * 0.044259 : 0.;
          qsimPerBX_.push_back( qtmp );
        }

        //get the conditions for this det id
        HGCSiliconDetId cellId(d.id());
        u_ = cellId.waferUV().first;
        v_ = cellId.waferUV().second;
        roc_    = sid2roc.getROCNumber(cellId);
        HGCalSiNoiseMap<HGCSiliconDetId>::SiCellOpCharacteristics siop 
          = scal_[i]->getSiCellOpCharacteristics(cellId);
        cce_       = siop.core.cce;
        thick_     = cellId.type();
        mipEqfC    = scal_[i]->getMipEqfC()[thick_];
        avgMipEqfC = avg_mipfC_[i][thick_];
        isSci_     = false;

        //additional info
        layer_ = cellId.layer();
        
        const HGCalDDDConstants &dddConst(i==0 ? dddConstCEE : dddConstCEH);
        const auto &xy(dddConst.locateCell(cellId.layer(), cellId.waferU(), cellId.waferV(), cellId.cellU(), cellId.cellV(), true, true));
        radius_ = sqrt(std::pow(xy.first, 2) + std::pow(xy.second, 2));  //in cm
        zside=cellId.zside();
        z_ = zside*dddConst.waferZ(layer_,true);

        // get tdcOnset from gain
        if (useTDCOnsetAuto_) {
          tdcOnset_fC_[i] = scal_[i]->getTDCOnsetAuto(gain);
        }
      }

      //Sci-specific
      //maybe here we could also add a ROC ?
      if(i==2) {

        HGCScintillatorDetId scId(d.id());
        u_=scId.iphi();
        v_=scId.ieta();

        //simulated "charge" (in reality this is in MIP units)
        qsim_ = simEexists ? simE[key] *1.0e+6 * sci_keV2MIP_ : 0.; // keV to mip
        qsimInBX_ = simEexists ? simEinBX[key] * 1.0e6 * sci_keV2MIP_ : 0.;
        qsimPreBX_ = simEexists ? simEpreBX[key] * 1.0e6 * sci_keV2MIP_ : 0.;
        qsimPostBX_ = simEexists ? simEpostBX[key] * 1.0e6 * sci_keV2MIP_ : 0.;
        qsimPerBX_.clear();
        for (std::vector<double>::iterator itBX = simEperBX[key].begin(); itBX != simEperBX[key].end(); ++itBX) {
          float qtmp = simEexists ? (*itBX) * 1.0e6 * sci_keV2MIP_ : 0.;
          qsimPerBX_.push_back( qtmp );
        }

        //get the conditions for this det id
        //signal scaled by tile and sipm area + dose
        HGCScintillatorDetId cellId(d.id());
        GlobalPoint global = geoList[2]->getPosition(cellId);
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
      mipsimInBX_ = qsimInBX_/mipEqfC;
      mipsimPreBX_ = qsimPreBX_/mipEqfC;
      mipsimPostBX_ = qsimPostBX_/mipEqfC;
      mipsimPerBX_.clear();
      for (std::vector<float>::iterator itqsim = qsimPerBX_.begin(); itqsim != qsimPerBX_.end(); ++itqsim) {
        mipsimPerBX_.push_back( (*itqsim)/mipEqfC );
      }
      avgmiprec_  = qrec_/avgMipEqfC;
      avgmipsim_  = qsim_/avgMipEqfC;

      //additional info
      //eta_    = TMath::ATanH(z_/sqrt(radius_*radius_+z_*z_));
      eta_    = atanh(z_/sqrt(radius_*radius_+z_*z_));

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
      if(!onlyROCTree_) tree_->Fill();

      //increment energy for the ROC (Si only)
      if(i==2) continue;
      rocKey_t k(layer_,(z_>0),u_,v_,roc_);
      if(rocs.find(k)==rocs.end()) rocs[k]=rocSummary_t(0.,0.);
      rocs[k].first += 1;
      rocs[k].second += miprec_/cce_;
    }
  }
  
  //save the summary of energy deposits in ROCs
  for(auto r : rocs ){
    layer_=std::get<0>(r.first) ;
    side_=std::get<1>(r.first) ;
    u_=std::get<2>(r.first) ;
    v_=std::get<3>(r.first) ;
    roc_=std::get<4>(r.first) ;
    nhits_=r.second.first;
    miprec_=r.second.second;
    rocTree_->Fill();
  }



}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCDigiTester);
