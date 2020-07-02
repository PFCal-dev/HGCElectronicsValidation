#ifndef _HGCOccupancyAnalyzer_h_
#define _HGCOccupancyAnalyzer_h_

#include "UserCode/HGCElectronicsValidation/interface/WaferOccupancyHisto.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerBase.h"  

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <string>

#include "TProfile.h"
#include "TH2F.h"

#include <unordered_map>


/**
   @class HGCOccupancyAnalyzer
   @short an EDAnalyzer to readout the digis in the event and steer the filling of the occupancy histograms
*/

class HGCOccupancyAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCOccupancyAnalyzer( const edm::ParameterSet& );
  ~HGCOccupancyAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  
  /**
     @short analyze a digi collection
   */
  void analyzeDigis(int ,edm::Handle<HGCalDigiCollection> &, const HGCalGeometry *);

  //histograms for the wafers
  std::map<WaferOccupancyHisto::WaferKey_t,WaferOccupancyHisto *> waferHistos_;

  //generator level information
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJets_;

  //geometry and digis to analyze
  std::string geoCEE_,geoCEH_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_; 

  std::vector<int> tdcHits_, toaHits_,adcHits_;
  TH1F *cellCount_;
  TProfile *tdcCountProf_,*toaCountProf_,*adcCountProf_;
  TH2F *adcHitsVsPU_;

  int nevts_;
  double adcThrMIP_,adcThrMIPbxm1_;

  //
  std::string doseMap_;
  std::map<int,HGCalSiNoiseMap *> noiseMaps_;
};
 

#endif
