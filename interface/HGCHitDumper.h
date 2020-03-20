#ifndef _HGCHitDumper_h_
#define _HGCHitDumper_h_

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
#include "DataFormats/JetReco/interface/GenJetCollection.h"

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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include <string>

#include "TTree.h"


/**
   @class HGCHitDumper
*/

class HGCHitDumper : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCHitDumper( const edm::ParameterSet& );
  ~HGCHitDumper();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  
  //geometry and hits
  std::string geo_;
  edm::EDGetTokenT<HGCalDigiCollection> digis_;
  edm::EDGetTokenT<HGCRecHitCollection> recHits_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJets_;
  
  TTree *t_;
  Int_t event_;
  Int_t ng_,layer_,waferL_,gain_,digiT_,digiE_;
  Bool_t isTDC_,hasTOA_,threshold_;
  Float_t en_,eta_,recT_,recE_;
};
 

#endif
