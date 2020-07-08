#ifndef waferoccupancyhisto_h
#define waferoccupancyhisto_h

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TString.h"
#include "TMath.h"

#include <exception>
#include <iostream>
#include <algorithm>

/**
   @class WaferOccupancyHisto
   @short book keeps a set of occupancy histograms for a wafer
*/
class WaferOccupancyHisto
{
public:

  typedef std::tuple<int,int,int,int> WaferKey_t;    
  typedef std::pair<int,int> UVKey_t;

  /**
     @short CTOR
  */
 WaferOccupancyHisto(WaferKey_t key) : ncells_(0), 
    nEvents_(0), 
    waferType_(-1),
    isInit_(false)
    { 

      int subDet=std::get<0>(key);
      int layer=std::get<1>(key);
      int u=std::get<2>(key);
      int v=std::get<3>(key);      
      myID_ = Form("%d_lay%d_%d_%d",subDet,layer,u,v);
    }

  /**
     increments pad counts
   */
  void addPad(int waferType) { 
    waferType_=waferType; 
    ncells_++; 
  }

  /**
     @short book all histograms
   */
  void bookHistos(edm::Service<TFileService> *fs) {

    if(isInit_ || fs==NULL) return;

    TFileDirectory mySubDir=(*fs)->mkdir(myID_.Data());

    TH1F *info = mySubDir.make<TH1F>(myID_+"_info","",2,0,2);
    info->SetBinContent(1,waferType_);
    info->SetBinContent(2,ncells_);
    
    //book histos for BX-1 and in-time bunch
    for(int i=-1; i<=0; i++){
      TString pfix(i==0 ? "" : "_bxm1");
      histos_["adc"+pfix] = mySubDir.make<TH1F>(myID_+"_adc"+pfix,";q [ADC];",100,0,100);
      histos_["adczs"+pfix] = mySubDir.make<TH1F>(myID_+"_adczs"+pfix,";q [ADC];",100,0,100);
      histos_["adcfull"+pfix] = mySubDir.make<TH1F>(myID_+"_adcfull"+pfix,";q [ADC];",100,0,500);
      
      adcCounts_[pfix]=0;
      histos_["counts"+pfix]    = mySubDir.make<TH1F>(myID_+"_counts"+pfix,";Counts above threshold;",ncells_+1,0,ncells_+1);

      adcCountsZS_[pfix]=0;
      histos_["countszs"+pfix]    = mySubDir.make<TH1F>(myID_+"_countszs"+pfix,";Counts above threshold;",ncells_+1,0,ncells_+1);

      toaCounts_[pfix]=0;
      histos_["toacounts"+pfix]  = mySubDir.make<TH1F>(myID_+"_toacounts"+pfix,";Cells in TOA;",ncells_+1,0,ncells_+1);

      tdcCounts_[pfix]=0;
      histos_["tdccounts"+pfix]  = mySubDir.make<TH1F>(myID_+"_tdccounts"+pfix,";Cells in TDC;",ncells_+1,0,ncells_+1);

      busyCounts_[pfix]=0;
      histos_["busycounts"+pfix] = mySubDir.make<TH1F>(myID_+"_busycounts"+pfix,";Busy cells;",ncells_+1,0,ncells_+1);
    }
    
    for(auto h : histos_) h.second->Sumw2();
    
    isInit_=true;
  }
  
  /**
     @short accumulate counts for a new hit
  */
  inline void count(uint32_t adc,bool isTOA,bool isTDC, bool isBusy, uint32_t thr, bool passZS, TString pfix="")
  {
    if(!isInit_) return;

    if(!isTDC){
      histos_["adc"+pfix]->Fill(adc);
      histos_["adcfull"+pfix]->Fill(adc);
      if(passZS)
        histos_["adczs"+pfix]->Fill(adc);
    }

    if(!isBusy) {
      
      //pads with valid energy measurement
      if( (!isTDC && adc>thr) || isTDC ) {
        
        adcCounts_[pfix]++;

        if(passZS || isTDC) adcCountsZS_[pfix]++;
        
        //+valid toa
        if(isTOA) toaCounts_[pfix]++;

        //+TDC
        if(isTDC && isTOA) tdcCounts_[pfix]++;
      }

    }
    else {

      busyCounts_[pfix]++;

    }
  }

  /**
     @short to be called once all hits have been counted
  */
  inline void analyze()
  {
    if(!isInit_) return;

    nEvents_++;
    
    for(int i=-1; i<=0; i++){
      TString pfix(i==0 ? "" : "_bxm1");
      histos_["counts"+pfix]->Fill(adcCounts_[pfix]);
      histos_["countszs"+pfix]->Fill(adcCountsZS_[pfix]);
      histos_["toacounts"+pfix]->Fill(toaCounts_[pfix]);
      histos_["tdccounts"+pfix]->Fill(tdcCounts_[pfix]);
      histos_["busycounts"+pfix]->Fill(busyCounts_[pfix]);
    }    
  }

  /**
     @short set all to 0
   */
  inline void resetCounters()
  {
    if(!isInit_) return;

    for(int i=-1; i<=0; i++){
      TString pfix(i==0 ? "" : "_bxm1");
      adcCounts_[pfix]=0;
      adcCountsZS_[pfix]=0;
      toaCounts_[pfix]=0;
      tdcCounts_[pfix]=0;
      busyCounts_[pfix]=0;
    }
  }


  /**
     @short normalize according to the number of events analyzed
   */
  inline void endJob() 
  {
    if(!isInit_ || nEvents_==0) return;

    //scale by the number of wafer equivalent
    float norm(nEvents_);    
    for(auto h : histos_) h.second->Scale(1./norm);
  }

  inline int getCells() { return ncells_; }
  inline int getCounts(TString k="",bool isZS=false) { return isZS ?  adcCountsZS_[k] : adcCounts_[k]; }

  /**
     @short DTOR
   */
  ~WaferOccupancyHisto() {}

 
 private:

  TString myID_;
  int ncells_,nEvents_,waferType_;
  bool isInit_;
  std::map<TString, TH1F *> histos_;
  std::map<TString,int> adcCounts_,adcCountsZS_,toaCounts_,tdcCounts_,busyCounts_;
};

#endif
