#ifndef waferoccupancyhisto_h
#define waferoccupancyhisto_h

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TMath.h"

#include <exception>
#include <iostream>
#include <algorithm>

namespace HGCalWafer{

  /**
     @struct a summary of the info for occupancy study of a pad
  */
  struct HitInfo_t {
    uint32_t adc,adcbxm1;
    float mip;
    bool isTOA,isTDC,isBusy;
    bool passThr,passLZSThr,passMZSThr,passTZSThr,passThrBxm1,passTightThrBxm1;
  };
  

  typedef std::tuple<int,int,int,int> WaferKey_t;    
  typedef std::pair<int,int> UVKey_t;

  /**
     @class WaferOccupancyHisto
     @short book keeps a set of occupancy histograms for a wafer
  */
  class WaferOccupancyHisto
  {
  public:

    /**
       @short CTOR
    */
  WaferOccupancyHisto(WaferKey_t key) : isInit_(false),
      nEvents_(0),
      ncells_(0),
      waferType_(-1),
      avgNoise_(0.),
      avgGain_(0.),
      avgThr_(0.),
      avgS_(0.),
      avgSN_(0.)
        { 
          vars_.push_back("");            //nominal threshold
          vars_.push_back("lzs");         //loose zero suppression
          vars_.push_back("mzs");         //loose zero suppression
          vars_.push_back("tzs");         //tight zero suppression

          int subDet=std::get<0>(key);
          int layer=std::get<1>(key);
          int u=std::get<2>(key);
          int v=std::get<3>(key);      
          myID_ = Form("%d_lay%d_%d_%d",subDet,layer,u,v);
	  std::cout << "[WaferOccupancyHisto] created : " << myID_ << std::endl;
        }

    /**
       increments pad counts
    */
    void addPad(int waferType,HGCalSiNoiseMap::SiCellOpCharacteristics &siop) { 
      waferType_=waferType; 

      avgNoise_ += siop.core.noise;
      avgGain_  += siop.core.gain;
      avgThr_   += siop.core.thrADC;
      avgS_     += siop.mipfC;
      if(siop.core.noise>0)
        avgSN_    += siop.mipfC/siop.core.noise;
      else
        avgSN_  += 0.;
      ncells_++; 
    }

    /**
       averages the Si Op Charateristics of the pads in the wafer
    */
    std::vector<float> getAveragedOpCharacteristics() {

      std::vector<float> avg(5,0.);

      if(ncells_>0){
        avg[0]=avgNoise_/ncells_;
        avg[1]=avgGain_/ncells_;
        avg[2]=avgThr_/ncells_;
        avg[3]=avgS_/ncells_;
        avg[4]=avgSN_/ncells_;
      }
    
      return avg;
    }


    /**
       @short book all histograms
    */
    void bookHistos(edm::Service<TFileService> *fs) {

      if(isInit_ || fs==NULL) {
	if (fs==NULL) std::cout << "WaferOccupancyHisto bookHistos: no TFileService, returning for:" <<  myID_ << std::endl;
	if (isInit_)  std::cout << "WaferOccupancyHisto bookHistos: already initialised, returning for:" <<  myID_ << std::endl;
	return;	}
      // std::cout << "WaferOccupancyHisto bookHistos: 1st call for :" <<  myID_ << std::endl;

      // create one directory per wafer
      TFileDirectory mySubDir=(*fs)->mkdir(myID_.Data());
    
      //book histos and counters (nominal threshold, loose and tight ZS)
      busyCounts_=0;
      histos_["busycounts"] = mySubDir.make<TH1F>(myID_+"_busycounts",";Busy cells;",ncells_+1,0,ncells_+1);    
      histos_["mip"]        = mySubDir.make<TH1F>(myID_+"_mip",";q [MIP];",200,0,100);
      for(auto &v : vars_) {
        
        //adc spectra
        histos_["adcext"+v]              = mySubDir.make<TH1F>(myID_+"_adcext"+v, ";q [ADC];",100,0,500);    
        histos_["adcexttoa"+v]           = mySubDir.make<TH1F>(myID_+"_adcexttoa"+v, ";q [ADC];",100,0,500);    
        histos_["adc"+v]                 = mySubDir.make<TH1F>(myID_+"_adc"+v,";q [ADC];",100,0,100);
        histos_["adcbxm1"+v]             = mySubDir.make<TH1F>(myID_+"_adcbxm1"+v,";q [ADC];",100,0,100);
        histos_["rejadc"+v]              = mySubDir.make<TH1F>(myID_+"_rejadc"+v,";q [ADC];",50,0,50);
        histos_["rejadcbxm1"+v]          = mySubDir.make<TH1F>(myID_+"_rejadcbxm1"+v,";q [ADC];",50,0,50);
        histos2D_["adcbxm1_vs_adc"+v]    = mySubDir.make<TH2F>(myID_+"_adcbxm1_vs_adc"+v,";q_{BX-1} [ADC];q_{BX} [ADC];",50,0,50,50,0,50);

        //counters
        adcCounts_[v]=0;
        adcCounts_[v+"bxm1"]=0;
        adcCounts_[v+"bxm1tight"]=0;
        adcCounts_[v+"bxm1notight"]=0;
        toaCounts_[v]=0;
        tdcCounts_[v]=0;
        histos_["counts"+v]          = mySubDir.make<TH1F>(myID_+"_counts"+v,";Cells;",ncells_+1,0,ncells_+1);      
        histos_["countsbxm1"+v]      = mySubDir.make<TH1F>(myID_+"_countsbxm1"+v,";Cells;",ncells_+1,0,ncells_+1);      
        histos_["countsbxm1notight"+v] = mySubDir.make<TH1F>(myID_+"_countsbxm1notight"+v,";Cells;",ncells_+1,0,ncells_+1);      
        histos_["countsbxm1tight"+v]   = mySubDir.make<TH1F>(myID_+"_countsbxm1tight"+v,";Cells;",ncells_+1,0,ncells_+1);      
        histos_["toacounts"+v]       = mySubDir.make<TH1F>(myID_+"_toacounts"+v,";Cells in TOA;",ncells_+1,0,ncells_+1);
        histos_["tdccounts"+v]       = mySubDir.make<TH1F>(myID_+"_tdccounts"+v,";Cells in TDC;",ncells_+1,0,ncells_+1);
      }
        
      for(auto h : histos_) h.second->Sumw2();
    
      isInit_=true;
    }
  
    /**
       @short accumulate counts for a new hit
    */
    inline void count(HitInfo_t &h)
    {
      if(!isInit_)  {
	std::cout << "WaferOccupancyHisto count: not initialised or no TFileService, returning" << std::endl;
	return;	}

      if(h.isBusy) {
        busyCounts_++;
        return;
      }

      for(auto &v : vars_) {

        //in-time bunch threshold
        bool passThr(true);
        if(v==""             && !h.passThr)    passThr=false;
        if(v.Contains("lzs") && !h.passLZSThr) passThr=false;
        if(v.Contains("mzs") && !h.passMZSThr) passThr=false;
        if(v.Contains("tzs") && !h.passTZSThr) passThr=false;

        //monitor rejected hits
        uint32_t adc(h.adc);
        if(!passThr){
          histos_["rejadc"+v]->Fill(adc);
          histos_["rejadcbxm1"+v]->Fill(h.adcbxm1);
          continue;
        }

        histos_["mip"]->Fill(h.mip);

        //adc spectra
        if(!h.isTDC) {
          histos_["adcext"+v]->Fill(adc);
          if(h.isTOA) histos_["adcexttoa"+v]->Fill(adc);
          histos_["adc"+v]->Fill(adc);
          histos_["adcbxm1"+v]->Fill(h.adcbxm1);
          histos2D_["adcbxm1_vs_adc"+v]->Fill(h.adcbxm1,adc);
        }

        //counters
        adcCounts_[v]++;

        if(h.isTOA) toaCounts_[v]++;
        if(h.isTDC) tdcCounts_[v]++;

        //energy threshold for bx-1
        if(h.passThrBxm1) {
          adcCounts_[v+"bxm1"]++;
        }
        
        //toa-based threshold for bx-1
        if(h.passTightThrBxm1) {
          adcCounts_[v+"bxm1tight"]++;

        }
        else {
          adcCounts_[v+"bxm1notight"]++;
        }
      }
    }

    /**
       @short to be called once all hits have been counted
    */
    inline void analyze()
    {
      if(!isInit_) return;

      nEvents_++;

      histos_["busycounts"]->Fill(busyCounts_);
      for(auto &v : vars_) {
        histos_["counts"+v]->Fill(adcCounts_[v]);
        histos_["countsbxm1"+v]->Fill(adcCounts_[v+"bxm1"]);
        histos_["countsbxm1tight"+v]->Fill(adcCounts_[v+"bxm1tight"]);
        histos_["countsbxm1notight"+v]->Fill(adcCounts_[v+"bxm1notight"]);
        histos_["toacounts"+v]->Fill(toaCounts_[v]);
        histos_["tdccounts"+v]->Fill(tdcCounts_[v]);
      }    
    }

    /**
       @short set all to 0
    */
    inline void resetCounters()
    {
      if(!isInit_) return;

      for(auto &v : vars_) {
        adcCounts_[v]=0;
        adcCounts_[v+"bxm1"]=0;
        adcCounts_[v+"bxm1tight"]=0;
        adcCounts_[v+"bxm1notight"]=0;
        toaCounts_[v]=0;
        tdcCounts_[v]=0;
      }
      busyCounts_=0;
      
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
    inline int getWaferType() { return waferType_; }
    inline int getCounts(TString v="") { return adcCounts_[v]; };

    /**
       @short DTOR
    */
    ~WaferOccupancyHisto() {}

 
  private:

    TString myID_;
    bool isInit_;
    int nEvents_;
    int ncells_,waferType_;
    float avgNoise_,avgGain_,avgThr_,avgS_,avgSN_;

    std::vector<TString> vars_;
    std::map<TString, TH1F *> histos_;
    std::map<TString, TH2F *> histos2D_;
    std::map<TString,int> adcCounts_,toaCounts_,tdcCounts_;
    int busyCounts_;
  };

}

#endif
