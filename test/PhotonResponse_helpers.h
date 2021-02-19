#ifndef __photonresponse_helpers_h__ 
#define __photonresponse_helpers_h__
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps; 
using rvec_i = const RVec<int>;
using rvec_f = const RVec<float>;

double dedx_wgt[]={0.0, 8.894541, 10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 
                   10.932882, 10.932882, 10.937907, 10.937907, 10.938169,10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
                   10.938169, 10.938169, 10.938169, 10.938169, 10.938169,10.938169, 10.938169, 10.938169, 32.332097, 51.574301, 
                   51.444192, 51.444192, 51.444192, 51.444192, 51.444192,51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
                   69.513118, 87.582044, 87.582044, 87.582044, 87.582044, 87.582044, 87.214571, 86.888309, 86.92952, 86.92952, 
                   86.92952};

/**
    @short applies dE/dx and CCE scaling to a list of hits
    The weights are used to scale by dE/dx in each layer and are directly retrieved from 
    [HGCalRecHit_cfi.py](https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py)
    The additional corrections above are stored in the ntuples per hit and named generically as CCE. 
    They must be applied to remove the losses induced by radiation 
    and also the variations in charge induced by the tile and SiPM areas for the scintillators.
*/
float scaledHit(const int &layer,const float &hit,const float &cce,bool applyDeDx=true, bool applyCCE=true)
{
    float newHit(hit);
    if (applyDeDx) newHit *= dedx_wgt[ layer ];
    if (applyCCE)  newHit /= (cce>0 ? cce : 1.);
    return newHit;
}

/**
*/
float swCompensatedHit(const float &hit, const int &thick, const bool &isSci,bool mpv=true)
{
     if(isSci)   return hit;
     if(hit<=0.) return hit;
    
     float R(1.0);
     if(mpv) {
         if(thick==0) R=0.08572457*TMath::Log(hit)+0.82158564;
         if(thick==1) R=0.07552747*TMath::Log(hit)+0.84078691;
         if(thick==2) R=0.06771299*TMath::Log(hit)+0.85585088;
     } else {
         R=0.05090596*TMath::Log(hit)+0.88952038;
     }

    return hit/R;
}
#endif
