## Deconvoluting radiation effects and gain setting at ROC level

The following intructions gather info and code snippets to handle the radiation-induced effects taken into account by the HGCAL digitisation in producing DIGI collections,
in order to produce derived collections siuch as ```recHits``` or ```trigger primitives```.
One can find in [this talk](https://indico.cern.ch/event/933714/contributions/3924245/) a comprehensive description of the model used to emulate the electronics;
in a nutshell, these instructions give guidance on how to get hold of the radiation scenario used in the DIGI step, and how to use it to present the amplitudes in units of MIP, taking consistenly into account the electronics gain.

It can be used to recompute the MIP from a digi correcting for the loss of charge collection efficiency and for the choice of the gain in the ROC.
The examples below are given for CE-E only.



### Configuring the radiation map tool

In the first place one needs to re-instantate the radiation map class which is used at digi level:

```
HGCalSiNoiseMap *rad_map=new HGCalSiNoiseMap;
rad_map->setDoseMap(iConfig.getParameter<std::string>("doseMap"),
                   iConfig.getParameter<uint32_t>("scaleByDoseAlgo"));
rad_map->setFluenceScaleFactor(iConfig.getParameter<double>("scaleByDoseFactor"));
rad_map->setIleakParam(iConfig.getParameter<edm::ParameterSet>("ileakParam").template getParameter<std::vector<double>>("ileakParam"));
rad_map->setCceParam(iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamFine"),
                     iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThin"),
                     iConfig.getParameter<edm::ParameterSet>("cceParams").template getParameter<std::vector<double>>("cceParamThick"));


//save also the TDC config parameters for later calibration
tdcOnsetfC_ = iConfig.getParameter<double>("tdcOnset");
tdcLSB_     = iConfig.getParameter<double>("tdcSaturation") / pow(2., iConfig.getParameter<uint32_t>("tdcNbits") );
```

The configuration parameters needed translate the operation conditions and must be aligned with the ones used at digitization level.
To use the default parameters one needs to add something like the following the configuration file.
Notice that the variable `integLumi` below needs to be set according to the sample being analyzed as the correction to the charge collection 
efficiency obtained will differ.

```
integLumi=3000.
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse,HGCAL_cceParams_toUse,hgceeDigitizer
process.foo = cms.EDAnalyzer("barAnalyzer",
                              ...
                              doseMap           = cms.string('SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt'),
                              scaleByDoseAlgo   = cms.uint32(0),
                              scaleByDoseFactor = cms.double(integLumi/3000.),
                              ileakParam        = HGCAL_ileakParam_toUse,
                              cceParams         = HGCAL_cceParams_toUse,
                              tdcNbits          = hgceeDigitizer.digiCfg.feCfg.tdcNbits,
                              tdcSaturation     = hgceeDigitizer.digiCfg.feCfg.tdcSaturation_fC,
                              tdcOnset          = hgceeDigitizer.digiCfg.feCfg.tdcOnset_fC,
                          )
```

### Producing uncalibrated hits (MIP units)

With the radiation map properly configured one can do the following when analyzing the digi collection:

```
//assign the geometry and tell the tool that the gain is automatically set to have the MIP close to 10ADC counts
rad_map->setGeometry(geometry, HGCalSiNoiseMap::AUTO, 10);
for(auto &hit : *digiColl)
    {
      if(hit.size()==0) continue;

      //digi information
      HGCSiliconDetId detId(hit.id());
      bool isTDC( hit.sample(itSample).mode() );
      bool isBusy( isTDC && rawData==0 );
      HGCalSiNoiseMap::GainRange_t gain( hit.gain() );
      double rawADC( double(hit.sample(itSample).data()) );


      //get the operation characterisics of this detId
      //and the nubmer of ADC counts corresponding to a MIP
      HGCalSiNoiseMap::SiCellOpCharacteristics siop=rad_map->getSiCellOpCharacteristics(detId);
      double mipADC=double(siop.mipADC);

      //number of MIPs
      double nmips( rawADC/mipADC );
      if(isTDC) {
                double adcLSB( 1./80.);
                if(gain==HGCalSiNoiseMap::q160fC) adcLSB=1./160.;
                if(gain==HGCalSiNoiseMap::q320fC) adcLSB=1./320.;
                 
                double charge( (std::floor(tdcOnsetfC_ / adcLSB) + 1.0) * adcLSB + (rawADC+0.5)*tdcLSB_ );
                nmips = rawCharge/double(mipADC);
      }
      if(isBusy) {
        nmips=0;
      }


      ...

}
```
