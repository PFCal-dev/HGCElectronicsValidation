import ROOT
import os
import optparse
import re
import pickle
import pandas as pd
import sys
import numpy as np

SIOPCOLS=['waferPreChoice','npads','waferSN','waferGain','waferThr']
DATACOLS=['counts',     'countslzs',     'countsmzs',    'countstzs',
          'toacounts',  'toacountslzs',  'toacountsmzs', 'toacountstzs',
          'countsbxm1', 'countsbxm1lzs', 'countsbxm1mzs', 'countsbxm1tzs', 
          'countsbxm1tight',    'countsbxm1notight',
          'countsbxm1tightlzs', 'countsbxm1notightlzs',
          'countsbxm1tightmzs', 'countsbxm1notightmzs',
          'countsbxm1tighttzs', 'countsbxm1notighttzs',
          'busycounts']
GEOMCOLS=['waferX','waferY','waferShape','waferRot','rho','phi','waferThickness']

def parseWaferGeometry(geomFile='data/geomnew_corrected_withmult_F_rotations_v11.1.txt'):

    """reads the ascii geometry file and converts it to a DataFrame with the geometry information"""

    df = pd.read_csv(geomFile, names=['layer','waferShape','waferThickness','waferX','waferY','waferRot','waferU','waferV'],skiprows=11,delim_whitespace=True)
    df=df[['layer','waferThickness','waferU','waferV','waferX','waferY','waferShape','waferRot']]
    df['rho']=np.sqrt(df['waferX']**2+df['waferY']**2)
    df['phi']=np.arctan2(df['waferY'],df['waferX'])*180./np.pi;
    return df


def fixExtremities(h,addOverflow=True,addUnderflow=True,norm=True):

    """increments the first and the last bin to show the under- and over-flows and normalizes the histogram"""

    if addUnderflow :
        fbin  = h.GetBinContent(0) + h.GetBinContent(1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(0)*h.GetBinError(0) + h.GetBinError(1)*h.GetBinError(1))
	h.SetBinContent(1,fbin)
	h.SetBinError(1,fbine)
	h.SetBinContent(0,0)
	h.SetBinError(0,0)
    if addOverflow:
        nbins = h.GetNbinsX();
	fbin  = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)*h.GetBinError(nbins)  + h.GetBinError(nbins+1)*h.GetBinError(nbins+1))
	h.SetBinContent(nbins,fbin)
	h.SetBinError(nbins,fbine)
	h.SetBinContent(nbins+1,0)
	h.SetBinError(nbins+1,0)
    if norm :
        n=h.Integral()
        if n>0:
            h.Scale(1./n)

def getQuantiles(p,q=[0.1,0.5,0.9]):

    """computes the quantiles for a plot"""

    nq=len(q)
    qval=np.array([0.]*nq)
    if p.Integral()>0:
        prob=np.array(q)
        p.GetQuantiles(nq,qval,prob)

    #subtract the 1/2 bin width=0.5
    return [max(qval[k]-0.5,0.) for k in range(nq)]

def getMean(p):

    """computes the quantiles for a plot"""

    #subtract the 1/2 bin width=0.5
    return [p.GetMean()-0.5,p.GetMeanError()]


def convertHistosToNumpy(plotColl):
    
    """converts a TH1 to a numpy array"""

    histos=[]
    for h in plotColl:
        y=np.array([h.GetBinContent(xbin+1) for xbin in range(h.GetNbinsX())])
        norm=h.Integral()
        if norm>0 : y /= norm
        histos.append(y)

    return histos


def analyzeOccupancyPlots(url,modulePos,siop,dirName='ana'):

    """reads the occupancy plots in a given analysis sub-directory and summarizes the quantiles"""

    occData=[]

    inF=ROOT.TFile.Open(url)
    d=inF.Get(dirName)

    #loop over directories
    for k in d.GetListOfKeys():

        sd=k.ReadObj()
        if not sd.InheritsFrom('TDirectory') : continue
        sdname=k.GetName()

        #parse sub-det, layer + u, v from name
        subdet,layer,u,v=[int(d) for d in re.findall(r'-?\d+', sdname)]

        #get wafer averaged Si operation
        mask     = (siop.section==subdet) & (siop.layer==layer) & (siop.waferU==u) & (siop.waferV==v)
        waf_siop = siop.loc[mask].iloc[0]

        #update layer index to be continuous
        if subdet==1:
            layer=layer+28
            
        #find module info (and discard if not available)
        mask=(modulePos.layer==layer) & (modulePos.waferU==u) & (modulePos.waferV==v)
        try:
            moduleGeom=modulePos.loc[mask].iloc[0]
        except:
            continue

        #build a summary based on quantiles
        qSummary =  [layer,u,v]
        qSummary += [moduleGeom[x] for x in GEOMCOLS]
        qSummary += [waf_siop[x] for x in SIOPCOLS]

        for c in DATACOLS:
            h=sd.Get(sdname+'_'+c)
            fixExtremities(h)
            qSummary += getQuantiles(h,q=[0.1,0.5,0.9])
            qSummary += getMean(h)

        occData.append( qSummary )

    return occData


def main():

    #get the sensor position map/wafer geometry etc
    cmssw=os.environ['CMSSW_BASE']
    modulePos = parseWaferGeometry(geomFile='%s/src/UserCode/HGCElectronicsValidation/data/geomnew_corrected.txt'%cmssw)

    url=sys.argv[1]

    #si-operation summary
    fIn=ROOT.TFile.Open(url)
    siop=ROOT.RDataFrame(fIn.Get('ana/data'))
    columns=['section','layer','waferU','waferV','waferPreChoice','npads']
    columns+=['waferNoise','waferGain','waferThr','waferS','waferSN']
    siop=siop.AsNumpy()
    siop=pd.DataFrame(data=siop, columns=columns)

    #summarize the occupancy for all wafers based on the CMSSW simulation
    occData=analyzeOccupancyPlots(url=url,modulePos=modulePos,siop=siop)        
    columns=['layer','waferU','waferV']
    columns+=GEOMCOLS
    columns+=SIOPCOLS
    for c in DATACOLS:
        if 'counts' in c :
            columns += [c+'_q10',c+'_q50',c+'_q90']
            columns += [c+'_mean',c+'_meanerr']
        else :
            columns +=[c]
    result=pd.DataFrame(data=np.array(occData), columns=columns)
    for c in columns:
        if c!='waferShape':
            result[c]=pd.to_numeric(result[c])

    result.rename(columns={'waferPreChoice':'waferPreChoiceCMSSW',
                           'waferThickness':'waferPreChoice'},
                  inplace=True)

    result.replace({'waferPreChoice':{120:0,200:1,300:2}},
                   inplace=True)

    result.to_hdf(url.replace('.root','.h5'), 
                  key='data', 
                  mode='w')

    print 'Summary with shape',result.shape
    print 'and columns',result.columns
    print result.head()
    print(result.dtypes)

if __name__ == "__main__":
    main()
