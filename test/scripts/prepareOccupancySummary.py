import ROOT
import os
import optparse
import re
import pickle
from sensorEquivalentMap import *

PLOTTITLES={
    'maxcounts'      : 'Hottest wafer equivalent occupancy',
    'busycounts'     : 'Busy cells',
    'toacounts'      : 'TOA gain regime',
    'tdccounts'      : 'TDC in progress',
    'hottestwafer0'  : 'Hottest wafer',
    'hottestwafer1'  : '1^{st} hottest wafer neighbor',
    'hottestwafer2'  : '2^{nd} hottest wafer neighbor',
    'hottestwafer3'  : '3^{rd} hottest wafer neighbor',
    'hottestwafer4'  : '4^{th} hottest wafer neighbor',
    'hottestwafer5'  : '5^{th} hottest wafer neighbor',
    'hottestwafer6'  : '6^{th} hottest wafer neighbor',
    'adc'            : 'Energy (ADC)',
    'adcfull'        : 'Energy (ADC)',
    'counts'         : 'Occupancy',
    'genmatchcounts' : 'Occupancy (best gen match)'
}

garbageList=[]

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


def getPlotsIn(url,dirName,title,pfix):

    """reads all the plots in a given analysis sub-directory"""

    layerPlots={}
    waferPlots={}

    inF=ROOT.TFile.Open(url)
    d=inF.Get(dirName)

    #loop over directories
    for k in d.GetListOfKeys():

        sd=k.ReadObj()
        sdname=k.GetName()
        
        #parse sub-det, layer + u, v from name
        vals = [int(d) for d in re.findall(r'-?\d+', sdname)]

        if sd.InheritsFrom('TH1') :
            layerKey=tuple(vals[0:2])
            pname=sdname.split('_')[-1]           
            if pname !='counts' and pname !='adc' : continue
            if not layerKey in layerPlots: layerPlots[layerKey]={}
            layerPlots[layerKey][pname]=sd.Clone(hname+pfix)
            fixExtremities(layerPlots[layerKey][pname])
            layerPlots[layerKey][pname].SetDirectory(0)
            layerPlots[layerKey][pname].SetTitle(title)
            layerPlots[layerKey][pname].SetLineWidth(2)

        else:            
            waferKey=tuple(vals[0:4])
            waferPlots[waferKey]={}
            for kk in sd.GetListOfKeys():
                hname=kk.GetName()
                pname=hname.replace(sdname+'_','')
                if 'genmatch' in pname: continue
                waferPlots[waferKey][pname]=kk.ReadObj().Clone(hname+pfix)
                fixExtremities(waferPlots[waferKey][pname])
                waferPlots[waferKey][pname].SetDirectory(0)
                waferPlots[waferKey][pname].SetTitle(title)
                waferPlots[waferKey][pname].SetLineWidth(2)
            
    return (layerPlots,waferPlots)

def getQuantiles(plotColl,q=[0.5,0.9]):

    momentSummary=[]
    for i in range(len(plotColl)):
        momentSummary.append({})
        for p in plotColl[i]:
            plotColl[i][p].SetBinContent(1,0) #FIXME in the analyzer
            qval=np.array([0.]*len(q))
            if plotColl[i][p].Integral()>0:
                prob=np.array(q)
                plotColl[i][p].GetQuantiles(len(q),qval,prob)
            momentSummary[i][p]=[qval[k] for k in range(len(q))]

    return momentSummary

def showPlotCollSummary(plotColl,extraText,pname,fitPeak=False,plotCDF=False):

    nPlots=len(plotColl[0])
    pnames=plotColl[0].keys()

    c=ROOT.TCanvas('c','c',600,600)
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetLogy()
    c.SetGridy()
    garbageList.append(c)

    langau=None
    if fitPeak:
        langau=ROOT.TF1('langau','gaus(0)+landau(3)',0,5)
        langau.SetParLimits(0,0.,20.)
        langau.FixParameter(1,0.)
        langau.SetParLimits(2,0.2,1.)
        langau.SetParLimits(3,0.,20.)
        langau.SetParLimits(4,0.2,1.5)
        langau.SetParLimits(5,0.5,2.)

    colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]    
    histos=[]
    for idx in range(nPlots):

        plot=pnames[idx]
        plotTitle=PLOTTITLES[plot]
                    
        drawOpt='hist' 
        fitParams=[]
        leg=ROOT.TLegend(0.6,0.86,0.9,0.68)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.SetFillStyle(0)
        for i in range(len(plotColl)):
            histos.append( plotColl[i][plot] )
            histos[-1].Draw(drawOpt)
            histos[-1].SetLineColor(colors[i])
            histos[-1].GetXaxis().SetLabelSize(0.045)
            histos[-1].GetXaxis().SetTitleSize(0.045)
            histos[-1].GetYaxis().SetLabelSize(0.045)
            histos[-1].GetYaxis().SetTitleSize(0.045)
            histos[-1].GetYaxis().SetTitleOffset(1.1)
            histos[-1].GetYaxis().SetTitle('PDF or CDF' if plotCDF else 'CDF')
            histos[-1].GetYaxis().SetRangeUser(1e-4,1)
            leg.AddEntry(histos[-1],histos[-1].GetTitle(),'l')

            if fitPeak and plot=='adc':
                histos[-1].Fit(langau,'MLRQ+','same',0.5,2)    
                fitParams.append([langau.GetParameter(ipar) for ipar in [4,5,2]])    

            if plotCDF:
                histos[-1].SetLineWidth(3)
                histos.append( plotColl[i][plot].GetCumulative(False) )
                histos[-1].SetLineColor(colors[i])
                histos[-1].SetLineStyle(2)
                histos[-1].Draw('histsame')

            drawOpt='histsame'

        leg.Draw()

        #print fit parameters
        if fitPeak and plot=='adc':
            fittex=ROOT.TLatex()
            fittex.SetTextFont(42)
            fittex.SetTextSize(0.04)
            fittex.SetNDC()
            fittex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
            for ipar,ilabel in [(0,'MPV'),(1,'#sigma'),(2,'#sigma_{noise}')]:                
                label='%s=%s'%(ilabel,'/'.join(['%3.2f'%fres[ipar] for fres in fitParams]))
                fittex.DrawLatex(0.95,0.8-ipar*0.05, label)

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.05)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
        tex.SetTextSize(0.035)
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        for i in range(len(extraText)): 
            tex.DrawLatex(0.95,0.65-i*0.045,extraText[i])

        titletex=ROOT.TLatex()
        titletex.SetTextFont(62)
        titletex.SetTextSize(0.045)
        titletex.SetNDC()
        titletex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        titletex.DrawLatex(0.95,0.9,plotTitle)

        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s_%s.%s'%(pname,plot,ext))


def main():

    #parse inputs
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)  
    parser.add_option('-o', '--out',            dest='output',        help='output directory [%default]',                 default='occ_plots', type='string')
    parser.add_option(      '--noSummary',      dest='noSummary',     help='skip creation of pck file [%default]',        default=False, action='store_true')    
    parser.add_option(      '--waferPlots',     dest='waferPlots',    help='show these wafer plots [%default]',             default=None, type='string')
    parser.add_option(      '--noGlobalPlots',  dest='noGlobalPlots', help='disable global plots [%default]',             default=False, action='store_true')
    parser.add_option(      '--noCountHistosSummary',  dest='noCountHistosSummary', help='disable saving the count histos [%default]',             
                            default=False, action='store_true')
    parser.add_option(      '--fitPeak',        dest='fitPeak',       help='enable gauss+noise fit to peak [%default]',  default=False, action='store_true')
    (opt, args) = parser.parse_args()

    #decode wafer plot list string 
    if opt.waferPlots:
        opt.waferPlots=[[int(y) for y in x.split(',')] for x in opt.waferPlots.split(':')]
        print opt.waferPlots

    #prepare output/root style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    os.system('mkdir -p '+opt.output)

    #get the sensor position map
    cmssw=os.environ['CMSSW_BASE']
    sensorPos=parseWaferPositionsMap(url='%s/src/UserCode/HGCElectronicsValidation/data/wafer_pos.dat'%cmssw)

    #define inputs
    procList=[ x.split(':') for x in args ]    

    #get the plots
    global_plots=[]
    plots=[]
    for i in range(len(args)):
        title,url,dirName=procList[i]
        layerPlots,waferPlots=getPlotsIn(url=url,dirName=dirName,title=title,pfix='_%d'%i)
        global_plots.append(layerPlots)
        plots.append(waferPlots)

    #compute the quantiles for the global plots
    globalQuantilesSummary={}
    for layerKey in global_plots[0]:
        plotColl=[ x[layerKey] for x in global_plots]
        globalQuantilesSummary[layerKey]=getQuantiles(plotColl=plotColl,q=[0.1,0.5,0.9])

        #plot if required
        if opt.noGlobalPlots : continue
        pname='summary_sd%d_lay%d'%layerKey
        extraText=[
            '%s layer %d'%('CEE' if layerKey[0]==0 else 'CEH', layerKey[1]),
        ]
        showPlotCollSummary(plotColl=plotColl,
                            extraText=extraText,
                            pname=os.path.join(opt.output,pname))
        
    #compute the quantiles for the wafer plots
    quantilesSummary={}
    for waferKey in plots[0]:
        if not waferKey in sensorPos : continue

        ncells,r,z,eta,phi=sensorPos[waferKey]
        if abs(eta)<1.5 or abs(eta)>3.0 : continue
        plotColl=[ x[waferKey] for x in plots]

        quantilesSummary[waferKey]=getQuantiles(plotColl=plotColl,q=[0.1,0.5,0.9])
        
        #skip if wafer plots were not required
        if not opt.waferPlots: continue
        waferUV=[waferKey[2],waferKey[3]]
        if not waferUV in opt.waferPlots: continue
            
        pname='sd%d_lay%d_%d_%d'%waferKey
        pname=pname.replace('-','m')
        extraText=[
            '%s layer %d'%('CEE' if waferKey[0]==0 else 'CEH', waferKey[1]),
            '(u,v)=(%s,%s)'%(waferKey[2],waferKey[3]),
            'R=%3.2f z=%3.2f'%(r,z),
            '#eta=%3.2f #phi=%3.2f'%(eta,phi)
        ]
        showPlotCollSummary(plotColl=plotColl,
                            extraText=extraText,
                            pname=os.path.join(opt.output,pname),
                            fitPeak=opt.fitPeak,
                            plotCDF=True)

    #save summary
    if not opt.noSummary:
        with open(os.path.join(opt.output,'summary.pck'),'w') as cache:
            pickle.dump(quantilesSummary,cache,pickle.HIGHEST_PROTOCOL)
            pickle.dump(globalQuantilesSummary,cache,pickle.HIGHEST_PROTOCOL)
            pickle.dump(sensorPos,cache,pickle.HIGHEST_PROTOCOL)
            
    #count histograms (for P. Bloch's FIFO studies)
    if not opt.noCountHistosSummary:
        countHistos={}
        maxCountHistos={}
        firstKey=plots[0]
        for waferKey in waferPlots:
            h=waferPlots[waferKey]['counts']
            y=np.array([h.GetBinContent(xbin+1) for xbin in range(h.GetNbinsX())])
            norm=h.Integral()
            if norm>0 : y /= norm
            countHistos[waferKey]=y

            hmax=waferPlots[waferKey]['maxcounts']
            ymax=np.array([hmax.GetBinContent(xbin+1) for xbin in range(hmax.GetNbinsX())])
            norm=hmax.Integral()
            if norm>0: ymax /=norm
            maxCountHistos[waferKey]=ymax

        with open(os.path.join(opt.output,'counts.pck'),'w') as cache:
            pickle.dump(countHistos,cache,pickle.HIGHEST_PROTOCOL)
            pickle.dump(maxCountHistos,cache,pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()
