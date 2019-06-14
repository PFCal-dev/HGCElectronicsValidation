import ROOT
import os
import optparse
import re
import pickle
from sensorEquivalentMap import *

PLOTTITLES={
    'maxcounts'     : 'Hottest wafer equivalent occupancy',
    'hottestwafer0' : 'Hottest wafer',
    'hottestwafer1' : '1^{st} hottest wafer neighbor',
    'hottestwafer2' : '2^{nd} hottest wafer neighbor',
    'hottestwafer3' : '3^{rd} hottest wafer neighbor',
    'hottestwafer4' : '4^{th} hottest wafer neighbor',
    'hottestwafer5' : '5^{th} hottest wafer neighbor',
    'hottestwafer6' : '6^{th} hottest wafer neighbor',
    'adc'           : 'Energy (MIP eq.)',
    'counts'        : 'Occupancy'
}

garbageList=[]

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
            if not layerKey in layerPlots: layerPlots[layerKey]={}
            layerPlots[layerKey][pname]=sd.Clone(hname+pfix)
            layerPlots[layerKey][pname].SetDirectory(0)
            layerPlots[layerKey][pname].SetTitle(title)
            layerPlots[layerKey][pname].SetLineWidth(2)

        else:            
            waferKey=tuple(vals[0:4])
            waferPlots[waferKey]={}
            for kk in sd.GetListOfKeys():
                hname=kk.GetName()
                pname=hname.replace(sdname+'_','')
                waferPlots[waferKey][pname]=kk.ReadObj().Clone(hname+pfix)
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

def showPlotCollSummary(plotColl,extraText,pname,fitPeak=False,nPerRow=2):

    nPlots=len(plotColl[0])
    pnames=plotColl[0].keys()
    nx=nPerRow
    ny=nPlots/nPerRow
    if nx*ny<nPlots : ny+=1
    dx=1./float(nPerRow)
    dy=1./float(ny)

    c=ROOT.TCanvas('c','c',800,ny*400)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)
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

    pads=[]
    colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]    
    for ix in range(0,nx):
        for iy in range(0,ny):
            c.cd()
            idx=len(pads)
            if idx>=len(pnames) : continue
            plot=pnames[idx]
            plotTitle=PLOTTITLES[plot]
            pads.append( ROOT.TPad(plot,plot,ix*dx,1-iy*dy,(ix+1)*dx,1-(iy+1)*dy) )            
            pads[-1].SetTopMargin(0.06)
            pads[-1].SetLeftMargin(0.12)
            pads[-1].SetRightMargin(0.02)
            pads[-1].SetBottomMargin(0.12)
            pads[-1].Draw()
            pads[-1].cd()
            pads[-1].SetLogy()
            
            drawOpt='hist' 
            fitParams=[]
            for i in range(len(plotColl)):
                plotColl[i][plot].Draw(drawOpt)

                if fitPeak and plot=='adc':
                    plotColl[i][plot].Fit(langau,'MLRQ+','same',0.5,2)    
                    fitParams.append([langau.GetParameter(ipar) for ipar in [4,5,2]])    

                plotColl[i][plot].SetLineColor(colors[i])
                plotColl[i][plot].GetXaxis().SetLabelSize(0.05)
                plotColl[i][plot].GetXaxis().SetTitleSize(0.05)
                plotColl[i][plot].GetYaxis().SetLabelSize(0.05)
                plotColl[i][plot].GetYaxis().SetTitleSize(0.05)
                drawOpt='histsame'

            #print fit parameters
            if fitPeak and plot=='adc':
                fittex=ROOT.TLatex()
                fittex.SetTextFont(42)
                fittex.SetTextSize(0.04)
                fittex.SetNDC()
                fittex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
                for ipar,ilabel in [(0,'MPV'),(1,'#sigma'),(2,'#sigma_{noise}')]:
                    print fitParams
                    label='%s=%s'%(ilabel,'/'.join(['%3.2f'%fres[ipar] for fres in fitParams]))
                    fittex.DrawLatex(0.95,0.8-ipar*0.05, label)

            if ix==0 and iy==0:
                leg=pads[-1].BuildLegend(0.6,0.86,0.9,0.68)
                leg.SetBorderSize(0)
                leg.SetTextFont(42)
                leg.SetTextSize(0.05)
                leg.SetFillStyle(0)

                tex=ROOT.TLatex()
                tex.SetTextFont(42)
                tex.SetTextSize(0.05)
                tex.SetNDC()
                tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
                tex.SetTextSize(0.04)
                for i in range(len(extraText)): 
                    tex.DrawLatex(0.15,0.9-i*0.05,extraText[i])

            titletex=ROOT.TLatex()
            titletex.SetTextFont(62)
            titletex.SetTextSize(0.05)
            titletex.SetNDC()
            titletex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
            titletex.DrawLatex(0.95,0.9,plotTitle)

    c.Modified()
    c.Update()
    c.SaveAs(pname+'.png')


def main():

    #parse inputs
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)  
    parser.add_option('-o', '--out',            dest='output',        help='output directory [%default]',                default='occ_plots', type='string')
    parser.add_option(      '--noSummary',      dest='noSummary',     help='skip creation of pck file [%default]',             default=False, action='store_true')    
    parser.add_option(      '--noWaferPlots',   dest='noWaferPlots',  help='disable wafer plots [%default]',             default=False, action='store_true')
    parser.add_option(      '--noGlobalPlots',  dest='noGlobalPlots', help='disable global plots [%default]',             default=False, action='store_true')
    parser.add_option(      '--fitPeak',        dest='fitPeak',       help='enable gauss+noise fit to peak [%default]',  default=False, action='store_true')
    (opt, args) = parser.parse_args()

    #prepare output/root style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    os.system('mkdir -p '+opt.output)

    #get the sensor position map
    cmssw=os.environ['CMSSW_BASE']
    sensorPos=parseWaferPositionsMap(url='%s/src/UserCode/HGCElectronicsValidation/test/scripts/wafer_pos.dat'%cmssw)

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
                            pname=os.path.join(opt.output,pname),
                            nPerRow=3)
        
    #compute the quantiles for the wafer plots
    quantilesSummary={}
    for waferKey in plots[0]:
        if not waferKey in sensorPos : continue

        ncells,r,z,eta,phi=sensorPos[waferKey]
        if abs(eta)<1.5 or abs(eta)>3.0 : continue
        plotColl=[ x[waferKey] for x in plots]

        quantilesSummary[waferKey]=getQuantiles(plotColl=plotColl,q=[0.1,0.5,0.9])
        
        if opt.noWaferPlots: continue
        pname='summary_sd%d_lay%d_%d_%d'%waferKey
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
                            fitPeak=opt.fitPeak
                            )

    #save summary
    if opt.noSummary: return
    with open(os.path.join(opt.output,'summary.pck'),'w') as cache:
        pickle.dump(quantilesSummary,cache,pickle.HIGHEST_PROTOCOL)
        pickle.dump(globalQuantilesSummary,cache,pickle.HIGHEST_PROTOCOL)
        pickle.dump(sensorPos,cache,pickle.HIGHEST_PROTOCOL)
            

if __name__ == "__main__":
    main()
