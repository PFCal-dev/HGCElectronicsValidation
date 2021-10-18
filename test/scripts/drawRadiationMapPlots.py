import ROOT
import sys
import os

colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
markers=[20,24,21,25,22,26]

def drawHeader(title=None):

    """ a lazy header for the plots """

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(ROOT.kHAlignCenter+ROOT.kVAlignCenter)
    txt.SetTextSize(0.035)
    if title:
        txt.DrawLatex(0.5,0.97,title)

def makeSciPlotsFrom(fIn,outName):

    """receives a TDirectoryFile with the plots from an analyzer and saves them in png/pdf"""

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)

    #2d maps
    plots=['tilecount',
           'doseMap',
           'fluenceMap',
           'scaleByDoseMap',
           'scaleByTileAreaMap',
           'scaleByDoseAreaMap',
           'noiseByFluenceMap', 
           'expNoiseMap',
           'probNoiseAboveHalfMip',
           'signalToNoiseFlatAreaMap',
           'signalToNoiseDoseMap',
           'signalToNoiseAreaMap',
           'signalToNoiseDoseAreaMap',
           'signalToNoiseDoseAreaSipmMap',
           'saturationMap']

    for p in plots:
        
        c.Clear()
        c.SetTopMargin(0.05)
        c.SetRightMargin(0.12)
        h=fIn.Get('plotter/{}'.format(p))
        h.SetTitle('')
        h.Draw('colz')
        h.GetZaxis().SetTitleOffset(-0.3)
        drawHeader(h.GetTitle())
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs(outName+'/sipontile_%s.%s'%(p,ext))



def makePlotsFrom(key,outName):

    """receives a TDirectoryFile with the plots from an analyzer and saves them in png/pdf"""

    plots    = ['sn','cce','noise','ileak','fluence','gain','mippeak']
    rangeMin = [ 0,   0,    0,      0,      0,        0.5,   3,      ]
    rangeMax = [ 12,  1.1,  0.7,    15,     10e15,    3.5,   20,     ]
    
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)

    tag=key.GetName()
    tag=tag.replace('plotter_','')

    for n, p in enumerate(plots):
        for d in [8,9]:

            #2d map
            c.Clear()
            c.SetTopMargin(0.05)
            c.SetRightMargin(0.12)
            h=key.Get('d%d_%s'%(d,p))                        
            h.Draw('colz')
            h.GetZaxis().SetTitleOffset(-0.3)
            drawHeader(h.GetTitle())
            c.Modified()
            c.Update()
            for ext in ['png','pdf']:
                c.SaveAs(outName+'/%s_%s_%s.%s'%(p,d,tag,ext))

            if 'fluencevssn' in p or 'fluencecount' in p : continue

            #per layer overview
            c.SetRightMargin(0.03)
            c.SetTopMargin(0.3)
            layers=range(1,29) if d==8 else range(1,23)            
            grs=[]
            mg=ROOT.TMultiGraph()
            xtitle,ytitle='',''
            for l in layers:

                pname='d%d_layer%d_%s'%(d,l,p)
                h=key.Get(pname)
                cellh=key.Get('d%d_layer%d_ncells'%(d,l))
                try:
                    cellh.Integral()
                except:
                    print(f'Could not get count for layer={l} det={d}')
                    continue

                xtitle=h.GetXaxis().GetTitle()
                ytitle=h.GetYaxis().GetTitle()

                #use only points where there were some sensors
                grs.append( ROOT.TGraphErrors( ) )
                grs[-1].SetLineColor(colors[(l-1)%len(colors)])
                grs[-1].SetMarkerColor(grs[-1].GetLineColor())
                grs[-1].SetMarkerStyle(markers[(l-1)%len(markers)])
                grs[-1].SetFillStyle(0)
                grs[-1].SetFillColor(0)
                grs[-1].SetName(pname)
                grs[-1].SetTitle('Layer %d'%l)
                for xbin in range(cellh.GetNbinsX()):
                    if cellh.GetBinContent(xbin+1)==0 : continue
                    np=grs[-1].GetN()
                    val=h.GetBinContent(xbin+1)
                    valUnc=h.GetBinError(xbin+1)
                    grs[-1].SetPoint(np,h.GetXaxis().GetBinCenter(xbin+1),val)
                    grs[-1].SetPointError(np,0,valUnc)
                mg.Add(grs[-1],'p')

            mg.Draw('ap')

            mg.GetXaxis().SetTitle(xtitle)
            mg.GetYaxis().SetTitle(ytitle)
            mg.GetYaxis().SetTitleOffset(1.2)
            mg.GetYaxis().SetRangeUser(rangeMin[n],rangeMax[n])
            leg=c.BuildLegend(0.12,0.72,0.95,0.92)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.03)
            leg.SetNColumns(5)
            leg.SetBorderSize(0)
            drawHeader()
            c.Modified()
            c.Update()
            for ext in ['png','pdf']:
                c.SaveAs(outName+'/%s_%s_perlayer_%s.%s'%(p,d,tag,ext))

def main():

    url='dosemap_output.root'
    if len(sys.argv)>1 :
        url=str(sys.argv[1])
    isSci=False
    if len(sys.argv)>2 and sys.argv[2]=='sci':
        isSci=True
    outName=url[:-5]
    if not os.path.isdir(outName):
        os.mkdir(outName)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    fIn=ROOT.TFile.Open(url)
    if isSci:
        makeSciPlotsFrom(fIn,outName)
    else:
        for key in fIn.GetListOfKeys(): 
            makePlotsFrom(key.ReadObj(),outName)

if __name__ == "__main__":
    main()
