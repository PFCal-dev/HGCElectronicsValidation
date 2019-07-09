import sys
import ROOT

def getPlotsFrom(t,url):

    histos={}
    fIn=ROOT.TFile.Open(url)
    for k in fIn.GetListOfKeys():
        kname=k.GetName()
        histos[kname]=k.ReadObj()
        histos[kname].SetDirectory(0)
        if 'emf' in kname: histos[kname].GetYaxis().SetTitle('<em energy fraction>')
        if 'hmf' in kname: histos[kname].GetYaxis().SetTitle('<had. energy fraction>')
        histos[kname].SetName(kname+'_'+t)
        histos[kname].SetTitle(t)
    fIn.Close()

    return histos

histos=[]
pList=[x.split('=') for x in sys.argv[1:]]
histos=[getPlotsFrom(t,url) for t,url in pList]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
colors=[ROOT.kBlack,     ROOT.kGray, ROOT.kMagenta,   ROOT.kMagenta+2,  ROOT.kMagenta-9, ROOT.kRed+1, ROOT.kAzure+7,   ROOT.kBlue-7]    
markers=[24,20,25,21,26,22,32,23]
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)
c.cd()
nprocs=len(histos)
for plot in histos[0].keys():

    frame=histos[0][plot].Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    maxY=1e-3

    isTProfile = True if frame.InheritsFrom('TProfile') else False    

    leg=ROOT.TLegend(0.15,0.95,0.96,0.89)
    leg.SetNColumns(nprocs)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    for i in range(nprocs):

        h=histos[i][plot]
        drawOpt='histsame'
        if isTProfile:
            h.SetMarkerStyle(markers[i])
            h.SetMarkerColor(colors[i])
            drawOpt='e1same'
        else:
            h.Scale(1./h.Integral())
        h.SetLineWidth(2)
        h.SetLineColor(colors[i])
        h.Draw(drawOpt)
        maxY=max(maxY,h.GetMaximum())
        leg.AddEntry(h,h.GetTitle(),'ep' if isTProfile else 'l')


    frame.GetYaxis().SetTitleOffset(0.8)
    frame.GetYaxis().SetRangeUser(1e-5,1.2*maxY)
    leg.Draw()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')        
    c.SetLogy(True if 'pt' in plot else False)
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(plot,ext))
