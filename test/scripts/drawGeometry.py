import ROOT
import sys

def drawLayer(gr):
    """draws the graph of sensors of a layer"""

    gr.Draw('apcol')
    gr.GetXaxis().SetTitle('x [cm]')
    gr.GetYaxis().SetTitle('y [cm]')

    tkns=gr.GetName().split('_')
    title='%s, layer %s'%(tkns[0],tkns[1].replace('lay',''))

    l=ROOT.TLine()
    l.DrawLine(0,-150,0,150)
    l.DrawLine(-150,0,150,0)


    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.045)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.12,0.92,title)
    tex.SetTextSize(0.04)
        
ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)    
c=ROOT.TCanvas('c','c',600,600) 

c.SetTheta(90)
c.SetPhi(0)

fIn=ROOT.TFile.Open(sys.argv[1])
for key in fIn.Get('ana').GetListOfKeys():
    drawLayer(key.ReadObj())
    c.Modified()
    c.Update()
    c.SaveAs(key.GetName()+'.png')
