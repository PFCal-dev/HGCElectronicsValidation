import ROOT
import os

def cmsHeader():
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{preliminary}')
    tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignBottom)
    tex.DrawLatex(0.97,0.96,'37.5 fb^{-1} (13 TeV)')  


baseDir='/eos/cms/store/cmst3/user/psilva/HGCal/Occupancies'


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetRightMargin(0.02)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)


histos=[]
for f,title,ci in [ #('FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_noNoise_20191010.root','noOOT,noNoise',ROOT.kGray),
                    #('FlatRandomPtGunProducer_NeutrinoGun_v11_aged_noOOT_20191010.root','noNoise',ROOT.kGreen+1),
                    ('FlatRandomPtGunProducer_NeutrinoGun_v11_aged_pu200_20191010.root','aged',ROOT.kAzure+1)]:
    fIn=ROOT.TFile.Open(os.path.join(baseDir,f))
    histos.append(fIn.Get('ana/adchitsvspu').Clone(title.replace(',','_')) )
    histos[-1].SetTitle(title)
    histos[-1].SetDirectory(0)
    histos[-1].SetLineColor(ci)
    fIn.Close()

for i in range(len(histos)):
    histos[i].Draw('box' if i==0 else 'boxsame')

leg=c.BuildLegend(0.15,0.94,0.4,0.8)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)

cmsHeader()
c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('adchitsvspu.'+ext)
