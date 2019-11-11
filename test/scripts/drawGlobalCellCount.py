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

baseDir='/eos/cms/store/cmst3/user/psilva/HGCal/Occupancies/4Nov/'

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetRightMargin(0.02)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)


histos=[]
for f,title,ci in [ ('FlatRandomPtGunProducer_NeutrinoGun_v11_noNoise_aged_20191101.root','aged, no noise',ROOT.kRed+1),
                    ('FlatRandomPtGunProducer_NeutrinoGun_v11_startup_20191104.root','startup',ROOT.kGreen+1),
                    ('ttbar_ak8GenJetsNoNu_100_1_ttbar_v11_aged_biased_20191101.root','aged',ROOT.kAzure+1),
                ]:
    fIn=ROOT.TFile.Open(os.path.join(baseDir,f))
    histos.append(fIn.Get('ana/adchitsvspu').Clone(title.replace(',','_')) )
    histos[-1].SetTitle(title)
    histos[-1].SetDirectory(0)
    histos[-1].SetLineColor(ci)
    histos[-1].SetFillColor(ci)
    histos[-1].SetFillStyle(1001)

    profs=[]
    drawOpt='e1'
    c.SetLogy()
    leg=ROOT.TLegend(0.65,0.94,0.95,0.8)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    c.SetGridy()
    for t in ['cell','adc','toa','tdc']:
        profs.append( fIn.Get('ana/%scount'%t) )
        profs[-1].SetTitle(t)
        #profs[-1].GetYaxis().SetRangeUser(1,5e4)
        ci=ROOT.kGray+1
        if t=='adc' : ci=ROOT.kBlack
        if t=='toa' : ci=ROOT.kGreen+1
        if t=='tdc' : ci=ROOT.kAzure+1
        profs[-1].SetLineColor(ci)
        profs[-1].SetLineWidth(2)
        profs[-1].Draw(drawOpt)
        drawOpt='e1same'
        leg.AddEntry(profs[-1],profs[-1].GetTitle(),'ep')

    leg.SetHeader('#it{%s}'%title)
    leg.Draw()
    cmsHeader()
    c.Modified()
    c.Update()

    cname=title.replace(' ','')
    cname=cname.replace(',','')
    for ext in ['png','pdf']:
        c.SaveAs('%s_count.%s'%(cname,ext))
    fIn.Close()

c.SetLogy(False)
leg=ROOT.TLegend(0.15,0.94,0.4,0.8)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
for i in range(len(histos)):
    histos[i].Draw('box' if i==0 else 'boxsame')
    leg.AddEntry(histos[i],histos[i].GetTitle(),'f')
leg.Draw()

cmsHeader()
c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('adchitsvspu.'+ext)

