import ROOT
import os
import optparse
import re
import pickle

from prepareOccupancySummary import PLOTTITLES
from sensorEquivalentMap import *

def getRatio(num,den):

    """takes the ratio of two occupancy graphs"""

    gr=num.Clone('%s_ratio_%s'%(num.GetName(),den.GetName()))
    gr.SetMarkerStyle(num.GetMarkerStyle())
    gr.SetMarkerColor(num.GetMarkerColor())
    gr.SetFillStyle(num.GetFillStyle())
    gr.SetFillColor(num.GetFillColor())
    gr.SetLineStyle(num.GetLineStyle())
    gr.SetLineColor(num.GetLineColor())
    gr.SetTitle(num.GetTitle())
    gr.Set(0)

    x,y_num,y_den=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for ipt in range(num.GetN()):
        num.GetPoint(ipt,x,y_num)
        ey_num=num.GetErrorYhigh(ipt)
        den.GetPoint(ipt,x,y_den)
        ey_den=den.GetErrorYhigh(ipt)
        if float(y_den)<=0 : continue

        ratio=float(y_num)/float(y_den)
        ratio90=float(y_num+ey_num)/float(y_den+ey_den)
        gr.SetPoint(ipt,x,ratio)
        gr.SetPointError(ipt,0,0,0,ratio90-ratio)

    gr.Sort()
    return gr


def showWaferMomentSummary(momentSummary,sensorPos,proc,outdir):

    """show the wafer data in a layer-by-layer representation"""
    
    #loop over each sub-detector layer
    subdets=set([x[0] for x in sensorPos])
    for sd in subdets:
        layers=set( [x[1] for x in sensorPos if x[0]==sd] )
        for lay in layers:
                
            for d in momentSummary:
                
                layerKey=(sd,lay)
                
                uvzlist=[]
                labels=[]
                for waferKey in momentSummary[d]:
                    isd,ilay,iu,iv=waferKey
                    if isd!=sd or ilay!=lay :continue

                    ncells,r,z,eta,phi=sensorPos[waferKey]               
                    occ=[float(x)/float(ncells) for x in momentSummary[d][waferKey]]

                    uvzlist.append( [iu,iv,occ[0]] )
                    labels.append( r'$%d^{+%d}_{-%d}$'%(100*occ[1],100*(occ[2]-occ[1]),100*(occ[1]-occ[0]) ) )

                if len(uvzlist)==0: continue
                extraText=[ proc,
                            PLOTTITLES[d],
                            '%s layer %d'%('CEE' if sd==0 else 'CEH', lay)
                            ]
                drawSensorEquivalentMap(uvzlist=uvzlist,
                                        labels=labels,
                                        outname=os.path.join(outdir,'%s_sd%d_lay%d'%(d,sd,lay)),
                                        extraText=extraText,
                                        cmapName='Wistia',
                                        zran=[0,1],
                                        labelSize=14)
                


def dumpCSVOccupancy(momentSummary,sensorPos,outdir):

    """dumps a CSV file with the summary of occupancies"""

    import csv
    fOut=open(os.path.join(outdir,'occupancy_summary.dat'),'w')
    csv_writer = csv.writer(fOut, delimiter=',')
    
    for waferKey,counts in sorted(momentSummary.items(), key=lambda x: x[0]):

        isd,ilay,iu,iv=waferKey
        ncells,r,z,eta,phi=sensorPos[waferKey]               
        occ=[float(x)/float(ncells) for x in counts]
        csv_writer.writerow( [isd,ilay,1 if ncells>400 else 0,iu,iv]+occ )

    fOut.close()

def getRadialSummary(momentSummary,sensorPos,proc):

    """builds the radial summary for low/high density wafers"""
    
    radialProf={}

    #loop over each sub-detector layer
    subdets=set([x[0] for x in sensorPos])
    for sd in subdets:
        layers=set( [x[1] for x in sensorPos if x[0]==sd] )
        for lay in layers:
            
            layKey=(sd,lay)
            radialProf[layKey]=[ROOT.TGraphAsymmErrors(),ROOT.TGraphAsymmErrors()]
            radialProf[layKey][0].SetName('ld_%d_%d'%(sd,lay))
            radialProf[layKey][1].SetName('hd_%d_%d'%(sd,lay))
            radialProf[layKey][0].SetTitle(proc+' (LD)')
            radialProf[layKey][1].SetTitle(proc+' (HD)')


            for waferKey in momentSummary:
                isd,ilay,iu,iv=waferKey
                if isd!=sd or ilay!=lay :continue

                ncells,r,z,eta,phi=sensorPos[waferKey]
                occ=[float(x)/float(ncells) for x in momentSummary[waferKey]]
                densIdx=1 if ncells>400 else 0
                ipt=radialProf[layKey][densIdx].GetN()
                radialProf[layKey][densIdx].SetPoint(ipt,r,occ[1])
                radialProf[layKey][densIdx].SetPointError(ipt,0,0,0,occ[2]-occ[1])
            
            for gr in radialProf[layKey]:
                gr.Sort()

    return radialProf



def main():

    #parse inputs
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--out',   dest='output',  help='output directory [%default]',  default='occ_plots', type='string')    
    parser.add_option('-p', '--plots', dest='plots',  help='plots (csv) [%default]',       default='counts', type='string')    
    parser.add_option(      '--noWaferPlots',   dest='noWaferPlots',  help='disable wafer plots [%default]',             default=False, action='store_true')
    (opt, args) = parser.parse_args()

    #generic root cfg
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #define inputs
    procList=[ x.split(':') for x in args ]    

    #show only these
    plotList=opt.plots.split(',')

    sensorPos=None
    iarg=0
    radialProfiles=[]
    for proc,cache,idx in procList:        
        idx=int(idx)

        #prepare output
        iarg+=1
        plotout='%s/proc%d'%(opt.output,iarg)
        os.system('mkdir -p '+plotout)

        #read summary from pickle file
        with open(cache,'r') as cache:
            momentSummary=pickle.load(cache)
            pickle.load(cache)
            sensorPos=pickle.load(cache)

        #select only the dists/process of interest
        selMomentSummary={}
        for dist in plotList:
            selMomentSummary[dist]={}
            for wafer in momentSummary:
                selMomentSummary[dist][wafer]=momentSummary[wafer][idx][dist]

        if 'counts' in selMomentSummary:
            countsSummary=selMomentSummary['counts']
            dumpCSVOccupancy(countsSummary,sensorPos,plotout)
            radialProfiles.append( getRadialSummary(countsSummary,sensorPos,proc) )

        #show summary
        if opt.noWaferPlots: continue
        showWaferMomentSummary(selMomentSummary,sensorPos,proc,plotout)

    showRadialProfiles(radialProfiles,opt.output)

def showRadialProfiles(radialProfiles,outdir):

    nprocs=len(radialProfiles)
    colors=[ROOT.kBlack,     ROOT.kGray, ROOT.kMagenta,   ROOT.kMagenta+2,  ROOT.kMagenta-9, ROOT.kRed+1, ROOT.kAzure+7,   ROOT.kBlue-7]    
    markers=[24,20,25,21,26,22,32,23]
    fills=[3001,3325,3352,3305]
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)
    c.cd()
    p1=ROOT.TPad('p1','p1',0,0 if nprocs==1 else 0.4,1,1)
    p1.SetTopMargin(0.06)
    p1.SetLeftMargin(0.12)
    p1.SetRightMargin(0.02)
    p1.SetBottomMargin(0.1 if nprocs==1 else 0.02)
    p1.SetLogy()
    p1.Draw()

    p2=None
    if nprocs>1:
        c.cd()
        p2=ROOT.TPad('p2','p2',0,0,1,0.4)
        p2.SetTopMargin(0.01)
        p2.SetLeftMargin(0.12)
        p2.SetRightMargin(0.02)
        p2.SetBottomMargin(0.18)
        p2.Draw()
        
    c.cd()


    nprocs=len(radialProfiles)
    for layKey in radialProfiles[0]:

        p1.cd()
        p1.Clear()

        dyleg=0.06 if nprocs==1 else 0.1
        leg=ROOT.TLegend(0.6,0.9,0.95,0.9-dyleg*nprocs)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035 if nprocs==1 else 0.05)
        leg.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)

        mg=ROOT.TMultiGraph()
        ratios=[]
        for iproc in range(nprocs):       

            for igr in range(2):
                ci=colors[2*iproc+igr]
                ms=markers[2*iproc+igr]
                fs=fills[iproc]
                gr=radialProfiles[iproc][layKey][igr]
                gr.SetLineColor(ci)
                gr.SetFillColor(ci)
                gr.SetFillStyle(fs)
                gr.SetMarkerColor(ci)
                gr.SetMarkerStyle(ms)
                mg.Add(gr)
                leg.AddEntry(gr,gr.GetTitle(),'pf')

                if iproc>0:
                    ratios.append( getRatio(radialProfiles[iproc][layKey][igr],radialProfiles[0][layKey][igr]) )
                    
        mg.Draw('a3p')
        mg.GetYaxis().SetTitle('Occupancy')
        mg.GetXaxis().SetTitle('R [cm]')
        mg.GetYaxis().SetMoreLogLabels()
        mg.GetYaxis().SetRangeUser(0.5e-2,1)
        if nprocs>1:
            mg.GetYaxis().SetTitleSize(0.05)
            mg.GetYaxis().SetLabelSize(0.05)
            mg.GetXaxis().SetTitleSize(0)
            mg.GetXaxis().SetLabelSize(0)

        leg.Draw()
        
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.05 if nprocs==1 else 0.06)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')        
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        tex.SetTextSize(0.04 if nprocs==1 else 0.05)
        tex.DrawLatex(0.97,0.97,'%s layer %d'%('CEE' if layKey[0]==0 else 'CEH',layKey[1]))

        #add ratios if available
        if p2:
            p2.cd()
            p2.SetGridy()
            p2.Clear()
            mg_ratio=ROOT.TMultiGraph()
            for g in ratios:
                mg_ratio.Add(g)
            mg_ratio.Draw('a3p')
            mg_ratio.GetYaxis().SetTitle('Ratio')
            mg_ratio.GetXaxis().SetTitle('R [cm]')
            mg_ratio.GetYaxis().SetTitleOffset(0.75)
            mg_ratio.GetXaxis().SetTitleOffset(0.9)
            mg_ratio.GetYaxis().SetRangeUser(0.6,2.8)        
            mg_ratio.GetYaxis().SetTitleSize(0.08)
            mg_ratio.GetYaxis().SetNdivisions(5)
            mg_ratio.GetYaxis().SetLabelSize(0.08)
            mg_ratio.GetXaxis().SetTitleSize(0.08)
            mg_ratio.GetXaxis().SetLabelSize(0.08)

        #synchronize x-axis ratio
        mg.GetXaxis().SetRangeUser( mg_ratio.GetXaxis().GetXmin(), mg_ratio.GetXaxis().GetXmax() )

        
        c.cd()
        c.Modified()
        c.Update()
        c.SaveAs(os.path.join(outdir,'radialprofile_%d_%d.png'%layKey))

    p1.Delete()
    if p2: p2.Delete()

if __name__ == "__main__":
    main()
