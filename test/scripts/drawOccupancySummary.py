import ROOT
import os
import optparse
import re
import pickle

from prepareOccupancySummary import PLOTTITLES
from sensorEquivalentMap import *
maxPhiSector={'CEH':2*ROOT.TMath.Pi()/3.,
              'CEE':ROOT.TMath.Pi()/3.}

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

    x,y_num,x_den,y_den=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for ipt in range(num.GetN()):
        num.GetPoint(ipt,x,y_num)
        y_numval=float(y_num)
        ey_num=num.GetErrorYhigh(ipt)
        y_num90=y_num+ey_num
        ey_num=num.GetErrorYhigh(ipt)

        #find closest
        closestj=-1
        closestdx=1e21
        for jpt in range(den.GetN()):
            den.GetPoint(jpt,x_den,y_den)         
            dx=abs(float(x_den)-float(x))
            if dx>closestdx: continue
            closestj=jpt
            closestdx=dx
        if closestj<0 : continue
        if closestdx>0.1 : continue

        jpt=closestj
        den.GetPoint(jpt,x_den,y_den)
        y_denval=float(y_den)
        ey_den=den.GetErrorYhigh(jpt)
        y_den90=y_denval+ey_den

        if y_denval==0:
            if y_numval==0:
                ratio=1
            else:
                continue
        else:
            ratio=y_numval/y_denval

        if y_den90==0:
            if y_num90==0:
                ratioUnc=0
            else:
                ratioUnc=1
        else:
            ratioUnc=abs(y_num90/y_den90-ratio)

        npt=gr.GetN()
        gr.SetPoint(npt,x,ratio)
        gr.SetPointError(npt,0,0,0,ratioUnc)

    gr.Sort()
    return gr


def showWaferMomentSummary(momentSummary,sensorPos,proc,outdir,tag,showUVLabeled=False):

    """show the wafer data in a layer-by-layer representation"""

    xbias=0.5 #the binning introduces a bias in the quantile
    
    #loop over each sub-detector layer
    subdets=set([x[0] for x in sensorPos])
    for sd in subdets:
        layers=set( [x[1] for x in sensorPos if x[0]==sd] )
        for lay in layers:
            layerKey=(sd,lay)
                
            uvzlist,sectoreq_uvzlist=[],[]
            labels,sectoreq_labels,sectoreq_uvlabels=[],[],[]
            for waferKey in momentSummary:
                isd,ilay,iu,iv=waferKey
                if isd!=sd or ilay!=lay :continue

                ncells,waferType,r,z,eta,phi,xpos,ypos=sensorPos[waferKey]               
                occ=[float(max(x-xbias,0))/float(ncells) for x in momentSummary[waferKey]]

                #uvzlist.append( [iu,iv,occ[0]] )
                if (phi>=0 or (iu<0 and iv>=0)) and phi<=maxPhiSector[sd]:
                    sectoreq_uvzlist.append( [xpos,ypos,occ[0]] )
                    sectoreq_labels.append( r'$%d^{+%d}_{-%d}$'%(round(100*occ[1]),
                                                                 round(100*(occ[2]-occ[1])),
                                                                 round(100*(occ[1]-occ[0])) ) )
                    sectoreq_uvlabels.append( r'(%d,%d)'%(iu,iv) )

                uvzlist.append( [xpos,ypos,occ[0]] )
                labels.append( '%d'%round(100*occ[1]) )

            if len(uvzlist)==0: continue
            extraText=[ proc, PLOTTITLES[tag], '%s layer %d'%(sd, lay)]
            drawSensorEquivalentMap(uvzlist=uvzlist,
                                    labels=labels,
                                    outname=os.path.join(outdir,'%s_%s_lay%d'%(tag,sd,lay)),
                                    extraText=extraText,
                                    cmapName='Wistia',
                                    zran=[0,1],
                                    labelSize=14)
            drawSensorEquivalentMap(uvzlist=sectoreq_uvzlist,
                                    labels=sectoreq_labels,
                                    outname=os.path.join(outdir,'%s_%s_lay%d_sectoreq'%(tag,sd,lay)),
                                    extraText=extraText,
                                    cmapName='Wistia',
                                    zran=[0,1],
                                    labelSize=14)
            if not showUVLabeled: continue
            drawSensorEquivalentMap(uvzlist=sectoreq_uvzlist,
                                    labels=sectoreq_uvlabels,
                                    outname=os.path.join(outdir,'%s_%s_lay%d_sectorequv'%(tag,sd,lay)),
                                    extraText=extraText,
                                    cmapName='Wistia',
                                    zran=[0,1],
                                    labelSize=14)
                


def dumpCSVOccupancy(momentSummary,sensorPos,outdir):

    """dumps a CSV file with the summary of occupancies"""

    xbias=0.5

    import csv
    fOut=open(os.path.join(outdir,'occupancy_summary.dat'),'w')
    csv_writer = csv.writer(fOut, delimiter=',')
    
    for waferKey,counts in sorted(momentSummary.items(), key=lambda x: x[0]):

        isd,ilay,iu,iv=waferKey
        ncells,waferType,r,z,eta,phi,xpos,ypos=sensorPos[waferKey] 
        occ=[float(max(x-xbias,0))/float(ncells) for x in counts]
        csv_writer.writerow( [isd,ilay,1 if waferType==0 else 0,ncells,xpos,ypos,iu,iv]+occ )

    fOut.close()

def getRadialSummary(momentSummary,sensorPos,proc):

    """builds the radial summary for low/high density wafers"""
    
    radialProf={}

    xbias=0.5

    #loop over each sub-detector layer
    subdets=set([x[0] for x in sensorPos])
    for sd in subdets:
        layers=set( [x[1] for x in sensorPos if x[0]==sd] )
        for lay in layers:
            
            layKey=(sd,lay)
            radialProf[layKey]=[ROOT.TGraphAsymmErrors(),ROOT.TGraphAsymmErrors()]
            radialProf[layKey][0].SetName('ld_%s_%d'%(sd,lay))
            radialProf[layKey][1].SetName('hd_%s_%d'%(sd,lay))
            radialProf[layKey][0].SetTitle(proc+' (LD)')
            radialProf[layKey][1].SetTitle(proc+' (HD)')


            for waferKey in momentSummary:
                isd,ilay,iu,iv=waferKey
                if isd!=sd or ilay!=lay :continue

                ncells,waferType,r,z,eta,phi,xpos,ypos=sensorPos[waferKey]
                if phi<0 or phi>maxPhiSector[isd]:
                    continue

                occ=[float(max(x-xbias,0))/float(ncells) for x in momentSummary[waferKey]]
                densIdx=1 if waferType==0 else 0
                ipt=radialProf[layKey][densIdx].GetN()
                radialProf[layKey][densIdx].SetPoint(ipt,r,occ[1])
                radialProf[layKey][densIdx].SetPointError(ipt,0,0,0,occ[2]-occ[1])
            
            for gr in radialProf[layKey]:
                gr.Sort()

    return radialProf


def getWaferUVSummary(momentSummary,sensorPos,proc,waferLongProfiles):

    """builds the longitudinal profile summary for a particular wafer"""
    
    waferUVSummary={}

    xbias=0.5

    #loop over each sub-detector layer
    subdets=set([x[0] for x in sensorPos])
    for sd in subdets:
        for waferKey in momentSummary:
            isd,ilay,iu,iv=waferKey

            if isd != sd: continue
            if not [iu,iv] in waferLongProfiles: continue
            
            sumKey=(isd,iu,iv)
            if not sumKey in waferUVSummary:
                waferUVSummary[sumKey]=ROOT.TGraphAsymmErrors()
                name='%s_%d_%d'%sumKey
                name=name.replace('-','m')
                waferUVSummary[sumKey].SetName(name)
                waferUVSummary[sumKey].SetTitle(proc)

            ncells,waferType,r,z,eta,phi,xpos,ypos=sensorPos[waferKey]
            occ=[float(max(x-xbias,0))/float(ncells) for x in momentSummary[waferKey]]
            ipt=waferUVSummary[sumKey].GetN()
            waferUVSummary[sumKey].SetPoint(ipt,z,occ[1])
            waferUVSummary[sumKey].SetPointError(ipt,0,0,0,occ[2]-occ[1])
            
           
    for sumKey in waferUVSummary:
        waferUVSummary[sumKey].Sort()
        

    return waferUVSummary

def showSummaryProfiles(profColl,outdir,tag,profType='radial',xtitle='R [cm]', ytitle='Occupancy', yran=(0.5e-2,1),yratioran=(0.92,1.08)):

    """ final grand summary plot """

    nprocs=len(profColl)
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
    p1.SetTopMargin(0.07)
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


    nprocs=len(profColl)
    for collKey in profColl[0]:

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
        x,y=ROOT.Double(0),ROOT.Double(0)
        xmin,xmax=1e21,-1e21
        for iproc in range(nprocs):       

            if type(profColl[iproc][collKey])==list:
                for igr in range(len(profColl[iproc][collKey])):
                    ci=colors[2*iproc+igr]
                    ms=markers[2*iproc+igr]
                    fs=fills[iproc]
                    gr=profColl[iproc][collKey][igr]

                    gr.GetPoint(0,x,y)
                    xmin=min(float(x),xmin)
                    gr.GetPoint(gr.GetN()-1,x,y)
                    xmax=max(float(x),xmax)

                    gr.SetLineColor(ci)
                    gr.SetFillColor(ci)
                    gr.SetFillStyle(fs)
                    gr.SetMarkerColor(ci)
                    gr.SetMarkerStyle(ms)
                    mg.Add(gr)
                    leg.AddEntry(gr,gr.GetTitle(),'pf')

                    if iproc>0:
                        ratios.append( getRatio(gr,profColl[0][collKey][igr]) )

            else:
                ci=colors[2*iproc]
                ms=markers[2*iproc]
                fs=fills[iproc]
                gr=profColl[iproc][collKey]
                gr.SetLineColor(ci)
                gr.SetFillColor(ci)
                gr.SetFillStyle(fs)
                gr.SetMarkerColor(ci)
                gr.SetMarkerStyle(ms)
                mg.Add(gr)
                leg.AddEntry(gr,gr.GetTitle(),'pf')

                if iproc>0:
                    ratios.append( getRatio(gr,profColl[0][collKey]) )

        frame=ROOT.TH1F('frame',';%s;Fraction of channels;'%xtitle,1,xmin,xmax)
        frame.SetDirectory(0)
        rframe=frame.Clone('rframe')
        rframe.SetDirectory(0)
        rframe.GetYaxis().SetTitle('Ratio to ref.')

        #draw            
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(0.8)
        frame.GetYaxis().SetRangeUser(yran[0],yran[1])
        frame.GetYaxis().SetTitle()
        mg.GetYaxis().SetMoreLogLabels()
        if nprocs>1:
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetYaxis().SetLabelSize(0.05)
            frame.GetXaxis().SetTitleSize(0)
            frame.GetXaxis().SetLabelSize(0)

        mg.Draw('3p')
        leg.Draw()
        
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.05 if nprocs==1 else 0.06)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.955,'#bf{CMS} #it{simulation preliminary}')        
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        tex.SetTextSize(0.04 if nprocs==1 else 0.05)
        if len(collKey)==2:
            tex.DrawLatex(0.97,0.965,'%s layer %d'%(collKey[0],collKey[1]))
        else:
            tex.DrawLatex(0.97,0.965,'%s (u,v)=(%d,%d)'%(collKey[0],collKey[1],collKey[2]))

        #add ratios if available
        if p2:
            p2.cd()
            p2.SetGridy()
            p2.Clear()
            rframe.GetYaxis().GetCenterTitle()
            rframe.GetYaxis().SetTitleOffset(0.75)
            rframe.GetXaxis().SetTitleOffset(0.9)
            rframe.GetYaxis().SetRangeUser(yratioran[0],yratioran[1])
            rframe.GetYaxis().SetTitleSize(0.08)
            rframe.GetYaxis().SetNdivisions(5)
            rframe.GetYaxis().SetLabelSize(0.08)
            rframe.GetXaxis().SetTitleSize(0.08)
            rframe.GetXaxis().SetLabelSize(0.08)
            rframe.Draw()

            mg_ratio=ROOT.TMultiGraph()
            for g in ratios:
                mg_ratio.Add(g)
            mg_ratio.Draw('p')
        
        c.cd()
        c.Modified()
        c.Update()        
        outNameFromKey='_'.join([str(x) for x in collKey])
        outNameFromKey=outNameFromKey.replace('-','m')
        c.SaveAs(os.path.join(outdir,'%sprofile_%s_%s.png'%(profType,tag,outNameFromKey)))
        rframe.Delete()

    try:
        frame.Delete()
        p1.Delete()
        if p2:
            rframe.Delete()
            p2.Delete()
    except Exception as e:
        print e


def main():

    #parse inputs
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--out',   dest='output',  help='output directory [%default]',  default='occ_plots', type='string')    
    parser.add_option('-p', '--plots', dest='plots',  help='plots (csv) [%default]',        default='counts', type='string')    
    parser.add_option('--yran',        dest='yran',  help='y range [%default]',             default='5e-3,1', type='string')    
    parser.add_option('--yratioran',   dest='yratioran',  help='y ratio range [%default]',  default='0.92,1.08', type='string')    
    parser.add_option(      '--waferLongProfiles',   dest='waferLongProfiles',    help='show these wafer longitudinal profiles [%default]',  default='0,3:0,5:0,10', type='string')
    parser.add_option(      '--noWaferPlots',        dest='noWaferPlots',         help='disable wafer plots [%default]',                     default=False,          action='store_true')
    (opt, args) = parser.parse_args()

    yran=[float(x) for x in opt.yran.split(',')]
    yratioran=[float(x) for x in opt.yratioran.split(',')]

    #decode wafer plot list string                                                                                                                                                                                                                                                                
    if opt.waferLongProfiles:
        opt.waferLongProfiles=[[int(y) for y in x.split(',')] for x in opt.waferLongProfiles.split(':')]
        print opt.waferLongProfiles

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
    radialProfiles={}
    longProfiles={}
    for proc,cache,idx in procList:        
        idx=int(idx)

        #prepare output
        iarg+=1
        plotout='%s/proc%d'%(opt.output,iarg)
        os.system('rm -rf {0} && mkdir -p {0}'.format(plotout))

        #read summary from pickle file
        with open(cache,'r') as cache:
            momentSummary=pickle.load(cache)
            pickle.load(cache)
            sensorPos=pickle.load(cache)

        #select only the dists/process of interest
        for dist in plotList:
            selMomentSummary={}
            for wafer in momentSummary:
                selMomentSummary[wafer]=momentSummary[wafer][idx][dist]
            if not 'counts' in dist: 
                continue
            if not dist in radialProfiles:
                radialProfiles[dist]=[]
                longProfiles[dist]=[]
            radialProfiles[dist].append( getRadialSummary(selMomentSummary,sensorPos,proc) )
            longProfiles[dist].append( getWaferUVSummary(selMomentSummary,sensorPos,proc,opt.waferLongProfiles) )

            if dist=='counts':         
                dumpCSVOccupancy(selMomentSummary,sensorPos,plotout)
                if not opt.noWaferPlots:
                    showWaferMomentSummary(selMomentSummary,sensorPos,proc,plotout,tag=dist,showUVLabeled=True if idx==1 else False)

    for profColl,profType,xtitle,ytitle in [(radialProfiles, 'radial', 'R [cm]', 'Occupancy'), 
                                            (longProfiles,   'long',   'z[cm]',  'Occupancy')]:
        for dist in profColl:
            if len(profColl[dist])==0 : continue
            showSummaryProfiles(profColl[dist],opt.output,tag=dist,profType=profType,xtitle=xtitle,ytitle=ytitle,yran=yran,yratioran=yratioran)



if __name__ == "__main__":
    main()
