import ROOT

from drawRadiationMapPlots import drawHeader

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

qe2fc=1.60217646E-4;
ileakParams = [0.993,-42.668]
noise_params={80 :[636., 15.6, 0.0328],
              160:[1045., 8.74, 0.0685],
              320:[1915., 2.79, 0.0878]}
cces=[]
signal=[]
ileak=[]
noise={80:[],160:[],320:[]}

SENSORDATA=[('dd-FZ 120#mum', [1.5e+15, -3.00394e-17, 0.318083],  120.*67.,  50, 0.52*(120.e-4), 1, ROOT.kBlack), 
            ('dd-FZ 200#mum', [1.5e+15, -3.09878e-16, 0.211207],  200.*70.,  65, 1.18*(200.e-4), 1, ROOT.kAzure-2),
            ('dd-FZ 300#mum', [1.5e+15, -7.96539e-16, 0.251751],  300.*73.,  45, 1.18*(300.e-4), 1, ROOT.kRed+1),
            ('epi 100#mum',   [3.5e+15, -9.73872e-19, 0.263812],  100.*73.,  50, 0.52*(100.e-4), 2, ROOT.kBlack),
            #('epi 200#mum',   [1.5e+15, -3.09878e-16, 0.211207],  200.*70.,  65, 1.18*(300.e-4), 2, ROOT.kAzure-2),
            #('epi 300#mum',   [6e+14,   -7.96539e-16, 0.251751],  300.*73.,  45, 1.18*(300.e-4), 2, ROOT.kRed+1),
    ]
noiseColors=[ROOT.kBlack,ROOT.kAzure-2,ROOT.kRed+1,ROOT.kBlack]
noiseLines=[1,1,1,2]

for title,params,mipeq,cap,vol,ls,ci in SENSORDATA:

    i=len(cces)
    cces.append( ROOT.TF1('cce%d'%i,'[3]*(x <[0] ? 1. + [1]*x : (1-[2]*log(x))+ ([1]*[0]+[2]*log([0])))',1e13,1e17) )
    cces[-1].SetTitle(title)
    for ip in range(len(params)):
        cces[-1].SetParameter(ip,params[ip])    
    cces[-1].SetParameter(3,1.)
    cces[-1].SetLineStyle(ls)
    cces[-1].SetLineWidth(2)
    cces[-1].SetLineColor(ci)
    
    ileak.append( ROOT.TF1('ileak%d'%i,'exp([0] *log(x) + [1]) * [2]',1e13,1e17) )
    for ip in range(len(ileakParams)):
        ileak[-1].SetParameter(ip,ileakParams[ip])                                        
    ileak[-1].SetParameter(2,vol*1e6)
    ileak[-1].SetTitle(title)
    ileak[-1].SetLineStyle(ls)
    ileak[-1].SetLineWidth(2)
    ileak[-1].SetLineColor(ci)

    icce=len(cces)
    signal.append( cces[-1].Clone('sig%d'%icce) )
    signal[-1].SetParameter(3,mipeq*qe2fc)
    
    for key in noise_params: 
        noise[key].append( ROOT.TF1('noiseq%d%d'%(key,icce),
                                    '%f*sqrt( pow(840.*sqrt(exp([0]*log(x)+[1])*[2]),2)+1.25*pow([3]+[4]*[6]+[5]*pow([6],2),2) )'%qe2fc,
                                    1e13,1e17) )
        for ip in range(len(ileakParams)):
            noise[key][-1].SetParameter(ip,ileakParams[ip])
        noise[key][-1].SetParameter(len(ileakParams),vol*1e6)
        for ip in range(len(noise_params[key])):
            noise[key][-1].SetParameter(ip+1+len(ileakParams),noise_params[key][ip])
        noise[key][-1].SetParameter(1+len(ileakParams)+len(noise_params[key]),cap)
        noise[key][-1].SetTitle(title)
        noise[key][-1].SetLineWidth(2)
        noise[key][-1].SetLineColor(noiseColors[len(noise[key])-1])
        noise[key][-1].SetLineStyle(noiseLines[len(noise[key])-1])



c=ROOT.TCanvas('c','c',600,600)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
c.SetGridy()
c.SetGridx()
c.SetLogx()

frame=ROOT.TH1F('frame',';Fluence [n_{eq}/cm^{2}];',1,1e13,1e17)
frame.GetXaxis().SetMoreLogLabels()
frame.GetYaxis().SetTitleOffset(1.1)
frame.GetXaxis().SetTitleOffset(1.3)

for gr,title,name,yran in [(cces,'CCE','CCE',(0,1)),
                           (ileak,'I_{leak} [#muA]','Ileak',(1e-1,100)),
                           (signal,'Signal [fC]','Signal',(0,5)),                           
                           (noise[80],'Noise (80fC) [fC]','Noise80',(0,5)),
                           (noise[160],'Noise (160fC) [fC]','Noise160',(0,5)),
                           (noise[320],'Noise (320fC) [fC]','Noise320',(0,5))
                       ]:
    
    c.SetLogy(False)
    if name=='Ileak' or 'Noise' in name or 'Signal' in name:
        leg=ROOT.TLegend(0.15,0.67,0.4,0.95)
        if name=='Ileak': 
            c.SetLogy()
    else:
        leg=ROOT.TLegend(0.15,0.12,0.4,0.4)

    frame.GetYaxis().SetTitle(title)
    frame.GetYaxis().SetRangeUser(yran[0],yran[1])
    frame.Draw()


    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    for i in range(len(gr)):
        gr[i].Draw('same')
        leg.AddEntry(gr[i],gr[i].GetTitle(),'l')
    leg.Draw()
    drawHeader()

    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(name,ext))

