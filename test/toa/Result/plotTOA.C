#ifndef plotTOA_C
#define plotTOA_C

// Custom headers
#include "../Macros/tdrstyle.C"
#include "../Macros/CMS_lumi.C"
#include "../Macros/TreeReader.h"
// ROOT headers
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TF1.h"
// C headers
#include <numeric>

// Global parameters
constexpr double v68 = 0.68268949;

// Functions
float divide(const float&, const float&);
void processEvent(std::map<float, std::map<std::string, TH2D>>&, std::map<float, TEfficiency>&, TreeReader&);
void drawHitToA(const std::vector<float>&, const std::vector<float>&, const std::vector<float>&, const std::string&);
void drawProj(std::map<std::string, TH2D>&, std::vector<float>&, const bool&, const float&, const std::string&);
void drawMeanVal(const std::map<float, std::vector<float>>&, const int&, const std::string&);
void drawProjVal(const std::map<float, std::vector<float>>&, const int&, const std::string&);
void drawEff(TEfficiency&, std::array<float,8>&, const float&, const std::string&);
void parmEffVal(const std::map<float, std::array<float,8>>&, const int&, const std::string&);
std::vector<float> fixSizeHighestDensity(const std::vector<float>&, std::vector<float>, const unsigned int& minNhits=3, const float& deltaT=0.210, const float& timeWidthBy=0.5);

bool isPion(false);
bool isPhoton(false);
bool isAged_2p8(false);
bool doEne(true);
const double maxRng = doEne ? 120 : 50;


const std::array<std::string, 4> getCut(const std::string& type)
{
  std::array<std::string, 4> res;
  // Create eta label
  auto d = type.substr(type.find("eta")+3);
  std::replace(d.begin(), d.end(), 'p', '.');
  res[0] = Form("#eta = %g", std::stof(d.substr(0, d.find("_"))));
  // Create hit cut label
  res[1] = "N^{rec}_{hits} #geq 3";
  // Create particle label
  if (type.find("Photon")!=std::string::npos) {
    res[2] = "#gamma";
    res[3] = "#rho < 3 cm";
  }
  else if (type.find("K0L")!=std::string::npos) {
    res[2] = "K^{0}_{L}";
    //res[3] = "#rho < [3, 8, 16] cm";
    res[3] = "#rho < 3 cm";
  }
  else if (type.find("Pion")!=std::string::npos) {
    res[2] = "#pi";
    res[3] = "";
  }
  if (type.find("InRad2p0")!=std::string::npos)
    res[3] = "#rho < 2 cm";
  else if (type.find("NoRad")!=std::string::npos)
    res[3] = "";
  return res;
}


void plotTOA(const std::string& tp="SinglePionGun_eta1p8_13_1_0_pre4_D99_startup", const size_t& opt=1)//_3iab
{
  // Extract ntuple
  const auto inputFile =  "Output/TREE_" + std::string(opt==2 ? "InRad2p0_" : (opt==1 ? "InRadius_" : "")) + tp + ".root";
  TreeReader reader(inputFile, "toa", false);

  if (tp.find("Pion")!=std::string::npos) isPion = true;
  if (tp.find("p22")!=std::string::npos || tp.find("hoton")!=std::string::npos) isPhoton = true;
  if (tp.find("iab")!=std::string::npos && tp.find("eta2p8")!=std::string::npos) isAged_2p8 = true;
  const auto type = (opt==2 ? "InRad2p0/" : (opt==1 ? "InRadius/" : "NoRad/")) + tp;

  // Define binning
  std::vector<float> minSimChgV({12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 60, 65, 70, 75, 80, 85, 90, 95, 100});
  //std::vector<float> minSimChgV({12, 22, 32, 42, 53, 65, 75, 85, 95, 100});
  //std::vector<float> minSimChgV({12});

  // Define histogram
  std::map<float, std::map<std::string, TH2D>> h2DM;
  std::map<float, TEfficiency> effM;
  const std::vector<std::string> varList = {"toaAtPV_corr"};//, "toaAt0_corr"};//, "toaDelta_avgQ"};

  // Define data container
  std::vector<std::map<std::string, float> > hitInfo;
  for (const auto& v : minSimChgV) {
    const auto lbl = doEne ? ";Energy [GeV];Efficiency" : ";p_{T} [GeV];Efficiency";
    effM[v] = TEfficiency(Form("eff_%.0f", v*10), lbl, floor(maxRng/1.0), 0, maxRng);
    for (const auto& s : varList) {
      const auto lbl = doEne ? ";Energy [GeV];time [ps]" : ";p_{T} [GeV];time [ps]";
	  h2DM[v][s] = TH2D(Form("h2D_%.0f_%s", v*10, s.c_str()), lbl,  floor(maxRng/1.0), 0, maxRng, 25E4, -5E3, 20E3);
	  h2DM[v][s].Sumw2();
    }
  }
  
  // Extract events from tree
  const auto& nentries = reader.getEntries();
  std::cout << "[INFO] Processing events from " << nentries << " entries" << std::endl;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    // Get ntuple entry
    reader.setEntry(jentry);
    if(!(jentry % (nentries/10))) std::cout << "Processed " << jentry << " entrries out of " << nentries << std::endl;
    processEvent(h2DM, effM, reader);
  }

  // Create output directories
  for (auto c : {"Event", "Histogram", "Mean", "Res", "Summary", "Efficiency", "Parametrize"})
    for (auto d : {"png", "pdf", "C"})
      gSystem->mkdir(("Plot/"+type+"/"+c+"/"+d+"/").c_str(), kTRUE);

  // Draw histograms
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetFillColor(kWhite);

  std::map<float, std::vector<float>> hValV, mValV;
  for (auto& h : h2DM) {
    //drawProj(h.second, mValV[h.first], true, h.first, type);
    drawProj(h.second, hValV[h.first], false, h.first, type);
  }
  //drawMeanVal(mValV, 0, type);
  for (int i=0; i<4; i++)
    drawProjVal(hValV, i, type);
 
  std::map<float, std::array<float,8>> effValV;
  for (auto& e : effM)
    drawEff(e.second, effValV[e.first], e.first, type);
  for (int i=0; i<3; i++)
    parmEffVal(effValV, i, type);
};


void processEvent(std::map<float, std::map<std::string, TH2D>>& h2DM, std::map<float, TEfficiency>& effM, TreeReader& reader)
{
  // Extract variables  
  std::map<std::string, float> genVarM;
  const auto varName = doEne ? "genergy" : "gpt";
  for (const auto& v : {varName})
    genVarM[v] = reader.getFloat(v);
  std::map<std::string, std::vector<float> > varM;
  for (const auto& v : {"qsim", "qrec", "toarecAtPV"})
    varM[v] = reader.getVFloat(v);
  
  // Require at least 3 hits
  const auto& nHits = varM.at("toarecAtPV").size();
  if (nHits<3) return;

  // Calculate information
  std::map<float, std::vector<float>> toarecV, qV;
  for (size_t i=0; i<nHits; i++) {
    // Process hits with toa
    const auto toarec = varM.at("toarecAtPV")[i];
    if (toarec<=-50) continue;

    const auto qrec = varM.at("qrec")[i];
    const auto qsim = varM.at("qsim")[i];
    
    for (const auto& h : h2DM) {
      const auto& v = h.first;
      if (qrec<=v) break;
      toarecV[v].push_back(toarec);
      qV[v].push_back(qrec);
    }
  }

  for (auto& hh : h2DM) {
    const auto& v = hh.first;
    
    // Require at least 3 selected hits
    const bool hasRecoToA = toarecV[v].size()>=3;

    // Fill efficiency
    effM[v].Fill(hasRecoToA, genVarM.at(varName));
    
    // Fill containers
    if (!hasRecoToA) continue;
    
    std::map<std::string, float> res, n;
    const auto toa_corr = fixSizeHighestDensity(toarecV[v], qV[v]);
    res["toaAtPV_corr"] = toa_corr[0];

    for (auto& h : hh.second)
      if (res.at(h.first)>-50)
	    h.second.Fill(genVarM.at(varName), res[h.first]*1000);
  }
};


float divide(const float& num, const float& den)
{
  if (den==0) return -99;
  return num/den;
};


template<class T>
void formatHist(T& h)
{
  h.SetMarkerSize(2);
  h.SetFillStyle(0);
  h.SetMarkerStyle(20);
  h.SetMarkerColor(kBlack);
  h.SetLineColor(kBlack);
  h.SetTitle("");
  h.GetXaxis()->CenterTitle(true);
  h.GetXaxis()->SetLabelSize(0.04);
  h.GetXaxis()->SetTitleOffset(0.81);
  h.GetXaxis()->SetTitleSize(0.063);
  h.GetXaxis()->SetTitleFont(42);
  h.GetYaxis()->CenterTitle(true);
  h.GetYaxis()->SetLabelSize(0.040);
  h.GetYaxis()->SetTitleSize(0.063);
  h.GetYaxis()->SetTitleOffset(0.85);
  h.GetYaxis()->SetTitleFont(42);
  h.GetYaxis()->SetNdivisions(509);
};


void drawHitToA(const std::vector<float>& toaV, const std::vector<float>& wV, const std::vector<float>& parV, const std::string& type)
{ 
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetFillColor(kWhite);
  
  // Fill histogram
  TH1D h("hHitToA", "", 175E2, -10, 25);
  h.Sumw2();
  for (size_t i=0; i<toaV.size(); i++)
    h.Fill(toaV[i], wV[i]);
  const auto& mean = h.GetMean();
  
  // Create canvas
  const std::string hName = h.GetName();
  TCanvas canvas(("c"+hName).c_str(), "", 1000, 800);
  canvas.cd();

  formatHist(h);
    
  h.GetXaxis()->SetTitle("ToA [ns]");
  h.GetYaxis()->SetTitle("q-Weighted Hits");
  h.GetXaxis()->SetRangeUser(-0.4, 1.0);
  
  // Draw projection
  h.Draw("P");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  TLine lin1(parV[0], 0, parV[0], 1000);
  lin1.SetLineStyle(1); lin1.SetLineWidth(3); lin1.SetLineColor(kBlue); lin1.Draw("SAME");
  TLine lin2(mean, 0, mean, 1000);
  lin2.SetLineStyle(1); lin2.SetLineWidth(3);  lin2.SetLineColor(kRed); lin2.Draw("SAME");
  TLine lin3(parV[1], 0, parV[1], 1600);
  lin3.SetLineStyle(7); lin3.SetLineWidth(3);  lin3.SetLineColor(kBlack); lin3.Draw("SAME");
  TLine lin4(parV[2], 0, parV[2], 1600);
  lin4.SetLineStyle(7); lin4.SetLineWidth(3);  lin4.SetLineColor(kBlack); lin4.Draw("SAME");
  
  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.85, Form("%s, %s, %s", lbl[2].c_str(), lbl[0].c_str(), lbl[3].c_str()));
  
  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Event/"+c+"/"+hName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing result in: " << ("Plot/"+type+"/Event/root/"+hName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Event/root/"+hName+".root").c_str(), "RECREATE");
  file.cd();
  h.Write("data");
  file.Close();
};


void drawHist(TH1D hh, const float& enMin, const float& enMax, const float& v, const std::string& type)
{
  const double qThr[3] = {(1-v68)/2., 0.5, (1+v68)/2.};
  double qVal[3] = {-1, -1, -1};
  hh.GetQuantiles(3, qVal, qThr);
  const auto rms = hh.GetRMS();
  
  auto& h = *dynamic_cast<TH1D*>(hh.Rebin(80));
  
  // Create canvas
  const std::string hName = h.GetName();
  TCanvas canvas(("c"+hName).c_str(), "", 1000, 800);
  canvas.cd();

  formatHist(h);
    
  h.GetXaxis()->SetTitle("ToA [ps]");
  h.GetYaxis()->SetTitle("Events");
  h.GetXaxis()->SetRangeUser(-50, 200);
  
  // Draw projection
  h.Draw("P");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  TLine lin1(qVal[1], 0, qVal[1], maxRng);
  lin1.SetLineStyle(1); lin1.SetLineWidth(3);  lin1.SetLineColor(kBlue); lin1.Draw("SAME");
  TLine lin2(qVal[0], 0, qVal[0], 100);
  lin2.SetLineStyle(7); lin2.SetLineWidth(3);  lin2.SetLineColor(kBlack); lin2.Draw("SAME");
  TLine lin3(qVal[2], 0, qVal[2], 100);
  lin3.SetLineStyle(7); lin3.SetLineWidth(3);  lin3.SetLineColor(kBlack); lin3.Draw("SAME");
  TLine lin4(qVal[1]-rms, 0, qVal[1]-rms, 150);
  TLine lin5(qVal[1]+rms, 0, qVal[1]+rms, 150);
  
  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.55, 0.85, Form("%s, %s, %s", lbl[2].c_str(), lbl[0].c_str(), lbl[1].c_str()));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.55, 0.79, lbl[3].c_str());
  const auto vName = doEne ? "E" : "p_{T}";
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.55, 0.73, Form("%g #leq %s < %g GeV", enMin, vName, enMax));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.55, 0.67, Form("q^{rec} > %g fC", v));
  
  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Histogram/"+c+"/"+hName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing result in: " << ("Plot/"+type+"/Histogram/root/"+hName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Histogram/root/"+hName+".root").c_str(), "RECREATE");
  file.cd();
  h.Write("data");
  file.Close();
};


void drawProj(std::map<std::string, TH2D>& h2DM, std::vector<float>& fVal, const bool& doMean, const float& v, const std::string& type)
{
  // Create canvas
  const std::string hName = h2DM.begin()->second.GetName();
  TCanvas canvas(("c"+hName).c_str(), "", 1000, 800);
  canvas.cd();

  std::map<std::string, std::unique_ptr<TF1>> f;
  std::map<std::string, double> minXM;

  // Create projection
  std::map<std::string, TGraphAsymmErrors> grM;
  std::vector<int> colorV({kBlack, kRed, kBlue, kOrange, kGreen-2, kViolet, kAzure-7, kGreen+2});
  std::vector<int> markerV({22, 23, 24, 22, 23, 24});
  size_t i(0);
  for (auto& h : h2DM) {
    if (h.first.rfind("toa", 0)!=0) continue;
    const auto& h2D = h.second;
    auto& gr = grM[h.first];
    gr = TGraphAsymmErrors(h2D.GetNbinsX());
    const auto maxN = gr.GetN();
    for (int i=1; i<=maxN; i++) {
      const auto xVal = h2D.GetXaxis()->GetBinCenter(i);
      const auto xWdt = h2D.GetXaxis()->GetBinWidth(i)/2.;
      const auto hp = std::unique_ptr<TH1D>(h2D.ProjectionY(Form("%s_E_%.0f", h.second.GetName(), xVal*100), i, i, "eo"));
      if (!hp) throw std::logic_error("[ERROR] Invalid histogram pointer!");
      if (!doMean && i==floor(0.5*maxN) && h.first=="toaAtPV_corr") drawHist(*hp, xVal-xWdt, xVal+xWdt, v, type);
      const auto n = hp->GetEntries();
      if (n<=0) continue;
      const double qThr[3] = {(1-v68)/2., 0.5, (1+v68)/2.};
      double qVal[3] = {-1, -1, -1};
      hp->GetQuantiles(3, qVal, qThr);
      const auto err = (qVal[2]-qVal[0])/2.;
      // Compute values
      const auto yVal = (doMean ? qVal[1] : err);
      const auto yErr = (doMean ? err/sqrt(n) : err/sqrt(2*n));
      if (!doMean && yVal<5.) continue;
      // Fill graph
      gr.SetPoint(i-1, xVal, yVal);
      gr.SetPointError(i-1, xWdt, xWdt, yErr, yErr);
    }

    // Find graph range
    int minI(999), maxI(-999);
    double minX(1E6), maxX(-1E6);
    for (int i=0; i<maxN; i++) {
      const auto x = gr.GetX()[i];
      const auto ex = gr.GetErrorX(i);
      if (gr.GetErrorY(i) > 0 && (i>=(maxN-1) || gr.GetErrorY(i+1) > 0)) {
        minX = std::min(minX, x-ex);
        maxX = std::max(maxX, x+ex);
        minI = std::min(minI, i);
        maxI = std::max(maxI, i);
      }
    }
    minXM[h.first] = minX;

    // Compute constact term
    double minR(0), maxR(0), tot1W(0), tot2W(0);
    for (int i=minI; i<maxI; i++) {
      const auto ey = gr.GetErrorY(i);
      if (ey<=0.) continue;
      const auto y = gr.GetY()[i];
      const auto w = 1./(ey*ey);
      if (i < floor(minI + maxI*0.15)) {
        maxR += y*w;
        tot2W += w;
      }
      else if (i >= floor(maxI*0.9)) {
        minR += y*w;
        tot1W += w;
      }
    }
    if (tot2W > 0) maxR /= tot2W;
    if (tot1W > 0) minR /= tot1W;
    const bool useSQRT = (maxR > minR*1.05);

    const bool isPion = type.rfind("Pion")!=std::string::npos;
    if (h.first=="toaAtPV_corr" && !doMean) {
		// f[0]: A , f[1]: B , f[2]: C , f[3]: MinE
		// f[4]: A , f[5]: B , f[6]: C , f[7]: MinE
        fVal.resize(8);
        for (size_t i=0; i<4; i++) fVal[i] = 0;
        for (size_t i=0; i<4; i++) fVal[i+4] = 0;
		fVal[3] = minX;
		fVal[7] = h2D.GetXaxis()->GetBinWidth(1)/2.;
        // Perform fit
        if (useSQRT) {
          f[h.first].reset(new TF1("f", "sqrt([0]*[0]/x/x + [1]*[1])", minX, maxX));
          if (isPhoton) {
            f[h.first]->SetParameters(40, minR);
            f[h.first]->SetParLimits(0, 0, 1000);
            f[h.first]->FixParameter(1, minR);
          }
          else {
            f[h.first]->SetParameters(300, minR);
            f[h.first]->SetParLimits(0, 5, 5000);
            f[h.first]->SetParLimits(1, minR/3., minR*3.);
          }
          if (isPion)
            f[h.first]->FixParameter(1, minR);
          gr.Fit(f[h.first].get(), "MRNQP", "", minX, maxX);
          for (size_t i=0; i<2; i++) fVal[i+1] = f[h.first]->GetParameter(i);
          for (size_t i=0; i<2; i++) fVal[i+5] = f[h.first]->GetParError(i);
          if (isPhoton || isPion)
            fVal[6] = 1./sqrt(tot1W);
        }
        else {
		  f[h.first].reset(new TF1("f", "[0]", minX, maxX));
          f[h.first]->FixParameter(0, minR);
          fVal[2] = f[h.first]->GetParameter(0);
          fVal[6] =1./sqrt(tot1W);
		}
	    f[h.first]->SetLineColor(kBlack);
	    f[h.first]->SetLineWidth(2);
    }
    else if (h.first=="toaAtPV_corr" && doMean) {
	  fVal.resize(1);
	  fVal[0] = minR;
    }
    
    formatHist(gr);
    
    gr.GetXaxis()->SetTitle(h2D.GetXaxis()->GetTitle());
    gr.GetYaxis()->SetTitle(doMean ? "ToA Median [ps]" : "ToA Resolution [ps]");
    gr.GetXaxis()->SetLimits(1, maxRng);
    if (doMean)
      gr.GetYaxis()->SetRangeUser(-20, 200);
    else if (grM.size() > 1)
      gr.GetYaxis()->SetRangeUser(-20, 400);
    else
      gr.GetYaxis()->SetRangeUser(-5, 150);
    gr.SetMarkerColor(colorV[i]);
    gr.SetLineColor(colorV[i]);
    gr.SetMarkerStyle(markerV[i]);
    i += 1;
  }
  
  // Draw projection
  grM.begin()->second.Draw("AP");
  for (auto& g : grM) {
    g.second.Draw("SAMEP");
    if (f[g.first]) f[g.first]->Draw("SAME");
  }
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  TLegend leg(0.6, 0.8, 0.95, 0.9);
  if (grM.size() > 1) {
    std::map<std::string, std::string> labelM({
        {"toaAtPV_corr", "Truncated Mean at PV"}, {"toaAt0_corr", "Truncated Mean"}, {"toaDelta_avgQ", "Q-weighted |rec-sim|"}
      });
    for (auto& g : grM)
      leg.AddEntry(&g.second, labelM[g.first].c_str(), "p");
    leg.Draw("SAME");
  }

  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  const auto iniX = grM.size() > 1 ? 0.30 : 0.25;
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(iniX, 0.85, (lbl[2]+", "+lbl[0]+", "+lbl[1]).c_str());
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(iniX, 0.79, lbl[3].c_str());
  const auto iniP = grM.size() > 1 ? 0.44 : 0.85;
  for (auto& g : grM)
    if (f[g.first]) {
      tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.63, iniP-0.00, Form("q^{rec} > %g fC", v));
      tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.63, iniP-0.06, Form(doEne ? "Min. Ene.: %g GeV" : "Min. p_{T}: %g GeV", minXM[g.first]));
      tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.63, iniP-0.12, Form("A: %.2f ps*GeV^{0.5}", fVal[0]));
      tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.63, iniP-0.18, Form("B: %.2f ps*GeV", fVal[1]));
      tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.63, iniP-0.24, Form("C: %.2f ps", fVal[2]));
    }

  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs((("Plot/"+type+(doMean?"/Mean":"/Res")+"/"+c+"/")+hName+(doMean?"_mean":"_res")+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing result in: " << ("Plot/"+type+(doMean?"/Mean":"/Res")+"/root/"+hName+(doMean?"_mean":"_res")+".root") << std::endl;
  TFile file(("Plot/"+type+(doMean?"/Mean":"/Res")+"/root/"+hName+(doMean?"_mean":"_res")+".root").c_str(), "RECREATE");
  file.cd();
  size_t k(0);
  for (auto& g : grM) {
    if (f[g.first]) f[g.first]->Write(Form("fit_%lu", k));
    g.second.Write(Form("data_%lu", k));
	k += 1;
  }  
  file.Close();
};


void drawMeanVal(const std::map<float, std::vector<float>>& hValV, const int& i, const std::string& type)
{
  // Create canvas
  const std::string eName = (i==0 ? "meanVal_C" : "meanVal_XXX");
  TCanvas canvas(("c"+eName).c_str(), "", 1000, 800);
  canvas.cd();

  size_t j(0);
  TGraphAsymmErrors gr(hValV.size());
  for (const auto& v : hValV)
    gr.SetPoint(j++, v.first, v.second[i]);

  formatHist(gr);
  gr.GetYaxis()->SetTitle(i==0 ? "Bias ToA [ps]" : "NOT DEFINED");
  gr.GetXaxis()->SetTitle("q^{rec} threshold [fC]");
  
  double min(999), max(-999);
  for (int i=0; i<gr.GetN(); i++) {
    const auto y = gr.GetY()[i];
    const auto e = gr.GetErrorY(i);
    min = std::min(min, y-e);
    max = std::max(max, y+e);
  }
  gr.GetYaxis()->SetRangeUser(min-0.1*(max-min)/0.75, max+0.15*(max-min)/0.75);
  
  // Draw efficiency
  gr.Draw("AP");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);
  
  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.20, 0.85, (lbl[2]+", "+lbl[0]+", "+lbl[1]).c_str());
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.20, 0.79, lbl[3].c_str());

  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Summary/"+c+"/"+eName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing result in: " << ("Plot/"+type+"/Summary/root/"+eName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Summary/root/"+eName+".root").c_str(), "RECREATE");
  file.cd();
  gr.Write("data");
  file.Close();
};


void drawProjVal(const std::map<float, std::vector<float>>& hValV, const int& i, const std::string& type)
{
  // Create canvas
  const std::string eName = (i==0 ? "resVal_A" : (i==1 ? "resVal_B" : (i==2 ? "resVal_C" : "resVal_MinE")));
  TCanvas canvas(("c"+eName).c_str(), "", 1000, 800);
  canvas.cd();

  size_t j(0);
  TGraphAsymmErrors gr(hValV.size());
  for (const auto& v : hValV) {
    gr.SetPoint(j, v.first, v.second[i]);
    gr.SetPointError(j++, 0, 0, v.second[i+4], v.second[i+4]);
  }

  formatHist(gr);
  gr.GetYaxis()->SetTitle(i==0 ? "Stochastic #sigma_{t} [ps*GeV^{0.5}]" : (i==1 ? "Noise #sigma_{t} [ps*GeV]" : (i==2 ? "Constant #sigma_{t} [ps]" : "Minimum Energy [GeV]")));
  gr.GetXaxis()->SetTitle("q^{rec} threshold [fC]");
  
  double minX(999), maxX(-999), minY(999), maxY(-999), minYp(999), maxYp(-999);
  for (int j=0; j<gr.GetN(); j++) {
    const auto x = gr.GetX()[j];
    const auto y = gr.GetY()[j];
    const auto ey = gr.GetErrorY(j);
    if (ey>0 && (!isPion || i!=2 || y<26)) {
      minX = std::min(minX, x-1);
      maxX = std::max(maxX, x+1);
    }
    minY = std::min(minY, y-ey);
    maxY = std::max(maxY, y+ey);
    minYp = std::min(minYp, y);
    maxYp = std::max(maxYp, y);
  }
  gr.GetYaxis()->SetRangeUser(minY-0.1*(maxY-minY)/0.65, maxY+0.25*(maxY-minY)/0.65);
  
  // Draw efficiency
  gr.Draw("AP");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.20, 0.85, (lbl[2]+", "+lbl[0]+", "+lbl[3]+" , "+lbl[1]).c_str());
  tex.SetTextSize(0.045); tex.SetTextFont(62);

  // Perform fit
  std::unique_ptr<TF1> f;
  if (i==3) {
	f.reset(new TF1("f", "[0] + [1]*x + [2]*x*x", minX, maxX));
    const double m = (maxYp - minYp) / (maxX - minX);
    f->SetParameters(minYp - m*minX, m, 0);
    gr.Fit(f.get(), "MRN", "", minX, maxX);
    tex.SetTextSize(0.045); tex.SetTextFont(62);
	tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q + %.3fE-4*q*q", f->GetParameter(0), (f->GetParameter(1)<0?"-":"+"), abs(f->GetParameter(1)*100), f->GetParameter(2)*10000));
  }
  else if (abs(maxY-minY)<0.001) {
    f.reset(new TF1("f", "[0]", minX, maxX));
    f->FixParameter(0, minY);
    tex.SetTextSize(0.045); tex.SetTextFont(62);
    tex.DrawLatex(0.20, 0.78, Form("%.3f", f->GetParameter(0)));
  }
  else if (isPhoton) {
    f.reset(new TF1("f", "[0] + [1]*x", minX, maxX));
    f->SetParameters(0.1, 0.1);
    gr.Fit(f.get(), "MRN", "", minX, maxX);
    tex.SetTextSize(0.045); tex.SetTextFont(62);
    tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q", f->GetParameter(0), (f->GetParameter(1)<0?"-":"+"), abs(f->GetParameter(1)*100)));
  }
  else {
    f.reset(new TF1("f", "[0] + [1]*x + [2]*x*x", minX, maxX));
    f->SetParameters(0.1, 0.1, 0.1);
    gr.Fit(f.get(), "MRN", "", minX, maxX);
    tex.SetTextSize(0.045); tex.SetTextFont(62);
    tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q %s %.3fE-4*q*q", f->GetParameter(0), (f->GetParameter(1)<0?"-":"+"), abs(f->GetParameter(1)*100), (f->GetParameter(2)<0?"-":"+"), abs(f->GetParameter(2)*10000)));
  }
  f->SetLineColor(kBlack);
  f->SetLineWidth(2);
  f->Draw("SAME");
  
  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Parametrize/"+c+"/"+eName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing fit result in: " << ("Plot/"+type+"/Parametrize/root/"+eName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Parametrize/root/"+eName+".root").c_str(), "RECREATE");
  file.cd();
  f->Write("fit");
  gr.Write("data");
  file.Close();
};


void drawEff(TEfficiency& eff, std::array<float,8>& fVal, const float& v, const std::string& type)
{
  // Create canvas
  const std::string eName = Form("eff_%.0f", v*10);
  TCanvas canvas(("c"+eName).c_str(), "", 1000, 800);
  canvas.cd();

  eff.Draw(); 
  gPad->Update();
  const auto& p = eff.GetPaintedGraph();
  if (!p) throw std::logic_error("[ERROR] Painted graph not available for "+eName);
  auto& gr = *p;

  formatHist(gr);
  gr.GetYaxis()->SetTitle("Efficiency");
  gr.GetXaxis()->SetLimits(0, maxRng);
  gr.GetYaxis()->SetRangeUser(0, 1.2);
  
  // Draw efficiency
  gr.Draw("AP");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  // Find x range
  double minX(1E6), maxX(-1E6), minY(1E6), maxY(-1E6);
  for (int i=0; i<gr.GetN(); i++) {
    const auto x = gr.GetX()[i];
    const auto ex = gr.GetErrorX(i);
    const auto y = gr.GetY()[i];
    if (gr.GetErrorY(i) > 0) {
      minX = std::min(minX, x-ex);
      maxX = std::max(maxX, x+ex);
      minY = std::min(minY, y);
      maxY = std::max(maxY, y);
    }
  }

  const auto useErf = (maxY-minY) > 0.05*maxY;
  std::unique_ptr<TF1> f;
  if (useErf) {
    f.reset(new TF1("f", "0.5*[0]*(1 + TMath::Erf( ((x - [1])/[2] + 1)*TMath::ErfInverse(0.9) ))", minX, maxX));
    if (isPhoton && doEne) {
      f->SetParameters(1.0, 10, 4);
      f->SetParLimits(2, 0.1, 100);
    }
    else if (isPhoton && !doEne)
      f->SetParameters(0.9, 3, 1);
    else if (doEne) {
      f->SetParameters(0.93, 15, 40);
      f->SetParLimits(0, 0, 1.01);
      if (isPion)
        f->SetParLimits(2, 0.1, 100);
    }
    else
      f->SetParameters(0.9, 3, 1);
    if (isPhoton)
      f->FixParameter(0, 1.0);
    gr.Fit(f.get(), "MRNQP", "", minX, maxX);
  }
  else {
    f.reset(new TF1("f", "[0]", minX, maxX));
    if (isPhoton)
      f->FixParameter(0, 1);
    else
      gr.Fit(f.get(), "FRNQ", "", minX, maxX);
  }
  f->SetLineColor(kBlack);
  f->SetLineWidth(2);
  f->Draw("SAME");
  
  for (size_t i=0; i<6; i++) fVal[i] = 0;
  if (useErf) {
    for (size_t i=0; i<3; i++) fVal[i+0] = std::max(0., f->GetParameter(i));
    for (size_t i=0; i<3; i++) fVal[i+3] = f->GetParError(i);
  }
  else {
    fVal[0] = f->GetParameter(0);
    fVal[1] = minX;
  }

  std::array<TLine, 3> lin;
  lin[0] = TLine(0, fVal[0], maxX, fVal[0]);
  //lin[1] = TLine(fVal[1]-2*fVal[2], 0, fVal[1]-2*fVal[2], 1);
  lin[1] = TLine(fVal[1], 0, fVal[1], fVal[0]);
  for (auto& l : lin) {
    l.SetLineStyle(7);
    l.SetLineColor(kRed);
    l.SetLineWidth(2);
    l.Draw("SAME");
  }
  
  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.20, 0.85, (lbl[2]+", "+lbl[0]+", "+lbl[3]).c_str());
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.46, lbl[1].c_str());
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.40, Form("q^{rec} > %g fC", v));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.34, Form(doEne ? "Min. Ene.: %g GeV" : "Min. p_{T}: %g GeV", minX));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.28, Form("Max. Eff.: %.2f", fVal[0]));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.22, Form("Eff. Thr.: %.2f GeV", fVal[1]));
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.64, 0.16, Form("Eff. Wdt.: %.2f GeV", fVal[2]));

  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Efficiency/"+c+"/"+eName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing result in: " << ("Plot/"+type+"/Efficiency/root/"+eName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Efficiency/root/"+eName+".root").c_str(), "RECREATE");
  file.cd();
  f->Write("fit");
  gr.Write("data");
  file.Close();
};


void parmEffVal(const std::map<float, std::array<float,8>>& effValV, const int& i, const std::string& type)
{
  // Create canvas
  const std::string eName = (i==0 ? "effVal_MaxEff" : (i==1 ? "effVal_Threshold" : "effVal_Width"));
  TCanvas canvas(("c"+eName).c_str(), "", 1000, 800);
  canvas.cd();

  size_t j(0);
  TGraphAsymmErrors gr(effValV.size());
  for (const auto& v : effValV) {
    gr.SetPoint(j, v.first, v.second[i]);
    gr.SetPointError(j++, 0, 0,  v.second[i+3], v.second[i+3]);
  }

  formatHist(gr);
  gr.GetYaxis()->SetTitle(i==0 ? "Maximum Efficiency" : (i==1 ? "Efficiency Threshold [GeV]" : "Efficiency Width [GeV]"));
  gr.GetXaxis()->SetTitle("q^{rec} threshold [fC]");
  if (i==0) gr.GetYaxis()->SetTitleOffset(1.0);

  bool doLine(false);
  const int minN(0), maxN(gr.GetN()-1);
  double minX(999), maxX(-999), minY(999), maxY(-999);
  for (int j=0; j<gr.GetN(); j++) {
    const auto y = gr.GetY()[j];
    const auto ey = gr.GetErrorY(j);
    minY = std::min(minY, y-ey);
    maxY = std::max(maxY, y+ey);
    if (ey>0 && (i!=2 || y<100)) {
      const auto x = gr.GetX()[j] + (gr.GetX()[j] - gr.GetX()[j-1])/2.;
      minX = std::min(minX, gr.GetX()[j] - (j>minN ? (gr.GetX()[j] - gr.GetX()[j-1])/2. : 1));
      maxX = std::max(maxX, gr.GetX()[j] + (j<maxN ? (gr.GetX()[j+1] - gr.GetX()[j])/2. : 1));
    }
    else doLine = true;
  }
  if (minX > maxX) { minX = gr.GetX()[minN]-1; maxX = gr.GetX()[maxN]+1; }
  gr.GetYaxis()->SetRangeUser(minY-0.1*(maxY-minY)/0.65, maxY+0.25*(maxY-minY)/0.65);

  // Draw efficiency
  gr.Draw("AP");
  CMS_lumi(&canvas, 33, "", "", false, -1, false);

  // Add text
  TLatex tex; tex.SetNDC();
  tex.SetTextSize(0.045); tex.SetTextFont(61); tex.DrawLatex(0.16, 0.94, "CMS");
  tex.SetTextSize(0.045); tex.SetTextFont(52); tex.DrawLatex(0.24, 0.94, "Simulation Preliminary");
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.73, 0.94, "#bf{Phase-2 HGCal}");
  const auto lbl = getCut(type);
  tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(0.20, 0.85, (lbl[2]+", "+lbl[0]+", "+lbl[3]+" , "+lbl[1]).c_str());

  // Fit
  std::vector<std::unique_ptr<TF1>> fV;
  if (isPhoton) {
    if (i==0) {
      fV.emplace_back(new TF1("f", "[0]", minX, maxX));
      gr.Fit(fV[0].get(), "FRNWQ", "", minX, maxX);
      tex.SetTextSize(0.045); tex.SetTextFont(62);
      tex.DrawLatex(0.20, 0.78, Form("%.3f", fV[0]->GetParameter(0)));
    }
    else if (doLine) {
      fV.emplace_back(new TF1("f", "[0]", 0, minX));
      fV[0]->SetParameters(minY, 0);
      gr.Fit(fV[0].get(), "FRNWQ", "", 0, minX);
      fV.emplace_back(new TF1("f2", "[0] + [1]*x + [2]*x*x", minX, maxX));
      fV[1]->SetParameters(minY, 0, 0);
      gr.Fit(fV[1].get(), "MRNQP", "", minX, maxX);
      tex.SetTextSize(0.039); tex.SetTextFont(62);
      tex.DrawLatex(0.20, 0.78, Form("q>%g ? (%.3f %s %.3fE-2*q %s %.3fE-4*q*q) : %.1f", minX, fV[1]->GetParameter(0), (fV[1]->GetParameter(1)<0?"-":"+"), abs(fV[1]->GetParameter(1)*100), (fV[1]->GetParameter(2)<0?"-":"+"), abs(fV[1]->GetParameter(2)*100*100), fV[0]->GetParameter(0)));
    }
    else {
      fV.emplace_back(new TF1("f", "[0] + [1]*x + [2]*x*x", minX, maxX));
      fV[0]->SetParameters(minY, 0.1, 0.1);
      gr.Fit(fV[0].get(), "MRNQP", "", minX, maxX);
      tex.SetTextSize(0.045); tex.SetTextFont(62);
      tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q %s %.3fE-4*q*q", fV[0]->GetParameter(0), (fV[0]->GetParameter(1)<0?"-":"+"), abs(fV[0]->GetParameter(1)*100), (fV[0]->GetParameter(2)<0?"-":"+"), abs(fV[0]->GetParameter(2)*100*100)));
    }
  }
  else {//if (type.find("eta2p8")!=std::string::npos || type.find("eta2p5")!=std::string::npos) {
	fV.emplace_back(new TF1("f", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", minX, maxX));
    fV[0]->SetParameters(minY, 0.1, 0.1, 0.0);
    gr.Fit(fV[0].get(), "MRNQP", "", minX, maxX);
    tex.SetTextSize(0.039); tex.SetTextFont(62);
    tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q %s %.3fE-4*q*q %s %.3fE-6*q*q*q", fV[0]->GetParameter(0), (fV[0]->GetParameter(1)<0?"-":"+"), abs(fV[0]->GetParameter(1)*100), (fV[0]->GetParameter(2)<0?"-":"+"), abs(fV[0]->GetParameter(2)*100*100), (fV[0]->GetParameter(3)<0?"-":"+"), abs(fV[0]->GetParameter(3)*100*100*100)));
  }
  /*
  else {
	fV.emplace_back(new TF1("f", "[0] + [1]*x + [2]*x*x", minX, maxX));
    fV[0]->SetParameters(minY, 0.1, 0.1);
    gr.Fit(fV[0].get(), "MRNQP", "", minX, maxX);
    tex.SetTextSize(0.045); tex.SetTextFont(62);
    tex.DrawLatex(0.20, 0.78, Form("%.3f %s %.3fE-2*q %s %.3fE-4*q*q", fV[0]->GetParameter(0), (fV[0]->GetParameter(1)<0?"-":"+"), abs(fV[0]->GetParameter(1)*100), (fV[0]->GetParameter(2)<0?"-":"+"), abs(fV[0]->GetParameter(2)*100*100)));
  }
  */

  for (auto& f : fV) {
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
    f->Draw("SAME");
  }

  // Save canvas
  for (auto c : {"png", "pdf", "C"})
    canvas.SaveAs(("Plot/"+type+"/Parametrize/"+c+"/"+eName+"."+c).c_str());
  canvas.Close();

  std::cout << "Storing fit result in: " << ("Plot/"+type+"/Parametrize/root/"+eName+".root") << std::endl;
  TFile file(("Plot/"+type+"/Parametrize/root/"+eName+".root").c_str(), "RECREATE");
  file.cd();
  fV[0]->Write("fit");
  for (size_t i=1; i<fV.size(); i++) fV[i]->Write(Form("fit_%lu", i));
  gr.Write("data");
  file.Close();
};


std::vector<size_t> decrease_sorted_indices(const std::vector<float>& v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v (decreasing order)
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
};
  

//time-interval based on that ~210ps wide and with the highest number of hits
//extension valid in high PU of taking smallest interval with (order of)68% of hits
std::vector<float> fixSizeHighestDensity(const std::vector<float>& time, std::vector<float> weight,
					 const unsigned int& minNhits, const float& deltaT, const float& timeWidthBy)
{
  if (time.size() < minNhits)
    return std::vector<float>({-99., -1., -99});

  if (weight.empty())
    weight.resize(time.size(), 1.);

  std::vector<float> t(time.size(), 0.);
  std::vector<float> w(time.size(), 0.);
  std::vector<size_t> sortedIndex = decrease_sorted_indices(time);
  for (std::size_t i = 0; i < sortedIndex.size(); ++i) {
    t[i] = time[sortedIndex[i]];
    w[i] = weight[sortedIndex[i]];
  }

  int max_elements = 0;
  int start_el = 0;
  int end_el = 0;
  float timeW = 0.f;
  float tolerance = 0.05f;

  for (auto start = t.begin(); start != t.end(); ++start) {
    const auto startRef = *start;
    int c = count_if(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
    if (c > max_elements) {
      max_elements = c;
      auto last_el = find_if_not(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
      auto valTostartDiff = *(--last_el) - startRef;
      if (std::abs(deltaT - valTostartDiff) < tolerance) {
        tolerance = std::abs(deltaT - valTostartDiff);
      }
      start_el = distance(t.begin(), start);
      end_el = distance(t.begin(), last_el);
      timeW = valTostartDiff;
    }
  }

  // further adjust time width around the chosen one based on the hits density
  // proved to improve the resolution: get as many hits as possible provided they are close in time
  float HalfTimeDiff = timeW * timeWidthBy;
  float sum = 0.;
  float num = 0;
  int totSize = t.size();
  
  for (int ij = 0; ij <= start_el; ++ij) {
    if (t[ij] > (t[start_el] - HalfTimeDiff)) {
      for (int kl = ij; kl < totSize; ++kl) {
        if (t[kl] < (t[end_el] + HalfTimeDiff)) {
          sum += t[kl] * w[kl];
          num += w[kl];
        } else
          break;
      }
      break;
    }
  }

  if (num == 0)
    return std::vector<float>({-99., -1., -99});
  
  return std::vector<float>({sum / num, float(1. / sqrt(num)), float(totSize), t[start_el] - HalfTimeDiff, t[end_el] + HalfTimeDiff});
};


#endif // #ifndef plotTOA_C
