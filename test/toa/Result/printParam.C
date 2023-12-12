#ifndef printParam_C
#define printParam_C

// Custom headers
// ROOT headers
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
// C headers
#include <iostream>
#include <fstream>
#include <dirent.h>


void printParam(const std::string& inDir="./Plot/InRadius/")
{
  const std::map<std::string, std::string> parV({{"K0L", "SingleK0L"}, {"Photon", "SinglePhoton"}, {"Pion", "SinglePion"}});
  const std::map<std::string, std::string> varV({{"effVal_MaxEff", "eff_M"}, {"effVal_Threshold", "eff_O"}, {"effVal_Width", "eff_W"},
                                                 {"resVal_A", "res_A"}, {"resVal_B", "res_B"}, {"resVal_C", "res_C"}, {"resVal_MinE", "res_minE"}});

  // find all directories
  std::map<bool, std::map<std::string, std::map<float, std::string>>> pathM;
  const auto& dir = opendir(inDir.c_str());
  auto entry = readdir(dir);
  while (entry!=NULL) {
    std::string d_name = entry->d_name;
    if (entry->d_type==DT_DIR && d_name!="." && d_name!="..") {
      const bool isAged = (d_name.rfind("_3iab")!=std::string::npos);
      if (!isAged && d_name.rfind("_startup")==std::string::npos)
        throw std::logic_error("Invalid aged condition type for "+d_name);
      std::string par("");
      for (const auto& p : parV)
        if (d_name.find(p.first)!=std::string::npos) { par = p.second; break; }
      if (par=="")
        throw std::logic_error("Invalid particle type for "+d_name);
      auto d = d_name.substr(d_name.find("eta")+3);
      std::replace(d.begin(), d.end(), 'p', '.');
      const auto eta = std::stof(d.substr(0, d.find("_")));
      pathM[isAged][par][eta] = inDir + d_name + "/Parametrize/root/";
    }
    entry = readdir(dir);
  }
  closedir(dir);

  // create file
  std::ofstream myfile;
  myfile.open(inDir+"param.py");
  myfile << "def _params(q=0, isStartUp=False):" << std::endl;
  myfile << "\tparams={}" << std::endl;
  for (const auto& a : pathM) {
    if (!a.first)
      myfile << "\tif isStartUp:" << std::endl;
    else
      myfile << "\telse:" << std::endl;
    myfile << "\t\tparams={" << std::endl;
    for (const auto& p : a.second) {
      myfile << "\t\t\t'" << p.first << "':{" << std::endl;
      for (const auto& v : varV) {
        myfile << "\t\t\t\t'" << v.second << "':{" << std::endl;
        for (const auto& e : p.second) {
          const auto f_name = e.second+v.first+".root";
          TFile file(f_name.c_str(), "READ");
          if (!file.IsOpen() || file.IsZombie())
            throw std::logic_error("Failed to open file: "+f_name);
          const auto& f0 = static_cast<TF1*>(file.Get("fit"));
          const auto& f = file.Get("fit_1") ? static_cast<TF1*>(file.Get("fit_1")) : f0;
          if (!f)
            throw std::logic_error("Fit not found in: "+f_name);
          std::string func("");
          if (f->GetNpar()>0)
            func = Form("%.4f", f->GetParameter(0));
          for (int i=1; i<f->GetNpar(); i++) {
            std::string q("q");
            for (int j=1; j<i; j++) q += "*q";
            const int n(2*i);
            func += Form(" %s %.4fE-%d*%s", (f->GetParameter(i)<0?"-":"+"), abs(f->GetParameter(i)*pow(10, n)), n, q.c_str());
          }
          if (f0!=f)
            func = Form("(%s) if q>%.1f else %.1f", func.c_str(), f->GetXmin(), abs(f0->GetParameter(0))); 
          if (v.second=="res_minE" || v.second=="res_B")
            func = "max("+func+", 0)";
          else if (v.second=="eff_W")
            func = "max("+func+", 1E-4)";
          file.Close();
          myfile << "\t\t\t\t\t" << Form("%.1f", e.first) << ":" << func << "," << std::endl;
        }
        myfile << "\t\t\t\t}," << std::endl;
      }
      myfile << "\t\t\t}," << std::endl;
    }
    myfile << "\t\t}" << std::endl;
  }
  myfile << "\treturn params" << std::endl;
  myfile.close();
}


#endif // #ifndef printParam_C
