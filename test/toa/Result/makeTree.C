#ifndef makeTree_C
#define makeTree_C

// Custom headers
#include "../Macros/TreeReader.h"
// ROOT headers
#include "TTree.h"
#include "TSystem.h"
// C headers
#include <map>


constexpr double c_cm_ns = 29.9792; //speed of light in cm / ns


//void makeTree(const std::string& inputFile="../Files/SinglePhotonGun_eta1p8_13_1_0_pre4_D99_startup.root", const size_t& opt=1, const std::string& inputNTuple="ana/hits")
void makeTree(const std::string& dir="/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/SingleK0LGun_eta1p8_13_1_0_pre4_D99_startup_15ns", const size_t& opt=1, const std::string& inputNTuple="ana/hits")
{ 
  // Extract ntuple

  std::vector<std::string> files;
  if(dir.find(".root")!=std::string::npos) {
    files.push_back(dir);
  }
  else{
    void *dirp = gSystem->OpenDirectory(dir.c_str());
    Char_t *afile;
    while ((afile = const_cast<Char_t *>(gSystem->GetDirEntry(dirp)))) {
      std::string url(dir + "/" + afile);
      if(url.find(".root")==std::string::npos) continue;
      files.push_back( url );
    }
  }
  
  TreeReader reader(files, inputNTuple, false);
  const auto& nentries = reader.getEntries();
  const bool isPhoton(dir.find("p22")!=std::string::npos || dir.find("hoton")!=std::string::npos);
  const bool isPion(dir.find("Pion")!=std::string::npos);

  // Extract event number
  std::vector<int> keys;
  std::map<int, std::vector<Long64_t> > evtMap;
  std::cout << "[INFO] First round: processing events from " << nentries << " entries" << std::endl;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    reader.setEntry(jentry);
    if(!(jentry % (nentries/10))) std::cout << "Processed " << jentry << " entrries out of " << nentries << std::endl;
    const auto iEvent = reader.getInt("event");
    auto ele = evtMap.find(iEvent);
    if (ele==evtMap.end()) {
      keys.push_back(iEvent);
      evtMap[iEvent].emplace_back(jentry);
    }
    else { ele->second.emplace_back(jentry); }
  }
  
  // Define output file
  std::string o1 = std::string(gSystem->BaseName(dir.c_str())) + ".root";
  gSystem->mkdir("Output");
  const auto outputFile = "Output/TREE_" + std::string(opt==2 ? "InRad2p0_" : (opt==1 ? "InRadius_" : "")) + o1;
  auto output = std::unique_ptr<TFile>(new TFile(outputFile.c_str(), "recreate"));
  if (!output || !output->IsOpen() || output->IsZombie()) throw std::logic_error("[ERROR] Failed to create output file");
  output->cd();

  // Define output tree
  TTree outTree("toa", "toa");
  int iEvent(-1);
  outTree.Branch("event", &iEvent);
  std::map<std::string, float> genVarM;
  for (const auto& v : {"gpt", "genergy"})
    outTree.Branch(v, &(genVarM[v]));
  std::map<std::string, std::vector<float> > varM;
  for (const auto& v : {"qsim", "qrec", "toarecAtPV"})
    outTree.Branch(v, &(varM[v]));
  
  // Loop over events
  const auto nEntries = keys.size();
  std::cout << "[INFO] Processing events from " << nEntries << " entries" << std::endl;
  for (size_t i=0; i<nEntries; i++) {
    if(!(i % (nEntries/10))) std::cout << "Processed " << i << " entrries out of " << nEntries << std::endl;
    // Initialize information
    iEvent = -1;
    for (auto& v : genVarM) v.second = -99;
    for (auto& v : varM) v.second.clear();
    // Fill information
    for (const auto& jentry : evtMap.at(keys[i])) {
      // Get ntuple entry
      reader.setEntry(jentry);
      // Pre-select hits
      // Exclude scintillators
      if (reader.getInt("isSci")!=0) continue;
      // Select hits that cross calorimeter
      if (reader.getInt("crossCalo")!=1) continue;
      // Select within 0.4 jet radius
      if (reader.getInt("inShower")!=1) continue;
      // Require in-time BX info
      if (reader.getFloat("toasim")<-50) continue;
      // Require photons to only be within CE-E
      if (isPhoton && reader.getInt("layer")>26) continue;
      // Require positive eta side
      if (reader.getFloat("geta") < 0) continue;
      // Matched to LC
      //if (reader.getInt("matchedToLC")!=1) continue;
      // Require to be within radius
      if (opt>0) {
        const auto aZ = std::abs(reader.getFloat("z"));
        //const auto thr = opt==2 ? 2. : (opt==1 ? ((aZ < 362.18) ? 3. : ((aZ < 437.87) ? 8. : 16.)) : 0);
        const auto thr = opt==2 ? 2. : (opt==1 ? 3. : 0);
        if (!isPion && reader.getFloat("gdradius")>=thr) continue;
      }
      // Process event
      if (iEvent<0) {
	    iEvent = reader.getInt("event");
	    for (auto& v : genVarM) v.second = reader.getFloat(v.first.c_str());
      }
      const auto beta = 1;//reader.getFloat("gbeta");
      const auto gvz = reader.getFloat("gvz");
      const auto gvt = reader.getFloat("gvt");
      const auto r = reader.getFloat("radius");
      const auto z = reader.getFloat("z");
      const auto dt_PV = sqrt(r*r + (z - gvz)*(z - gvz))/beta/c_cm_ns + gvt;
      for (auto& v : varM) {
	    if (v.first=="toarecAtPV") v.second.push_back(reader.getFloat("toarec") - dt_PV);
	    else v.second.push_back(reader.getFloat(v.first));
      }
    }
    // Store event
    outTree.Fill();
  }

  // Store tree
  outTree.Write();

  // Close output file
  output->Close();
};


#endif // #ifndef makeTree_C
