#ifndef TreeReader_h
#define TreeReader_h

#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <typeinfo>


class TreeReader
{
 public:
  TreeReader(){};

  TreeReader(const std::vector<std::string>& paths, const std::string& treeName, const bool& setAllBranches=true)
  {
    init(paths, treeName, setAllBranches);
  };

  TreeReader(const std::string& path, const std::string& treeName, const bool& setAllBranches=true)
  {
    init({path}, treeName, setAllBranches);
  };

  ~TreeReader()
  {
  };

  // getters
  TTree* getTree() { return reader_->GetTree(); };
  const Long64_t&  getEntries() const { return nEntries_;  };
  const Int_t&     getTreeNumber() const { return treeIndex_; };

  const UChar_t&   getUChar  (const std::string& s) const { return getValue(objUChar_,   s); };
  const UShort_t&  getUShort (const std::string& s) const { return getValue(objUShort_,  s); };
  const UInt_t&    getUInt   (const std::string& s) const { return getValue(objUInt_,    s); };
  const ULong64_t& getULong64(const std::string& s) const { return getValue(objULong64_, s); };
  const bool&      getBool   (const std::string& s) const { return getValue(objBool_,    s); };
  const char&      getChar   (const std::string& s) const { return getValue(objChar_,    s); };
  const short&     getShort  (const std::string& s) const { return getValue(objShort_,   s); };
  const int&       getInt    (const std::string& s) const { return getValue(objInt_,     s); };
  const Long64_t&  getLong64 (const std::string& s) const { return getValue(objLong64_,  s); };
  const float&     getFloat  (const std::string& s) const { return getValue(objFloat_,   s); };
  const double&    getDouble (const std::string& s) const { return getValue(objDouble_,  s); };

  const TTreeReaderArray<UChar_t>&   getAUChar  (const std::string& s) const { return getArray(arrUChar_,   s); };
  const TTreeReaderArray<UShort_t>&  getAUShort (const std::string& s) const { return getArray(arrUShort_,  s); };
  const TTreeReaderArray<UInt_t>&    getAUInt   (const std::string& s) const { return getArray(arrUInt_,    s); };
  const TTreeReaderArray<ULong64_t>& getAULong64(const std::string& s) const { return getArray(arrULong64_, s); };
  const TTreeReaderArray<bool>&      getABool   (const std::string& s) const { return getArray(arrBool_,    s); };
  const TTreeReaderArray<char>&      getAChar   (const std::string& s) const { return getArray(arrChar_,    s); };
  const TTreeReaderArray<short>&     getAShort  (const std::string& s) const { return getArray(arrShort_,   s); };
  const TTreeReaderArray<int>&       getAInt    (const std::string& s) const { return getArray(arrInt_,     s); };
  const TTreeReaderArray<Long64_t>&  getALong64 (const std::string& s) const { return getArray(arrLong64_,  s); };
  const TTreeReaderArray<float>&     getAFloat  (const std::string& s) const { return getArray(arrFloat_,   s); };
  const TTreeReaderArray<double>&    getADouble (const std::string& s) const { return getArray(arrDouble_,  s); };

  const std::vector<UChar_t>&   getVUChar  (const std::string& s) const { return getValue(vecUChar_,   s); };
  const std::vector<UShort_t>&  getVUShort (const std::string& s) const { return getValue(vecUShort_,  s); };
  const std::vector<UInt_t>&    getVUInt   (const std::string& s) const { return getValue(vecUInt_,    s); };
  const std::vector<ULong64_t>& getVULong64(const std::string& s) const { return getValue(vecULong64_, s); };
  const std::vector<bool>&      getVBool   (const std::string& s) const { return getValue(vecBool_,    s); };
  const std::vector<char>&      getVChar   (const std::string& s) const { return getValue(vecChar_,    s); };
  const std::vector<short>&     getVShort  (const std::string& s) const { return getValue(vecShort_,   s); };
  const std::vector<int>&       getVInt    (const std::string& s) const { return getValue(vecInt_,     s); };
  const std::vector<Long64_t>&  getVLong64 (const std::string& s) const { return getValue(vecLong64_,  s); };
  const std::vector<float>&     getVFloat  (const std::string& s) const { return getValue(vecFloat_,   s); };
  const std::vector<double>&    getVDouble (const std::string& s) const { return getValue(vecDouble_,  s); };

  const std::vector<std::vector<UChar_t> >& getVVUChar (const std::string& s) const { return getValue(vecVUChar_, s); };
  const std::vector<std::vector<int> >&     getVVInt   (const std::string& s) const { return getValue(vecVInt_,   s); };
  const std::vector<std::vector<float> >&   getVVFloat (const std::string& s) const { return getValue(vecVFloat_, s); };

  // setters
  void setChain(const std::vector<std::string>& paths, const std::string& treeName)
  {
    chain_.reset(new TChain(treeName.c_str()));
    for (const auto& path : paths) {
      const auto& file = TFile::Open(path.c_str(), "READ");
      if (!file || !file->IsOpen() || file->IsZombie() || !chain_->Add(path.c_str())) {
        throw std::runtime_error("[ERROR] File "+path+" failed to open!");
      }
      file->Close();
      chain_->AddFile(path.c_str());
    }
  };

  void setEntry(const Long64_t& i)
  { 
    if (i<0 || i>=nEntries_) { throw std::runtime_error(Form("[ERROR] Invalid index: %lld", i)); }
    if (i==index_) return;
    index_ = i;
    loadTreeEntry(reader_.get());
    if (chain_) treeIndex_ = chain_->GetTreeNumber();
  };

 private:
  void init(const std::vector<std::string>& paths, const std::string& treeName, const bool& setAllBranches)
  {
    setChain(paths, treeName);
    setTreeReader();
    if (!reader_) throw std::runtime_error("[ERROR] TTreeReader is null!");
    if (setAllBranches) setBranches();
  };
  
  // getters
  template<typename T>
  const T& getValue(const std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<T> > >& m, const std::string& s) const
  {
    const auto& p = m.find(s);
    if (p!=m.end()) return *p->second->Get();
    return *loadBranch(m, s)->Get();
  };

  template<typename T>
  const TTreeReaderArray<T>& getArray(const std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<T> > >& m, const std::string& s) const
  {
    const auto& p = m.find(s);
    if (p!=m.end()) return *p->second;
    return *loadBranch(m, s);
  };

  void checkValue(ROOT::Internal::TTreeReaderValueBase* value) const
  {
    if (value==NULL) {
      throw std::runtime_error("[ERROR] Value pointer is null");
    }
    else if (value->GetSetupStatus() < 0) {
      throw std::runtime_error(Form("[ERROR] Setup status %d when setting up reader for %s", value->GetSetupStatus(), value->GetBranchName()));
    }
    else if (value->GetReadStatus() > 1) {
      throw std::runtime_error(Form("[ERROR] Read status %d when setting up reader for %s", value->GetReadStatus(), value->GetBranchName()));
    }
  };

  void checkValues() const
  { 
    for (const auto& o : objUChar_  ) { checkValue(o.second.get()); }
    for (const auto& o : objUShort_ ) { checkValue(o.second.get()); }
    for (const auto& o : objUInt_   ) { checkValue(o.second.get()); }
    for (const auto& o : objULong64_) { checkValue(o.second.get()); }
    for (const auto& o : objBool_   ) { checkValue(o.second.get()); }
    for (const auto& o : objChar_   ) { checkValue(o.second.get()); }
    for (const auto& o : objShort_  ) { checkValue(o.second.get()); }
    for (const auto& o : objInt_    ) { checkValue(o.second.get()); }
    for (const auto& o : objLong64_ ) { checkValue(o.second.get()); }
    for (const auto& o : objFloat_  ) { checkValue(o.second.get()); }
    for (const auto& o : objDouble_ ) { checkValue(o.second.get()); }
    for (const auto& o : arrUChar_  ) { checkValue(o.second.get()); }
    for (const auto& o : arrUShort_ ) { checkValue(o.second.get()); }
    for (const auto& o : arrUInt_   ) { checkValue(o.second.get()); }
    for (const auto& o : arrULong64_) { checkValue(o.second.get()); }
    for (const auto& o : arrBool_   ) { checkValue(o.second.get()); }
    for (const auto& o : arrChar_   ) { checkValue(o.second.get()); }
    for (const auto& o : arrShort_  ) { checkValue(o.second.get()); }
    for (const auto& o : arrInt_    ) { checkValue(o.second.get()); }
    for (const auto& o : arrLong64_ ) { checkValue(o.second.get()); }
    for (const auto& o : arrFloat_  ) { checkValue(o.second.get()); }
    for (const auto& o : arrDouble_ ) { checkValue(o.second.get()); }
    for (const auto& o : vecUChar_  ) { checkValue(o.second.get()); }
    for (const auto& o : vecUShort_ ) { checkValue(o.second.get()); }
    for (const auto& o : vecUInt_   ) { checkValue(o.second.get()); }
    for (const auto& o : vecULong64_) { checkValue(o.second.get()); }
    for (const auto& o : vecBool_   ) { checkValue(o.second.get()); }
    for (const auto& o : vecChar_   ) { checkValue(o.second.get()); }
    for (const auto& o : vecShort_  ) { checkValue(o.second.get()); }
    for (const auto& o : vecInt_    ) { checkValue(o.second.get()); }
    for (const auto& o : vecLong64_ ) { checkValue(o.second.get()); }
    for (const auto& o : vecFloat_  ) { checkValue(o.second.get()); }
    for (const auto& o : vecDouble_ ) { checkValue(o.second.get()); }
    for (const auto& o : vecVUChar_ ) { checkValue(o.second.get()); }
    for (const auto& o : vecVInt_   ) { checkValue(o.second.get()); }
    for (const auto& o : vecVFloat_ ) { checkValue(o.second.get()); }
  };

  void loadTreeEntry(TTreeReader* r) const
  {
    if (!r) return;
    const auto status = r->SetEntry(index_);
    if (status!=TTreeReader::kEntryValid) {
      std::string msg = "";
      if      (status==TTreeReader::kEntryNotLoaded) { msg = "no entry has been loaded yet"; }
      else if (status==TTreeReader::kEntryNoTree) { msg = "the tree does not exist"; }
      else if (status==TTreeReader::kEntryNotFound) { msg = "the tree entry number does not exist"; }
      else if (status==TTreeReader::kEntryChainSetupError) { msg = "problem in accessing a chain element, e.g. file without the tree"; }
      else if (status==TTreeReader::kEntryChainFileError) { msg = "problem in opening a chain's file"; }
      else if (status==TTreeReader::kEntryDictionaryError) { msg = "problem reading dictionary info from tree"; }
      throw std::runtime_error("[ERROR] Invalid entry: "+msg);
    }
  };

  template<typename T>
  const std::unique_ptr<T>& loadBranch(const std::unordered_map<std::string, std::unique_ptr<T> >& m, const std::string& s) const
  {
    reader_->Restart();
    auto& p = (*const_cast<std::unordered_map<std::string, std::unique_ptr<T> >* >(&m))[s];
    p.reset(new T(*reader_, s.c_str()));
    loadTreeEntry(reader_.get());
    if (checkValue_) checkValue(p.get());
    return p;
  };

  // setters
  void setTreeReader()
  {
    if (reader_ || !chain_) return;
    reader_.reset(new TTreeReader(chain_.get()));
    if (!reader_->GetTree()) {
      throw std::runtime_error(Form("[ERROR] Failed to open chain %s !", chain_->GetName()));
    }
    index_ = -1;
    nEntries_ = reader_->GetEntries(true);
  };

  void setBranches()
  {
    // loop over input tree branches
    const auto& branchList = reader_->GetTree()->GetListOfBranches();
    for (int b=0; b<branchList->GetEntries(); b++) {
      const auto& branch = branchList->At(b);
      const auto& brName = branch->GetName();
      std::string brType = branch->GetTitle();

      // TEMPORAL FIX TO STUPID BUG
      if (std::string(brName).rfind("HPMuon",0)==0) continue;

      if (brType.rfind("]/")!=std::string::npos) { brType = "A"+brType.substr(brType.rfind("/")+1); }
      else if (brType.rfind("/")!=std::string::npos) { brType = brType.substr(brType.rfind("/")+1); }
      else { brType = reader_->GetTree()->GetBranch(brName)->GetClassName(); }

      if      (brType=="b" ) { objUChar_[brName].reset(new TTreeReaderValue<UChar_t>(*reader_, brName)); }
      else if (brType=="s" ) { objUShort_[brName].reset(new TTreeReaderValue<UShort_t>(*reader_, brName)); }
      else if (brType=="i" ) { objUInt_[brName].reset(new TTreeReaderValue<UInt_t>(*reader_, brName)); }
      else if (brType=="l" ) { objULong64_[brName].reset(new TTreeReaderValue<ULong64_t>(*reader_, brName)); }
      else if (brType=="O" ) { objBool_[brName].reset(new TTreeReaderValue<bool>(*reader_, brName)); }
      else if (brType=="B" ) { objChar_[brName].reset(new TTreeReaderValue<char>(*reader_, brName)); }
      else if (brType=="S" ) { objShort_[brName].reset(new TTreeReaderValue<short>(*reader_, brName)); }
      else if (brType=="I" ) { objInt_[brName].reset(new TTreeReaderValue<int>(*reader_, brName)); }
      else if (brType=="L" ) { objLong64_[brName].reset(new TTreeReaderValue<Long64_t>(*reader_, brName)); }
      else if (brType=="F" ) { objFloat_[brName].reset(new TTreeReaderValue<float>(*reader_, brName)); }
      else if (brType=="D" ) { objDouble_[brName].reset(new TTreeReaderValue<double>(*reader_, brName)); }
      else if (brType=="Ab") { arrUChar_[brName].reset(new TTreeReaderArray<UChar_t>(*reader_, brName)); }
      else if (brType=="As") { arrUShort_[brName].reset(new TTreeReaderArray<UShort_t>(*reader_, brName)); }
      else if (brType=="Ai") { arrUInt_[brName].reset(new TTreeReaderArray<UInt_t>(*reader_, brName)); }
      else if (brType=="Al") { arrULong64_[brName].reset(new TTreeReaderArray<ULong64_t>(*reader_, brName)); }
      else if (brType=="AO") { arrBool_[brName].reset(new TTreeReaderArray<bool>(*reader_, brName)); }
      else if (brType=="AB") { arrChar_[brName].reset(new TTreeReaderArray<char>(*reader_, brName)); }
      else if (brType=="AS") { arrShort_[brName].reset(new TTreeReaderArray<short>(*reader_, brName)); }
      else if (brType=="AI") { arrInt_[brName].reset(new TTreeReaderArray<int>(*reader_, brName)); }
      else if (brType=="AL") { arrLong64_[brName].reset(new TTreeReaderArray<Long64_t>(*reader_, brName)); }
      else if (brType=="AF") { arrFloat_[brName].reset(new TTreeReaderArray<float>(*reader_, brName)); }
      else if (brType=="AD") { arrDouble_[brName].reset(new TTreeReaderArray<double>(*reader_, brName)); }
      else if (brType=="vector<unsigned char>"    ) { vecUChar_[brName].reset(new TTreeReaderValue<std::vector<UChar_t> >(*reader_, brName)); }
      else if (brType=="vector<unsigned short>"   ) { vecUShort_[brName].reset(new TTreeReaderValue<std::vector<UShort_t> >(*reader_, brName)); }
      else if (brType=="vector<unsigned int>"     ) { vecUInt_[brName].reset(new TTreeReaderValue<std::vector<UInt_t> >(*reader_, brName)); }
      else if (brType=="vector<unsigned long int>") { vecULong64_[brName].reset(new TTreeReaderValue<std::vector<ULong64_t> >(*reader_, brName)); }
      else if (brType=="vector<bool>"             ) { vecBool_[brName].reset(new TTreeReaderValue<std::vector<bool> >(*reader_, brName)); }
      else if (brType=="vector<char>"             ) { vecChar_[brName].reset(new TTreeReaderValue<std::vector<char> >(*reader_, brName)); }
      else if (brType=="vector<short>"            ) { vecShort_[brName].reset(new TTreeReaderValue<std::vector<short> >(*reader_, brName)); }
      else if (brType=="vector<int>"              ) { vecInt_[brName].reset(new TTreeReaderValue<std::vector<int> >(*reader_, brName)); }
      else if (brType=="vector<long int>"         ) { vecLong64_[brName].reset(new TTreeReaderValue<std::vector<Long64_t> >(*reader_, brName)); }
      else if (brType=="vector<float>"            ) { vecFloat_[brName].reset(new TTreeReaderValue<std::vector<float> >(*reader_, brName)); }
      else if (brType=="vector<double>"           ) { vecDouble_[brName].reset(new TTreeReaderValue<std::vector<double> >(*reader_, brName)); }
      else if (brType=="vector<vector<unsigned char> >") { vecVUChar_[brName].reset(new TTreeReaderValue<std::vector<std::vector<UChar_t> > >(*reader_, brName)); }
      else if (brType=="vector<vector<int> >") { vecVInt_[brName].reset(new TTreeReaderValue<std::vector<std::vector<int> > >(*reader_, brName)); }
      else if (brType=="vector<vector<float> >") { vecVFloat_[brName].reset(new TTreeReaderValue<std::vector<std::vector<float> > >(*reader_, brName)); }
    }
  };

  std::unique_ptr<TChain> chain_;
  std::unique_ptr<TTreeReader> reader_;
  bool checkValue_;
  Long64_t index_, nEntries_;
  Int_t treeIndex_;
  
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UChar_t  > > > objUChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UShort_t > > > objUShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t   > > > objUInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<ULong64_t> > > objULong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<bool     > > > objBool_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Char_t   > > > objChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<short    > > > objShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<int      > > > objInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Long64_t > > > objLong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<float    > > > objFloat_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<double   > > > objDouble_;

  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UChar_t  > > > arrUChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UShort_t > > > arrUShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UInt_t   > > > arrUInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<ULong64_t> > > arrULong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<bool     > > > arrBool_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Char_t   > > > arrChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<short    > > > arrShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<int      > > > arrInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Long64_t > > > arrLong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<float    > > > arrFloat_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<double   > > > arrDouble_;

  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<UChar_t  > > > > vecUChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<UShort_t > > > > vecUShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<UInt_t   > > > > vecUInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<ULong64_t> > > > vecULong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<bool     > > > > vecBool_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<Char_t   > > > > vecChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<short    > > > > vecShort_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<int      > > > > vecInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<Long64_t > > > > vecLong64_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<float    > > > > vecFloat_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<double   > > > > vecDouble_;

  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<std::vector<UChar_t> > > > > vecVUChar_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<std::vector<int    > > > > > vecVInt_;
  std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<std::vector<float  > > > > > vecVFloat_;
};


#endif // ifndef TreeReader_h
