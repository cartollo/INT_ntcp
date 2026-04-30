#include "myntcpfunctions.h"



void bookHisto(TFile *outrootfile, const   map<int, PatientData> &sample){
   TDirectory* dvhdir = outrootfile->GetDirectory("dvhplots");
  if (!dvhdir) dvhdir = outrootfile->mkdir("dvhplots");
  dvhdir->cd();

  return;
}

void PrintSampleLine(int idx, const   map<int, PatientData> &sample){
  cout<<"Print sample at index: "<<idx<<endl;  
  if(sample.count(idx)==0 ){
    throw std::runtime_error("PrintSampleLine: idx not present in sample" + idx);
    return;
  }
  cout<<sample.at(idx).id<<" "<<sample.at(idx).volume<<" "<<sample.at(idx).max_dose_plan<<" "<<sample.at(idx).max_dose_st<<" "<<sample.at(idx).min_dose_st<<" "<<sample.at(idx).mean_dose_st<<endl;
  for(const double& value : sample.at(idx).dvhmap)
    cout<<value<<" ";
  cout<<endl; 
}

vector<string> splitCsvLine(const string& line, const TString delimiter) {
  vector<string> out;
  std::stringstream ss(line);
  string item;
  char dummy=delimiter[0];
  while (getline(ss, item, dummy)) 
    out.push_back(item);
  return out;
}

string trim(const string& s) {
  size_t start = 0;
  while (start < s.size() && isspace(static_cast<unsigned char>(s[start])))
    ++start;
  size_t end = s.size();
  while (end > start && isspace(static_cast<unsigned char>(s[end - 1])))
    --end;
  return s.substr(start, end - start);
}

 int loadDvhFile(const string& filename,   map<int, PatientData> &sample) {
  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("loadDvhCsv: cannot open file " + filename);
  }
  string line;
  // ---------- prima riga: header ----------
  if (!std::getline(in, line)){
    throw std::runtime_error("loadDvhCsv: empty file " + filename);
    return 1;
  }

  // regex per id tipo RQ12345-6
  std::regex id_pattern(R"(ML-\d{4})");
  // ---------- righe dati ----------
  while (std::getline(in, line)) {
    PatientData patient;
    if (line.empty())
      continue;
    TString delimiter(" ");
    vector<string> parts = splitCsvLine(line, delimiter);
    if (parts.size() < 8){
      // throw std::runtime_error("loadDvhCsv: check your dvh file, there is a line with less than 8 elements: " + filename);
      continue;
    }

    string readed = trim(parts[0]);
    if (!std::regex_search(readed, id_pattern)) {
      throw std::runtime_error("loadDvhCsv: check your dvh file, there is a non valid id: " +readed+ " " +filename);
      continue;
    }

    int id=std::stoi(readed.substr(readed.size() - 4));
    readed = trim(parts[2]);
    patient.volume=std::stod(readed);
    readed = trim(parts[3]);
    patient.max_dose_plan=std::stod(readed);
    readed = trim(parts[4]);
    patient.max_dose_st=std::stod(readed);
    readed = trim(parts[5]);
    patient.min_dose_st=std::stod(readed);
    readed = trim(parts[6]);
    patient.mean_dose_st=std::stod(readed);
    for (size_t i = 7; i < parts.size(); ++i) {
      if(parts[i].size()==0)
        break;
      readed = trim(parts[i]);
      patient.dvhmap.push_back(std::stod(readed));
    }
  sample.insert(std::pair{id,patient});
  }
  cout<<"loadDvhCsv done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;
  return 0;
}

int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname){

  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("loadMetaFile: cannot open file " + filename);
  }
  string line;
  // ---------- prima riga: header ----------
  if (!std::getline(in, line))
    throw std::runtime_error("loadMetaFile: empty file " + filename);
    TString delimiter(",");
    vector<string> parts = splitCsvLine(line, delimiter);
    
    int tgtpos=-1;
    for(int i=0;i<parts.size();i++){
      if(tgtname.CompareTo(parts[i])==0){
        tgtpos=i;
        break;
      }
    }
    if(tgtpos==-1) {
      throw std::runtime_error("loadMetaFile: target not found " + tgtname);
      return 1;
    }


  // regex per id tipo RQ12345-6
  std::regex id_pattern(R"(ML-\d{4})");
  // ---------- righe dati ----------
  int addedinmap=0;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    TString delimiter(",");
    vector<string> parts = splitCsvLine(line, delimiter);
    if (parts.size() < tgtpos)
      continue;

    string readed = trim(parts[0]);
    if (!std::regex_search(readed, id_pattern)) {
      throw std::runtime_error("loadDvhCsv: check your dvh file, there is a non valid id: " +readed+ " " +filename);
      continue;
    }
    auto it=sample.find(std::stoi(readed.substr(readed.size() - 4)));
    if(it!=sample.end()){
      it->second.tgt_acutegitox=std::stoi(trim(parts[tgtpos]));
      addedinmap++;
    }
  }
  if(addedinmap!=sample.size())
    cout<<"WARNING: mismatch between addedinmap="<<addedinmap<<"  and  sample.size()="<<sample.size()<<endl;
  else
    cout<<"loadMetaFile done, target aquired target="<<tgtname<<"  overall sample size="<<sample.size()<<endl;

  return 0;
}
