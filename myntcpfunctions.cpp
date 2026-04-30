#include "myntcpfunctions.h"



void bookHisto(TFile *outrootfile, const vector<PatientData> &sample){
   TDirectory* dvhdir = outrootfile->GetDirectory("dvhplots");
  if (!dvhdir) dvhdir = outrootfile->mkdir("dvhplots");
  dvhdir->cd();

  return;
}

void PrintSampleLine(int idx, const vector<PatientData> &sample){
  cout<<"Print sample at index: "<<idx<<endl;  
  if(idx<0 || idx>=sample.size()){
    throw std::runtime_error("PrintSampleLine: wrong idx" + idx);
    return;
  }
  cout<<sample.at(idx).id<<" "<<sample.at(idx).volume<<" "<<sample.at(idx).max_dose_plan<<" "<<sample.at(idx).max_dose_st<<" "<<sample.at(idx).min_dose_st<<" "<<sample.at(idx).mean_dose_st<<endl;
  for(const double& value : sample.at(idx).dvhmap)
    cout<<value<<" ";
  cout<<endl; 
}

vector<string> splitCsvLine(const string& line) {
  vector<string> out;
  std::stringstream ss(line);
  string item;
  while (getline(ss, item, ' ')) 
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

 void loadDvhCsv(const string& filename, vector<PatientData> &sample) {
  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("loadDvhCsv: cannot open file " + filename);
  }
  string line;
  // ---------- prima riga: header ----------
  if (!std::getline(in, line))
    throw std::runtime_error("loadDvhCsv: empty file " + filename);

  // regex per id tipo RQ12345-6
  std::regex id_pattern(R"(ML-\d{4})");
  // ---------- righe dati ----------
  while (std::getline(in, line)) {
    PatientData patient;
    if (line.empty())
      continue;
    vector<string> parts = splitCsvLine(line);
    if (parts.size() < 8){
      // throw std::runtime_error("loadDvhCsv: check your dvh file, there is a line with less than 8 elements: " + filename);
      continue;
    }

    string readed = trim(parts[0]);
    if (!std::regex_search(readed, id_pattern)) {
      throw std::runtime_error("loadDvhCsv: check your dvh file, there is a non valid id: " +readed+ " " +filename);
      continue;
    }

    patient.id=std::stoi(readed.substr(readed.size() - 4));;
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
  sample.push_back(patient);
  }
  cout<<"loadDvhCsv done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;
  return;
}
