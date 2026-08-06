#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMultiGraph.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define debug 0

struct ParameterResult
{
  double value = 0.0;
  double error = 0.0;
};

struct BootstrapResult
{
  int seed = -1;
  int status = -999;
  int cov_matrix_status = -999;
  double loglikelihood_min_value = 0.0;
  map<string, ParameterResult> parameters;
  map<int, vector<double>> ptidscore;
  string source_file;
};

enum class OobpleParseResult {
    NotMatched,
    Parsed,
    Error
};


//toglie gli spazi bianchi prima e dopo
string Trim(const string &text){
  const size_t first = text.find_first_not_of(" \t\r\n");
  if (first == string::npos)
    return "";
  const size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

bool ExtractProp2Dose(const std::string& line, int& prop2dose) {
    static const std::regex pattern(R"(_prop2dose:([01])_)");
    std::smatch match;
    if (!std::regex_search(line, match, pattern))
        return false;
    prop2dose = std::stoi(match[1].str());
    cout<<"found prop2dose="<<prop2dose<<endl;
    return true;
}

void ExtractUsedDoseVar(const std::string& line, int& useddosevar) {
    static const std::regex pattern(R"(_selectededoses4volume:([0-9]+)_)");
    std::smatch match;
    if (!std::regex_search(line, match, pattern)){
      useddosevar=-1;
      cout<<"_selectededoses4volume not found, useddosevar set to eud, useddosevar="<<useddosevar<<endl;
      return;
    }
    useddosevar = std::stoi(match[1].str());
    cout<<"found _selectededoses4volume; useddosevar="<<useddosevar<<endl;
    return ;
}

bool ExtractClinicalFactors(const std::string& line, int& clinicalfactors) {
    static const std::regex pattern(R"(_clinicalfactors:([0-9]+)_)");
    std::smatch match;
    if (!std::regex_search(line, match, pattern)){
      return false;
    }
    clinicalfactors = std::stoi(match[1].str());
    cout<<"found clinicalfactors="<<clinicalfactors<<endl;
    return true;
}

bool ParseSeed(const string &line, int &seed){
  static const regex pattern(R"(^\s*seed\s*=\s*(-?\d+)\s*$)");
  smatch match;
  if (!regex_match(line, match, pattern))
    return false;
  seed = stoi(match[1].str());
  return true;
}

bool ParseStatusLine(const string &line,int &status,int &cov_matrix_status,double &loglikelihood_min_value){
  static const regex pattern(
      R"(^\s*status\s*=\s*(-?\d+)\s+CovMatrixStatus\s*=\s*(-?\d+)\s+loglikehood_min_value\s*=\s*([-+0-9.eE]+)\s*$)");
  smatch match;
  if (!regex_match(line, match, pattern))
    return false;
  status = stoi(match[1].str());
  cov_matrix_status = stoi(match[2].str());
  loglikelihood_min_value = stod(match[3].str());
  return true;
}

bool ParseParameterLine(const string &line,string &parameter_name,double &value,double &error){
  static const regex pattern(R"(^\s*([^:]+?)\s*:\s*([-+0-9.eE]+)\s*\+\-\s*([-+0-9.eE]+)\s*$)");
  smatch match;
  if (!regex_match(line, match, pattern))
    return false;
  parameter_name = Trim(match[1].str());
  value = stod(match[2].str());
  error = stod(match[3].str());
  return true;
}

OobpleParseResult ParseOobpleLine(const std::string &line, std::istream &input, int &oobplesize, vector<pair<int, double>> &id_prediction){
  static const std::regex size_pattern(R"(^\s*oobplesize\s*=\s*([0-9]+)\s*$)");
  static const std::regex prediction_pattern(R"(^\s*patient_id\s*:\s*(-?[0-9]+)\s+prediction\s*:\s*([-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\s*$)");
  std::smatch match;
  // La riga corrente non è l'inizio di un blocco OOB.
  if (!std::regex_match(line, match, size_pattern))
    return OobpleParseResult::NotMatched;
  try{
    oobplesize = std::stoi(match[1].str());
  }catch(const std::exception &e){
    std::cerr<<"ERROR: invalid oobplesize in line: "<<line<<endl;
    return OobpleParseResult::Error;
  }

  id_prediction.clear();
  id_prediction.reserve(static_cast<std::size_t>(oobplesize));
  for (int i=0;i<oobplesize;++i){
    std::string prediction_line;
    if (!std::getline(input, prediction_line)){
      std::cerr<<"ERROR: unexpected end of file while reading OOB predictions. Expected"<<oobplesize<<" predictions, found only "<<i<<endl;
      return OobpleParseResult::Error;
    }
    prediction_line = Trim(prediction_line);
    std::smatch prediction_match;
    if (!std::regex_match(prediction_line,prediction_match,prediction_pattern)){
      std::cerr<<"ERROR: invalid OOB prediction line:"<<endl<<prediction_line<<endl<<"Expected format:"<<endl<<"patient_id: <int>  prediction: <double>"<<endl;
      return OobpleParseResult::Error;
    }

    try{
      const int patient_id=std::stoi(prediction_match[1].str());
      const double prediction=std::stod(prediction_match[2].str());
      id_prediction.emplace_back(patient_id, prediction);
    }catch(const std::exception &e){
      std::cerr<<"ERROR: invalid numeric value in OOB prediction line:"<<endl<<prediction_line <<endl;
      return OobpleParseResult::Error;
    }
  }
  return OobpleParseResult::Parsed;
}

bool ReadBootstrapFile(const string &filename,string &reference_header,bool &reference_header_initialized,vector<BootstrapResult> &results, vector<string> &reference_parameter_names,bool &reference_parameters_initialized, int &prop2dose, int &useddosevar, int &clinicalfactors){
  ifstream input(filename);

  if (!input.is_open()){
    cerr << "ERROR: cannot open " << filename << '\n';
    return false;
  }

  string line;
  if (!getline(input, line)){
    cerr << "ERROR: empty file " << filename << '\n';
    return false;
  }

  const string current_header = line;

  if (!reference_header_initialized){
    reference_header = current_header;
    ExtractProp2Dose(line, prop2dose);
    ExtractUsedDoseVar(line, useddosevar);
    ExtractClinicalFactors(line, clinicalfactors);
    reference_header_initialized = true;
  } else if (current_header != reference_header){
    cout<<"ERROR: first line differs in file:\n"<<"  "<<filename<<"Reference line:"<<endl;
    cout<<reference_header<<endl;
    cout<<"Current line:"<<endl;
    cout<< current_header<<endl;
    cout<<"this file will be skipped, but please, take a look at it"<<endl;
    return false;
  }
  
  BootstrapResult current_result;
  bool inside_block = false;
  bool status_found = false;

  //definizione funzione lambda per leggere i blocchi!
  auto close_current_block = [&]() -> bool{
    if (!inside_block)
      return true;

    if (!status_found){
      cerr<<"WARNING: missing status line for seed "<<current_result.seed<<" in "<<filename<<endl;
      inside_block = false;
      return true;
    }
    
    vector<string> current_parameter_names;

    for (const auto &[name, result] : current_result.parameters)
      current_parameter_names.push_back(name);

    if (!reference_parameters_initialized){
      reference_parameter_names = current_parameter_names;
      reference_parameters_initialized = true;
    }else if (current_parameter_names != reference_parameter_names){
      cerr << "ERROR: variable list differs for seed "<<current_result.seed<<" in "<<filename<<endl;
      cerr << "Expected variables:"<<endl;
      for (const auto &name : reference_parameter_names)
        cerr<<"  "<<name<<endl;
      cerr << "Found variables:"<<endl;
      for (const auto &name : current_parameter_names)
        cerr<<"  "<<name<<endl;
      inside_block = false;
      return true;
    }
    results.push_back(current_result);
    return true;
  };

  while (getline(input, line)){
    line = Trim(line);
    if (line.empty())
      continue;
    int seed = -1;

    if (ParseSeed(line, seed)){
      if (!close_current_block())
        return false;
      current_result = BootstrapResult();
      current_result.seed = seed;
      current_result.source_file = filename;
      inside_block = true;
      status_found = false;
      continue;
    }

    if (!inside_block){
      cerr<<"WARNING: ignoring line outside a block in "<<filename<< ": "<<line<<endl;
      continue;
    }

    int status = 0;
    int cov_matrix_status = 0;
    double loglikelihood = 0.0;
    if (ParseStatusLine(line,status,cov_matrix_status,loglikelihood)){
      current_result.status = status;
      current_result.cov_matrix_status = cov_matrix_status;
      current_result.loglikelihood_min_value = loglikelihood;
      status_found = true;
      continue;
    }
    string parameter_name;
    double value = 0.0;
    double error = 0.0;
    if (ParseParameterLine(line, parameter_name, value, error)){
      if (current_result.parameters.find(parameter_name) != current_result.parameters.end()){
        cerr<<"ERROR: duplicated parameter "<<parameter_name<<" for seed "<<current_result.seed<<" in "<<filename<<endl;
        return false;
      }
      current_result.parameters[parameter_name] = {value, error};
      continue;
    }
    int oobplesize=0;
    vector<pair<int, double>> id_prediction;
    const OobpleParseResult oobple_result=ParseOobpleLine(line,input,oobplesize,id_prediction);
    if (oobple_result == OobpleParseResult::Error)
      return false;
    if(oobple_result==OobpleParseResult::Parsed){
      for(const auto& [patient_id, prediction] : id_prediction) {
        current_result.ptidscore[patient_id].push_back(prediction);
      }
      continue;
    }

    cerr<<"WARNING: unrecognized line in "<<filename<<": "<<line<<endl;
  }

  if (!close_current_block())
    return false;

  return true;
}

void PlotBootstrapResults(const string directory_name=".",const string prefix = "bootstrapout",const string root_filename = "bootstrapresults.root"){
  namespace fs = filesystem;
  int prop2dose, useddosevar, clinicalfactors;

  if (!fs::exists(directory_name)){
    cerr << "ERROR: directory does not exist: "<< directory_name << '\n';
    return;
  }

  const regex filename_pattern("^" + prefix + R"(_([0-9]+)\.txt$)");
  vector<pair<int, string>> input_files;

  if(debug)
    cout<<"start loop on folder"<<endl;

  //loop su cartella per caricare tutti gli input file validi
  for (const auto &entry : fs::directory_iterator(directory_name)){
    if (!entry.is_regular_file())
      continue;
    const string basename = entry.path().filename().string();
    smatch match;
    if (!regex_match(basename, match, filename_pattern))
      continue;
    const int file_number = stoi(match[1].str());
    input_files.push_back({file_number,entry.path().string()});
  }

  if (input_files.empty()){
    cerr<<"ERROR: no files matching "<<prefix<<"_numero.txt found in "<<directory_name<<endl;
    return;
  }

  //lambda function
  sort(input_files.begin(),input_files.end(),[](const auto &first, const auto &second){return first.first < second.first;});

  cout<<"Trovati"<< input_files.size()<< "input files:"<<endl;
  for (const auto &[number, filename] : input_files)
    cout<<"  "<<filename<<endl;

  vector<BootstrapResult> results;
  string reference_header;
  bool reference_header_initialized = false;

  vector<string> parameter_names;
  bool reference_parameters_initialized = false;

  for (const auto &[number, filename] : input_files){
    if (!ReadBootstrapFile(filename,reference_header,reference_header_initialized,results,parameter_names,reference_parameters_initialized, prop2dose, useddosevar, clinicalfactors)){
      cerr << "ERROR while reading " << filename << '\n';
      return;
    }
  }

  if (results.empty()){
    cerr << "ERROR: no bootstrap blocks found.\n";
    return;
  }

  sort(results.begin(),results.end(),[](const BootstrapResult &first,const BootstrapResult &second){return first.seed<second.seed;});

  //faccio checks
  set<int> seeds;
  for (const BootstrapResult &result : results){
    if (!seeds.insert(result.seed).second){
      cerr<<"ERROR: duplicated seed "<<result.seed<<" found in file "<<result.source_file<<endl;
      return;
    }
  }

  cout<<"Header check passed"<<endl;
  cout<<"Total cases="<<results.size()<<endl;
  cout<<"Number of variables: "<<parameter_names.size()<<endl;

  for (const string &name : parameter_names)
    cout<<name<<endl;

  TFile output_file(root_filename.c_str(), "RECREATE");

  if (output_file.IsZombie()){
    cerr<<"ERROR: cannot create "<<root_filename<<endl;
    return;
  }

  //let's do parameters plots  
  for (const string &parameter_name : parameter_names){
    TMultiGraph *multigr=new TMultiGraph(("multi_"+parameter_name).c_str(), Form("%s;Seed;%s",parameter_name.c_str(),parameter_name.c_str()));
    TGraphErrors *valid_graph=new TGraphErrors();
    TGraphErrors *invalid_graph=new TGraphErrors();
    valid_graph->SetName(("graph_" + parameter_name + "_valid").c_str());
    invalid_graph->SetName(("graph_" + parameter_name + "_invalid").c_str());
    int valid_point=0, invalid_point=0;
    for (const BootstrapResult &result : results){
      const ParameterResult &parameter=result.parameters.at(parameter_name);
      const bool valid_fit=(result.status == 0 && result.cov_matrix_status == 3);
      if (valid_fit){
        valid_graph->SetPoint(valid_point,result.seed,parameter.value);
        valid_graph->SetPointError(valid_point,0.0,parameter.error);
        ++valid_point;
      }else{
        invalid_graph->SetPoint(invalid_point,result.seed,parameter.value);
        invalid_graph->SetPointError(invalid_point,0.0,parameter.error);
        ++invalid_point;
      }
    }
    valid_graph->SetMarkerStyle(20);
    valid_graph->SetMarkerSize(0.9);
    valid_graph->SetMarkerColor(kBlue + 1);
    valid_graph->SetLineColor(kBlue + 1);
    invalid_graph->SetMarkerStyle(24);
    invalid_graph->SetMarkerSize(1.0);
    invalid_graph->SetMarkerColor(kRed + 1);
    invalid_graph->SetLineColor(kRed + 1);

    multigr->Add(valid_graph);
    multigr->Add(invalid_graph);
    multigr->Write();
  }

  //other plots
  TMultiGraph *likeminimum=new TMultiGraph("multi_likehood_minimum", "seeds;loglikehood minumim value");
  TGraphErrors *lkmin_val=new TGraphErrors();
  lkmin_val->SetMarkerStyle(20);
  lkmin_val->SetMarkerColor(kBlue + 1);
  lkmin_val->SetLineColor(kBlue + 1);
  
  TGraphErrors *lkmin_inv=new TGraphErrors();
  lkmin_inv->SetMarkerStyle(24);
  lkmin_inv->SetMarkerColor(kRed + 1);
  lkmin_inv->SetLineColor(kRed + 1);

  TH1D* h;
  h=new TH1D("status", "status should be 0;status;counts",12,-1.5,10.5);  
  h=new TH1D("CovMatrixStatus", "CovMatrixStatus shuld be 3;CovMatrixStatus;counts",12,-1.5,10.5);  
  int valid_point = 0;
  int invalid_point = 0;
  double maxxvalue= (useddosevar==-1) ? 70 : 0.5;
  TH2D* noclfactors_th2d=new TH2D("noclfactors",Form("noclinical ntcp curve heatmap;%s;counts", useddosevar == -1 ? "EUD" : "Volume"),100,0,maxxvalue, 200, 0., 0.);
  TH2D* clfactor0_th2d=new TH2D("clfactor_0",Form("noclinical ntcp curve heatmap;%s;counts", useddosevar == -1 ? "EUD" : "Volume"),100,0,maxxvalue, 200, 0., 0.);
  TH2D* clfactor1_th2d=new TH2D("clfactor_1",Form("noclinical ntcp curve heatmap;%s;counts", useddosevar == -1 ? "EUD" : "Volume"),100,0,maxxvalue, 200, 0., 0.);
  TF1 noclinicaltf1("noclinicaltf1","1./(1.+exp(-[0]-[1]*x))",0.,maxxvalue);
  TF1 clfactor_0tf1("clfactor_0tf1",(prop2dose==1) ? "1./(1.+exp(-[0]-([1]+[2])*x))" : "1./(1.+exp(-[0]-[1]*x-[2]))",0.,maxxvalue);
  TF1 clfactor_1tf1("clfactor_1tf1",(prop2dose==1) ? "1./(1.+exp(-[0]-([1]+[2])*x))" : "1./(1.+exp(-[0]-[1]*x-[2]))",0.,maxxvalue);
  TDirectory *pzscoredir =output_file.mkdir("patients_scores");
  //loop to fill stuff
  for (const BootstrapResult &result : results){
    if ((result.status==0 && result.cov_matrix_status==3)){
      lkmin_val->SetPoint(valid_point,result.seed,result.loglikelihood_min_value);
      ++valid_point;
      noclinicaltf1.SetParameters(result.parameters.at("beta_zero").value, result.parameters.at("beta_eud_a").value);
      if(clinicalfactors>0)
      clfactor_0tf1.SetParameters(result.parameters.at("beta_zero").value, result.parameters.at("beta_eud_a").value, result.parameters.at("clinical_factor_0").value);
      if(clinicalfactors>1)
      clfactor_1tf1.SetParameters(result.parameters.at("beta_zero").value, result.parameters.at("beta_eud_a").value, result.parameters.at("clinical_factor_1").value);
      for(int i=1;i<=noclfactors_th2d->GetXaxis()->GetNbins();i++){
        double xval=noclfactors_th2d->GetXaxis()->GetBinCenter(i);
        noclfactors_th2d->Fill(xval,noclinicaltf1.Eval(xval));
        if(clinicalfactors>0)
        clfactor0_th2d->Fill(xval,clfactor_0tf1.Eval(xval));
        if(clinicalfactors>1)
        clfactor1_th2d->Fill(xval,clfactor_1tf1.Eval(xval));
      }
      //patients scores
      TH1D *h = nullptr;
      for(const auto &pzt:result.ptidscore){
        pzscoredir->GetObject(Form("pzscore_%i",pzt.first),h);
        if(!h){
          h=new TH1D(Form("pzscore_%i",pzt.first), Form("patient %i OOB NTCP score;bootstrap repetitions;NTCP score",pzt.first), 120,-0.1,1.1);
          h->SetDirectory(pzscoredir);
        }
        for(const auto &score:pzt.second)
          h->Fill(score);
      }
    }else{
      lkmin_inv->SetPoint(invalid_point,result.seed,result.loglikelihood_min_value);
      ++invalid_point;
    }
    (dynamic_cast<TH1D*>(gDirectory->Get("status")))->Fill(result.status);    
    (dynamic_cast<TH1D*>(gDirectory->Get("CovMatrixStatus")))->Fill(result.cov_matrix_status);    
  }

  //OOB calibration curve
    std::map<int, int> pzid_tox;
    std::ifstream pzidtoxfile("pazient_tox.txt");
    int key, value;
    while(pzidtoxfile>>key>>value)
      pzid_tox[key]=value;
    pzidtoxfile.close();

  //write stuff
  likeminimum->Add(lkmin_val);
  likeminimum->Add(lkmin_inv);
  likeminimum->Write();
  output_file.Write();
  cout<<"program ended, Created: "<<root_filename<<endl;

  return;

}