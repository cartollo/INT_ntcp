#include "myntcpfunctions.h"

extern int debug;

void bookHisto(TFile *outrootfile, const vector<double> &alfabeta, const vector<double> &nvalue4eud, double eqd2binwidth){

  if(debug)
    cout<<"start bookHisto"<<endl;

  TDirectory* dvhdir = outrootfile->GetDirectory("dvhplots");
  if (!dvhdir) dvhdir = outrootfile->mkdir("dvhplots");
  TH1D* h;
  TH2D* h2;

  //map properties
  h=new TH1D("sample_volumes", "sample volumes;Volume;patients [cm3];Number of patients;Number of patients",300, 0., 3000);
  h=new TH1D("max_dose_plan", "max dose plan;max dose plan [gy];Number of patients",100, 50., 100);
  h=new TH1D("max_dose_st", "max dose st;max dose st [gy];Number of patients",200, 0., 100);
  h=new TH1D("min_dose_st", "min dose st;min dose st [gy];Number of patients",20, 0., 10);
  h=new TH1D("mean_dose_st", "mean dose st;min dose st [gymax_dose_st];Number of patients",100, 0., 50);
  h=new TH1D("dose_ptv", "Dose to the PTV;ptv dose [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("nfraction", "Total number of fractions;Number of fractions;Number of patients",101, -0.5, 100.5);
  h=new TH1D("dose_per_fraction", "mean dose st;min dose st [gymax_dose_st];Number of patients",42, -0.1, 4.1);
  h=new TH1D("mean_dose_rectum", "mean dose st;min dose st [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("tgt_acutegitox", "tgt_acutegitox;tgt_acutegitox dose st [gymax_dose_st];Number of patients",2, -0.5, 1.5);
  
  //normalized dvh fofr targets:
  for(auto const &asub:alfabeta){
    h=new TH1D(Form("dvhnorm_tgt_%.3f_Yes",asub), "mean dvh plot for patients with tox;normalized dose [gy];mean volume",100/eqd2binwidth, 0., 100);
    h=new TH1D(Form("dvhnorm_tgt_%.3f_No",asub), "mean dvh plot for patients without tox;normalized dose [gy];mean volume",100/eqd2binwidth, 0., 100);
    h=new TH1D(Form("dvhnorm_tgt_%.3f_Diff",asub), "mean dvh plot for difference between patients w/o tox;normalized dose [gy];mean volume",100/eqd2binwidth, 0., 100);
    
    for(auto const &n:nvalue4eud)
      h=new TH1D(Form("eud_n_%.3f_asub%.3f",n,asub), "eud value;eud;Number of patients",100, 0., 100);    
  }
    
  if(debug)
    cout<<"bookHisto done"<<endl;

  return;
}


void fillHisto(TFile *outrootfile,  map<int, PatientData> &sample, const vector<double> &alfabeta, const vector<double> &nvalue4eud, double eqd2binwidth){

  if(debug)
    cout<<"start fillHisto"<<endl;

  gDirectory->cd("dvhplots");
  map<double,vector<double>> dvhnorm_tgtyes;
  map<double,vector<double>> dvhnorm_tgtno;
  int tgtyescounter=0, tgtnocounter=0;
  for(auto const &asub:alfabeta){
    dvhnorm_tgtyes[asub].resize(100/eqd2binwidth,0); //ensure that this length is the same of Form("dvhnorm_tgt_%.3f_Yes",asub)
    dvhnorm_tgtno[asub].resize(100/eqd2binwidth,0);
  }
  for(auto &paziente : sample){
    //dvh plots
    gDirectory->mkdir(Form("pz_%i",paziente.second.id));
    gDirectory->cd(Form("pz_%i",paziente.second.id));
    TH1D* hdvh=new TH1D(Form("pzdvh_%i",paziente.second.id), "DVH plot;Dose [gy];Volume",paziente.second.dvhmap.size(), 0., paziente.second.dvhmap.size());
    bool isboring=true; //si dice vingardium monotonic
    for(int i=0;i<paziente.second.dvhmap.size();i++){
      hdvh->SetBinContent(i+1,paziente.second.dvhmap.at(i));
      if(i>0){
        if(paziente.second.dvhmap.at(i)>paziente.second.dvhmap.at(i-1)){
          cout<<"WARNING in fillHisto: the dvh function of patient id="<<paziente.second.id<<"is not monotonous"<<endl;
          throw std::runtime_error("patient dhv not monotonous, patient id: " + paziente.second.id);         
          isboring=false;
          break;
        }
      }
    }
    if(!isboring){
      hdvh->SetName(Form("pzdvh_%i_notmonotonous",paziente.second.id));
      paziente.second.status=1;
    }    

    //scaled dose plots
    map <double, vector<TH1D*>> heqd2vec, hdvhnormvec;
    for(auto const &asub:alfabeta){
      TH1D* heqd2=new TH1D(Form("pzeqd2_%i_%.2f",paziente.second.id, asub), "eqd2 plot;bin;equivalent dose value [gy]",paziente.second.eqd2map.at(asub).size(), 0., paziente.second.eqd2map.at(asub).size());
      heqd2vec[asub].push_back(heqd2);
      TH1D* hdvhnorm=new TH1D(Form("pzdvhnorm_%i_%.2f",paziente.second.id, asub), "Normalized DVH plot;Equivalent Dose [gy];Volume",paziente.second.dvhnormmap.at(asub).size(), 0., paziente.second.dvhnormmap.at(asub).size()*eqd2binwidth);
      hdvhnormvec[asub].push_back(hdvhnorm);
      for(int i=0;i<paziente.second.eqd2map.at(asub).size();i++)
        heqd2->SetBinContent(i+1,paziente.second.eqd2map.at(asub).at(i));
      for(int i=0;i<paziente.second.dvhnormmap.at(asub).size();i++){
        hdvhnorm->SetBinContent(i+1,paziente.second.dvhnormmap.at(asub).at(i));
      }
    }
    
    gDirectory->cd(".."); 

    (dynamic_cast<TH1D*>(gDirectory->Get("../sample_volumes")))->Fill(paziente.second.volume);
    (dynamic_cast<TH1D*>(gDirectory->Get("../max_dose_plan")))->Fill(paziente.second.max_dose_plan);
    (dynamic_cast<TH1D*>(gDirectory->Get("../max_dose_st")))->Fill(paziente.second.max_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../min_dose_st")))->Fill(paziente.second.min_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../mean_dose_st")))->Fill(paziente.second.mean_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dose_ptv")))->Fill(paziente.second.dose_ptv);
    (dynamic_cast<TH1D*>(gDirectory->Get("../nfraction")))->Fill(paziente.second.nfraction);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dose_per_fraction")))->Fill(paziente.second.dose_per_fraction);
    (dynamic_cast<TH1D*>(gDirectory->Get("../mean_dose_rectum")))->Fill(paziente.second.mean_dose_rectum);
    (dynamic_cast<TH1D*>(gDirectory->Get("../tgt_acutegitox")))->Fill(paziente.second.tgt_acutegitox);
    
    if(paziente.second.tgt_acutegitox>0) tgtyescounter++; else tgtnocounter++;
    
    for(auto const &asub:alfabeta){
      for(int i=0;i<paziente.second.dvhnormmap.at(asub).size();i++){
        if(paziente.second.tgt_acutegitox>0)
        dvhnorm_tgtyes[asub].at(i)+=paziente.second.dvhnormmap.at(asub).at(i);
        else
        dvhnorm_tgtno[asub].at(i)+=paziente.second.dvhnormmap.at(asub).at(i);
      }
      
      for(auto const &n:nvalue4eud)
        (dynamic_cast<TH1D*>(gDirectory->Get(Form("../eud_n_%.3f_asub%.3f",n,asub))))->Fill(paziente.second.eudmap.at(make_pair(n,asub)));

    }
  }//end of loop on patient

  for(auto const &asub:alfabeta){
    for(int i=0;i<dvhnorm_tgtyes.at(asub).size();i++){
      // cout<<"suba="<<asub<<" i="<<i<<" dvhnorm_tgtyes.at(asub).at(i)/sample.size()="<<dvhnorm_tgtyes.at(asub).at(i)/sample.size()<<endl;
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhnorm_tgt_%.3f_Yes",asub))))->SetBinContent(i+1,dvhnorm_tgtyes.at(asub).at(i)/tgtyescounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhnorm_tgt_%.3f_No",asub))))->SetBinContent(i+1,dvhnorm_tgtno.at(asub).at(i)/tgtnocounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhnorm_tgt_%.3f_Diff",asub))))->SetBinContent(i+1,dvhnorm_tgtyes.at(asub).at(i)/tgtyescounter-dvhnorm_tgtno.at(asub).at(i)/tgtnocounter);
    }
  }

  gDirectory->cd(".."); 

  if(debug)
    cout<<"fillHisto done"<<endl;

  return;
}

double fitSigmoidal(TGraph* graph, int parnum, TFitResultPtr &fitresults){

  if(debug)
    cout<<"start fitSigmoidal"<<endl;

  if(parnum<2)
    throw std::runtime_error("fitSigmoidal: parnum too small" + parnum);
    
    TF1 *sigmoid;
  if(parnum==2)
    sigmoid=new TF1("sigmoid", "1./(1.+exp(-[0]-[1]*x))", 0, 100);
  else //per ora solo caso semplice con due parametri liberi ed una sola variabile X=EUD
    throw std::runtime_error("fitSigmoidal: parnum too large" + parnum);

  fitresults=graph->Fit(sigmoid, "EMQ+","",0,100);

  if(debug)
    cout<<"fitSigmoidal done"<<endl;

  // return fitresults.Chi2();
  return 2;
}


int CreateNtcpSigmoidal(map<int, PatientData> &sample, TString tgtname, const vector<double> &alfabeta, const vector<double> &nvalue4eud){

  gDirectory->mkdir("ntpc_linear");
  gDirectory->cd("ntpc_linear");

  map <string, TGraph*> ntcp_linear_models;
  TGraph *gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_volume_"+tgtname);
  ntcp_linear_models["volume"]=gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_max_dose_plan_"+tgtname);
  ntcp_linear_models["max_dose_plan"]=gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_max_dose_st_"+tgtname);
  ntcp_linear_models["max_dose_st"]=gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_min_dose_st_"+tgtname);
  ntcp_linear_models["min_dose_st"]=gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_mean_dose_st_"+tgtname);
  ntcp_linear_models["mean_dose_st"]=gr;
  
  for(auto const &paziente:sample){
    for(auto const &comb:paziente.second.eudmap){
      gr=new TGraph(sample.size());
      gr->SetName(Form("ntcp_linear_eud_n%f_alphabeta%.3f_",comb.first.first, comb.first.second)+tgtname);
      ntcp_linear_models[Form("ntcp_linear_eud_n%f_alphabeta%.3f_",comb.first.first, comb.first.second)]=gr;
    }
  }

  int index=0;
  for(auto &paziente : sample){
    ntcp_linear_models["volume"]->SetPoint(index, paziente.second.volume, paziente.second.tgt_acutegitox);
    ntcp_linear_models["max_dose_plan"]->SetPoint(index, paziente.second.max_dose_plan, paziente.second.tgt_acutegitox);
    ntcp_linear_models["max_dose_st"]->SetPoint(index, paziente.second.max_dose_st, paziente.second.tgt_acutegitox);
    ntcp_linear_models["min_dose_st"]->SetPoint(index, paziente.second.min_dose_st, paziente.second.tgt_acutegitox);
    ntcp_linear_models["mean_dose_st"]->SetPoint(index, paziente.second.mean_dose_st, paziente.second.tgt_acutegitox);
    for(auto const &comb:paziente.second.eudmap)
      ntcp_linear_models[Form("ntcp_linear_eud_n%f_alphabeta%.3f_",comb.first.first, comb.first.second)]->SetPoint(index, comb.second,paziente.second.tgt_acutegitox);
    index++;
  }

  for(auto &graph:ntcp_linear_models){
    graph.second->Sort();
    TFitResultPtr fitresults;
    fitSigmoidal(graph.second,2,fitresults);
    graph.second->SetMarkerStyle(20);
    graph.second->SetMarkerColor(2);
    graph.second->SetLineWidth(0);
    graph.second->SetLineColor(0);
    graph.second->SetDrawOption("AP*");
    graph.second->Write();
  }

  gDirectory->cd("..");

  return 0;
}

void PrintSampleLine(int idx, const   map<int, PatientData> &sample){
  cout<<"Print sample at index: "<<idx<<endl;  
  if(sample.count(idx)==0 ){
    throw std::runtime_error("PrintSampleLine: idx not present in sample" + idx);
    return;
  }
  cout<<sample.at(idx).id<<" "<<sample.at(idx).volume<<" "<<sample.at(idx).max_dose_plan<<" "<<sample.at(idx).max_dose_st<<" "<<sample.at(idx).min_dose_st<<" "<<sample.at(idx).mean_dose_st<<" "<<sample.at(idx).dose_ptv<<" "<<sample.at(idx).nfraction<<" "<<sample.at(idx).dose_per_fraction<<" "<<sample.at(idx).mean_dose_rectum<<endl;
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

 int loadDvhFile(const string& filename,   map<int, PatientData> &sample, const vector<double> &alfabeta) {
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
    patient.id=id;
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
  patient.status=-1;
  for(auto const asub:alfabeta){
    vector<double> tmpvec(patient.dvhmap.size(),0);
    patient.eqd2map[asub]=tmpvec;
  }
  sample.insert(std::pair{id,patient});
  }

  if(debug)
    cout<<"loadDvhCsv done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;
  
    return 0;
}

double CalculateEudFromScratch(PatientData &paziente, double eqd2binwidth, double alfabeta, double nvalue){

  if(debug)
    cout<<"start CalculateEudFromScratch with paziente.id="<<paziente.id<<" eqd2binwidth="<<eqd2binwidth<<endl;

    vector<double> eqd2map;
    for(int dose=0;dose<paziente.dvhmap.size();dose++){
      if(paziente.dvhmap.at(dose)!=0){
        eqd2map.push_back(((double)dose)*(alfabeta+((double)dose)/paziente.nfraction)/(alfabeta+2.));
      }else
        break;
    }     

    int eqd2index=0;
    double eud=0;
    vector<double> dvhnormmap(eqd2map.size()/eqd2binwidth+0.5,0);
    for(int i=0;i<dvhnormmap.size();i++){
      double dose=eqd2binwidth*(i+1);
      if(eqd2map.back()<dose){
        dvhnormmap.at(i)=0;
        break;
      }
      for(int k=eqd2index;k<eqd2map.size();k++){
        if(eqd2map.at(k)>dose){
          break;
        }
        eqd2index=k;
      }
      dvhnormmap.at(i)=paziente.dvhmap.at(eqd2index);     
      eud+=pow(dose,nvalue)*dvhnormmap.at(i)/dvhnormmap.at(0);
    }
    eud=pow(eud,1./nvalue);
    
  if(debug)
    cout<<"CalculateEudFromScratch done; eud="<<eud<<endl;

    return eud;
}

int evaluateEqdEud(map<int, PatientData> &sample, const vector<double> &alfabeta, const vector<double> &nvalue4eud, double eqd2binwidth){
  
  if(debug)
  cout<<"start evaluateEqdEud"<<endl;
  
  for(auto &paziente : sample){
    for(auto const &asub:alfabeta){
      //fill eqd2map
      for(int dose=0;dose<paziente.second.eqd2map.at(asub).size();dose++){
        if(paziente.second.dvhmap.at(dose)==0){ //non dovrebbe succedere
          throw std::runtime_error(Form("Error in evaluateEqdEud: dvhmap contains a 0 value dose=%i",dose));
          paziente.second.eqd2map.at(asub).resize(dose+1);
          paziente.second.eqd2map.at(asub).shrink_to_fit();
          break;
        }
        paziente.second.eqd2map.at(asub).at(dose)=((double)dose)*(asub+((double)dose)/paziente.second.nfraction)/(asub+2.);
      }
  
      //calculate dvhnormmap
      vector<double> tmpmap(paziente.second.eqd2map.at(asub).size()/eqd2binwidth+0.5,0);
      paziente.second.dvhnormmap[asub]=tmpmap;
      int eqd2index=0;
      for(int i=0;i<paziente.second.dvhnormmap[asub].size();i++){ 
        double dose=eqd2binwidth*(i+1);
        if(paziente.second.eqd2map.at(asub).back()<dose){
          paziente.second.dvhnormmap[asub].at(i)=0;
          break;
        }
        for(int k=eqd2index;k<paziente.second.eqd2map.at(asub).size();k++){
          if(paziente.second.eqd2map.at(asub).at(k)>dose){
            break;
          }
          eqd2index=k;
        }
        paziente.second.dvhnormmap[asub].at(i)=paziente.second.dvhmap.at(eqd2index);
      
        for(auto const &n:nvalue4eud){
          paziente.second.eudmap[make_pair(n,asub)]+=pow(dose,n)*paziente.second.dvhnormmap[asub].at(i)/paziente.second.dvhnormmap[asub].at(0);
        }

      }
      
      for(auto const &n:nvalue4eud){
        paziente.second.eudmap.at(make_pair(n,asub))=pow(paziente.second.eudmap.at(make_pair(n,asub)),1./n);
        // //check consistency with function (it works!)
        // double tmpeud=CalculateEudFromScratch(paziente.second, eqd2binwidth, asub, n);
        // if((fabs(tmpeud)-paziente.second.eudmap.at(make_pair(n,asub)))>0.1)
        //   cout<<"inconsistency: tmpeud="<<tmpeud<<"  paziente.second.eudmap.at(make_pair(n,asub)))="<<paziente.second.eudmap.at(make_pair(n,asub))<<endl;
      }

    } //end of loop on alfabeta
  }//end of loop on patients
  

  if(debug)
    cout<<"evaluateEqdEud done"<<endl;
  
  return 0;
}

int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname){

  if(debug)
    cout<<"start loadMetaFile"<<endl;

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
    
    int tgtpos=-1, doseptv_pos=-1, nfrac_pos=-1, meandoserectum_pos=-1, doseperfraction_pos=-1;
    TString doseptv_name("dose to the prostate PTV (Gy)");
    TString nfrac_name("fraction number");
    TString meandoserectum_name("mean dose to the rectum (Gy)");
    TString doseperfraction_name("dose/fraction");
    for(int i=0;i<parts.size();i++){
      if(tgtname.CompareTo(parts[i])==0)
        tgtpos=i;
      if(doseptv_name.CompareTo(parts[i])==0)
        doseptv_pos=i;
      if(nfrac_name.CompareTo(parts[i])==0)
        nfrac_pos=i;
      if(meandoserectum_name.CompareTo(parts[i])==0)
        meandoserectum_pos=i;
      if(doseperfraction_name.CompareTo(parts[i])==0)
        doseperfraction_pos=i;
    }
    if(tgtpos==-1 || doseptv_pos==-1 || nfrac_pos==-1 || meandoserectum_pos==-1 || doseperfraction_pos==-1) {
      throw std::runtime_error(Form("loadMetaFile: error, some variable not found: tgtpos=%i  doseptv_pos=%i  nfrac_pos=%i  meandoserectum_pos=%i  doseperfraction_pos=%i", tgtpos, doseptv_pos,  nfrac_pos,  meandoserectum_pos, doseperfraction_pos));
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
    if (parts.size() < tgtpos || parts.size() < doseptv_pos|| parts.size() < nfrac_pos|| parts.size() < meandoserectum_pos|| parts.size() < doseperfraction_pos)
      continue;

    string readed = trim(parts[0]);
    if (!std::regex_search(readed, id_pattern)) {
      throw std::runtime_error("loadDvhCsv: check your dvh file, there is a non valid id: " +readed+ " " +filename);
      continue;
    }
    auto it=sample.find(std::stoi(readed.substr(readed.size() - 4)));//trovo il paziente in mappa
    if(it!=sample.end()){
      it->second.tgt_acutegitox=std::stoi(trim(parts[tgtpos]));
      it->second.dose_ptv=std::stod(trim(parts[doseptv_pos]));
      it->second.nfraction=std::stoi(trim(parts[nfrac_pos]));
      it->second.dose_per_fraction=std::stod(trim(parts[doseperfraction_pos]));
      it->second.mean_dose_rectum=std::stod(trim(parts[meandoserectum_pos]));
      addedinmap++;
    }
  }
  
  if(addedinmap!=sample.size())
    cout<<"WARNING: mismatch between addedinmap="<<addedinmap<<"  and  sample.size()="<<sample.size()<<endl;
  if(debug)
    cout<<"loadMetaFile done, target aquired target="<<tgtname<<"  overall sample size="<<sample.size()<<endl;

  return 0;
}
