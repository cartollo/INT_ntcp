#include "myntcpfunctions.h"

extern int debug;
void bookHisto(TFile *outrootfile, const map<int, PatientData> &sample, const globalstuff &glbstuff, const bool multivariate){

  if(debug)
    cout<<"start bookHisto"<<endl;

  TDirectory* dvhdir = outrootfile->GetDirectory("dvhplots");
  if (!dvhdir) dvhdir = outrootfile->mkdir("dvhplots");
  TH1D* h;
  TH2D* h2;

  //map properties 
  h=new TH1D("sample_volumesAll", "sample volumes;Volume;Number of patients",300, 0., 3000);
  h=new TH1D("sample_volumesYes", "sample volumes;Volume;Number of patients",300, 0., 3000);
  h=new TH1D("sample_volumesNo", "sample volumes;Volume;Number of patients",300, 0., 3000);
  h=new TH1D("max_dose_planAll", "max dose plan;max dose plan [gy];Number of patients",100, 50., 100);
  h=new TH1D("max_dose_planYes", "max dose plan;max dose plan [gy];Number of patients",100, 50., 100);
  h=new TH1D("max_dose_planNo", "max dose plan;max dose plan [gy];Number of patients",100, 50., 100);
  h=new TH1D("max_dose_stAll", "max dose st;max dose st [gy];Number of patients",200, 0., 100);
  h=new TH1D("max_dose_stYes", "max dose st;max dose st [gy];Number of patients",200, 0., 100);
  h=new TH1D("max_dose_stNo", "max dose st;max dose st [gy];Number of patients",200, 0., 100);
  h=new TH1D("min_dose_stAll", "min dose st;min dose st [gy];Number of patients",20, 0., 10);
  h=new TH1D("min_dose_stYes", "min dose st;min dose st [gy];Number of patients",20, 0., 10);
  h=new TH1D("min_dose_stNo", "min dose st;min dose st [gy];Number of patients",20, 0., 10);
  h=new TH1D("mean_dose_stAll", "mean dose st;min dose st [gymax_dose_st];Number of patients",100, 0., 50);
  h=new TH1D("mean_dose_stYes", "mean dose st;min dose st [gymax_dose_st];Number of patients",100, 0., 50);
  h=new TH1D("mean_dose_stNo", "mean dose st;min dose st [gymax_dose_st];Number of patients",100, 0., 50);
  h=new TH1D("dose_ptvAll", "Dose to the PTV;ptv dose [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("dose_ptvYes", "Dose to the PTV;ptv dose [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("dose_ptvNo", "Dose to the PTV;ptv dose [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("nfractionAll", "Total number of fractions;Number of fractions;Number of patients",101, -0.5, 100.5);
  h=new TH1D("nfractionYes", "Total number of fractions;Number of fractions;Number of patients",101, -0.5, 100.5);
  h=new TH1D("nfractionNo", "Total number of fractions;Number of fractions;Number of patients",101, -0.5, 100.5);
  h=new TH1D("dose_per_fractionAll", "mean dose st;min dose st [gymax_dose_st];Number of patients",42, -0.1, 4.1);
  h=new TH1D("dose_per_fractionYes", "mean dose st;min dose st [gymax_dose_st];Number of patients",42, -0.1, 4.1);
  h=new TH1D("dose_per_fractionNo", "mean dose st;min dose st [gymax_dose_st];Number of patients",42, -0.1, 4.1);
  h=new TH1D("mean_dose_rectumAll", "mean dose st;min dose st [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("mean_dose_rectumYes", "mean dose st;min dose st [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("mean_dose_rectumNo", "mean dose st;min dose st [gymax_dose_st];Number of patients",101, -0.5, 100.5);
  h=new TH1D("tgt_acutegitox", "tgt_acutegitox;tgt_acutegitox dose st [gymax_dose_st];Number of patients",2, -0.5, 1.5);
  if(glbstuff.clusternum>0){
    h=new TH1D("cluster", "cluster;cluster;Number of patients",glbstuff.clusternum, -0.5, glbstuff.clusternum-0.5);
    h2=new TH2D("cluster_vs_toxicity", "cluster vs toxicity;cluster;toxicity",glbstuff.clusternum, -0.5, glbstuff.clusternum-0.5, 2,-0.5,1.5);
  }
  
  //normalized dvh fofr targets:
  for(auto const &asub:glbstuff.alfabeta){
    h=new TH1D(Form("dvhdiffnorm_tgt_%.3f_Yes",asub), "mean dvh plot for patients with tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhdiffnorm_tgt_%.3f_No",asub), "mean dvh plot for patients without tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhdiffnorm_tgt_%.3f_Diff",asub), "mean dvh plot for difference between patients w/o tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhdiffnorm_tgt_%.3f_All",asub), "mean dvh plot for difference between patients w/o tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhcumnorm_tgt_%.3f_Yes",asub), "mean dvh plot for patients with tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhcumnorm_tgt_%.3f_No",asub), "mean dvh plot for patients without tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhcumnorm_tgt_%.3f_Diff",asub), "mean dvh plot for difference between patients w/o tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    h=new TH1D(Form("dvhcumnorm_tgt_%.3f_All",asub), "mean dvh plot for difference between patients w/o tox;normalized dose [gy];mean volume",glbstuff.maxbin/glbstuff.eqd2binwidth, 0., glbstuff.maxbin);
    
    for(auto const &n:glbstuff.nvalue4eud){
      h=new TH1D(Form("eud_n_%.3f_asub%.3f_All",n,asub), "eud value;eud;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin);    
      h=new TH1D(Form("eud_n_%.3f_asub%.3f_Yes",n,asub), "eud value;eud;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin);    
      h=new TH1D(Form("eud_n_%.3f_asub%.3f_No",n,asub), "eud value;eud;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin);    
    }
  }

  h=new TH1D("dctdiffnorm_tgt_Yes", "mean dct plot for patients with tox;index;value",glbstuff.maxbin, 0., glbstuff.maxbin);
  h=new TH1D("dctdiffnorm_tgt_No", "mean dct plot for patients without tox;index;value",glbstuff.maxbin, 0., glbstuff.maxbin);
  h=new TH1D("dctdiffnorm_tgt_Diff", "mean dct plot for difference between patients w/o tox;index;value",glbstuff.maxbin, 0., glbstuff.maxbin);
  h=new TH1D("dctdiffnorm_tgt_All", "mean dct plot for difference between patients w/o tox;index;value",glbstuff.maxbin, 0., glbstuff.maxbin);
  
  auto paziente=sample.begin();
  for(int i=0;i<paziente->second.clinical_factor.size();i++)
    h2=new TH2D(Form("clinicalfactor_%i_vs_tox",i), "clinical factors vs toxicity;clinical factor;toxycity",6, -1.5, 4.5, 2,-0.5,1.5);

  if(debug)
    cout<<"bookHisto done"<<endl;

  return;
}



void computeDCT(const vector<double>& x, vector<double>& c){
  c.clear();
  c.resize(x.size());
  for (int k = 0; k < x.size(); k++) {
      double sum = 0.0;
      for (int n = 0; n < x.size(); ++n)
          sum += x[n] * cos(M_PI / x.size() * (n + 0.5) * k);
      // normalizzazione ortonormale
      c[k] = ((k == 0) ? sqrt(1.0 / x.size()) : sqrt(2.0 / x.size())) * sum;
  }

    return;
}

// const vector<double> &glbstuff.alfabeta, const vector<double> &glbstuff.nvalue4eud, double glbstuff.eqd2binwidth
void fillHisto(map<int, PatientData> &sample,const globalstuff &glbstuff){

  if(debug)
    cout<<"start fillHisto"<<endl;

  gDirectory->cd("dvhplots");
  map<double,vector<double>> dvhdiffnorm_tgtyes, dvhdiffnorm_tgtno, dvhcumnorm_tgtyes, dvhcumnorm_tgtno;
  int tgtyescounter=0, tgtnocounter=0;
  for(auto const &asub:glbstuff.alfabeta){
    dvhdiffnorm_tgtyes[asub].resize(glbstuff.maxbin/glbstuff.eqd2binwidth,0); //ensure that this length is the same of Form("dvhdiffnorm_tgt_%.3f_Yes",asub)
    dvhdiffnorm_tgtno[asub].resize(glbstuff.maxbin/glbstuff.eqd2binwidth,0);
    dvhcumnorm_tgtyes[asub].resize(glbstuff.maxbin/glbstuff.eqd2binwidth,0); //ensure that this length is the same of Form("dvhcumnorm_tgt_%.3f_Yes",asub)
    dvhcumnorm_tgtno[asub].resize(glbstuff.maxbin/glbstuff.eqd2binwidth,0);
  }

  if(debug>5)
    cout<<"fillHisto: start first loop on sample"<<endl;
  
    vector<double> dctdiff_tgtyes(glbstuff.maxbin,0), dctdiff_tgtno(glbstuff.maxbin,0);
  for(auto &paziente : sample){
    //dvh plots
    gDirectory->mkdir(Form("pz_%i",paziente.second.id));
    gDirectory->cd(Form("pz_%i",paziente.second.id));
    TH1D* hdvhcum=new TH1D(Form("pzdvhcum_%i",paziente.second.id), "DVH cumulative plot;Dose [gy];Volume",paziente.second.dvhmapcum.size(), 0., paziente.second.dvhmapcum.size());
    TH1D* hdvhdiff=new TH1D(Form("pzdvhdiff_%i",paziente.second.id), "DVH differential plot;Dose [gy];Volume",paziente.second.dvhmapdiff.size(), 0., paziente.second.dvhmapdiff.size());
    TH1D* hdctdiff=new TH1D(Form("pzdctdiff_%i",paziente.second.id), "DCT of DVH differential plot;index;value",paziente.second.dvhdctdiff.size(), 0., paziente.second.dvhdctdiff.size());
    bool isboring=true; //si dice vingardium monotonic
    for(int i=0;i<paziente.second.dvhmapcum.size();i++){
      hdvhcum->SetBinContent(i+1,paziente.second.dvhmapcum.at(i));
      if(i>0){
        if(paziente.second.dvhmapcum.at(i)>paziente.second.dvhmapcum.at(i-1)){
          cout<<"WARNING in fillHisto: the dvh function of patient id="<<paziente.second.id<<"is not monotonous"<<endl;
          throw std::runtime_error("patient dhv not monotonous, patient id: " + paziente.second.id);         
          isboring=false;
          break;
        }
      }
    }
    
    if(debug>5)
      cout<<"fillHisto; first sample loop: loop on paziente.second.dvhmapcum.size() done"<<endl;
    
    for(int i=0;i<paziente.second.dvhmapdiff.size();i++)
    hdvhdiff->SetBinContent(i+1,paziente.second.dvhmapdiff.at(i));
    for(int i=0;i<paziente.second.dvhdctdiff.size();i++)
    hdctdiff->SetBinContent(i+1,paziente.second.dvhdctdiff.at(i));
    if(!isboring){
      hdvhcum->SetName(Form("pzdvhcum_%i_notmonotonous",paziente.second.id));
      paziente.second.status=1;
    }    
    
    if(debug>5)
      cout<<"fillHisto; first sample loop: loop on paziente.second.dvhmapdiff.size() done"<<endl;
  
    //scaled dose plots
    for(auto const &asub:glbstuff.alfabeta){
      TH1D* hdvhcumnorm=new TH1D(Form("pzdvhcumnorm_%i_%.2f",paziente.second.id, asub), "Normalized DVH cumulative plot;Equivalent Dose [gy];Volume",paziente.second.dvhcumnormmap.at(asub).size(), 0., paziente.second.dvhcumnormmap.at(asub).size()*glbstuff.eqd2binwidth);
      for(int i=0;i<paziente.second.dvhcumnormmap.at(asub).size();i++)
        hdvhcumnorm->SetBinContent(i+1,paziente.second.dvhcumnormmap.at(asub).at(i));
      
      TH1D* hdvhdiffnorm=new TH1D(Form("pzdvhdiffnorm_%i_%.2f",paziente.second.id, asub), "Normalized DVH differential plot;Equivalent Dose [gy];Volume",paziente.second.dvhdiffnormmap.at(asub).size(), 0., paziente.second.dvhdiffnormmap.at(asub).size()*glbstuff.eqd2binwidth);
      for(int i=0;i<paziente.second.dvhdiffnormmap.at(asub).size();i++)
        hdvhdiffnorm->SetBinContent(i+1,paziente.second.dvhdiffnormmap.at(asub).at(i));      
    }

    if(debug>5)
    cout<<"fillHisto; first sample loop: loop on scaled dose plots done"<<endl;
    
    gDirectory->cd(".."); 
    
    (dynamic_cast<TH1D*>(gDirectory->Get("../sample_volumesAll")))->Fill(paziente.second.volume);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../sample_volumesYes" : "../sample_volumesNo")))->Fill(paziente.second.volume);
    (dynamic_cast<TH1D*>(gDirectory->Get("../max_dose_planAll")))->Fill(paziente.second.max_dose_plan);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../max_dose_planYes" : "../max_dose_planNo")))->Fill(paziente.second.max_dose_plan);
    (dynamic_cast<TH1D*>(gDirectory->Get("../max_dose_stAll")))->Fill(paziente.second.max_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../max_dose_stYes":"../max_dose_stNo")))->Fill(paziente.second.max_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../min_dose_stAll")))->Fill(paziente.second.min_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../min_dose_stYes":"../min_dose_stNo")))->Fill(paziente.second.min_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../mean_dose_stAll")))->Fill(paziente.second.mean_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../mean_dose_stYes" : "../mean_dose_stNo")))->Fill(paziente.second.mean_dose_st);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dose_ptvAll")))->Fill(paziente.second.dose_ptv);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../dose_ptvYes" : "../dose_ptvNo")))->Fill(paziente.second.dose_ptv);
    (dynamic_cast<TH1D*>(gDirectory->Get("../nfractionAll")))->Fill(paziente.second.nfraction);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../nfractionYes" : "../nfractionNo")))->Fill(paziente.second.nfraction);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dose_per_fractionAll")))->Fill(paziente.second.dose_per_fraction);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../dose_per_fractionYes": "../dose_per_fractionNo")))->Fill(paziente.second.dose_per_fraction);
    (dynamic_cast<TH1D*>(gDirectory->Get("../mean_dose_rectumAll")))->Fill(paziente.second.mean_dose_rectum);
    (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? "../mean_dose_rectumYes" : "../mean_dose_rectumNo")))->Fill(paziente.second.mean_dose_rectum);
    (dynamic_cast<TH1D*>(gDirectory->Get("../tgt_acutegitox")))->Fill(paziente.second.tgt_acutegitox);
    if(glbstuff.clusternum>0){
      (dynamic_cast<TH1D*>(gDirectory->Get("../cluster")))->Fill(paziente.second.cluster);
      (dynamic_cast<TH2D*>(gDirectory->Get("../cluster_vs_toxicity")))->Fill(paziente.second.cluster, paziente.second.tgt_acutegitox);
    }
    for(int i=0;i<paziente.second.clinical_factor.size();i++)
      (dynamic_cast<TH2D*>(gDirectory->Get(Form("../clinicalfactor_%i_vs_tox",i))))->Fill(paziente.second.clinical_factor.at(i),paziente.second.tgt_acutegitox);
    
    if(debug>5)
      cout<<"fillHisto; first sample loop: fill sample_volumes etc plots done"<<endl;
    
    if(paziente.second.tgt_acutegitox>0) tgtyescounter++; else tgtnocounter++;
    
    for(auto const &asub:glbstuff.alfabeta){
      for(int i=0;i<paziente.second.dvhdiffnormmap.at(asub).size();i++){
        if(paziente.second.tgt_acutegitox>0)
          dvhdiffnorm_tgtyes[asub].at(i)+=paziente.second.dvhdiffnormmap.at(asub).at(i);
        else
          dvhdiffnorm_tgtno[asub].at(i)+=paziente.second.dvhdiffnormmap.at(asub).at(i);
      }

      for(int i=0;i<paziente.second.dvhcumnormmap.at(asub).size();i++){
        if(paziente.second.tgt_acutegitox>0)
          dvhcumnorm_tgtyes[asub].at(i)+=paziente.second.dvhcumnormmap.at(asub).at(i);
        else
          dvhcumnorm_tgtno[asub].at(i)+=paziente.second.dvhcumnormmap.at(asub).at(i);
      }
      
      for(auto const &n:glbstuff.nvalue4eud){
        (dynamic_cast<TH1D*>(gDirectory->Get(Form("../eud_n_%.3f_asub%.3f_All",n,asub))))->Fill(paziente.second.eudmap.at(make_pair(n,asub)));
        (dynamic_cast<TH1D*>(gDirectory->Get(paziente.second.tgt_acutegitox>0 ? Form("../eud_n_%.3f_asub%.3f_Yes",n,asub) : Form("../eud_n_%.3f_asub%.3f_No",n,asub)  )))->Fill(paziente.second.eudmap.at(make_pair(n,asub)));
      }
    }

    if(debug>5)
      cout<<"fillHisto; first sample loop: fill dct stuff"<<endl;    

      if(paziente.second.tgt_acutegitox>0)
      for(int i=0;i<paziente.second.dvhdctdiff.size();i++)
        dctdiff_tgtyes.at(i)+=paziente.second.dvhdctdiff.at(i);
    else
      for(int i=0;i<paziente.second.dvhdctdiff.size();i++)
        dctdiff_tgtno.at(i)+=paziente.second.dvhdctdiff.at(i);
  }//end of loop on patient

  if(debug>5)
  cout<<"fillHisto: start fill dctdiffnorm_tgt_*"<<endl;
  
  for(int i=0;i<glbstuff.maxbin;i++){
    (dynamic_cast<TH1D*>(gDirectory->Get("../dctdiffnorm_tgt_All")))->SetBinContent(i+1,(dctdiff_tgtyes.at(i)+dctdiff_tgtno.at(i))/(tgtnocounter+tgtyescounter));
    (dynamic_cast<TH1D*>(gDirectory->Get("../dctdiffnorm_tgt_Yes")))->SetBinContent(i+1,dctdiff_tgtyes.at(i)/tgtyescounter);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dctdiffnorm_tgt_No")))->SetBinContent(i+1,dctdiff_tgtno.at(i)/tgtnocounter);
    (dynamic_cast<TH1D*>(gDirectory->Get("../dctdiffnorm_tgt_Diff")))->SetBinContent(i+1,dctdiff_tgtyes.at(i)/tgtyescounter-dctdiff_tgtno.at(i)/tgtnocounter);
  }
  
  if(debug>5)
    cout<<"fillHisto: start last loop on alfabeta and dvhdiffnorm_tgtyes.at(asub).size()"<<endl;

  
  for(auto const &asub:glbstuff.alfabeta){
    for(int i=0;i<dvhdiffnorm_tgtyes.at(asub).size();i++){
      // cout<<"suba="<<asub<<" i="<<i<<" dvhdiffnorm_tgtyes.at(asub).at(i)/sample.size()="<<dvhdiffnorm_tgtyes.at(asub).at(i)/sample.size()<<endl;
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhdiffnorm_tgt_%.3f_All",asub))))->SetBinContent(i+1,(dvhdiffnorm_tgtyes.at(asub).at(i)+dvhdiffnorm_tgtno.at(asub).at(i))/(tgtnocounter+tgtyescounter));
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhdiffnorm_tgt_%.3f_Yes",asub))))->SetBinContent(i+1,dvhdiffnorm_tgtyes.at(asub).at(i)/tgtyescounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhdiffnorm_tgt_%.3f_No",asub))))->SetBinContent(i+1,dvhdiffnorm_tgtno.at(asub).at(i)/tgtnocounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhdiffnorm_tgt_%.3f_Diff",asub))))->SetBinContent(i+1,dvhdiffnorm_tgtyes.at(asub).at(i)/tgtyescounter-dvhdiffnorm_tgtno.at(asub).at(i)/tgtnocounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhcumnorm_tgt_%.3f_All",asub))))->SetBinContent(i+1,(dvhcumnorm_tgtyes.at(asub).at(i)+dvhdiffnorm_tgtno.at(asub).at(i))/(tgtnocounter+tgtyescounter));
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhcumnorm_tgt_%.3f_Yes",asub))))->SetBinContent(i+1,dvhcumnorm_tgtyes.at(asub).at(i)/tgtyescounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhcumnorm_tgt_%.3f_No",asub))))->SetBinContent(i+1,dvhcumnorm_tgtno.at(asub).at(i)/tgtnocounter);
      (dynamic_cast<TH1D*>(gDirectory->Get(Form("../dvhcumnorm_tgt_%.3f_Diff",asub))))->SetBinContent(i+1,dvhcumnorm_tgtyes.at(asub).at(i)/tgtyescounter-dvhcumnorm_tgtno.at(asub).at(i)/tgtnocounter);
    }
  }

  gDirectory->cd(".."); 

  if(debug)
    cout<<"fillHisto done"<<endl;

  return;
}

void PostLoopAnalysis(map<int, PatientData> &sample, const globalstuff &glbstuff){

  if(debug)
    cout<<"start PostLoopAnalysis"<<endl;
  //WARNING: THE PEAK POSITION IS SET BY HAND!!!!!!!!!!!!!! HERE
  //tentative to fit the peak at "high" dose in dvhcumnorm_tgt_%.3f_Diff ...and try to do stuff...
  TF1 mygaus=TF1("mygaus","-gaus", 0., glbstuff.maxbin);
  for(auto const &asub:glbstuff.alfabeta){
    mygaus.SetParameter(0,-(dynamic_cast<TH1D*>(gDirectory->Get(Form("dvhcumnorm_tgt_%.3f_Diff",asub))))->GetBinContent((dynamic_cast<TH1D*>(gDirectory->Get(Form("dvhcumnorm_tgt_%.3f_Diff",asub))))->FindBin(30)));
    mygaus.SetParameter(1,20);
    (dynamic_cast<TH1D*>(gDirectory->Get(Form("dvhcumnorm_tgt_%.3f_Diff",asub))))->Fit("mygaus","QBM+","",20,40);
    if(mygaus.GetParameter(1)<30 && mygaus.GetParameter(1)>10){
      TH1D* hAll=new TH1D(Form("dvhcumnorm_maxdiff_%.3f_All",asub), "Volume for dvhcumnnorm maxdiff dose value;Volume;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin*10);
      TH1D* hYes=new TH1D(Form("dvhcumnorm_maxdiff_%.3f_Yes",asub), "Volume for dvhcumnorm maxdiff dose value;Volume;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin*10);
      TH1D* hNo=new TH1D(Form("dvhcumnorm_maxdiff_%.3f_No",asub), "Volume for dvhcumnorm maxdiff dose value;Volume;Number of patients",glbstuff.maxbin, 0., glbstuff.maxbin*10);
      TVectorD xv(sample.size()), yv(sample.size());
      int counter=0;
      for(auto const &paziente:sample){
        int binpos=(dynamic_cast<TH1D*>(gDirectory->Get(Form("dvhcumnorm_tgt_%.3f_Diff",asub))))->FindBin(mygaus.GetParameter(1));
        if(paziente.second.dvhcumnormmap.at(asub).size()<=binpos)
          continue;
        xv(counter)=paziente.second.dvhcumnormmap.at(asub).at(binpos);
        yv(counter)=paziente.second.tgt_acutegitox;
        hAll->Fill(xv(counter));
        if(yv(counter)>0)
          hYes->Fill(xv(counter));
        else
          hNo->Fill(xv(counter));
        counter++;
      }
      CreateNtcpSigmoidalSingle(Form("ntcp_linear_maxdiff_dvhcumnorm_%.3f_",asub),xv,yv);     
    }else{
      cout<<"PostLoopAnalysis: fit on "<<Form("dvhcumnorm_tgt_%.3f_Diff",asub)<<" failed, skip the postloop analysis on this case"<<endl;
    }
  }
  
return;
}

//select the algorithm with the best auc, cosidering only tracking status==0 and covmatrixstatus==3 (accurate)
void ChooseBestFit(globalstuff &glbstuff){

  if(debug)
    cout<<"start ChooseBestFit"<<endl;

  glbstuff.bestvalue=-1;
  double bestauc=-1;
  for(const auto &v:glbstuff.fitresults){
    if(v.second.at(0)==0 && v.second.at(1)==3){ //status==0 and covariance matrix=3: "accurate"
      if(v.second.at(7)>bestauc){
        glbstuff.bestvalue=v.first;
        bestauc=v.second.at(7);
      }
    }
  }

  if(glbstuff.bestvalue>=0)
    cout<<"ChooseBestFit done, bestvalue="<<glbstuff.bestvalue<<"  best algorithm="<<glbstuff.fitalgo.at(glbstuff.bestvalue).first<<"/"<<glbstuff.fitalgo.at(glbstuff.bestvalue).second<<endl;
  else
    cout<<"WARNING: ChooseBestFit: glbstuff.bestvalue cannot be set"<<endl;

  return;
}

double fitSigmoidal(TGraph* graph, int parnum, int functype){

  if(debug)
    cout<<"start fitSigmoidal on "<<graph->GetName()<<"  graph->GetN()="<<graph->GetN()<<endl;

  if(parnum<2)
    throw std::runtime_error("fitSigmoidal: parnum too small" + parnum);
    
    TF1 *sigmoid;
  if(parnum==2){
    if(functype==0)
      sigmoid=new TF1("sigmoid", "1./(1.+exp(-[0]-[1]*x))", 0, 1000);
    else if(functype==1)
      sigmoid=new TF1("signoexp", "1./(1.+[0]/(x^[1]))", 0, 1000);
  }
  else //per ora solo caso semplice con due parametri liberi ed una sola variabile X=EUD
    throw std::runtime_error("fitSigmoidal: parnum too large" + parnum);

  TFitResultPtr fitres=graph->Fit(sigmoid, "EMSQ+","",0,1000);

  if(debug)
    cout<<"fitSigmoidal done"<<endl;

  if(!fitres->IsEmpty())
    return fitres->Chi2();

  return -1;
}

int CreateNtcpSigmoidalSingle(TString pltname, const TVectorD &vx, const TVectorD &vy){
  
  if(debug)
    cout<<"start CreateNtcpSigmoidalSingle on "<<pltname<<endl;

  TDirectory *dir=gDirectory->GetDirectory("ntcp_linear");
  if(dir){
    gDirectory->cd("ntcp_linear");
  }else{
    cout<<"WARNING in CreateNtcpSigmoidalSingle: ntcp_linear not present... I'll create a new one"<<endl;
    gDirectory->mkdir("ntcp_linear");
    gDirectory->cd("ntcp_linear");
  }
  
  TGraph *gr=new TGraph(vx, vy);
  TH1D* h;
  gr->SetName(pltname);
  fitSigmoidal(gr,2,1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(2);
  gr->SetLineWidth(0);
  gr->SetLineColor(0);
  gr->SetDrawOption("AP*");
  gr->Write();
  // CreateHistoFromTgraph(gr,h); //to be fixed

  gDirectory->cd("..");
  
  return 0;
}

int CreateNtcpSigmoidalPredefined(map<int, PatientData> &sample, TString tgtname){

  if(debug)
    cout<<"start CreateNtcpSigmoidalPredefined"<<endl;

  gDirectory->mkdir("ntcp_linear");
  gDirectory->cd("ntcp_linear");

  map <string, pair<TGraph*,TH1D*>> ntcp_linear_models;
  TGraph *gr;
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_volume_"+tgtname);
  ntcp_linear_models["volume"]=make_pair(gr, nullptr);
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_max_dose_plan_"+tgtname);
  ntcp_linear_models["max_dose_plan"]=make_pair(gr, nullptr);
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_max_dose_st_"+tgtname);
  ntcp_linear_models["max_dose_st"]=make_pair(gr, nullptr);
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_min_dose_st_"+tgtname);
  ntcp_linear_models["min_dose_st"]=make_pair(gr, nullptr);
  gr=new TGraph(sample.size());
  gr->SetName("ntcp_linear_mean_dose_st_"+tgtname);
  ntcp_linear_models["mean_dose_st"]=make_pair(gr, nullptr);
  TH1D* h_chi2=new TH1D("fittedchi2", "reduced chi2;Reduced chi2;functions",100, 0., 50);  
  
  for(auto const &paziente:sample){
    for(auto const &comb:paziente.second.eudmap){
      gr=new TGraph(sample.size());
      gr->SetName(Form("ntcp_linear_eud_n%.2f_alphabeta%.3f_",comb.first.first, comb.first.second)+tgtname);
      ntcp_linear_models[Form("ntcp_linear_eud_n%.2f_alphabeta%.3f_",comb.first.first, comb.first.second)]=make_pair(gr, nullptr);
    }
  }

  int index=0;
  for(auto &paziente : sample){
    ntcp_linear_models.at("volume").first->SetPoint(index, paziente.second.volume, paziente.second.tgt_acutegitox);
    ntcp_linear_models.at("max_dose_plan").first->SetPoint(index, paziente.second.max_dose_plan, paziente.second.tgt_acutegitox);
    ntcp_linear_models.at("max_dose_st").first->SetPoint(index, paziente.second.max_dose_st, paziente.second.tgt_acutegitox);
    ntcp_linear_models.at("min_dose_st").first->SetPoint(index, paziente.second.min_dose_st, paziente.second.tgt_acutegitox);
    ntcp_linear_models.at("mean_dose_st").first->SetPoint(index, paziente.second.mean_dose_st, paziente.second.tgt_acutegitox);
    for(auto const &comb:paziente.second.eudmap)
      ntcp_linear_models.at(Form("ntcp_linear_eud_n%.2f_alphabeta%.3f_",comb.first.first, comb.first.second)).first->SetPoint(index, comb.second,paziente.second.tgt_acutegitox);
    index++;
  }

  for(auto &graph:ntcp_linear_models){
    graph.second.first->Sort();
    double chi2=fitSigmoidal(graph.second.first,2,1);
    h_chi2->Fill(chi2);
    graph.second.first->SetMarkerStyle(20);
    graph.second.first->SetMarkerColor(2);
    graph.second.first->SetLineWidth(0);
    graph.second.first->SetLineColor(0);
    graph.second.first->SetDrawOption("AP*");
    graph.second.first->Write();
    // CreateHistoFromTgraph(graph.second.first,graph.second.second); //TODO: to be fixed
  }

  gDirectory->cd("..");

  if(debug)
    cout<<"CreateNtcpSigmoidalPredefined done"<<endl;  

  return 0;
}

//TODO: To be fixed
void CreateHistoFromTgraph(TGraph *gr, TH1D *h){

  if(debug)
    cout<<"start CreateHistoFromTgraph on "<<gr->GetName()<<endl;    

  double grxmin, grxmax, grymin, grymax;
  gr->ComputeRange(grxmin, grymin, grxmax, grymax);
  if(grymax!=1 || grymin!=0){
    cout<<"Error in CreateHistoFromTgraph: grymax or grymin are not 1 or 0, probably no target case!  grymax="<<grymax<<"  grymin="<<grymin<<endl;
    cout<<"gr->GetName()="<<gr->GetName()<<endl;
  }
  // int hbinnum=(int)((grxmax-grxmin)/2);
  int hbinnum=20;
  TString histtitle(gr->GetName());
  histtitle+="_hist";
  h=new TH1D(histtitle,gr->GetTitle(),hbinnum,grxmin, grxmax);
  vector<pair<double, int>> valuenum(hbinnum,make_pair(0,0)); //first=bin value, second=counter
  for(int i=0;i<gr->GetN();i++){
    int bin=(int)((gr->GetPointX(i)-grxmin)/hbinnum);
    if(bin<hbinnum){
      valuenum.at(bin).first+= gr->GetPointY(i);
      valuenum.at(bin).second++;
    }
  }
  //set the histo
  double valbefore=0;
  int indexbefore=0;
  for(int i=0;i<valuenum.size();i++){
    if(valuenum.at(i).second==0){
      double valafter=valbefore;
      int indexafter=i+1;
      for(int k=i+1;k<valuenum.size();k++){
        if(valuenum.at(k).second!=0){
          valafter=valuenum.at(k).first/valuenum.at(k).second;
          indexafter=k;
          break;
        }
      }
      double weightedvalue=(valbefore/(i-indexbefore)+valafter/(indexafter-i))/(indexafter-indexbefore);
      valuenum.at(i).first=weightedvalue;
      for(int k=i+1;k<valuenum.size() && k<=indexafter;k++){
        valuenum.at(i).first=weightedvalue;
        i=k;
      }
     }else{
      valuenum.at(i).first=valuenum.at(i).first/valuenum.at(i).second;
      valbefore=valuenum.at(i).first;
      indexbefore=i;
    }
  }
  //fill the histo
  for(int i=0;i<hbinnum;i++)
    h->SetBinContent(i+1,valuenum.at(i).first);

  if(debug)
    cout<<"CreateHistoFromTgraph done"<<endl;    

  return;
}

void PrintSampleLine(int idx, const   map<int, PatientData> &sample){
  cout<<"Print sample at index: "<<idx<<endl;  
  if(sample.count(idx)==0 ){
    throw std::runtime_error("PrintSampleLine: idx not present in sample" + idx);
    return;
  }
  cout<<sample.at(idx).id<<" "<<sample.at(idx).volume<<" "<<sample.at(idx).max_dose_plan<<" "<<sample.at(idx).max_dose_st<<" "<<sample.at(idx).min_dose_st<<" "<<sample.at(idx).mean_dose_st<<" "<<sample.at(idx).dose_ptv<<" "<<sample.at(idx).nfraction<<" "<<sample.at(idx).dose_per_fraction<<" "<<sample.at(idx).mean_dose_rectum<<endl;
  for(const double& value : sample.at(idx).dvhmapcum)
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

int loadSyntheticFile(const string& filename,   map<int, PatientData> &sample){

  if(debug)
    cout<<"start loadSyntheticFile on "<<filename<<endl;

  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("loadSyntheticFile: cannot open file " + filename);
  }
  string line;

  if (!std::getline(in, line)){ //throw away first header line
    throw std::runtime_error("loadSyntheticFile: empty file " + filename);
    return 1;
  }

  while (std::getline(in, line)) {

    if(debug>5)
      cout<<"readed new line: "<<line<<endl;

      PatientData patient;
    if (line.empty())
      continue;
    TString delimiter(" ");
    vector<string> parts = splitCsvLine(line, delimiter);

    if(debug>5)
      cout<<"parts.size()= "<<parts.size()<<endl;

    if (parts.size() < 6){
      // throw std::runtime_error("loadDvhFile: check your dvh file, there is a line with less than 8 elements: " + filename);
      continue;
    }

    patient.id=std::stoi(trim(parts[0]));
    patient.clinical_factor.push_back(std::stoi(trim(parts[1])));
    patient.nfraction=std::stoi(trim(parts[2]));
    patient.tgt_acutegitox=std::stoi(trim(parts[3]));
    patient.operation=std::stoi(trim(parts[4]));

    for (size_t i = 5; i < parts.size(); i+=2) {
      if(parts[i].size()==0)
        break;
      patient.dvhmapcum.push_back(std::stod(trim(parts[i]))/100.);
    }

    patient.status=-1;
    patient.dvhmapdiff.push_back(0);
    for(int i=1; i<patient.dvhmapcum.size();i++)
      patient.dvhmapdiff.push_back(patient.dvhmapcum.at(i-1)-patient.dvhmapcum.at(i));
    computeDCT(patient.dvhmapdiff, patient.dvhdctdiff);
    sample.insert(std::pair{patient.id,patient});
}

  if(debug)
    cout<<"loadSyntheticFile done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;
  
    return 0;

}

 int loadDvhFile(const string& filename,   map<int, PatientData> &sample) {

  if(debug)
    cout<<"start loadDvhFile on "<<filename<<endl;

  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("loadDvhFile: cannot open file " + filename);
  }
  string line;
  // ---------- prima riga: header ----------
  if (!std::getline(in, line)){
    throw std::runtime_error("loadDvhFile: empty file " + filename);
    return 1;
  }

  // regex per id tipo RQ12345-6
  std::regex id_pattern(R"(ML-\d{4})");
  // ---------- righe dati ----------
  while (std::getline(in, line)) {

    if(debug>5)
      cout<<"readed new line: "<<line<<endl;

      PatientData patient;
    if (line.empty())
      continue;
    TString delimiter(" ");
    vector<string> parts = splitCsvLine(line, delimiter);

    if(debug>5)
      cout<<"parts.size()= "<<parts.size()<<endl;

    if (parts.size() < 8){
      // throw std::runtime_error("loadDvhFile: check your dvh file, there is a line with less than 8 elements: " + filename);
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
      patient.dvhmapcum.push_back(std::stod(readed)/patient.volume);
    }
    patient.status=-1;
    patient.dvhmapdiff.push_back(0);
    double sum=0;
    for(int i=1; i<patient.dvhmapcum.size();i++){
      patient.dvhmapdiff.push_back(patient.dvhmapcum.at(i-1)-patient.dvhmapcum.at(i));
      sum+=patient.dvhmapdiff.back();
    }
    //normalize dvhmapdiff
    for(int i=0;i<patient.dvhmapdiff.size();i++)
      patient.dvhmapdiff.at(i)=patient.dvhmapdiff.at(i)/sum;
    
    computeDCT(patient.dvhmapdiff, patient.dvhdctdiff);

    patient.dose_ptv=-1;
    patient.nfraction=-1;
    patient.dose_per_fraction=-1;
    patient.mean_dose_rectum=-1; 
    patient.tgt_acutegitox=-1;
    patient.operation=-1; 
    patient.prostatectomy=-1;
    patient.appendectomy=-1;
    patient.age=-1;
    patient.bmi=-1;
    patient.smoke=-1;
    patient.alcohol=-1;
    patient.diabetes=-1;
    patient.diverticulitis=-1;
    patient.psa=-1;
    patient.tumour_risk=-1;
    patient.lymph_irr=-1;
    patient.semivesicle_irr=-1;
    patient.mb_risk=-1;
    patient.cluster=-99;
    sample.insert(std::pair{id,patient});
}

  if(debug)
    cout<<"loadDvhFile done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;

    return 0;
}


int CheckSampleSamrectConsistency(const map<int, PatientData> &sample, const map<int, PatientData> &samrect){
  if(sample.size()!=samrect.size()){
    cout<<"ERROR: CheckSampleSamrectConsistency:: sample.size()="<<sample.size()<<"  samrect.size()="<<samrect.size()<<endl;
    return 1;
  }
  for(const auto &paziente:sample){
    if(samrect.at(paziente.second.id).tgt_acutegitox!=paziente.second.tgt_acutegitox || samrect.at(paziente.second.id).powp!=paziente.second.powp || samrect.at(paziente.second.id).mb_risk!=paziente.second.mb_risk ||samrect.at(paziente.second.id).cluster!=paziente.second.cluster || samrect.at(paziente.second.id).mean_dose_rectum!=paziente.second.mean_dose_rectum  ){
      cout<<"ERROR in CheckSampleSamrectConsistency"<<endl;
      return 1;
    }
  }

  return 0;
}



int optimizeLikehood(map<int, PatientData> &sample, globalstuff &glbstuff, const int fitalgindex,  map<int, PatientData> &samrect, const pair<int,double> fixedpar) {
  
  if(debug)
    cout<<"start optimizeLikehood, fitalgindex="<<fitalgindex<<endl;

  //lambda
  auto lamdalikehoodFull = [&](const double* par) {return functorLikehoodFull(sample, par);};
  auto lamdalikehoodAlfabdone = [&](const double* par) {return functorLikehoodAlfabdone(sample, par);};
  auto lamdalikehoodFullClinical_0 = [&](const double* par) {return functorLikehoodFullClinical_0(sample, par);};
  auto lamdalikehoodAlfabdoneClinical_0 = [&](const double* par) {return functorLikehoodAlfabdoneClinical_0(sample, par);};
  auto lamdalikehoodFullClinical_1 = [&](const double* par) {return functorLikehoodFullClinical_1(sample, par);};
  auto lamdalikehoodAlfabdoneClinical_1 = [&](const double* par) {return functorLikehoodAlfabdoneClinical_1(sample, par);};
  auto lamdalikehoodAlfabdone2DvhClinical_1 = [&](const double* par) {return functorLikehoodAlfabdone2DvhClinical_1(sample, samrect, par);};  

  auto lambdalikehoodDose4Vol = [&](const double* par) {return functorLikehoodDose4Vol(sample, par);};  
  auto lambdalikehoodDose4VolClinical = [&](const double* par) {return functorLikehoodDose4VolClinical(sample, par);};  



  // ROOT::Math::Functor fpFunctor= (glbstuff.clinicalfactors.size()==0) ? ((glbstuff.alfabdone < 0) ?  ROOT::Math::Functor(lamdalikehoodFull, 4) : ROOT::Math::Functor(lamdalikehoodAlfabdone, 3)) : ((glbstuff.alfabdone < 0) ? ROOT::Math::Functor(lamdalikehoodFullClinical_0, 5) : ROOT::Math::Functor(lamdalikehoodAlfabdoneClinical_0, 4 ));
  ROOT::Math::Functor fpFunctor;
  if(glbstuff.usedosevar==-1){
    if(glbstuff.clinicalfactors == 0 && glbstuff.alfabdone<0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodFull"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodFull, 4);
    }else if(glbstuff.clinicalfactors==0 && glbstuff.alfabdone>=0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodAlfabdone"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdone, 3);
    }else if(glbstuff.clinicalfactors==1 && glbstuff.alfabdone<0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodFullClinical_0"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodFullClinical_0, 5);
    }else if(glbstuff.clinicalfactors==1 && glbstuff.alfabdone>=0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodAlfabdoneClinical_0"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdoneClinical_0, 4);
    }else if(glbstuff.clinicalfactors==2 && glbstuff.alfabdone<0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodFullClinical_1"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodFullClinical_1, 6);
    }else if(glbstuff.clinicalfactors==2 && glbstuff.alfabdone>=0 && glbstuff.twodvh==0){
      cout<<"optimizeLikehood: functor=lamdalikehoodAlfabdoneClinical_1"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdoneClinical_1, 5);
    }else if(glbstuff.clinicalfactors==2 && glbstuff.alfabdone>=0 && glbstuff.twodvh>0){
      cout<<"optimizeLikehood: functor=lamdalikehoodAlfabdone2DvhClinical_1"<<endl;
      fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdone2DvhClinical_1, 6);
    }
  }else{
    if(glbstuff.clinicalfactors==2 && glbstuff.twodvh==0)
      fpFunctor = ROOT::Math::Functor(lambdalikehoodDose4VolClinical, 6); 
    else
      fpFunctor = ROOT::Math::Functor(lambdalikehoodDose4Vol, 4); 

    cout<<"optimizeLikehood: functor=lambdalikehoodDose4Vol, usedosevar="<<glbstuff.usedosevar<<"  doses4volume.at(glbstuff.usedosevar)="<<glbstuff.doses4volume.at(glbstuff.usedosevar)<<endl;
  }

  std::unique_ptr<ROOT::Math::Minimizer> fpMinimizer(ROOT::Math::Factory::CreateMinimizer(glbstuff.fitalgo.at(fitalgindex).first, glbstuff.fitalgo.at(fitalgindex).second));

  fpMinimizer->SetFunction(fpFunctor);
  fpMinimizer->SetMaxFunctionCalls(1000000);
  fpMinimizer->SetMaxIterations(100000);
  fpMinimizer->SetTolerance(1e-3);//to be checked
  fpMinimizer->SetStrategy(2);
  fpMinimizer->SetErrorDef(0.5); // I minimize a -log likehoood, so the error def should be 0.5
  if(debug)
    fpMinimizer->SetPrintLevel(2);
  
  std::map<int, decltype(glbstuff.fitpars.begin())> ordered;

  if(debug)
    cout<<"optimizeLikehood: set variables"<<endl;

  // ROOT/Minuit2 uses the variable registration order as the FCN par[] order.
  // Therefore variables must be registered sorted by index.
  for (auto it = glbstuff.fitpars.begin(); it != glbstuff.fitpars.end();++it)
      ordered[it->second.first] = it;

  for (const auto& p : ordered){
      const auto& name = p.second->first;
      const auto& lim  = p.second->second.second;
      if(fixedpar.first!=p.first){
        if(lim[1]>0){
          if(debug)
            cout<<"fpMinimizer->setlimtedvariable index="<<p.first<<"  name="<<name<<endl;
            fpMinimizer->SetLimitedVariable(p.first, name.c_str(), lim[0], lim[1], lim[2],lim[3]);
          }else{
            if(debug)
              cout<<"fpMinimizer->setfixed variableindex="<<p.first<<"  name="<<name<<endl;
          fpMinimizer->SetFixedVariable(p.first, name.c_str(), lim[0]);
        }
        glbstuff.nameindex[p.first]=name.c_str();
      }else{
        if(debug)
          cout<<"fpMinimizer->setFixedVariable index="<<p.first<<"  name="<<name<<endl;        
        fpMinimizer->SetFixedVariable(p.first, name.c_str(), fixedpar.second);
        glbstuff.nameindex[p.first]=name.c_str();
      }
  }  

  if(debug)
    cout<<"optimizeLikehood: before minimization"<<endl;  

  fpMinimizer->Minimize();
  if(debug)
    cout<<"optimizeLikehood: minimization done"<<endl; 
// using Minuit2 you have as status: (https://root-forum.cern.ch/t/is-fit-validity-or-minimizer-status-more-important/30637)
// status = 1 : Covariance was made pos defined
// status = 2 : Hesse is invalid
// status = 3 : Edm is above max
// status = 4 : Reached call limit
// status = 5 : Covariance is not positive defined   
  Int_t status=fpMinimizer->Status();
  // if(fixedpar.first<0 || (fixedpar.first>=0 && )){ //no parameter was fixed, this is a minimization of the full model
  if(fixedpar.first<0){ //no parameter was fixed, this is a minimization of the full model
    cout<<endl<<endl<<"optimizeLikehood: minimization done with "<<glbstuff.fitalgo.at(fitalgindex).first<<"/"<<glbstuff.fitalgo.at(fitalgindex).second<<endl;
    if(status==0){
      double par[10]; //10 per essere sicuri di prendere tutti i parametri
      for(int i=0;i<fpMinimizer->NDim();i++){
        par[i]=fpMinimizer->X()[i];
      }
      for(auto &paziente : sample){//fill scores
        paziente.second.optlike_ntcpscore[fitalgindex]= (glbstuff.twodvh==0) ? EvalScoreSelector(glbstuff, paziente.second, par) : EvalScoreAlfabdone2DvhClinical_1(paziente.second, samrect.at(paziente.second.id), par)  ;
      }
    }else{
      cout<<"minimization failed"<<endl;  
      cout<<"print correlations:"<<endl;
      for(int i=0;i<fpMinimizer->NFree();++i)
        for(int j=0;j<fpMinimizer->NFree();++j)
          std::cout << i << " " << j << " " << fpMinimizer->CovMatrix(i,j)/sqrt(fpMinimizer->CovMatrix(i,i)*fpMinimizer->CovMatrix(j,j)) << std::endl;  
    }
    cout<<"status="<<status<<" CovMatrixStatus="<<fpMinimizer->CovMatrixStatus()<<"  Edm="<<fpMinimizer->Edm()<<"  degree of freedom:"<<fpMinimizer->NFree()<<"  fpMinimizer->NDim()="<<fpMinimizer->NDim()<<"   -loglikehood minimum value:"<<fpMinimizer->MinValue()<<"  AIC="<<2*fpMinimizer->NFree()+2*fpMinimizer->MinValue()<<"  AICc_reduced="<<2*fpMinimizer->NFree()+2*fpMinimizer->MinValue()+2*fpMinimizer->NFree()*(fpMinimizer->NFree()+1)/(sample.size()-fpMinimizer->NFree()-1)<<"  deviance/dof="<<2*fpMinimizer->MinValue()/(sample.size()-fpMinimizer->NFree())<<endl;

    for(int i=0;i<fpMinimizer->NDim();i++){
      if(glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size()==0){
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].push_back(fpMinimizer->X()[i]);
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].push_back(fpMinimizer->Errors()[i]);
        TString tmp_title=glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second+"_likelihood_"+fpMinimizer->VariableName(i);
        TH1D* h=new TH1D(tmp_title, "variable value;likelihood value",200, fpMinimizer->X()[i]-3*fpMinimizer->Errors()[i], fpMinimizer->X()[i]+3.*fpMinimizer->Errors()[i]);   
      }else{
        cout<<"WARNING:optimizeLikehood: glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size() not =0, the fitted values has already filled... this is strange, please take a look; fitalgindex="<<fitalgindex<<"  size="<<glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size()<<endl;
        cout<<"I'll overwrite the new values:"<<endl;
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].at(0)=fpMinimizer->X()[i];
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].at(1)=fpMinimizer->Errors()[i];
      }
      cout<<fpMinimizer->VariableName(i)<<":  "<<fpMinimizer->X()[i]<<" +- "<<fpMinimizer->Errors()[i]<<endl;
    }
    glbstuff.fitresults[fitalgindex].insert(glbstuff.fitresults[fitalgindex].end(), {(double)status, (double)fpMinimizer->CovMatrixStatus(), fpMinimizer->Edm(), (double)fpMinimizer->NFree(), fpMinimizer->MinValue(), 2*fpMinimizer->NFree()+2*fpMinimizer->MinValue(), 2*fpMinimizer->MinValue()/(sample.size()-fpMinimizer->NFree())});    
  }else{ //fixed par just for the calculation of Likehood Ratio Test (LRT) TODO: LRT case for dose4volume case need to be done
    if(debug)
      cout<<"minimization done for LRT for variable="<<fpMinimizer->VariableName(fixedpar.first)<<" with status="<<status<<endl;
    if(status==0){
      if(glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].size()==2){
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].push_back(fpMinimizer->MinValue()); //append -loglikehood value for fixed par
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].push_back( ROOT::Math::chisquared_cdf_c(2.* (fpMinimizer->MinValue()-glbstuff.fitresults[fitalgindex].at(4)),1) ); //append pvalue
        cout<<"LRT for "<<fpMinimizer->VariableName(fixedpar.first)<<"  H_reduced w/o variable likehood minimum="<<glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].at(2)<<" H_full likehood minimum="<<glbstuff.fitresults[fitalgindex].at(4)<<", pvalue from LRT="<<glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].at(3)<<endl;
      }else{
        cout<<"WARNING:optimizeLikehood:something wrong happend, I'm trying to do the LR test, but glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)] has a wrong size: "<<glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].size()<<"  fitalgindex="<<fitalgindex<<"  fpMinimizer->VariableName(fixedpar.first)="<<fpMinimizer->VariableName(fixedpar.first)<<" I'll overwrite everything, but be aware that I don't know what is happening"<<endl;
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].at(2)=fpMinimizer->MinValue(); //append -loglikehood value for fixed par
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(fixedpar.first)].at(3)= ROOT::Math::chisquared_cdf_c(2.* (fpMinimizer->MinValue()-glbstuff.fitresults[fitalgindex].at(4)),1); //append pvalue        
      }
    }
  }

  if(debug)
    cout<<"optimizeLikehood done, fitalgindex="<<fitalgindex<<endl;

  if(status==0){
    vector<double> cov(fpMinimizer->NDim() * fpMinimizer->NDim());
    fpMinimizer->GetCovMatrix(cov.data());
    optlike_fill(sample, glbstuff,fitalgindex, samrect,cov); // TODO: need to be fixed
    pair<double,double>aucprecres=optlike_aucROC(sample, glbstuff, fitalgindex);
    cout<<"fitalgo with "<<glbstuff.fitalgo.at(fitalgindex).first<<" and "<<glbstuff.fitalgo.at(fitalgindex).second<<" done. AUC="<<aucprecres.first<<"  average_precision="<<aucprecres.second<<endl;
    SetAucAvgPrec(fitalgindex, aucprecres, glbstuff); 
    PlotCalibrationCurveQuantilesAndHLtest(sample, glbstuff, fitalgindex, 3);    
    if(fixedpar.first>0){
      for(auto &par:glbstuff.fitpars)
      if(par.second.second.at(1)!=0)//skip fixed pars
          optimizeLikehood(sample, glbstuff, fitalgindex, samrect, make_pair(par.second.first,0.)); //for Likelihood ratio test, dove parametri in più sono messi a 0 pef fare contfronto. 
        // for(auto &par:fitpars) //TODO: to be developed
        //   optimizeLikehood(sample, glbstuff, i, make_pair(par.second.first,par.second)); //for Likelihood profile, dove parametri in più sono messi al minimo per fare likelihood profile.
    }
  }

  return status;
}

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=alfabeta
double functorLikehoodFull(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-=std::log(EvalScoreLikehoodFull(paziente.second, par));
  return eval;
}

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=alfabeta
double functorLikehoodAlfabdone(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-=std::log(EvalScoreLikehoodAlfabdone(paziente.second, par)); 
  return eval;
}

//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=alfabeta, 4=clinical factor at 0
double functorLikehoodFullClinical_0(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-=std::log(EvalScoreLikehoodFullClinical_0(paziente.second, par));
    return eval;
}

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=clinical factor at 0
double functorLikehoodAlfabdoneClinical_0(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-= std::log(EvalScoreLikehoodAlfabdoneClinical_0(paziente.second, par));
  return eval;
}
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=alfabeta, 4=clinical factor at 0
double functorLikehoodFullClinical_1(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-=std::log(EvalScoreLikehoodFullClinical_1(paziente.second, par));
    return eval;
}

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=clinical factor at 0
double functorLikehoodAlfabdoneClinical_1(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-= std::log(EvalScoreLikehoodAlfabdoneClinical_1(paziente.second, par));
  return eval;
}

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=clinical factor at 0
double functorLikehoodAlfabdone2DvhClinical_1(const map<int, PatientData> &sample, const map<int, PatientData> &samrect , const double* par){
  double eval=0.;
  for(const auto &pazienta : sample){
    auto it = samrect.find(pazienta.first);
    if(it!=samrect.end())
      eval-= std::log(EvalScoreLikehoodAlfabdone2DvhClinical_1(pazienta.second, samrect.at(pazienta.second.id), par));
    else
      cout<<"functorLikehoodAlfabdone2DvhClinical_1: no matching found for sample.second.id="<<pazienta.second.id<<endl;
  }
  return eval;
}

//likehood:
//par0=alfabeta, par1=dose value
double functorLikehoodDose4Vol(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample){
    if(paziente.second.dvhcumnormmap.at(par[2]).size()>par[3])
      eval-= std::log( EvalScoreLikehoodDose4Vol(paziente.second, par));
  }
  return eval;
}

double functorLikehoodDose4VolClinical(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample){
    if(paziente.second.dvhcumnormmap.at(par[2]).size()>par[3])
      eval-= std::log( EvalScoreLikehoodDose4VolClinical(paziente.second, par));
  }
  return eval;
}


//evaluate eud given an alfabeta and an nvalue
double CalculateEudFromScratch(const PatientData &paziente, double alfabeta, double nvalue){

  if(debug>5)
    cout<<"start CalculateEudFromScratch with paziente.id="<<paziente.id<<"  alfabeta="<<alfabeta<<"  nvalue="<<nvalue<<endl;

    double eud=0;
    for(int i=0;i<paziente.dvhmapdiff.size();i++)
      eud+=pow(((double)i)*(alfabeta+((double)i)/paziente.nfraction)/(alfabeta+2.),1./nvalue)*paziente.dvhmapdiff.at(i);
    eud=pow(eud,nvalue);
    
  if(debug>5)
    cout<<"CalculateEudFromScratch done; eud="<<eud<<"  paziente.nfraction="<<paziente.nfraction<<endl;

    return eud;
}

double CalculateEudEqdAlreadyDone(const PatientData &paziente, double nvalue){

  if(debug)
    cout<<"start  CalculateEudEqdAlreadyDone with paziente.id="<<paziente.id<<endl;

  double eud=0;
  for(int i=0;i<paziente.dvhmapdiff.size();i++)
    eud+=  pow((double)i,1./nvalue)*paziente.dvhmapdiff.at(i);
  eud=pow(eud,nvalue);
    
  if(debug)
    cout<<"CalculateEudEqdAlreadyDone done; eud="<<eud<<endl;

    return eud;
}

//fill paziente.second.dvhdiffnormmap, paziente.second.dvhcumnormmap, paziente.second.eudmap if alfabdone==false
void evaluateEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff){
  
  if(debug)
  cout<<"start evaluateEqdEud"<<endl;
  
  for(auto &paziente : sample){
    for(auto const &asub:glbstuff.alfabeta){      

      //calculate dvhdiffnormmap
      paziente.second.dvhdiffnormmap[asub].clear();
      paziente.second.dvhdiffnormmap[asub].resize(paziente.second.dvhmapdiff.size()/glbstuff.eqd2binwidth+0.5);
      paziente.second.dvhcumnormmap[asub].clear();
      paziente.second.dvhcumnormmap[asub].resize(paziente.second.dvhmapcum.size()/glbstuff.eqd2binwidth+0.5);
      
      vector<pair<double, double>> tmpdiffvec; //eqd dose, volume
      for(int i=0;i<paziente.second.dvhmapdiff.size();i++)
        tmpdiffvec.push_back(make_pair(EqdDose(paziente.second, asub, (double)i),paziente.second.dvhmapdiff.at(i)));//eqd
        
        int afterdiffidx=0;
        for(int i=0;i<paziente.second.dvhdiffnormmap[asub].size();i++){ 
          double dose=glbstuff.eqd2binwidth*(i);
          for(int k=afterdiffidx;k<tmpdiffvec.size();k++){
            if(tmpdiffvec.at(k).first>dose)
            break;
            afterdiffidx=k;
          }
          if(afterdiffidx>0)
            paziente.second.dvhdiffnormmap[asub].at(i)=(tmpdiffvec.at(afterdiffidx).second+tmpdiffvec.at(afterdiffidx-1).second)/2.;
          else
            paziente.second.dvhdiffnormmap[asub].at(i)=tmpdiffvec.at(afterdiffidx).second;
          
          for(auto const &n:glbstuff.nvalue4eud)
            paziente.second.eudmap[make_pair(n,asub)]+=pow(EqdDose(paziente.second, asub, dose),1./n)*paziente.second.dvhdiffnormmap[asub].at(i);
        }
        
        vector<pair<double, double>> tmpcumvec; //eqd dose, volume
        for(int i=0;i<paziente.second.dvhmapcum.size();i++)
          tmpcumvec.push_back(make_pair(EqdDose(paziente.second, asub, (double)i),paziente.second.dvhmapcum.at(i)));//eqd
        int aftercumidx=0;
        for(int i=0;i<paziente.second.dvhcumnormmap.at(asub).size();i++){ 
          double dose=glbstuff.eqd2binwidth*(i);
          for(int k=aftercumidx;k<tmpcumvec.size();k++){
            if(tmpcumvec.at(k).first>dose)
              break;
            aftercumidx=k;
          }
          if(aftercumidx>0)
            paziente.second.dvhcumnormmap.at(asub).at(i)=(tmpcumvec.at(aftercumidx).second+tmpcumvec.at(aftercumidx-1).second)/2.;
          else
            paziente.second.dvhcumnormmap.at(asub).at(i)=tmpcumvec.at(aftercumidx).second;
        }

      for(auto const &n:glbstuff.nvalue4eud)
        paziente.second.eudmap.at(make_pair(n,asub))=pow(paziente.second.eudmap.at(make_pair(n,asub)),n);
    } //end of loop on alfabeta
  }//end of loop on patients

  if(debug)
    cout<<"evaluateEqdEud done"<<endl;
  
  return;
}

//fill paziente.second.dvhdiffnormmap, paziente.second.dvhcumnormmap, paziente.second.eudmap if alfabdone==false
void FillEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff){
  
  if(debug)
  cout<<"start FillEqdEud"<<endl;
  
  for(auto &paziente : sample){

    //calculate dvhdiffnormmap
    paziente.second.dvhdiffnormmap[glbstuff.alfabdone].clear();
    paziente.second.dvhdiffnormmap[glbstuff.alfabdone].resize(paziente.second.dvhmapdiff.size()/glbstuff.eqd2binwidth+0.5);
    paziente.second.dvhcumnormmap[glbstuff.alfabdone].clear();
    paziente.second.dvhcumnormmap[glbstuff.alfabdone].resize(paziente.second.dvhmapcum.size()/glbstuff.eqd2binwidth+0.5);
          
      int afterdiffidx=0;
      for(int i=0;i<paziente.second.dvhdiffnormmap[glbstuff.alfabdone].size();i++){ 
        double dose=glbstuff.eqd2binwidth*(i);
        for(int k=afterdiffidx;k<paziente.second.dvhmapdiff.size();k++){
          if(k>dose)
          break;
          afterdiffidx=k;
        }
        if(afterdiffidx>0)
          paziente.second.dvhdiffnormmap[glbstuff.alfabdone].at(i)=(paziente.second.dvhmapdiff.at(afterdiffidx)+paziente.second.dvhmapdiff.at(afterdiffidx-1))/2.;
        else
          paziente.second.dvhdiffnormmap[glbstuff.alfabdone].at(i)=paziente.second.dvhmapdiff.at(afterdiffidx);
        
        for(auto const &n:glbstuff.nvalue4eud)
          paziente.second.eudmap[make_pair(n,glbstuff.alfabdone)]+=pow(dose,1./n)*paziente.second.dvhdiffnormmap[glbstuff.alfabdone].at(i);
      }
      
      int aftercumidx=0;
      for(int i=0;i<paziente.second.dvhcumnormmap.at(glbstuff.alfabdone).size();i++){ 
        double dose=glbstuff.eqd2binwidth*(i);
        for(int k=aftercumidx;k<paziente.second.dvhmapcum.size();k++){
          if(k>dose)
            break;
          aftercumidx=k;
        }
        if(aftercumidx>0)
          paziente.second.dvhcumnormmap.at(glbstuff.alfabdone).at(i)=(paziente.second.dvhmapcum.at(aftercumidx)+paziente.second.dvhmapcum.at(aftercumidx-1))/2.;
        else
          paziente.second.dvhcumnormmap.at(glbstuff.alfabdone).at(i)=paziente.second.dvhmapcum.at(aftercumidx);
      }

    for(auto const &n:glbstuff.nvalue4eud)
      paziente.second.eudmap.at(make_pair(n,glbstuff.alfabdone))=pow(paziente.second.eudmap.at(make_pair(n,glbstuff.alfabdone)),n);
  }//end of loop on patients

  if(debug)
    cout<<"FillEqdEud done"<<endl;
  
  return;
}

void optlike_fill(map<int, PatientData> &sample, const globalstuff &glbstuff, int fitalgindex, map<int, PatientData> &samrect, vector<double> &cov){

  if(debug)
    cout<<"start optlike_fill"<<endl;
  TString varname=(glbstuff.usedosevar==-1) ? "eud":"dosevol";
  TString name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+ "_"+varname+ "_optlike_All";
  TH1D *hall=new TH1D(name, varname+" value;"+varname+";Number of patients",100, 0., 100);        
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+varname+"_optlike_No";
  TH1D *hno=new TH1D(name, varname+" value;"+varname+";Number of patients",100, 0., 100);        
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"eud_optlike_Yes";
  TH1D *hyes=new TH1D(name, varname+" value;"+varname+";Number of patients",100, 0., 100);        
  
  gDirectory->cd("ntcp_linear");
  TGraph* gr_eud_vs_tox=new TGraph();
  // vector<TGraph*> gr_clfac_eud_vs_tox;
  // for(int i=0;i<glbstuff.clusternum;i++){
  //   TGraph *gr=new TGraph();
  //   name=Form("gr_clfac_eud_vs_tox_%i",i);
  //   gr->SetName(name);
  //   gr->SetTitle("eud vs toxicity optimized NTCP ;Eud;Toxicity");
  //   gr_clfac_eud_vs_tox.push_back(gr);
  // }
  TGraph* gr_eud_vs_pred=new TGraph();
  TGraph2D* gr2d_eud_vs_tox;
  if(glbstuff.twodvh)
    gr2d_eud_vs_tox=new TGraph2D();
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_ntcp_linear_"+varname+"_vs_tox";
  gr_eud_vs_tox->SetName(name);
  gr_eud_vs_tox->SetTitle("eud vs toxicity optimized NTCP ;"+varname+";toxicity");
  if(glbstuff.twodvh){ //non previsto per volume
    name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_2d_ntcp_linear_"+varname+"_vs_tox";
    gr2d_eud_vs_tox->SetName(name);
    gr2d_eud_vs_tox->SetTitle(varname+" vs toxicity optimized NTCP ;"+varname+" A; "+varname+" B;toxicity");
  }
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_ntcp_linear_"+varname+"_vs_pred";
  gr_eud_vs_pred->SetName(name);
  gr_eud_vs_pred->SetTitle(varname+" vs NTCP prediction with optimized NTCP ;"+varname+";prediction");
  vector<TGraph*> gr_eudwithtox_vs_tox; //mi aspetto che numerazione cluster vada da 0 in su
  if(glbstuff.clusternum>0){
    for(int i=0;i<glbstuff.clusternum;i++){
      TGraph* gr=new TGraph((dynamic_cast<TH1D*>(gDirectory->Get("../cluster")))->GetBinContent(i+1));
      name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+Form("_ntcp_linear_"+varname+"_vs_tox_withclinicalfactors%i",i);
      gr->SetName(name);
      gr->SetTitle(Form(varname+" vs toxicity with cluster%i;eud;toxicity",i));
      gr_eudwithtox_vs_tox.push_back(gr);
    }
  }

  if(debug)
    cout<<"optlike_fill: fill plots"<<endl;

  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_ntcp_linear_score_vs_tox";
  int index=0;
  vector<int> clsindex(glbstuff.clusternum,0);
  vector<int> gr_eudwithtox_index(gr_eudwithtox_vs_tox.size(),0);
  for(auto &paziente : sample){
    paziente.second.optlike_eud[fitalgindex]= (glbstuff.usedosevar==-1) ?
    ( (glbstuff.alfabdone<0) ? CalculateEudFromScratch(paziente.second,glbstuff.fittedpar.at(fitalgindex).at("alfabeta").at(0) , glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0)) : CalculateEudEqdAlreadyDone(paziente.second, glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0)) ) :
    ((paziente.second.dvhcumnormmap.at(glbstuff.fittedpar.at(fitalgindex).at("alfabeta").at(0)).size()>glbstuff.fittedpar.at(fitalgindex).at("volume").at(0)) ? 
      paziente.second.dvhcumnormmap.at(glbstuff.fittedpar.at(fitalgindex).at("alfabeta").at(0)).at(glbstuff.fittedpar.at(fitalgindex).at("volume").at(0) ) : -1);

    if(glbstuff.twodvh)
      samrect.at(paziente.second.id).optlike_eud[fitalgindex]= (glbstuff.alfabdone<0) ? CalculateEudFromScratch(samrect.at(paziente.second.id),glbstuff.fittedpar.at(fitalgindex).at("alfabeta").at(0) , glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0)) : CalculateEudEqdAlreadyDone(samrect.at(paziente.second.id), glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0));    
    hall->Fill(paziente.second.optlike_eud.at(fitalgindex));
    if(paziente.second.tgt_acutegitox>0.5)
      hyes->Fill(paziente.second.optlike_eud.at(fitalgindex));
    else
      hno->Fill(paziente.second.optlike_eud.at(fitalgindex));
      
    gr_eud_vs_tox->SetPoint(index, paziente.second.optlike_eud.at(fitalgindex), paziente.second.tgt_acutegitox);
    if(glbstuff.twodvh)
      gr2d_eud_vs_tox->SetPoint(index, paziente.second.optlike_eud.at(fitalgindex), samrect.at(paziente.second.id).optlike_eud.at(fitalgindex), paziente.second.tgt_acutegitox);
    gr_eud_vs_pred->SetPoint(index, paziente.second.optlike_eud.at(fitalgindex), paziente.second.optlike_ntcpscore.at(fitalgindex));
    if(glbstuff.clusternum>0){
      gr_eudwithtox_vs_tox.at(paziente.second.cluster)->SetPoint(gr_eudwithtox_index.at(paziente.second.cluster), paziente.second.optlike_eud.at(fitalgindex), paziente.second.tgt_acutegitox+paziente.second.cluster*0.05);
      gr_eudwithtox_index.at(paziente.second.cluster)++;
    }
    index++;
  }
  gr_eud_vs_tox->Sort();
  gr_eud_vs_pred->Sort();
  for(auto &gr:gr_eudwithtox_vs_tox)
  gr->Sort();
  
  TF1::EAddToList::kNo;
  TF1* sigmoidbestnoclinical_0=nullptr,*sigmoidbestnoclinical_1=nullptr,*sigmoidbestnoclinical_2=nullptr;
  sigmoidbestnoclinical_0=new TF1(Form("sigmoidbestnoclinical_0_%i",glbstuff.fitpars.at("beta_zero").first ),"1./(1.+exp(-[0]-[1]*x))", 0, 100);
  sigmoidbestnoclinical_0->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
  sigmoidbestnoclinical_0->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
  gr_eud_vs_tox->Fit(sigmoidbestnoclinical_0, "BS+","",0,100);
  TGraphErrors *tgrer_sigmoidbestnoclinical_0, *tgrer_sigmoidbestclinical_1, *tgrer_sigmoidbestclinical_2;
  vector<int> covindexxmakeband;
  if(glbstuff.fitpars.at("beta_zero").second.at(1)>0)
  covindexxmakeband.push_back(glbstuff.fitpars.at("beta_zero").first);
  if(glbstuff.fitpars.at("beta_eud_a").second.at(1)>0)
    covindexxmakeband.push_back(glbstuff.fitpars.at("beta_eud_a").first);

  if(covindexxmakeband.size()>0){
    tgrer_sigmoidbestnoclinical_0=MakeBandFromMinimizer(sigmoidbestnoclinical_0, covindexxmakeband, 2, cov, glbstuff,fitalgindex, 200, 1);
    tgrer_sigmoidbestnoclinical_0->SetLineColor(kGreen+3);
    tgrer_sigmoidbestnoclinical_0->SetMarkerColor(kGreen+3);
  }

  if(glbstuff.clinicalfactors>0 && glbstuff.clusternum==3){
    gr_eudwithtox_vs_tox.at(0)->SetMarkerColor(kGreen+3);
    gr_eudwithtox_vs_tox.at(0)->SetMarkerStyle(20);
    sigmoidbestnoclinical_1=new TF1("sigmoidbestnoclinical_1",(glbstuff.prop2dose==1) ? "1./(1.+exp(-[0]-([1]+[2])*x))" : "1./(1.+exp(-[0]-[1]*x-[2]))", 0, 100);
    sigmoidbestnoclinical_1->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
    sigmoidbestnoclinical_1->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
    sigmoidbestnoclinical_1->FixParameter(2, glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0));
    gr_eud_vs_tox->Fit(sigmoidbestnoclinical_1, "BS+","",0,100);
    if(glbstuff.fitpars.at("clinical_factor_0").second.at(1)>0)
    covindexxmakeband.push_back(glbstuff.fitpars.at("clinical_factor_0").first);
    if(covindexxmakeband.size()>0){
      tgrer_sigmoidbestclinical_1=MakeBandFromMinimizer(sigmoidbestnoclinical_1, covindexxmakeband, 3, cov, glbstuff,fitalgindex, 210, 1);
      tgrer_sigmoidbestclinical_1->SetLineColor(kBlue);
      tgrer_sigmoidbestclinical_1->SetMarkerColor(kBlue);
      // tgrer_sigmoidbestclinical_1->Draw("LCE same");
    }
    gr_eudwithtox_vs_tox.at(1)->SetMarkerColor(kBlue);
    gr_eudwithtox_vs_tox.at(1)->SetMarkerStyle(21);
    sigmoidbestnoclinical_2=new TF1("sigmoidbestnoclinical_2",(glbstuff.prop2dose==1) ? "1./(1.+exp(-[0]-([1]+[2])*x))" :  "1./(1.+exp(-[0]-[1]*x-[2]))", 0, 100);
    sigmoidbestnoclinical_2->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
    sigmoidbestnoclinical_2->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
    sigmoidbestnoclinical_2->FixParameter(2, (glbstuff.clinicalfactors==2) ? glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_1").at(0) : glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0)*2.);
    gr_eud_vs_tox->Fit(sigmoidbestnoclinical_2, "BS+","",0,100);    
    if(glbstuff.clinicalfactors==2)
      covindexxmakeband.at(2)=glbstuff.fitpars.at("clinical_factor_1").first;
    if(covindexxmakeband.size()>0){
      tgrer_sigmoidbestclinical_2=MakeBandFromMinimizer(sigmoidbestnoclinical_2, covindexxmakeband, 3, cov, glbstuff,fitalgindex, 220, 1);
      tgrer_sigmoidbestclinical_2->SetLineColor(kRed);
      tgrer_sigmoidbestclinical_2->SetMarkerColor(kRed);
      // tgrer_sigmoidbestclinical_2->Draw("LCE same");
    }
    gr_eudwithtox_vs_tox.at(2)->SetMarkerColor(kRed);
    gr_eudwithtox_vs_tox.at(2)->SetMarkerStyle(22);
  }
  gr_eud_vs_tox->SetMarkerStyle(24);
  gr_eud_vs_tox->SetLineWidth(0);
  gr_eud_vs_tox->SetDrawOption("AP*");
  gr_eud_vs_tox->Write();
  TCanvas *canvas_ntcps_eud_vs_tox = new TCanvas("canvas_ntcps_eud_vs_tox", "ntcp curves", 800, 600);
  // gr_eud_vs_tox->GetListOfFunctions()->Clear();
  if(!(glbstuff.clinicalfactors>0 && glbstuff.clusternum==3)){
    gr_eud_vs_tox->Draw("AP*");
    gr_eud_vs_tox->GetXaxis()->SetLimits(0, (glbstuff.usedosevar==-1) ? 1:40);   // solo asse X
    gr_eud_vs_tox->SetMinimum(-0.2);
    gr_eud_vs_tox->SetMaximum(1.2);         
  }else{
    gr_eudwithtox_vs_tox.at(0)->Draw("AP");    
    gr_eudwithtox_vs_tox.at(0)->GetXaxis()->SetLimits(0, (glbstuff.usedosevar==-1) ? 1:40);   // solo asse X
    gr_eudwithtox_vs_tox.at(0)->SetMinimum(-0.2);
    gr_eudwithtox_vs_tox.at(0)->SetMaximum(1.2);         

  }
  sigmoidbestnoclinical_0->SetLineColor(kGreen+3);
  tgrer_sigmoidbestnoclinical_0->Draw("E same");
  sigmoidbestnoclinical_0->Draw("same");
  if(glbstuff.clinicalfactors>0 && glbstuff.clusternum==3){
    for(int i=1;i<glbstuff.clusternum;i++)
      gr_eudwithtox_vs_tox.at(i)->Draw("P same");
    tgrer_sigmoidbestclinical_1->Draw("E same");
    tgrer_sigmoidbestclinical_2->Draw("E same");
    sigmoidbestnoclinical_1->SetLineColor(kBlue);
    sigmoidbestnoclinical_1->Draw("same");
    sigmoidbestnoclinical_2->SetLineColor(kRed);
    sigmoidbestnoclinical_2->Draw("same");
  }
  canvas_ntcps_eud_vs_tox->Write();
  delete canvas_ntcps_eud_vs_tox;
  
  gr_eud_vs_pred->SetMarkerStyle(10);
  gr_eud_vs_pred->SetMarkerColor(3);
  gr_eud_vs_pred->SetLineWidth(2);
  gr_eud_vs_pred->SetLineColor(3);
  gr_eud_vs_pred->SetDrawOption("AP*");
  gr_eud_vs_pred->Write();

  //TODO: TO BE FIXED 
  // if(glbstuff.twodvh){
  //   TF2* sigmoid2dbestnoclinical_0=new TF2("sigmoid2dbestnoclinical_0","1./(1.+exp(-[0]-[1]*x-[2]*y))", 0, 1000, 0, 1000);
  //   sigmoid2dbestnoclinical_0->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
  //   sigmoid2dbestnoclinical_0->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
  //   sigmoid2dbestnoclinical_0->FixParameter(2, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_b").at(0));
  //   sigmoid2dbestnoclinical_0->Write();
  //   gr2d_eud_vs_tox->Fit(sigmoid2dbestnoclinical_0, "BS+","");
  //   if(glbstuff.clinicalfactors>0 && glbstuff.clusternum==3){
  //     TF2* sigmoid2dbestnoclinical_1=new TF2("sigmoid2dbestnoclinical_1","1./(1.+exp(-[0]-[1]*x-[2]*y-[3]))", 0, 1000, 0, 1000);
  //     sigmoid2dbestnoclinical_1->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
  //     sigmoid2dbestnoclinical_1->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
  //     sigmoid2dbestnoclinical_1->FixParameter(2, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_b").at(0));
  //     sigmoid2dbestnoclinical_1->FixParameter(3, glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0));
  //     sigmoid2dbestnoclinical_1->Write();
  //     gr2d_eud_vs_tox->Fit(sigmoid2dbestnoclinical_1, "BS+","");
  //     TF2* sigmoid2dbestnoclinical_2=new TF2("sigmoid2dbestnoclinical_2","1./(1.+exp(-[0]-[1]*x-[2]*y -[3]))", 0, 1000, 0, 1000);
  //     sigmoid2dbestnoclinical_2->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
  //     sigmoid2dbestnoclinical_2->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
  //     sigmoid2dbestnoclinical_2->FixParameter(2, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_b").at(0));
  //     sigmoid2dbestnoclinical_2->FixParameter(3, (glbstuff.clinicalfactors==2) ? glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_1").at(0) : glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0)*2.);
  //     sigmoid2dbestnoclinical_2->Write();
  //     gr2d_eud_vs_tox->Fit(sigmoid2dbestnoclinical_2, "BS+","");
  //     gr2d_eud_vs_tox->GetListOfFunctions()->Print();
  //   }
  //   gr2d_eud_vs_tox->SetMarkerStyle(20);
  //   gr2d_eud_vs_tox->SetMarkerColor(2);
  //   gr2d_eud_vs_tox->SetLineWidth(0);
  //   gr2d_eud_vs_tox->SetLineColor(0);
  //   gr2d_eud_vs_tox->SetDrawOption("AP*");
  //   gr2d_eud_vs_tox->Write();    
  // }

  if(debug)
    cout<<"optlike_fill: fit gr_eudwithtox_vs_tox"<<endl;

  if(glbstuff.clinicalfactors>0){ 
    for(int i=0;i<gr_eudwithtox_vs_tox.size();i++){
      TF1* sigmoidfromlikehood=new TF1("sigmoidfromlikehood",  (glbstuff.prop2dose==1) ? "1./(1.+exp(-[0]-([1]+[2])*x))" : "1./(1.+exp(-[0]-[1]*x-[2]))", 0, 100);
      sigmoidfromlikehood->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
      sigmoidfromlikehood->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud_a").at(0));
      sigmoidfromlikehood->FixParameter(2, (glbstuff.clinicalfactors==1) ? glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0)*i   :    
        (  (i==2) ? glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_1").at(0) :  (  (i==1) ? glbstuff.fittedpar.at(fitalgindex).at("clinical_factor_0").at(0) : 0   ) )  );
      gr_eudwithtox_vs_tox.at(i)->Fit(sigmoidfromlikehood,"BS+","",0,100);
      gr_eudwithtox_vs_tox.at(i)->SetMarkerStyle(20);
      gr_eudwithtox_vs_tox.at(i)->SetMarkerColor(2);
      gr_eudwithtox_vs_tox.at(i)->SetLineWidth(0);
      gr_eudwithtox_vs_tox.at(i)->SetLineColor(0);
      gr_eudwithtox_vs_tox.at(i)->SetDrawOption("AP*");
      gr_eudwithtox_vs_tox.at(i)->Write();
    }
  }  
  gDirectory->cd("..");

  if(debug)
    cout<<"optlike_fill done"<<endl;

  return;
}





pair<double,double> optlike_aucROC(const map<int, PatientData> &sample, const globalstuff &glbstuff, const int fitalgindex)
{
  map<double, vector<int>, std::greater<double>> data;

  int nPos = 0;
  int nNeg = 0;

  for(const auto &paziente:sample){
    if(paziente.second.tgt_acutegitox>0.5)
      nPos++;
    else
      nNeg++;     
    data[paziente.second.optlike_ntcpscore.at(fitalgindex)].push_back(paziente.second.tgt_acutegitox);
  }

  if (nPos == 0 || nNeg == 0)
      return make_pair(-1,-1);

  double tp = 0.0;
  double fp = 0.0;
  double prevTPR = 0.0;
  double prevFPR = 0.0;
  double auc = 0.0, avg_precision=0.; //AUC is the area under the curve given by the point, avg_precision is the actual area under each point as a step function
  int ip=0;
  TGraph* grROC = new TGraph();
  TString pltname="roccurve_"+glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second;
  grROC->SetName(pltname);
  grROC->SetPoint(ip++, 0.0, 0.0);
  for (const auto& kv : data) {
    for (int label : kv.second) {
      if (label == 1) 
        tp++;
      else 
        fp++;
    }
    double TPR = tp /nPos;
    double FPR = fp /nNeg;
    auc += (FPR - prevFPR) * (TPR + prevTPR) / 2.0; // area trapezoidale
    avg_precision+=  tp/(tp+fp)*(TPR-prevTPR);
    grROC->SetPoint(ip++, FPR, TPR);
    prevTPR = TPR;
    prevFPR = FPR;
  }
  grROC->SetTitle(Form("ROC curve for likehood optmized simple NTCP model with AUC=%f;False Positive Rate;True Positive Rate",auc));
  grROC->SetLineWidth(2);
  grROC->SetMarkerStyle(20);
  TF1* bis = new TF1("roc_bisector", "x", 0.0, 1.0);
  bis->SetLineColor(kRed);
  bis->SetLineStyle(2);
  bis->SetLineWidth(2);
  grROC->GetListOfFunctions()->Add(bis);
  grROC->Draw("ALP");
  grROC->Write();

  return make_pair(auc,avg_precision);
}


// it's only a likelihood slice, not a real likelihood profile in which each variable is minimized 
void DrawLikeHood(std::map<int, PatientData>& sample, const globalstuff& glbstuff){

  if(debug)
    cout<<"start DrawLikeHood"<<endl;

  if(!glbstuff.fittedpar.count(glbstuff.bestvalue)){
    std::cerr<<"DrawLikeHood: bestvalue not found"<<std::endl;
    return;
  }

  std::map<int, std::string> orderedNames;
  for(const auto &p : glbstuff.fitpars)
    orderedNames[p.second.first]=p.first; //key=index, value=fitted parameter name

  std::vector<double> pars(orderedNames.size());
  std::vector<double> errs(orderedNames.size());

  for(const auto &p : orderedNames){
    auto fitit=glbstuff.fittedpar.at(glbstuff.bestvalue).find(p.second);
    if(fitit==glbstuff.fittedpar.at(glbstuff.bestvalue).end()){
      cout<<"ERROR: DrawLikeHood: parameter of glbstuff.fitpars not found in glbstuff.fittedpar: "<<p.second<<endl;
      continue;  
    }

    pars[p.first]=fitit->second.at(0);
    errs[p.first]=fitit->second.at(1);
  }

  for(auto it1=orderedNames.begin();it1!=orderedNames.end();++it1){
    if (glbstuff.fitpars.at(it1->second).second.at(1)==0) //it's a fixed parameter
      continue;
    for(auto it2=std::next(it1);it2!=orderedNames.end();++it2){
      if (glbstuff.fitpars.at(it2->second).second.at(1)==0)
        continue;

      double ex=errs[it1->first];
      double ey=errs[it2->first];

      if(!std::isfinite(ex) || ex<=0.) ex=std::abs(pars[it1->first])*0.2+1.;
      if(!std::isfinite(ey) || ey<=0.) ey=std::abs(pars[it2->first])*0.2+1.;
      double xmin=std::max(pars[it1->first]-3.*ex, glbstuff.fitpars.at(it1->second).second.at(2));
      double xmax=std::min(pars[it1->first]+3.*ex, glbstuff.fitpars.at(it1->second).second.at(3));

      double ymin=std::max(pars[it2->first]-3.*ey, glbstuff.fitpars.at(it2->second).second.at(2));
      double ymax=std::min(pars[it2->first]+3.*ey, glbstuff.fitpars.at(it2->second).second.at(3));      

      TH2D* h2=new TH2D(("hlike_"+it1->second+"_vs_"+it2->second).c_str(),("Likelihood scan "+it1->second+" vs "+it2->second+";"+it1->second+";"+it2->second+";NLL").c_str(),100,xmin, xmax, 100, ymin, ymax);
      std::vector<double> trial=pars;

      for(int ix=1;ix<=h2->GetNbinsX();ix++){
        trial[it1->first]=h2->GetXaxis()->GetBinCenter(ix);

        for(int iy=1;iy<=h2->GetNbinsY();iy++){
          trial[it2->first]=h2->GetYaxis()->GetBinCenter(iy);

          double val=functorSelector(glbstuff, sample, trial.data());

          if(!std::isfinite(val)) val=9999;

          h2->SetBinContent(ix, iy, val);
        }
      }
    }
  }

return;
}


void PlotCalibrationCurveQuantilesAndHLtest(const std::map<int, PatientData>& sample, const globalstuff& glbstuff, int fitalgindex, int nbins = 4){

  if(debug)
    cout<<"start PlotCalibrationCurveQuantiles"<<endl;

  vector<pair<double,double>> data; // pred, obs
  for(const auto& paziente : sample){
      double pred = paziente.second.optlike_ntcpscore.at(fitalgindex);
      double obs  = paziente.second.tgt_acutegitox;

      if((pred > 1.0) || (pred < 0. && glbstuff.usedosevar==-1 ) ){
        cout<<"ERROR:PlotCalibrationCurveQuantiles: predicted ntcpscore is "<<pred<<endl;
        continue;
      }
      data.push_back({pred, obs});
  }

  std::sort(data.begin(), data.end(),[](const auto& a, const auto& b){return a.first < b.first;});
  int n = data.size();

  if(n == 0){
    cout<<"ERROR:PlotCalibrationCurveQuantiles: number of data is =0"<<endl;
    return;
  } 
  if(nbins > n) 
    nbins = n;
  double hosmerlemeshow_test=0.;
  std::vector<double> x, y, ex, ey;
  for(int b = 0; b < nbins; b++){
      int start = b * n / nbins;
      int end   = (b + 1) * n / nbins;
      double sum_pred = 0.0;
      double sum_obs  = 0.0;
      int count = end - start;
      if(count <= 0) continue;
      for(int i = start; i < end; i++){
          sum_pred += data[i].first;
          sum_obs  += data[i].second;
      }
      double obs_rate  = sum_obs  / count;
      ex.push_back(0.5*(data[end-1].first - data[start].first));
      x.push_back(sum_pred / count);
      y.push_back(obs_rate);
      ey.push_back(std::sqrt(obs_rate * (1.0 - obs_rate) / count));//TODO: errore efficienza si può migliorare usando formula modificagta, da trovare 

      double term_events= (sum_pred>0) ? (sum_obs-sum_pred)*(sum_obs-sum_pred)/sum_pred : 0.;
      double term_noevents= ((count-sum_pred) > 0 ) ? (sum_pred-sum_obs)*(sum_pred-sum_obs)/(count-sum_pred) : 0.;
      hosmerlemeshow_test+=term_events+term_noevents;
  }

  int ndof=nbins-2;
  if(ndof<=0){
    cout<<"ERROR:PlotCalibrationCurveQuantilesAndHLtest: cannot calculate Hosmer-Lemeshow test since ndof="<<ndof<<", please change the number of bins"<<endl;
  }else{
    double pvalue=ROOT::Math::chisquared_cdf_c(hosmerlemeshow_test, ndof);
    cout<<"Hosmer-Lemeshow test="<<hosmerlemeshow_test<<",  ndof="<<ndof<<",  pvalue="<<pvalue<<endl;
    if(pvalue<0.05)
      cout<<"WARNING:PlotCalibrationCurveQuantilesAndHLtest: Hosmer-Lemeshow test pvalue is suspiciously low, check the fit"<<endl;
  }

  TString pltname="calib_"+glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second;
  // TCanvas* c = new TCanvas(pltname, "Calibration quantiles", 800, 800);
  TGraphErrors* gr = new TGraphErrors(x.size(), x.data(), y.data(),ex.data(),ey.data());
  pltname="calib_"+glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second;
  gr->SetName(pltname);
  gr->SetTitle(";Mean predicted NTCP;Observed toxicity rate");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetLimits(0.0, 1.0);
  gr->GetYaxis()->SetRangeUser(0.0, 1.0);
  TF1* bis = new TF1("roc_bisector", "x", 0.0, 1.0);
  bis->SetLineColor(kRed);
  bis->SetLineStyle(2);
  bis->SetLineWidth(2);
  gr->GetListOfFunctions()->Add(bis);
  gr->Draw("APLE");
  gr->Write();
  // TLine* line = new TLine(0.0, 0.0, 1.0, 1.0);
  // line->SetLineStyle(2);
  // line->SetLineWidth(2);
  // line->Draw("same");
  // c->Write();
  // c->SetGrid();
  // c->SaveAs(outname.c_str());

  return;
}

void SetAucAvgPrec(int index, const pair<double,double> aucavgin, globalstuff& glbstuff){
  if(glbstuff.fitresults.at(index).size()==7){
    glbstuff.fitresults.at(index).push_back(aucavgin.first);
    glbstuff.fitresults.at(index).push_back(aucavgin.second);
    return;
  }
  if(glbstuff.fitresults.at(index).size()==9){
    cout<<"SetAucAvgPrec: WARNING; I'm overwriting the previous auc and avg values, be aware if this is correct or not... previous auc="<<glbstuff.fitresults.at(index).at(7)<<"  new auc value="<<aucavgin.first<<"  previous avg="<<glbstuff.fitresults.at(index).at(8)<<"  new avg="<<aucavgin.second<<endl;
    glbstuff.fitresults.at(index).at(7)=aucavgin.first;
    glbstuff.fitresults.at(index).at(8)=aucavgin.second;
    return;
  }

  cout<<"ERROR: SetAucAvgPrec: length of glbstuff.fitresults.at(index) is not correct, I don't know what happened, but I cannot save correctly auc and avg since this vector has a length of "<<glbstuff.fitresults.at(index).size()<<"  in principle if auc/avg are not stored yet it should be of 7, or 9 if they already has been saved"<<endl;

return;
}


int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname, const globalstuff &glbstuff){

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
    
    int tgtpos=-1, powp_pos, doseptv_pos=-1, nfrac_pos=-1, meandoserectum_pos=-1, doseperfraction_pos=-1, prostatectomy_pos=-1, appendectomy_pos=-1, age_pos=-1,bmi_pos=-1,smoke_pos=-1,alcohol_pos=-1,diabetes_pos=-1,diverticulitis_pos=-1, psa_pos=-1, tumour_risk_pos=-1, lymph_irr_pos=-1, semivesicle_irr_pos=-1,mb_risk_pos=-1, cluster_pos=-1;
    TString doseptv_name("dose to the prostate PTV (Gy)");
    TString powp_name("WPRT");
    TString nfrac_name("fraction number");
    TString meandoserectum_name("mean dose to the rectum (Gy)");
    TString doseperfraction_name("dose/fraction");
    TString prostatectomy_name("Prostatectomy");
    TString appendectomy_name("Appendectomy");
    TString age_name("Age");
    TString bmi_name("BMI");
    TString smoke_name("Smoke");
    TString alcohol_name("Alcohol");
    TString diabetes_name("Diabetes");
    TString diverticulitis_name("Diverticulitis");
    TString psa_name("PSA at diagnosis (ng/ml)");
    TString tumour_risk_name("tumour risk class");
    TString lymph_irr_name("lymph node irradiation");
    TString semivesicle_irr_name("semil vesicle irradiation");
    TString mb_risk_name("MB risk class");
    TString cluster_name((glbstuff.datatype==4) ? "cluster2" : "cluster");
    for(int i=0;i<parts.size();i++){
      if(tgtname.CompareTo(parts[i])==0)
        tgtpos=i;
      if(doseptv_name.CompareTo(parts[i])==0)
        doseptv_pos=i;
      if(powp_name.CompareTo(parts[i])==0)
        powp_pos=i;
      if(nfrac_name.CompareTo(parts[i])==0)
        nfrac_pos=i;
      if(meandoserectum_name.CompareTo(parts[i])==0)
        meandoserectum_pos=i;
      if(doseperfraction_name.CompareTo(parts[i])==0)
        doseperfraction_pos=i;
      if(prostatectomy_name.CompareTo(parts[i])==0)
        prostatectomy_pos=i;
      if(appendectomy_name.CompareTo(parts[i])==0)
        appendectomy_pos=i;
      if(age_name.CompareTo(parts[i])==0)
        age_pos=i;
      if(bmi_name.CompareTo(parts[i])==0)
        bmi_pos=i;
      if(smoke_name.CompareTo(parts[i])==0)
        smoke_pos=i;
      if(alcohol_name.CompareTo(parts[i])==0)
        alcohol_pos=i;
      if(diabetes_name.CompareTo(parts[i])==0)
        diabetes_pos=i;
      if(diverticulitis_name.CompareTo(parts[i])==0)
        diverticulitis_pos=i;
      if(psa_name.CompareTo(parts[i])==0)
        psa_pos=i;
      if(tumour_risk_name.CompareTo(parts[i])==0)
        tumour_risk_pos=i;
      if(lymph_irr_name.CompareTo(parts[i])==0)
        lymph_irr_pos=i;
      if(semivesicle_irr_name.CompareTo(parts[i])==0)
        semivesicle_irr_pos=i;
      if(mb_risk_name.CompareTo(parts[i])==0)
        mb_risk_pos=i;
      if(cluster_name.CompareTo(parts[i])==0)
        cluster_pos=i;
    }
    if(tgtpos==-1 || doseptv_pos==-1 || nfrac_pos==-1 || meandoserectum_pos==-1 || doseperfraction_pos==-1 || prostatectomy_pos==-1 ||appendectomy_pos==-1 || age_pos==-1 ||bmi_pos==-1 ||smoke_pos==-1 ||alcohol_pos==-1 ||diabetes_pos==-1 ||diverticulitis_pos==-1 || psa_pos==-1 || tumour_risk_pos==-1 || lymph_irr_pos==-1 || semivesicle_irr_pos==-1 ||mb_risk_pos==-1 || cluster_pos==-1 || powp_pos==-1) {
      throw std::runtime_error(Form("loadMetaFile: error, some variable not found: tgtpos=%i  powp_pos=%i  doseptv_pos=%i  nfrac_pos=%i  meandoserectum_pos=%i  doseperfraction_pos=%i , prostatectomy_pos=%i, appendectomy_pos=%i, age_pos=%i,bmi_pos=%i,smoke_pos=%i,alcohol_pos=%i,diabetes_pos=%i,diverticulitis_pos=%i, psa_pos=%i, tumour_risk_pos=%i, lymph_irr_pos=%i, semivesicle_irr_pos=%i,mb_risk_pos=%i, cluster_pos=%i", tgtpos, powp_pos, doseptv_pos,  nfrac_pos,  meandoserectum_pos, doseperfraction_pos, prostatectomy_pos, appendectomy_pos, age_pos,bmi_pos,smoke_pos,alcohol_pos,diabetes_pos,diverticulitis_pos, psa_pos, tumour_risk_pos, lymph_irr_pos, semivesicle_irr_pos,mb_risk_pos, cluster_pos));
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
    //stoi->int, stod->float
    auto it=sample.find(std::stoi(readed.substr(readed.size() - 4)));//trovo il paziente in mappa
    if(it!=sample.end()){
      it->second.tgt_acutegitox=std::stoi(trim(parts[tgtpos]));
      it->second.dose_ptv=std::stod(trim(parts[doseptv_pos]));
      it->second.nfraction=std::stoi(trim(parts[nfrac_pos]));
      it->second.dose_per_fraction=std::stod(trim(parts[doseperfraction_pos]));
      it->second.mean_dose_rectum=std::stod(trim(parts[meandoserectum_pos]));
      
      it->second.powp= trim(parts[powp_pos]) == "WPRT" ? 1 : 0;
      it->second.prostatectomy= trim(parts[prostatectomy_pos]) == "Prostatectomy" ? 1 : 0;
      it->second.appendectomy= trim(parts[appendectomy_pos]) == "Appendectomy" ? 1:0;
      it->second.age=std::stoi(trim(parts[age_pos]));
      it->second.bmi=std::stod(trim(parts[bmi_pos]));
      it->second.smoke= (trim(parts[bmi_pos])=="Smoker");
      it->second.alcohol= (trim(parts[bmi_pos])=="Alcohol use");
      it->second.diabetes= (trim(parts[diabetes_pos])=="Diabetes");
      it->second.diverticulitis= (trim(parts[diverticulitis_pos])=="Diverticulitis");
      it->second.psa= std::stod(trim(parts[psa_pos]));
      it->second.tumour_risk= trim(parts[tumour_risk_pos])=="oligometastatic" ? 4: trim(parts[tumour_risk_pos])=="high" ? 3 : trim(parts[tumour_risk_pos])=="intermediate" ? 2 : trim(parts[tumour_risk_pos])=="low" ? 1 : -1;
      it->second.lymph_irr=std::stoi(trim(parts[lymph_irr_pos]));
      it->second.semivesicle_irr=std::stoi(trim(parts[semivesicle_irr_pos]));
      it->second.mb_risk= trim(parts[mb_risk_pos])=="high.risk" ? 2: trim(parts[mb_risk_pos])=="moderate.risk" ? 1 : trim(parts[mb_risk_pos])=="low.risk" ? 0 : -1;
      it->second.cluster=std::stoi(trim(parts[cluster_pos]))-1; //-1 added so the cluster index start from 0 instead of 1
      addedinmap++;
    }
  }
  in.close();

  if(glbstuff.datatype==2 || glbstuff.datatype==4){//nanoport data, remove non nanoport patients
    for(auto it = sample.begin(); it != sample.end();){
      if(it->second.cluster <0 )
        it = sample.erase(it);
      else
        ++it;
    }
  }else if(glbstuff.datatype==3){//old article clustering with MB class risk
    for(auto it = sample.begin(); it != sample.end();){
      it->second.cluster=it->second.mb_risk;
      if(it->second.cluster <0 )
        it = sample.erase(it);
      else
        ++it;
    }
  }
  
  //wprt / port selection
  if(glbstuff.powptype!=-1){
    for(auto it = sample.begin(); it != sample.end();){
      if(it->second.powp!=glbstuff.powptype)
        it = sample.erase(it);
      else
        ++it;
    }
  }
  
  if(addedinmap!=sample.size())
    cout<<"WARNING: mismatch between addedinmap="<<addedinmap<<"  and  sample.size()="<<sample.size()<<endl;
  cout<<"loadMetaFile done, target aquired target="<<tgtname<<"  overall sample size="<<sample.size()<<endl;

  return 0;
}

int SetClusterAsClinicalFactor(map<int, PatientData> &sample, const globalstuff &glbstuff){
  for(auto &paziente : sample){
    if(paziente.second.cluster>=0 && paziente.second.clinical_factor.size()==0){
      if(glbstuff.clinicalfactors==1)
        paziente.second.clinical_factor.push_back(paziente.second.cluster); //valore cluster è valore di clinical factor
      if(glbstuff.clinicalfactors==2){//WARNING: SE SONO DUE VALORI SI PRESUPPONE SIANO DA CLUSTER E CHE SIANO NORMALIZZATI: cluster=0->(0,0) cluster=1->(1,0), cluster=2->(0,1)
        paziente.second.clinical_factor.push_back( (paziente.second.cluster==1) ? 1 : 0  );
        paziente.second.clinical_factor.push_back( (paziente.second.cluster==2) ? 1 : 0  );
      }
    }else{
      cout<<"ERROR in SetClusterAsClinicalFactor: paziente.cluster="<<paziente.second.cluster<<" paziente.clinical_factor.size()="<<paziente.second.clinical_factor.size()<<endl;
      return 1;
    }
  }
  
  return 0;
}

void fillGlobalStuff(globalstuff &glbstuff, double alfabdone, double eqd2binwidth, const vector<double> &nvalue4eud, const vector<double> &alfabeta, const map<string, pair<int,vector<double>>> &fitpars,   const vector<pair<string,string>> &fitalgo, int datatype, int clinicalfactors, int clusternum, int powptype, int twodvh, int prop2dose, const vector<double> &doses4volume, int usedosevar){
  glbstuff.alfabdone=alfabdone;
  glbstuff.eqd2binwidth=eqd2binwidth;
  glbstuff.nvalue4eud=nvalue4eud;
  glbstuff.alfabeta=alfabeta;
  glbstuff.fitpars=fitpars;
  glbstuff.fitalgo=fitalgo;
  glbstuff.datatype=datatype;
  glbstuff.clinicalfactors=clinicalfactors;
  glbstuff.clusternum=clusternum;
  glbstuff.maxbin= (datatype==1) ? 200:100;
  glbstuff.powptype=powptype;
  glbstuff.twodvh=twodvh;
  glbstuff.prop2dose=prop2dose;
  glbstuff.doses4volume=doses4volume;
  glbstuff.usedosevar=usedosevar;
  return;
}

// clsscale: 1 sigma; usa 1.96 per ~95%
TGraphErrors *MakeBandFromMinimizer(TF1 *f, vector<int> &cov_indices , const int npar, const vector<double> &cov, const globalstuff &glbstuff, int fitalgindex, int npoints, double clscale) {

  if(debug>3)
    cout<<"start MakeBandFromMinimizer with f="<<f->GetName()<<endl;
  
  int covndim = static_cast<int>(std::sqrt(cov.size()));
  if (covndim * covndim != static_cast<int>(cov.size())){
    throw std::runtime_error("La dimensione del vettore non corrisponde a una matrice quadrata.");
  }    
  for (int idx : cov_indices) {
    if (idx < 0 || idx >= covndim) {
      throw std::runtime_error("Indice parametro fuori range rispetto alla matrice di covarianza.");
    }
  }

  vector<double> cov2(npar*npar);
  for (int i = 0; i < npar; ++i) {
    for (int j = 0; j < npar; ++j) {
      cov2[i * npar + j] = cov[cov_indices[i] * covndim + cov_indices[j]];
    }
  }

  for (int i = 0; i < npar; i++)
    f->SetParameter(i, glbstuff.fittedpar.at(fitalgindex).at(glbstuff.nameindex.at(cov_indices[i])).at(0));
  if(std::strcmp(f->GetName(), "sigmoidbestnoclinical_2") == 0 && glbstuff.clinicalfactors==1)
    f->SetParameter(2,f->GetParameter(2)*2);
  

  TGraphErrors *band = new TGraphErrors(npoints);
  double xmin = f->GetXmin();
  double xmax = f->GetXmax();

  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    double x = xmin + (xmax - xmin) * ipoint / double(npoints - 1);
    double y = f->Eval(x);
    std::vector<double> grad(npar);
    for (int ipar = 0; ipar < npar; ++ipar) {
      double p0 = f->GetParameter(ipar);
      int ipar_global = cov_indices[ipar];
      double h = 0.1 * glbstuff.fittedpar.at(fitalgindex).at(glbstuff.nameindex.at(ipar_global)).at(1);
      if (h == 0.0)
          h = 1e-5 * std::max(1.0, std::abs(p0));
      f->SetParameter(ipar, p0 + h);
      if(std::strcmp(f->GetName(), "sigmoidbestnoclinical_2") == 0 && glbstuff.clinicalfactors==1)
        f->SetParameter(2,f->GetParameter(2)*2);
      double yp = f->Eval(x);
      f->SetParameter(ipar, p0 - h);
      if(std::strcmp(f->GetName(), "sigmoidbestnoclinical_2") == 0 && glbstuff.clinicalfactors==1)
        f->SetParameter(2,f->GetParameter(2)*2);
      double ym = f->Eval(x);
      f->SetParameter(ipar, p0);
      if(std::strcmp(f->GetName(), "sigmoidbestnoclinical_2") == 0 && glbstuff.clinicalfactors==1)
        f->SetParameter(2,f->GetParameter(2)*2);
      grad[ipar] = (yp - ym) / (2.0 * h);
    }
    double var = 0.0;
    for (int a = 0; a < npar; ++a) {
      for (int b = 0; b < npar; ++b) {
        var += grad[a] * cov2[a * npar + b] * grad[b];
      }
    }
    double ey = clscale * std::sqrt(std::max(0.0, var));
    band->SetPoint(ipoint, x, y);
    band->SetPointError(ipoint, 0.0, ey);
  }

  if(debug>3)
    cout<<"MakeBandFromMinimizer done"<<endl;

    return band;
}