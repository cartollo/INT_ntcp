#include "myntcpfunctions.h"

extern int debug;
void bookHisto(TFile *outrootfile,  const globalstuff &glbstuff, const bool multivariate){

  if(debug)
    cout<<"start bookHisto"<<endl;
  
  // if(multivariate)

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
        xv(counter)=paziente.second.dvhcumnormmap.at(asub).at( (dynamic_cast<TH1D*>(gDirectory->Get(Form("dvhcumnorm_tgt_%.3f_Diff",asub))))->FindBin(mygaus.GetParameter(1)));
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

  TDirectory *dir=gDirectory->GetDirectory("ntpc_linear");
  if(dir){
    gDirectory->cd("ntpc_linear");
  }else{
    cout<<"WARNING in CreateNtcpSigmoidalSingle: ntpc_linear not present... I'll create a new one"<<endl;
    gDirectory->mkdir("ntpc_linear");
    gDirectory->cd("ntpc_linear");
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

  gDirectory->mkdir("ntpc_linear");
  gDirectory->cd("ntpc_linear");

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
    for(int i=1; i<patient.dvhmapcum.size();i++)
      patient.dvhmapdiff.push_back(patient.dvhmapcum.at(i-1)-patient.dvhmapcum.at(i));
    computeDCT(patient.dvhmapdiff, patient.dvhdctdiff);
    sample.insert(std::pair{id,patient});
}

  if(debug)
    cout<<"loadDvhFile done, read from "<<filename<<" "<<sample.size()<<" elements"<<endl;
  
    return 0;
}



int optimizeLikehood(map<int, PatientData> &sample, globalstuff &glbstuff, const int fitalgindex, const pair<int,double> fixedpar){
  
  if(debug)
    cout<<"start optimizeLikehood, fitalgindex="<<fitalgindex<<endl;

  //lambda
  auto lamdalikehoodFull = [&](const double* par) {
      return functorLikehoodFull(sample, par);
  };
  auto lamdalikehoodAlfabdone = [&](const double* par) {
      return functorLikehoodAlfabdone(sample, par);
  };
  
  auto lamdalikehoodFullClinical_0 = [&](const double* par) {
    return functorLikehoodFullClinical_0(sample, par);
  };
  auto lamdalikehoodAlfabdoneClinical_0 = [&](const double* par) {
      return functorLikehoodAlfabdoneClinical_0(sample, par);
  };

  ROOT::Math::Functor fpFunctor;
  if(glbstuff.clinicalfactors.size()==0){
    if (glbstuff.alfabdone < 0) 
    fpFunctor = ROOT::Math::Functor(lamdalikehoodFull, 4);
    else
    fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdone, 3);
  }else if(glbstuff.clinicalfactors.size()==1){
    if (glbstuff.alfabdone < 0) 
        fpFunctor = ROOT::Math::Functor(lamdalikehoodFullClinical_0, 5);
    else
        fpFunctor = ROOT::Math::Functor(lamdalikehoodAlfabdoneClinical_0, 4 );
  }

  std::unique_ptr<ROOT::Math::Minimizer> fpMinimizer(ROOT::Math::Factory::CreateMinimizer(glbstuff.fitalgo.at(fitalgindex).first, glbstuff.fitalgo.at(fitalgindex).second));

  fpMinimizer->SetFunction(fpFunctor);
  fpMinimizer->SetMaxFunctionCalls(1000000);
  fpMinimizer->SetMaxIterations(100000);
  fpMinimizer->SetTolerance(1e-2);
  fpMinimizer->SetStrategy(2);
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
      if(fixedpar.first!=p.first)
        fpMinimizer->SetLimitedVariable(p.first, name.c_str(), lim[0], lim[1], lim[2],lim[3]);
      else
        fpMinimizer->SetFixedVariable(p.first, name.c_str(), fixedpar.second);

  }  

  if(debug)
    cout<<"optimizeLikehood: before minimization"<<endl;  

  fpMinimizer->Minimize();
  if(debug)
    cout<<"optimizeLikehood: minimization odne"<<endl;    
  Int_t status=fpMinimizer->Status();
  if(fixedpar.first<0){ //no paramter was fixed, this is a minimization of the full model
    cout<<endl<<endl<<"optimizeLikehood: minimization done with "<<glbstuff.fitalgo.at(fitalgindex).first<<"/"<<glbstuff.fitalgo.at(fitalgindex).second<<endl;
    if(status==0){
      for(auto &paziente : sample){//fill scores
        double par[10]; //10 per essere sicuri di prendere tutti i parametri
        for(int i=0;i<fpMinimizer->NDim();i++)
          par[i]=fpMinimizer->X()[i];
        paziente.second.optlike_ntcpscore[fitalgindex]=EvalScoreSelector(glbstuff, paziente.second, par);
      }
    }else{
      cout<<"minimization failed"<<endl;  
      cout<<"print correlations:"<<endl;
      for(int i=0;i<fpMinimizer->NFree();++i)
        for(int j=0;j<fpMinimizer->NFree();++j)
          std::cout << i << " " << j << " " << fpMinimizer->CovMatrix(i,j)/sqrt(fpMinimizer->CovMatrix(i,i)*fpMinimizer->CovMatrix(j,j)) << std::endl;  
    }
    cout<<"status="<<status<<"CovMatrixStatus="<<fpMinimizer->CovMatrixStatus()<<"  Edm="<<fpMinimizer->Edm()<<"  degree of freedom:"<<fpMinimizer->NFree()<<"   -loglikehood minimum value:"<<fpMinimizer->MinValue()<<"  AIC="<<2*fpMinimizer->NFree()+2*fpMinimizer->MinValue()<<"  deviance/dof="<<2*fpMinimizer->MinValue()/(sample.size()-fpMinimizer->NFree())<<endl;
    for(int i=0;i<fpMinimizer->NDim();i++){
      if(glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size()==0){
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].push_back(fpMinimizer->X()[i]);
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].push_back(fpMinimizer->Errors()[i]);
      }else{
        cout<<"WARNING:optimizeLikehood: glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size() not =0, the fitted values has already filled... this is strange, please take a look; fitalgindex="<<fitalgindex<<"  size="<<glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].size()<<endl;
        cout<<"I'll overwrite the new values:"<<endl;
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].at(0)=fpMinimizer->X()[i];
        glbstuff.fittedpar[fitalgindex][fpMinimizer->VariableName(i)].at(1)=fpMinimizer->Errors()[i];
      }
      cout<<fpMinimizer->VariableName(i)<<":  "<<fpMinimizer->X()[i]<<" +- "<<fpMinimizer->Errors()[i]<<endl;
    }
    glbstuff.fitresults[fitalgindex].insert(glbstuff.fitresults[fitalgindex].end(), {(double)status, (double)fpMinimizer->CovMatrixStatus(), fpMinimizer->Edm(), (double)fpMinimizer->NFree(), fpMinimizer->MinValue(), 2*fpMinimizer->NFree()+2*fpMinimizer->MinValue(), 2*fpMinimizer->MinValue()/(sample.size()-fpMinimizer->NFree())});    
  }else{ //fixed par just for the calculation of Likehood Ratio Test (LRT)
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


//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=alfabeta, 4=clinical factor at 0
double functorLikehoodFullClinical_0(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-=std::log(EvalScoreLikehoodFullClinical_0(paziente.second, par));
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

//likehood: 
//par: 0=[0] firstfitted value in tf1, 1=secondfitted value in tf1, 2=nvalue, 3=clinical factor at 0
double functorLikehoodAlfabdoneClinical_0(const map<int, PatientData> &sample, const double* par){
  double eval=0.;
  for(const auto &paziente : sample)
    eval-= std::log(EvalScoreLikehoodAlfabdoneClinical_0(paziente.second, par));
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
    cout<<"start CalculateEudEqdAlreadyDone with paziente.id="<<paziente.id<<endl;

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


void optlike_fill(map<int, PatientData> &sample, const globalstuff &glbstuff, int fitalgindex){

  if(debug)
    cout<<"start optlike_fill"<<endl;
  TString name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_eud_optlike_All";
  TH1D *hall=new TH1D(name, "eud value;eud;Number of patients",100, 0., 100);        
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"eud_optlike_No";
  TH1D *hno=new TH1D(name, "eud value;eud;Number of patients",100, 0., 100);        
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"eud_optlike_Yes";
  TH1D *hyes=new TH1D(name, "eud value;eud;Number of patients",100, 0., 100);        
  
  gDirectory->cd("ntpc_linear");
  TGraph* gr_eud_vs_tox=new TGraph(sample.size());
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_ntcp_linear_eud_vs_tox";
  gr_eud_vs_tox->SetName(name);
  gr_eud_vs_tox->SetTitle("eud vs toxicity;eud;toxicity");
  TGraph* gr_score_vs_tox=new TGraph(sample.size());
  name=glbstuff.fitalgo.at(fitalgindex).first+glbstuff.fitalgo.at(fitalgindex).second+"_ntcp_linear_score_vs_tox";
  gr_score_vs_tox->SetName(name);
  gr_eud_vs_tox->SetTitle("eud vs NTCP score ;eud;NTCP model score");
  int index=0;
  for(auto &paziente : sample){
    paziente.second.optlike_eud[fitalgindex]= (glbstuff.alfabdone<0) ? CalculateEudFromScratch(paziente.second,glbstuff.fittedpar.at(fitalgindex).at("alfabeta").at(0) , glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0)) : CalculateEudEqdAlreadyDone(paziente.second, glbstuff.fittedpar.at(fitalgindex).at("nvalue").at(0));
    
    hall->Fill(paziente.second.optlike_eud.at(fitalgindex));
    if(paziente.second.tgt_acutegitox>0.5)
      hyes->Fill(paziente.second.optlike_eud.at(fitalgindex));
    else
      hno->Fill(paziente.second.optlike_eud.at(fitalgindex));
    gr_eud_vs_tox->SetPoint(index, paziente.second.optlike_eud.at(fitalgindex), paziente.second.tgt_acutegitox);
    gr_score_vs_tox->SetPoint(index, paziente.second.optlike_ntcpscore.at(fitalgindex), paziente.second.tgt_acutegitox);
    index++;
  }
  gr_eud_vs_tox->Sort();

  TF1* sigmoidbest=new TF1("sigmoidbest", "1./(1.+exp(-[0]-[1]*x))", 0, 1000);
  sigmoidbest->FixParameter(0, glbstuff.fittedpar.at(fitalgindex).at("beta_zero").at(0));
  sigmoidbest->FixParameter(1, glbstuff.fittedpar.at(fitalgindex).at("beta_eud").at(0));
  gr_eud_vs_tox->Fit(sigmoidbest, "BS+","",0,100);
  gr_eud_vs_tox->SetMarkerStyle(20);
  gr_eud_vs_tox->SetMarkerColor(2);
  gr_eud_vs_tox->SetLineWidth(0);
  gr_eud_vs_tox->SetLineColor(0);
  gr_eud_vs_tox->SetDrawOption("AP*");
  gr_eud_vs_tox->Write();
  gr_score_vs_tox->SetMarkerStyle(20);
  gr_score_vs_tox->SetMarkerColor(2);
  gr_score_vs_tox->SetLineWidth(0);
  gr_score_vs_tox->SetLineColor(0);
  gr_score_vs_tox->SetDrawOption("AP*");
  gr_score_vs_tox->Write();
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
    auc += (FPR - prevFPR) * (TPR + prevTPR) / 2.0;
    avg_precision+=  tp/(tp+fp)*(TPR-prevTPR);
    grROC->SetPoint(ip++, FPR, TPR);
    prevTPR = TPR;
    prevFPR = FPR;
  }
  grROC->SetTitle(Form("ROC curve for likehood optmized simple NTCP model with AUC=%f;False Positive Rate;True Positive Rate",auc));
  grROC->SetLineWidth(2);
  grROC->SetMarkerStyle(20);
  grROC->Draw("ALP");
  grROC->Write();

  return make_pair(auc,avg_precision);
}



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
    for(auto it2=std::next(it1);it2!=orderedNames.end();++it2){

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


void PlotCalibrationCurveQuantilesAndHLtest(const std::map<int, PatientData>& sample, const globalstuff& glbstuff, int fitalgindex, int nbins = 10){

  if(debug)
    cout<<"start PlotCalibrationCurveQuantiles"<<endl;

  vector<pair<double,double>> data; // pred, obs
  for(const auto& paziente : sample){
      double pred = paziente.second.optlike_ntcpscore.at(fitalgindex);
      double obs  = paziente.second.tgt_acutegitox;

      if(pred < 0.0 || pred > 1.0){
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

  TString pltname="canvascalib_"+glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second;
  TCanvas* c = new TCanvas("c_calib_quant", "Calibration quantiles", 800, 800);
  TGraphErrors* gr = new TGraphErrors(x.size(), x.data(), y.data(),ex.data(),ey.data());
  pltname="calib_"+glbstuff.fitalgo.at(fitalgindex).first+"/"+glbstuff.fitalgo.at(fitalgindex).second;
  gr->SetName(pltname);
  gr->SetTitle(";Mean predicted NTCP;Observed toxicity rate");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetLimits(0.0, 1.0);
  gr->GetYaxis()->SetRangeUser(0.0, 1.0);
  gr->Write();
  gr->Draw("APLE");
  TLine* line = new TLine(0.0, 0.0, 1.0, 1.0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");
  c->Write();
  // c->SetGrid();
  // c->SaveAs(outname.c_str());
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


void fillGlobalStuff(globalstuff &glbstuff, double alfabdone, double eqd2binwidth, const vector<double> &nvalue4eud, const vector<double> &alfabeta, const map<string, pair<int,vector<double>>> &fitpars,   const vector<pair<string,string>> &fitalgo, int issynthetic, const vector<int> &clinicalfactors){
  glbstuff.alfabdone=alfabdone;
  glbstuff.eqd2binwidth=eqd2binwidth;
  glbstuff.nvalue4eud=nvalue4eud;
  glbstuff.alfabeta=alfabeta;
  glbstuff.fitpars=fitpars;
  glbstuff.fitalgo=fitalgo;
  glbstuff.issynthetic=issynthetic;
  glbstuff.clinicalfactors=clinicalfactors;
  glbstuff.maxbin= (issynthetic==true) ? 200:100;

  return;
}
