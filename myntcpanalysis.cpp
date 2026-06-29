#include "myntcpfunctions.h"

int debug = 0;

int main(int argc, char* argv[]) {
  
  TString dvhbfilename("");
  TString dvhafilename("");
  // dvhafilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv";
  // dvhafilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT_plus_mbprofile_pts.csv";
  // dvhafilename="DVH_PRIME_bowel_by_Boris_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv";
  dvhafilename="DVH_Bowel_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT_plus_mbprofile_pts.csv";
  TString metafilename("metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24_MODIFIEDYUNWITHCLUSTER.csv");
  TString outrootname("ntcp_outputs.root");
  TString txtappended("");
  TString tgtname("acute GI toxicity");
  int datatype=3; //0=not specified, 1=hiroc synthetic, 2=nanoport, 3=old article clustering with MB class risk
  int clusternum=3; //number of cluster considered, it is related to clinicalfactors, 
  int clinicalfactors=1; //0=no clinical factors, 1=only one value as clinical factor, 2= three values of cluster that actually are normalized
  // vector<int> clinicalfactors;
  double alfabdone=10; //if the dose are already normalized for fractions and alfa/beta, otherwise set to -1
  double eqd2binwidth=1.; //binwidth in gy di eqd2 normalizzato
  vector<double> alfabeta={0.1,0.5,1,2,4,10};
  if(alfabdone>0){
    alfabeta.clear();
    alfabeta.push_back(alfabdone);
  }
  vector<double> nvalue4eud={ 0.01,0.1,1,5,10,100}; //1->eud=mean dose, +inf->eud=max dose
  map<string, pair<int,vector<double>>> fitpars;
  vector<pair<string,string>> fitalgo;
  // fitalgo.push_back(make_pair("Minuit2","Simplex"));
  fitalgo.push_back(make_pair("Minuit2","Migrad"));
  fitalgo.push_back(make_pair("Minuit2","Hesse"));
  fitalgo.push_back(make_pair("Minuit2","Combined"));
  // fitalgo.push_back(make_pair("Minuit2","Scan"));
  // fitalgo.push_back(make_pair("GSLMultiMin","BFGS2"));
  // fitalgo.push_back(make_pair("Generic",""));
  // fitalgo.push_back(make_pair("GSLSimAn",""));
  
  //value: 0:number of parameter index, 1=initial value, 2=step, 3=lower, 4=upper (for setlimitedvariable)
  vector parlimits = {-2.0, 0.01, -50.0, 50.0};
  int parnum=0;
  fitpars["beta_zero"]=make_pair(parnum++, parlimits);
  parlimits={0.05, 0.01, -5.0, 5.0};
  fitpars["beta_eud"]=make_pair(parnum++, parlimits);
  parlimits={1.0, 0.01, 0.01, 100.0};
  fitpars["nvalue"]=make_pair(parnum++, parlimits);
  if(alfabdone<0){
    parlimits={0.5, 0.01, 0.01, 5.0};
    fitpars["alfabeta"]=make_pair(parnum++, parlimits);
  }
  for(int i=0;i<clinicalfactors;i++){
    parlimits={1.0, 0.01, -5., 5.0};
    fitpars[Form("clinical_factor_%i",i)]=make_pair(parnum++, parlimits);
  }

  if(datatype==1){
    dvhafilename="Synthetic_Data/clinical_nfrac_out_risk.csv";
    outrootname="ntcp_synthetic_out.root";
  }

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)         debug = atoi(argv[++i]);
    if (strcmp(argv[i], "-dvha") == 0)      dvhafilename = TString(argv[++i]);
    if (strcmp(argv[i], "-dvhb") == 0)      dvhbfilename = TString(argv[++i]);
    if (strcmp(argv[i], "-meta") == 0)      metafilename = TString(argv[++i]);
    if (strcmp(argv[i], "-out") == 0)       outrootname = TString(argv[++i]);
    if (strcmp(argv[i], "-txt") == 0)       txtappended = TString(argv[++i]);
    if (strcmp(argv[i], "-tgt") == 0)       tgtname = TString(argv[++i]);
  }

  if(txtappended.Length()>0){
    freopen(txtappended.Data(),"w",stdout);
    cout<<"dvha="<<dvhafilename.Data()<<"  dvhbfilename="<<dvhbfilename.Data()<<"  metafilename="<<metafilename.Data()<<"  outrootname="<<outrootname.Data()<<"  tgtname="<<tgtname.Data()<<endl;
    cout<<"datatype="<<datatype<<"  clinicalfactors="<<clinicalfactors<<"  alfabdone="<<alfabdone<<"  eqd2binwidth="<<eqd2binwidth<<endl;
    cout<<"alfabeta=(";
    for(int i=0;i<alfabeta.size();i++)
      cout<<alfabeta.at(i)<<", ";
    cout<<")"<<endl;
    cout<<"nvalue4eud=(";
    for(int i=0;i<nvalue4eud.size();i++)
      cout<<nvalue4eud.at(i)<<", ";
    cout<<")"<<endl;
    for(const auto &p : fitpars){
        cout<<p.first<<" "<< p.second.first;
        for (double x : p.second.second)
            cout<<", "<<x;
        cout << endl;
    }    
  }

  map<int, PatientData> sample, samrect;
  
  if(datatype==1){
    if(loadSyntheticFile(dvhafilename.Data(), sample))
      return 1;
  }else{ 
    if(loadDvhFile(dvhafilename.Data(), sample))
      return 1;
    if(loadMetaFile(metafilename.Data(), sample, tgtname, datatype))
      return 1;
    if(dvhbfilename.Length()>0){
      if(loadDvhFile(dvhbfilename.Data(), samrect))
        return 1;
    }
  }

  globalstuff glbstuff; 
  fillGlobalStuff(glbstuff, alfabdone, eqd2binwidth, nvalue4eud, alfabeta, fitpars, fitalgo, datatype, clinicalfactors, clusternum);

  if(clinicalfactors>0)
    SetClusterAsClinicalFactor(sample, glbstuff);

  if(alfabdone<0){
    evaluateEqdEud(sample, glbstuff);
    if(samrect.size()>0)
    evaluateEqdEud(samrect, glbstuff);
  }
  else{
    FillEqdEud(sample, glbstuff);
    if(samrect.size()>0)
      FillEqdEud(samrect,glbstuff);
  }

  TFile* outrootfile = new TFile(outrootname, "RECREATE");
  if (!outrootfile || outrootfile->IsZombie()) {
    cerr << "ERROR: cannot create output ROOT file " << outrootname.Data() << endl;
    return 1;
  }
 
  bookHisto(outrootfile, sample, glbstuff, samrect.size());
  fillHisto(sample,glbstuff);
  CreateNtcpSigmoidalPredefined(sample, tgtname);
  PostLoopAnalysis(sample, glbstuff);
  
  //le'ts find the best eud, n, alfabeta, n values minimizing a likehood function
  for(int i=0;i<fitalgo.size();i++){
    int optstatus=optimizeLikehood(sample, glbstuff, i);//0 because it's the first index for fitalgo which has only one enters, for the moment
    if(optstatus==0){
      optlike_fill(sample, glbstuff,i);
      pair<double,double>aucprecres=optlike_aucROC(sample, glbstuff, i);
      cout<<"fitalgo with "<<fitalgo.at(i).first<<" and "<<fitalgo.at(i).second<<" done. AUC="<<aucprecres.first<<"  average_precision="<<aucprecres.second<<endl;
      SetAucAvgPrec(i, aucprecres, glbstuff); 
      PlotCalibrationCurveQuantilesAndHLtest(sample, glbstuff, i, 5);    
      for(auto &par:fitpars)
        optimizeLikehood(sample, glbstuff, i, make_pair(par.second.first,0.)); //for Likelihood ratio test, dove parametri in più sono messi a 0 pef fare contfronto.
      // for(auto &par:fitpars) //TODO: to be developed
      //   optimizeLikehood(sample, glbstuff, i, make_pair(par.second.first,par.second)); //for Likelihood profile, dove parametri in più sono messi al minimo per fare likelihood profile.
    }
  }
  
  ChooseBestFit(glbstuff);
  // DrawLikeHood(sample, glbstuff);  //it took time to draw likehood for each point, disable it for code development purpouse

  outrootfile->Write();
  outrootfile->Close();
  cout<<"myntcpanalysis done: outputfile: "<<outrootfile->GetName()<<endl;

return 0;
}