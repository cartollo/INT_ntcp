#include "myntcpfunctions.h"

int debug = 0;

int main(int argc, char* argv[]) {
  
  TString dvhbfilename("");
  // dvhbfilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv";
  TString dvhafilename("DVH_PRIME_bowel_by_Boris_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv");
  TString metafilename("metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24.csv");
  TString outrootname("ntcp_outputs.root");
  TString txtappended(".txt");
  TString tgtname("acute GI toxicity");
  bool issynthetic=true;
  vector<int> clinicalfactors{0};
  double alfabdone=-1; //if the dose are already normalized for fractions and alfa/beta, otherwise set to -1
  double eqd2binwidth=1.; //binwidth in gy di eqd2 normalizzato
  vector<double> alfabeta={0.1,0.5,1,2,4,10};
  vector<double> nvalue4eud={ 0.01,0.1,1,5,10,100}; //1->eud=mean dose, +inf->eud=max dose
  map<string, pair<int,vector<double>>> fitpars;
  vector<pair<string,string>> fitalgo;
  // fitalgo.push_back(make_pair("Minuit2","Simplex"));
  fitalgo.push_back(make_pair("Minuit2","Migrad"));
  // fitalgo.push_back(make_pair("Minuit2","Hesse"));
  fitalgo.push_back(make_pair("Minuit2","Combined"));
  // fitalgo.push_back(make_pair("Minuit2","Scan"));
  // fitalgo.push_back(make_pair("GSLMultiMin","BFGS2"));
  // fitalgo.push_back(make_pair("Generic",""));
  // fitalgo.push_back(make_pair("GSLSimAn",""));
  
  //value: 0:number of parameter index, 1=initial value, 2=step, 3=lower, 4=upper (for setlimitedvariable)
  vector parlimits = {-2.0, 0.01, -50.0, 50.0};
  fitpars["beta_zero"]=make_pair(0, parlimits);
  parlimits={0.05, 0.01, -5.0, 5.0};
  fitpars["beta_eud"]=make_pair(1, parlimits);
  parlimits={1.0, 0.01, 0.01, 100.0};
  fitpars["nvalue"]=make_pair(2, parlimits);
  if(alfabdone<0){
    parlimits={0.5, 0.01, 0.01, 5.0};
    fitpars["alfabeta"]=make_pair(3, parlimits);
  }

  parlimits={1.0, 0.01, -5., 5.0};
  fitpars["clinical_factor_0"]=make_pair(4, parlimits);

  if(issynthetic){
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

  map<int, PatientData> sample, samrect;
  if(issynthetic){
    if(loadSyntheticFile(dvhafilename.Data(), sample))
      return 1;
  }else{
    if(loadDvhFile(dvhafilename.Data(), sample))
      return 1;
    if(loadMetaFile(metafilename.Data(), sample, tgtname))
      return 1;
  }
  if(dvhbfilename.Length()>0){
    if(loadDvhFile(dvhbfilename.Data(), samrect))
      return 1;
  }

  globalstuff glbstuff; 
  fillGlobalStuff(glbstuff, alfabdone, eqd2binwidth, nvalue4eud, alfabeta, fitpars, fitalgo, issynthetic, clinicalfactors);

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
 
  bookHisto(outrootfile, glbstuff, samrect.size());
  fillHisto(sample,glbstuff);
  CreateNtcpSigmoidalPredefined(sample, tgtname);
  PostLoopAnalysis(sample, glbstuff);
  
  //le'ts find the best eud, n, alfabeta, n values minimizing a likehood function
  for(int i=0;i<fitalgo.size();i++){
    int optstatus=optimizeLikehood(sample, glbstuff, i);//0 because it's the first index for fitalgo which has only one enters, for the moment
    if(optstatus==0){
      optlike_fill(sample, glbstuff,i);
      double optauc=optlike_aucROC(sample, i);
      cout<<"fitalgo with "<<fitalgo.at(i).first<<" and "<<fitalgo.at(i).second<<" done. AUC="<<optauc<<endl;
      SetAuc(i, optauc, glbstuff); 
    }
  }
  
  ChooseBestFit(glbstuff);
  DrawLikeHood(sample, glbstuff);

  outrootfile->Write();
  outrootfile->Close();
  cout<<"myntcpanalysis done: outputfile: "<<outrootfile->GetName()<<endl;

return 0;
}