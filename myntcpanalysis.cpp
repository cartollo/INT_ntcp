#include "myntcpfunctions.h"

int debug = 0;

int main(int argc, char* argv[]) {
  
  TString dvhbfilename("");
  dvhbfilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT_plus_mbprofile_pts.csv";
  TString dvhafilename("");
  // dvhafilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv";
  // dvhafilename="DVH_Rectum_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT_plus_mbprofile_pts.csv";
  // dvhafilename="DVH_PRIME_bowel_by_Boris_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv";
  dvhafilename="DVH_Bowel_MIM_1Gy_10ab_single_dose_converted_ML_prostate_WPRT_plus_mbprofile_pts.csv";
  Int_t powptype=1; //-1=no selection, 0= only prostate, 1= whole pelvis 
  TString metafilename("metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24_MODIFIEDYUNWITHCLUSTER.csv");
  TString outrootname("ntcp_outputs.root");
  TString txtappended("");
  TString tgtname("acute GI toxicity");
  int twodvh=0; //1= use both dvhb and dvha, otherwise only 0 WARNING: if twodvh==1, only clusterfactor 2 and clinicalfactors=2 and alfabdoneshould be set
  int datatype=3; //0=not specified, 1=hiroc synthetic, 2=nanopore by michele, 3=old article clustering with MB class risk
  int clusternum=0; //number of cluster considered, it is related to clinicalfactors, (per ora è 0 o 3) 
  int clinicalfactors=0; //0=no clinical factors, 1=only one value as clinical factor, 2= three values of cluster that actually are normalized
  // vector<int> clinicalfactors;
  int drawlikehood=1; //draw or not the likehood function (it took time)
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
  // fitalgo.push_back(make_pair("Minuit2","Hesse")); // tanto vengono uguali...
  // fitalgo.push_back(make_pair("Minuit2","Combined"));
  // fitalgo.push_back(make_pair("Minuit2","Scan"));
  // fitalgo.push_back(make_pair("GSLMultiMin","BFGS2"));
  // fitalgo.push_back(make_pair("Generic",""));
  // fitalgo.push_back(make_pair("GSLSimAn",""));
  
  //value: 0:number of parameter index, 1=initial value, 2=step, 3=lower, 4=upper (for setlimitedvariable)
  int parnum=0;
  vector parlimits = {-2.0, 0.01, -5.0, 5.0};
  fitpars["beta_zero"]=make_pair(parnum++, parlimits);
  parlimits={0.0, 0.01, -10., 10.0};
  fitpars["beta_eud_a"]=make_pair(parnum++, parlimits);
  parlimits={90, 0.01, 0, 150.0};
  fitpars["nvalue"]=make_pair(parnum++, parlimits);
  if(twodvh){
    parlimits={0.0, 0.01, -10., 10.0};
    fitpars["beta_eud_b"]=make_pair(parnum++, parlimits);
  }
  if(alfabdone<0){
    parlimits={0.5, 0.01, 0.01, 5.0};
    fitpars["alfabeta"]=make_pair(parnum++, parlimits);
  }
  for(int i=0;i<clinicalfactors;i++){
    if(i==0)
      parlimits={1.0, 0.01, -10, 10};
    else
      parlimits={5.0, 0.01, 0., 10.};
    fitpars[Form("clinical_factor_%i",i)]=make_pair(parnum++, parlimits);
  }

  if(datatype==1){
    dvhafilename="Synthetic_Data/clinical_nfrac_out_risk.csv";
    outrootname="ntcp_synthetic_out.root";
  }

//do some checks on dvhb:
if(twodvh>0 && (alfabdone<0 || dvhbfilename.Length()<4 || (datatype!=2 && datatype!=3))){
  cout<<"ERROR twodvh="<<twodvh<<"  alfabdone="<<alfabdone<<"  dvhbfilename="<<dvhbfilename.Data()<<" datatype="<<datatype<<endl<<"Check the parameters, recompile and relauch, asino"<<endl;
  return 1;  
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
    cout<<"datatype="<<datatype<<"  powptype="<<powptype<<"  clusternum="<<clusternum<<"  clinicalfactors="<<clinicalfactors<<"  alfabdone="<<alfabdone<<"  eqd2binwidth="<<eqd2binwidth<<endl;
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
    if(loadMetaFile(metafilename.Data(), sample, tgtname, datatype, powptype))
      return 1;
    if(twodvh){//double file with double sample!
      if(loadDvhFile(dvhbfilename.Data(), samrect))
        return 1;
      if(loadMetaFile(metafilename.Data(), samrect, tgtname, datatype, powptype))
        return 1;
      if(CheckSampleSamrectConsistency(sample, samrect))
        return 1;
    }
  }

  globalstuff glbstuff; 
  fillGlobalStuff(glbstuff, alfabdone, eqd2binwidth, nvalue4eud, alfabeta, fitpars, fitalgo, datatype, clinicalfactors, clusternum, powptype, twodvh);

  if(clinicalfactors>0)
    SetClusterAsClinicalFactor(sample, glbstuff);
  if(clinicalfactors>0 && twodvh)
    SetClusterAsClinicalFactor(samrect, glbstuff);

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
 
  bookHisto(outrootfile, sample, glbstuff, samrect.size());  //samrect non usato in realtà
  fillHisto(sample,glbstuff);
  CreateNtcpSigmoidalPredefined(sample, tgtname);
  PostLoopAnalysis(sample, glbstuff);
  
  //le'ts find the best eud, n, alfabeta, n values minimizing a likehood function
  for(int i=0;i<fitalgo.size();i++){
    int optstatus=optimizeLikehood(sample, glbstuff, i, samrect);//0 because it's the first index for fitalgo which has only one enters, for the moment
  }
  
  ChooseBestFit(glbstuff);
  if(drawlikehood)
    DrawLikeHood(sample, glbstuff);  //it took time to draw likehood for each point, disable it for code development purpouse

  outrootfile->Write();
  outrootfile->Close();
  cout<<"myntcpanalysis done: outputfile: "<<outrootfile->GetName()<<endl;

return 0;
}