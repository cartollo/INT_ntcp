#include "myntcpfunctions.h"

int debug = 0;

int main(int argc, char* argv[]) {
  
  TString dvhfilename("DVH_PRIME_bowel_by_Boris_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv");
  TString metafilename("metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24.csv");
  TString outrootname("ntcp_outputs.root");
  TString txtappended(".txt");
  TString tgtname("acute GI toxicity");
  double eqd2binwidth=0.5; //binwidth in gy di eqd2 normalizzato
  vector<double> alfabeta={1.,1.25,1.5,1.75,2.,2.5,3};
  vector<double> nvalue4eud={1,5,10,100}; //1->eud=mean dose, +inf->eud=max dose

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)         debug = atoi(argv[++i]);
    if (strcmp(argv[i], "-dvh") == 0)       dvhfilename = TString(argv[++i]);
    if (strcmp(argv[i], "-meta") == 0)      metafilename = TString(argv[++i]);
    if (strcmp(argv[i], "-out") == 0)       outrootname = TString(argv[++i]);
    if (strcmp(argv[i], "-txt") == 0)       txtappended = TString(argv[++i]);
    if (strcmp(argv[i], "-tgt") == 0)       tgtname = TString(argv[++i]);
  }

  map<int, PatientData> sample;
  if(loadDvhFile(dvhfilename.Data(), sample, alfabeta))
    return 1;
  if(loadMetaFile(metafilename.Data(), sample, tgtname))
    return 1;

  if(evaluateEqdEud(sample, alfabeta, nvalue4eud,eqd2binwidth))
    return 1;

  TFile* outrootfile = new TFile(outrootname, "RECREATE");
  if (!outrootfile || outrootfile->IsZombie()) {
    cerr << "ERROR: cannot create output ROOT file " << outrootname.Data() << endl;
    return 1;
  }
 
  bookHisto(outrootfile, alfabeta,nvalue4eud, eqd2binwidth);
  fillHisto(outrootfile, sample, alfabeta, nvalue4eud, eqd2binwidth);
  CreateNtcpSigmoidal(sample, tgtname, alfabeta, nvalue4eud);

  outrootfile->Write();
  outrootfile->Close();

return 0;
}