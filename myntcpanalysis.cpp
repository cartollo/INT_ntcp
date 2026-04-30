#include "myntcpfunctions.h"


int main(int argc, char* argv[]) {
  
  int debug = 0;
  TString dvhfilename("DVH_PRIME_bowel_by_Boris_1Gy_10ab_single_dose_converted_ML_prostate_WPRT.csv");
  TString metafilename("metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24.csv");
  TString outrootname("ntcp_outputs.root");
  TString txtappended(".txt");
  TString tgtname("acute GI toxicity");

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)         debug = atoi(argv[++i]);
    if (strcmp(argv[i], "-dvh") == 0)       dvhfilename = TString(argv[++i]);
    if (strcmp(argv[i], "-meta") == 0)      metafilename = TString(argv[++i]);
    if (strcmp(argv[i], "-out") == 0)       outrootname = TString(argv[++i]);
    if (strcmp(argv[i], "-txt") == 0)       txtappended = TString(argv[++i]);
    if (strcmp(argv[i], "-tgt") == 0)       tgtname = TString(argv[++i]);
  }

  map<int, PatientData> sample;
  if(loadDvhFile(dvhfilename.Data(), sample))
    return 1;
  if(loadMetaFile(metafilename.Data(), sample, tgtname))
    return 1;
return 0;
}