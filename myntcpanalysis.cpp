#include "myntcpfunctions.h"


int main(int argc, char* argv[]) {
  
  int debug = 0;
  TString dvhfilename("");
  TString metafilename("");
  TString outrootname("ntcp_outputs.root");
  TString txtappended(".txt");

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0)         debug = atoi(argv[++i]);
    if (strcmp(argv[i], "-dvh") == 0)       dvhfilename = TString(argv[++i]);
    if (strcmp(argv[i], "-meta") == 0)      metafilename = TString(argv[++i]);
    if (strcmp(argv[i], "-out") == 0)       outrootname = TString(argv[++i]);
    if (strcmp(argv[i], "-txt") == 0)       txtappended = TString(argv[++i]);
  }

  vector<PatientData> sample;
  loadDvhCsv(dvhfilename.Data(), sample);
  PrintSampleLine(0,sample);
  PrintSampleLine(3,sample);

return 0;
}