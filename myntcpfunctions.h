// std libraries
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cctype>
#include <regex>
#include <string>
#include <sstream>
#include <vector>
#include <map>


// ROOT
#include <TApplication.h>
#include <TString.h>
#include <TH2D.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include "TROOT.h"
#include <TCanvas.h>


#define maxdvhgy 77


using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;

struct PatientData {
int id;
double volume;
double max_dose_plan;
double max_dose_st;
double min_dose_st;
double mean_dose_st;
vector<double> dvhmap;
int tgt_acutegitox;
};



void bookHisto(TFile *outrootfile, const   map<int, PatientData> &sample);
void PrintSampleLine(int idx, const   map<int, PatientData> &sample);
int loadDvhFile(const string& filename,   map<int, PatientData> &sample);
int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname);
string trim(const string& s);
vector<string> splitCsvLine(const string& line, const TString delimiter);
