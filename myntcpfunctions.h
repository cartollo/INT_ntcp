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
#include <TGraph.h>


#define maxdvhgy 77


using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::cerr;

struct PatientData {
int id;

//from dsh file
double volume;
double max_dose_plan;
double max_dose_st;
double min_dose_st;
double mean_dose_st;
//from meta file
double dose_ptv;
int nfraction;
double dose_per_fraction;
double mean_dose_rectum; 
int tgt_acutegitox;
//other stuff
vector<double> dvhmap;
vector<double> eqdmap;
vector<double> eudmap;
int status; //-1=not set, 0=all ok, 1=not monotonous
};



void fillHisto(TFile *outrootfile, map<int, PatientData> &sample);
void bookHisto(TFile *outrootfile);
int fitNtcpLinearRegression(map<int, PatientData> &sample, TString tgtname);
void PrintSampleLine(int idx, const   map<int, PatientData> &sample);
int loadDvhFile(const string& filename,   map<int, PatientData> &sample);
int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname);
int evaluateEqdEud(PatientData &patient);
string trim(const string& s);
vector<string> splitCsvLine(const string& line, const TString delimiter);
