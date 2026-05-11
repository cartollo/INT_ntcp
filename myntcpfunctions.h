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
#include <TF1.h>
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
using std::make_pair;
using std::pair;

struct PatientData {
int id;

//from dsh file
double volume;
double max_dose_plan;
double max_dose_st;
double min_dose_st;
double mean_dose_st;
vector<double> dvhmap; //binwidth =1 gy, each value is the volume for a given dose (dose=index)

//from meta file
double dose_ptv;
int nfraction;
double dose_per_fraction;
double mean_dose_rectum; 
int tgt_acutegitox;

//calculated stuff
map<double, vector<double>> eqd2map; //equivalent dose in 2 gy fraction, each value in this vector is a dose value, key=alfabeta
map<double,vector<double>> dvhnormmap; //equivalent dvh normalized to the same scale for all the patients with eqd2map, key=alfabeta
map<pair<double,double>, double> eudmap; //eud value (value) for given n and alfa/beta (key=make_pair(n,alfabeta))
int status; //-1=not set, 0=all ok, 1=not monotonous
};



void fillHisto(TFile *outrootfile,  map<int, PatientData> &sample, const vector<double> &alfabeta, const vector<double> &nvalue4eud, double eqd2binwidth);
void bookHisto(TFile *outrootfile, const vector<double> &alfabeta, const vector<double> &nvalue4eud, double eqd2binwidth);
int CreateNtcpSigmoidal(map<int, PatientData> &sample, TString tgtname, const vector<double> &alfabeta, const vector<double> &nvalue4eud);
void PrintSampleLine(int idx, const   map<int, PatientData> &sample);
int loadDvhFile(const string& filename,   map<int, PatientData> &sample, const vector<double> &alfabeta);
int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname);
int evaluateEqdEud(map<int, PatientData> &sample, const vector<double> &alfabeta,const vector<double> &nvalue4eud, double eqd2binwidth);
double CalculateEudFromScratch(PatientData &paziente, double eqd2binwidth, double alfabeta, double nvalue);
double fitSigmoidal(TGraph* graph, int parnum, TFitResultPtr &fitresults);
string trim(const string& s);
vector<string> splitCsvLine(const string& line, const TString delimiter);
