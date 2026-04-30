// std libraries
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cctype>
#include <regex>
#include <string>
#include <sstream>
#include <vector>


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

struct PatientData {
int id;
double volume;
double max_dose_plan;
double max_dose_st;
double min_dose_st;
double mean_dose_st;
vector<double> dvhmap;
};



void bookHisto(TFile *outrootfile, const vector<PatientData> &sample);
void PrintSampleLine(int idx, const vector<PatientData> &sample);
void loadDvhCsv(const string& filename, vector<PatientData> &sample);
string trim(const string& s);
vector<string> splitCsvLine(const string& line);
