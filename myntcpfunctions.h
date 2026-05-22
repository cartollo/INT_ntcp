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
#include <algorithm>
#include <functional>
#include <cmath>


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
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TVectorD.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::cerr;
using std::make_pair;
using std::pair;
using std::tuple;

struct PatientData {
int id;

//from dsh file
double volume;
double max_dose_plan;
double max_dose_st;
double min_dose_st;
double mean_dose_st;
vector<double> dvhmapcum; //binwidth =1 gy, each value is the volume for a given dose (dose=index), cumulative
vector<double> dvhmapdiff; //binwidth =1 gy, each value is the volume for a given dose (dose=index), differential
vector<double> dvhdctdiff; //differential dct transform

//from meta file
double dose_ptv;
int nfraction;
double dose_per_fraction;
double mean_dose_rectum; 
int tgt_acutegitox;

//calculated stuff
// map<double, vector<double>> eqd2map; //equivalent dose in 2 gy fraction, each value in this vector is a dose value, key=alfabeta
map<double,vector<double>> dvhcumnormmap; //equivalent cumulative dvh normalized for alfabeta and fractions for all the patients with eqd2map, key=alfabeta
map<double,vector<double>> dvhdiffnormmap; //equivalent differential dvh normalized for alfabeta and fractions for all the patients with eqd2map, key=alfabeta
map<pair<double,double>, double> eudmap; //eud value (value) for given n and alfa/beta (key=make_pair(n,alfabeta))
double optlike_eud; //from optimizeLikehood
map<double,double> optlike_ntcpscore; //key: fitalgo index, value: ntcp score
int status; //-1=not set, 0=all ok, 1=not monotonous
};

struct globalstuff{
  //parameters:
  double alfabdone;
  double eqd2binwidth;
  vector<double> alfabeta;
  vector<double> nvalue4eud;

  //fitted stuff
  vector<pair<string,string>> fitalgo; //algorithm used by tminimizer for fitting
  map<int, pair<string,vector<double>>> fitpars; //ntcp model parameters key: number of parameters, value: name, value, step, lower, upper (for setlimitedvariable)
  map<int, vector<double>> fitresults; //key: fitalgo index, value: 0=status, 1=CovMatrixStatus, 2=edm, 3=dof, 4=likehood->minvalue, 5=aic, 6=dev/dof
  map<int, map<string,pair<double,double>>> fittedpar; //key: fitalgo index, value: fitted parameter name, value and error
  int bestvalue; //index of the best fit 
};


void fillHisto(map<int, PatientData> &sample,const globalstuff &glbstuff);
void PostLoopAnalysis(map<int, PatientData> &sample, const globalstuff &glbstuff);
void bookHisto(TFile *outrootfile, const globalstuff &glbstuff, const bool multivariate);
int CreateNtcpSigmoidalSingle(TString pltname, const TVectorD &vx, const TVectorD &vy);
int CreateNtcpSigmoidalPredefined(map<int, PatientData> &sample, TString tgtname);
void PrintSampleLine(int idx, const   map<int, PatientData> &sample);
int loadDvhFile(const string& filename,   map<int, PatientData> &sample);
int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname);
void evaluateEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff);
void FillEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff);
double CalculateEudFromScratch(const PatientData &paziente, double alfabeta, double nvalue);
double CalculateEudEqdAlreadyDone(const PatientData &paziente, double nvalue);
double functorLikehoodFull(const map<int, PatientData> &sample, const double* par);
double functorLikehoodAlfabdone(const map<int, PatientData> &sample, const double* par);
double fitSigmoidal(TGraph* graph, int parnum, int functype);
string trim(const string& s);
vector<string> splitCsvLine(const string& line, const TString delimiter);
void CreateHistoFromTgraph(TGraph *gr, TH1D *h);
int optimizeLikehood(const map<int, PatientData> &sample, globalstuff &glbstuff, const int fitalgindex);
void optlike_fill(map<int, PatientData> &sample, const globalstuff &glbstuff, int fitalgindex);
double optlike_aucROC(const map<int, PatientData> &sample, const int fitalgindex);
void computeDCT(const vector<double>& x, vector<double>& c);
void fillGlobalStuff(globalstuff &glbstuff, double alfabdone, double eqd2binwidth, const vector<double> &nvalue4eud, const vector<double> &alfabeta, const map<int, pair<string,vector<double>>> &fitpars,   const vector<pair<string,string>> &fitalgo);

inline double EqdDose(PatientData &paziente, double alfabeta, double dose){return dose*(alfabeta+dose/paziente.nfraction)/(alfabeta+2.);};