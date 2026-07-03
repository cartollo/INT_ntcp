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
#include <TF2.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include "TROOT.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <Math/DistFunc.h>
#include <Math/ProbFuncMathCore.h>


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
vector<double> dvhdctdiff; //differential dct transform (a tentative to use dct transform to ealuate dvh differences)

//from meta file
double dose_ptv;
int powp; //prostate only=0 or whole pelvis=1
int nfraction;
double dose_per_fraction;
double mean_dose_rectum; 
int tgt_acutegitox;
int operation; //at present, =TURP for synthetic data
vector<int> clinical_factor; 
int prostatectomy;
int appendectomy;

int age;
double bmi;
int smoke;
int alcohol;
int diabetes;
int diverticulitis;
double psa;
int tumour_risk;
int lymph_irr;
int semivesicle_irr;
int mb_risk;
int cluster; //microbiota cluster classification



//calculated stuff
// map<double, vector<double>> eqd2map; //equivalent dose in 2 gy fraction, each value in this vector is a dose value, key=alfabeta
map<double,vector<double>> dvhcumnormmap; //equivalent cumulative dvh normalized for alfabeta and fractions for all the patients with eqd2map, key=alfabeta
map<double,vector<double>> dvhdiffnormmap; //equivalent differential dvh normalized for alfabeta and fractions for all the patients with eqd2map, key=alfabeta
map<pair<double,double>, double> eudmap; //eud value (value) for given n and alfa/beta (key=make_pair(n,alfabeta))
map<int,double> optlike_eud; //from optimizeLikehood
map<int,double> optlike_ntcpscore; //key: fitalgo index, value: ntcp score
int status; //-1=not set, 0=all ok, 1=not monotonous
};

struct globalstuff{
  //parameters:
  double alfabdone;
  int clinicalfactors; //number of clinical factor index 
  int clusternum; //number of clusters 
  double eqd2binwidth;
  vector<double> alfabeta;
  vector<double> nvalue4eud;
  int maxbin;
  int datatype;
  int powptype;
  int twodvh;

  //fitted stuff
  vector<pair<string,string>> fitalgo; //algorithm used by tminimizer for fitting
  map<string, pair<int,vector<double>>> fitpars; //ntcp model parameters key: name, value: 0:number of parameter index, 1=initial value, 2=step, 3=lower, 4=upper (for setlimitedvariable)
  map<int, vector<double>> fitresults; //key: fitalgo index, value: 0=status, 1=CovMatrixStatus, 2=edm, 3=dof, 4=likehood->minvalue, 5=aic, 6=dev/dof, 7=auc, 8=avg_precision
  map<int, map<string,vector<double>>> fittedpar; //key: fitalgo index, value: fitted parameter name, second. 0=value 1=error, 2=likehoodH0->minvalue(H0=LH value w/o the parameter), 3=pvalue
  int bestvalue; //index of the best fit 
};


void fillHisto(map<int, PatientData> &sample,const globalstuff &glbstuff);
void PostLoopAnalysis(map<int, PatientData> &sample, const globalstuff &glbstuff);
void bookHisto(TFile *outrootfile, const map<int, PatientData> &sample, const globalstuff &glbstuff, const bool multivariate);
int CreateNtcpSigmoidalSingle(TString pltname, const TVectorD &vx, const TVectorD &vy);
int CreateNtcpSigmoidalPredefined(map<int, PatientData> &sample, TString tgtname);
void PrintSampleLine(int idx, const   map<int, PatientData> &sample);
int loadDvhFile(const string& filename,   map<int, PatientData> &sample);
int loadMetaFile(const string& filename,   map<int, PatientData> &sample, TString tgtname, const int datatype, const int powptype);
int loadSyntheticFile(const string& filename,   map<int, PatientData> &sample);
void evaluateEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff);
int CheckSampleSamrectConsistency(const map<int, PatientData> &sample, const map<int, PatientData> &samrect);
void FillEqdEud(map<int, PatientData> &sample, const globalstuff &glbstuff);
double CalculateEudFromScratch(const PatientData &paziente, double alfabeta, double nvalue);
double CalculateEudEqdAlreadyDone(const PatientData &paziente, double nvalue);
double fitSigmoidal(TGraph* graph, int parnum, int functype);
string trim(const string& s);
vector<string> splitCsvLine(const string& line, const TString delimiter);
void CreateHistoFromTgraph(TGraph *gr, TH1D *h);
int optimizeLikehood(map<int, PatientData> &sample, globalstuff &glbstuff, const int fitalgindex, map<int, PatientData> &samrect, const pair<int,double> fixedpar=make_pair(-1,-1) );
void optlike_fill(map<int, PatientData> &sample, const globalstuff &glbstuff, int fitalgindex, map<int, PatientData> &samrect);
pair<double,double> optlike_aucROC(const map<int, PatientData> &sample,const globalstuff &glbstuff, const int fitalgindex);
void computeDCT(const vector<double>& x, vector<double>& c);
void DrawLikeHood(std::map<int, PatientData>& sample, const globalstuff& glbstuff);
void fillGlobalStuff(globalstuff &glbstuff, double alfabdone, double eqd2binwidth, const vector<double> &nvalue4eud, const vector<double> &alfabeta, const map<string, pair<int,vector<double>>> &fitpars,   const vector<pair<string,string>> &fitalgo, int datatype, int clinicalfactors, int clusternum, int powptype, int twodvh);
void ChooseBestFit(globalstuff &glbstuff);
void PlotCalibrationCurveQuantilesAndHLtest(const std::map<int, PatientData>& sample, const globalstuff& glbstuff,int fitalgindex, int nbins);
int SetClusterAsClinicalFactor(map<int, PatientData> &sample, const globalstuff &glbstuff);
TGraphErrors *MakeBandFromMinimizer(TF1 *f, ROOT::Math::Minimizer *minimizer, int npoints, double clscale);

  
double functorLikehoodFull(const map<int, PatientData> &sample, const double* par);
double functorLikehoodAlfabdone(const map<int, PatientData> &sample, const double* par);
double functorLikehoodFullClinical_0(const map<int, PatientData> &sample, const double* par);
double functorLikehoodFullClinical_1(const map<int, PatientData> &sample, const double* par);
double functorLikehoodAlfabdoneClinical_0(const map<int, PatientData> &sample, const double* par);
double functorLikehoodAlfabdoneClinical_1(const map<int, PatientData> &sample, const double* par);
double functorLikehoodAlfabdone2DvhClinical_1(const map<int, PatientData> &sample, const map<int, PatientData> &samrect, const double* par);

void SetAucAvgPrec(int index, const pair<double,double> aucavgin, globalstuff& glbstuff);

inline double EqdDose(PatientData &paziente, double alfabeta, double dose){return dose*(alfabeta+dose/paziente.nfraction)/(alfabeta+2.);};

inline double EvalScoreFull(const PatientData& paziente, const double *par){return 1./(1.+exp(-par[0]-par[1]*CalculateEudFromScratch(paziente, par[3], par[2])));};
inline double EvalScoreLikehoodFull(const PatientData& paziente, const double *par){
  return (paziente.tgt_acutegitox<0.5) ? 1.-EvalScoreFull(paziente, par) : EvalScoreFull(paziente,par);};

  inline double EvalScoreAlfabdone(const PatientData& paziente, const double *par){ return 1./(1.+exp(-par[0]-par[1]*CalculateEudEqdAlreadyDone(paziente, par[2])));};
  inline double EvalScoreLikehoodAlfabdone(const PatientData& paziente, const double *par){
    return (paziente.tgt_acutegitox<0.5) ? 1.- EvalScoreAlfabdone(paziente, par): EvalScoreAlfabdone(paziente, par);};

inline double EvalScoreFullClinical_0(const PatientData& paziente, const double *par){return  1./(1.+exp(-par[0]-par[1]*CalculateEudFromScratch(paziente, par[3], par[2])-par[4]*paziente.clinical_factor[0]));};
inline double EvalScoreLikehoodFullClinical_0(const PatientData& paziente, const double *par){
  return (paziente.tgt_acutegitox<0.5) ? 1.- EvalScoreFullClinical_0(paziente, par): EvalScoreFullClinical_0(paziente, par);};

inline double EvalScoreAlfabdoneClinical_0(const PatientData& paziente, const double *par){return 1./(1.+exp(-par[0]-par[1]*CalculateEudEqdAlreadyDone(paziente, par[2])-par[3]*paziente.clinical_factor[0]));};
inline double EvalScoreLikehoodAlfabdoneClinical_0(const PatientData& paziente, const double *par){
  return (paziente.tgt_acutegitox<0.5) ? 1.- EvalScoreAlfabdoneClinical_0(paziente, par):EvalScoreAlfabdoneClinical_0(paziente, par) ;};

inline double EvalScoreFullClinical_1(const PatientData& paziente, const double *par){return  1./(1.+exp(-par[0]-par[1]*CalculateEudFromScratch(paziente, par[3], par[2])-par[4]*paziente.clinical_factor[0]-par[5]*paziente.clinical_factor[1]));};
inline double EvalScoreLikehoodFullClinical_1(const PatientData& paziente, const double *par){
  return (paziente.tgt_acutegitox<0.5) ? 1.- EvalScoreFullClinical_1(paziente, par): EvalScoreFullClinical_1(paziente, par);};

inline double EvalScoreAlfabdoneClinical_1(const PatientData& paziente, const double *par){return 1./(1.+exp(-par[0]-par[1]*CalculateEudEqdAlreadyDone(paziente, par[2])-par[3]*paziente.clinical_factor[0]-par[4]*paziente.clinical_factor[1]));};
inline double EvalScoreLikehoodAlfabdoneClinical_1(const PatientData& paziente, const double *par){
  return (paziente.tgt_acutegitox<0.5) ? 1.- EvalScoreAlfabdoneClinical_1(paziente, par):EvalScoreAlfabdoneClinical_1(paziente, par) ;};

inline double EvalScoreAlfabdone2DvhClinical_1(const PatientData& pazienta, const PatientData& pazientb, const double *par){return 1./(1.+exp(-par[0]-par[1]*CalculateEudEqdAlreadyDone(pazienta, par[2])-par[3]*CalculateEudEqdAlreadyDone(pazientb, par[2])-par[4]*pazienta.clinical_factor[0]-par[5]*pazienta.clinical_factor[1]));};
inline double EvalScoreLikehoodAlfabdone2DvhClinical_1(const PatientData& pazienta, const PatientData& pazientb, const double *par){
  return (pazienta.tgt_acutegitox<0.5) ? 1.- EvalScoreAlfabdone2DvhClinical_1(pazienta, pazientb, par):EvalScoreAlfabdone2DvhClinical_1(pazienta, pazientb, par) ;};


inline double EvalScoreSelector(const globalstuff& glbstuff, const PatientData& paziente, const double *par){
  if(glbstuff.clinicalfactors==0){
    return (glbstuff.alfabdone < 0) ?  EvalScoreFull(paziente, par) : EvalScoreAlfabdone(paziente, par);
  }else if(glbstuff.clinicalfactors==1){
    return (glbstuff.alfabdone < 0) ? EvalScoreFullClinical_0(paziente, par) : EvalScoreAlfabdoneClinical_0(paziente, par);
  }else if(glbstuff.clinicalfactors==2){
    return (glbstuff.alfabdone < 0) ? EvalScoreFullClinical_1(paziente, par) : EvalScoreAlfabdoneClinical_1(paziente, par);
  }
};

inline double functorSelector(const globalstuff& glbstuff, const map<int, PatientData> &sample, const double* par){
  if(glbstuff.clinicalfactors==0){
    return (glbstuff.alfabdone < 0) ?  functorLikehoodFull(sample, par) : functorLikehoodAlfabdone(sample, par);
  }else if(glbstuff.clinicalfactors==1){
    return (glbstuff.alfabdone < 0) ? functorLikehoodFullClinical_0(sample, par) : functorLikehoodAlfabdoneClinical_0(sample, par);
  }else if(glbstuff.clinicalfactors==2){
    return (glbstuff.alfabdone < 0) ? functorLikehoodFullClinical_1(sample, par) : functorLikehoodAlfabdoneClinical_1(sample, par);
  }
};
