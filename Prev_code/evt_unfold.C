// Only for generated new bin arrangement. Similar what did before unolding. 
//run : root -l 
//        .L evt_unfold.C
//genhist()    




//#include "CLHEP/Vector/TwoVector.h"
//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/LorentzVector.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"


#include <TRandom.h>

//#include "RooUnfold.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldResponse.h"
//#include "RooUnfHistoSvd.h"
#include "TLegend.h"
#include "TLatex.h"



//using namespace std;
//using namespace CLHEP;

const double mass_pi=0.139;

bool pedsub=true;
//bool pedsub=false;
const unsigned int njetsetmx=14; //two for reco jets (pfjet, calojet) and
                                //one for MC 
const unsigned int ntrksetmx=2; //Two for reco+gen track
const unsigned int nmcjtsetmx=8; //Number of MC jet combinations

const int njetptcomb= 3; //combinations of Pt thresholds
//const int njetptall = njetptmn; //;+ njetptcomb; //sum of all and excluding first one
const unsigned int njetetamn=2;
const int njetptall =8;
const int ptstrt=0;
const int etastrt=0; //0; //1; //0; //1;

//const int nvarx=1;
//const int varindex[nvarx+1]={15};

const int nvarx=4;
const int varindex[nvarx+1]={3,9,18,24};

// 10, 3,4,5,6,9
const int igenres=10;

//#define JERUP
//#define JERDN



const unsigned nvar=34;

const char* varname[nvar]={"y3anti", "y3ceanti", "y3cranti", "thrustc", "thrustce", "thrustcr",
			     "minorc", "minorce", "minorcr", "t3mass", "t3masse", "t3massr",
			     "h3mass", "h3masse", "h3massr", "y3c", "y3ce", "y3cr",
			     "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
			     "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
			     "sphericity", "cparameter", "p3byp12", "p3byp12c"};

// with eta < 1.3 data using pythia6
/*
  int strtiter[nvar][njetptall+4]={{0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {4,3,3,3,3,4,3,3},  // thrustc
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {6,5,4,3,3,6,4,4},  // minor 
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {8,8,8,9,9,8,7,7}, //t3mass
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {6,5,5,6,6,6,5,5}, //h3mass 
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {7,6,7,9,5,7,6,5},  // y3c
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {4,4,5,5,5,4,4,4}, //broadt
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {4,4,5,5,5,4,4,4}, //broadw 
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {5,4,4,4,4,5,4,3}, //ttmass
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {6,5,5,6,7,6,5,5}, //htmass
  {0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0},
  {6,4,3,3,3,5,4,3},    //sphericity 
  {6,5,4,3,3,6,4,4},    //cparameter
  {0,0,0,0,0,0,0,0},
  {4,4,4,4,3,4,3,3}}; //p3byp12c 
*/
int strtiter[4][3][8]={{{4,4,4,4,5,4,4,3},
			{3,4,5,5,4,5,5,7},
			{4,4,4,4,4,5,5,4}},// thrustc
		       {{4,4,4,4,4,4,4,4},
			{4,4,4,4,4,6,4,5},
			{5,6,6,5,6,6,6,5}},//t3mass
		       //{{6,5,8,4,4,4,4,4},
		       	//{6,5,8,4,4,4,4,4},
			//{6,5,8,4,4,4,4,4}},//y3c
		       {{4,7,7,4,5,5,5,4},
			{3,5,5,6,7,5,4,7},
			{6,6,7,7,6,5,8,5}},//broadt
		       {{4,5,4,5,4,5,5,5},
			{4,5,5,7,8,8,8,7},
			{5,5,4,5,4,5,5,4}}};//ttmass



const char* vartitle[nvar]={"Anti-Y_{23,C} ", "Anti-Y_{23,E} ", "Anti-Y_{23,R} ", 
			    "#tau_{_{#perp} _{   ,C}} ", "#tau_{_{#perp} _{   ,E}} ", "#tau_{_{#perp} _{   ,R}} ",
			    "T_{ m,C} ", "T_{ m,E} ", "T_{ m,R} ",
			    "#rho_{Tot,C} ", "#rho_{Tot,E} ", "#rho_{Tot,R} ",
			    "#rho_{H,C} ", "#rho_{H,E} ", "#rho_{H,R} ",
			    "Y_{23,C} ", "Y_{23,E} ", "Y_{23,R} ",
			    "B_{ T,C} ", "B_{ T,E} ", "B_{ T,R} ", 
			    "B_{ W,C} ", "B_{ W,E} ", "B_{ W,R} ",
			    "#rho^{T}_{Tot,C} ", "#rho^{T}_{Tot,E} ", "#rho^{T}_{Tot,R} ",
			    "#rho^{T}_{H,C} ", "#rho^{T}_{H,E} ", "#rho^{T}_{H,R} ",
			    "S_{_{#perp} _{   ,C}}", "C-parameter_{C}", 
			    "2#times P_{T3}/(P_{T2} + P_{T3})", "2#times P_{T3}/(P_{T2} + P_{T3})c"};



const char* jetsname[6]={ "Jets", "All particle","All Particle" "All particle: P_{T}>0.25", "All particle: P_{T}>0.50", "All particle: P_{T}>1"};
const int norder[6]={0,1,2,3};
const int nextraUnfold=0;//7;
int nAlgoorder[nextraUnfold+1]={0};//,1};//,2,3}; //10,3,4,5,6,9,7,8};

//const char* jetsname[njetsetmx+ntrksetmx] = {"Genjet", "PFjet"};
//const char* jetsname[njetsetmx+ntrksetmx] = {"Genjet", "Jetensmr", "Jetangsmr", "Jetsmr", "JetsmrUp", "JetsmrDown", "JetExtSmr", "JetESPlus", "JetESMinus", "JetExtSmrUp", "PFjet", "PFjetak7", "PFjetkt4", "Calojet","Track", "GenParticle"};

//only GenLevel files
const int ngenfiles=33;

//const char* genrootname[ngenfiles]={"/home/tanmay/QCD/13TeV/miniaod/rootfile/modify/rootfile_formodify_plot/pythia_weight_logfile.root"};//"/home/tanmay/QCD/13TeV/miniaod/merge-rootfile/pythia8_weightadd_25ns_v3.root"}; //home/tanmay/pythia8_flat_25_ns.root"}; 
//const char* genrootname[ngenfiles]={"/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Herwigpp/herwigpp.root"};//"/home/tanmay/QCD/13TeV/miniaod/merge-rootfile/pythia8_weightadd_25ns_v3.root"}; //home/tanmay/pythia8_flat_25_ns.root"}; 
//const char* genrootname[ngenfiles]={"/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_monash_30_default_13tev_30_flat.root","/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_mpioff.root", "/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_haddoff.root", "hist_BeamRemnants_reconnectRange_2.2.root", "hist_qcdevt_pythia8_monash_30_StringZ_bLund_1.2.root", "hist_qcdevt_pythia8_monash_30_StringZ_aLund_0.9.root", "hist_qcdevt_pythia8_monash_30_TimeShower_pTmin_0.6.root", "hist_qcdevt_pythia8_monash_30_TimeShower_pTminChgQ_0.6.root", "hist_qcdevt_pythia8_monash_30_StringPT_sigma_0.38.root", "hist_TimeShower_SpaceShower_alphaSvalue_12285_re.root", "hist_MultipartonInteractions_pT0Ref_2.5.root", "hist_MultipartonInteractions_expPow_2.2.root", "hist_TimeShower_SpaceShower_alphaSvalue_15015_re.root", "hist_qcdevt_pythia8_monash_30_TimeShower_alphaSvalue_12285.root", "hist_qcdevt_pythia8_monash_30_TimeShower_alphaSvalue_15015.root", "hist_qcdevt_pythia8_monash_30_SpaceShower_alphaSvalue_12285.root", "hist_qcdevt_pythia8_monash_30_SpaceShower_alphaSvalue_15015.root"};//"/home/tanmay/QCD/13TeV/miniaod/merge-rootfile/pythia8_weightadd_25ns_v3.root"}; //home/tanmay/pythia8_flat_25_ns.root"}; 
const char* genrootname[ngenfiles+1]={"/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_monash_30_default_13tev_30_flat.root", "/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_mpioff.root", "/home/tanmay/QCD/13TeV/miniaod/combined/Tune/hist_haddoff.root", "hist_BeamRemnants_reconnectRange_2.2.root", "hist_qcdevt_pythia8_monash_30_StringZ_bLund_1.2.root", "hist_qcdevt_pythia8_monash_30_StringZ_aLund_0.9.root", "hist_qcdevt_pythia8_monash_30_TimeShower_pTmin_0.6.root", "hist_qcdevt_pythia8_monash_30_TimeShower_pTminChgQ_0.6.root", "hist_qcdevt_pythia8_monash_30_StringPT_sigma_0.38.root", "hist_TimeShower_SpaceShower_alphaSvalue_12285_re.root", "hist_MultipartonInteractions_pT0Ref_2.5.root", "hist_MultipartonInteractions_expPow_2.2.root", "hist_TimeShower_SpaceShower_alphaSvalue_15015_re.root", "hist_qcdevt_pythia8_monash_30_TimeShower_alphaSvalue_12285.root", "hist_qcdevt_pythia8_monash_30_TimeShower_alphaSvalue_15015.root", "hist_qcdevt_pythia8_monash_30_SpaceShower_alphaSvalue_12285.root", "hist_qcdevt_pythia8_monash_30_SpaceShower_alphaSvalue_15015.root", "hist_qcdevt_pythia8_monash_30_isroff.root", "hist_qcdevt_pythia8_monash_30_fsroff.root", "hist_qcdevt_pythia8_monash_30_remnantsoff.root", "hist_qcdevt_pythia8_monash_30_MultipartonInteractions_pT0Ref_23.root", "hist_qcdevt_pythia8_monash_30_MultipartonInteractions_pT0Ref_2.root", "hist_qcdevt_pythia8_monash_30_MultipartonInteractions_ecmPow_23.root", "hist_qcdevt_pythia8_monash_30_MultipartonInteractions_ecmPow_19.root", "hist_qcdevt_pythia8_monash_30_MultipartonInteractions_expPow_17.root", "hist_qcdevt_pythia8_monash_30_down.root", "hist_qcdevt_pythia8_monash_30_up.root", "hist_qcdevt_pythia8_monash_30_CUETP8M1.root", "pt205.root", "py8_diff_cross.root", "tocheck_pt_inned_generated_sample_pythia_weight_qcd_evt_logfile_generated.root", "pythia_ptbin_weight_logfile_generated_monash.root", "pythis8_flat_CU_tune_modified.root"}; 
//const char* genmodelname[ngenfiles]={"gen_herwigpp"}; 
const char* genmodelname[ngenfiles+1]={"default","mpioff", "hadoff", "BeamRemnants_reconnectRange_22", "StringZ_bLund_12", "StringZ_aLund_9", "TimeShower_pTmin_6", "TimeShower_pTminChgQ_6", "StringPT_sigma_38", "TimeShower_SpaceShower_alphaSvalue_12285", "MultipartonInteractions_pT0Ref_2.5", "MultipartonInteractions_expPow_2.2", "TimeShower_SpaceShower_alphaSvalue_15015", "TimeShower_alphaSvalue_12285", "TimeShower_alphaSvalue_15015", "SpaceShower_alphaSvalue_12285", "SpaceShower_alphaSvalue_15015", "isr_off", "fsr_off", "remnants_off", "MultipartonInteractions_pT0Ref_23", "MultipartonInteractions_pT0Ref_2", "MultipartonInteractions_ecmPow_23", "MultipartonInteractions_ecmPow_19", "MultipartonInteractions_expPow_17", "CUETP8M1_dn", "CUETP8M1_up", "CUETP8M1", "pt8205", "pt8_official", "tocheck_py8_generated_ptin_sample", "py8_monash_ptbin_sample", "py8_205_flat_official_modified"}; 

const int nfiles=4;

TLatex l;
const char* rootname[nfiles]={"/home/tanmay/ttmppythia8/data_2015all.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Pythia8/v1/evt_wpythia_weight_logfile.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Madgh/v1/only_typ0to4_madgh_weight_logfile.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Herwigpp/herwigpp.root"};
const char* modelname[nfiles]={"Data", "Pythia8", "Madgh", "Herwigpp"};
//const char* rootname[nfiles]={"/home/tanmay/ttmppythia8/data_2015all.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Madgh/v1/only_typ0to4_madgh_weight_logfile.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Pythia8/v1/evt_wpythia_weight_logfile.root", "/home/tanmay/QCD/13TeV/miniaod/combined/rootfile/Herwigpp/herwigpp.root"};
//const char* modelname[nfiles]={"Data", "Madgh", "Pythia8", "Herwigpp"};


const int indexx=1; ///5; // 4; //sample to use for response matrix

TFile* fileInput[ngenfiles]; // index is lm 
int nxmod2 = 0;
double xmod2[200];

const int nitermn=2; ///4; // 3; // 4; // 2;
const int nitermx=50; ///5; //20; // 5; //10;

TH1D* h_mcgeninput[ngenfiles]={0};
TH1D* h_recoinput[nfiles]={0};
TH1D* h_recoinputwobkgsub[nfiles]={0};

TH1D* h_unfoldedbayes[nfiles]={0};
TH1D* h_unfoldedsvd[nfiles]={0};

TH2D* h2d_mccorrel=0;
TH1D* h1d_mcreco=0;


TH1D* tmphist2=0;
TH1D* tmphist2x=0;
TH1D* tmphist1=0;
TH1D* tmphist1x=0;

TH1D* h_tmpgen=0;
TH1D* rathist1[nfiles][2] ={0};
TH1D* rathist2[nfiles][2] ={0};
TH1D* rathist3[nfiles][2] ={0};
TH1D* rathist4[nfiles][2] ={0};
TH1D* rathist5[nfiles][2] ={0};
TH1D* rathist6[nfiles][2] ={0};
TH1D* rathist7[nfiles][2] ={0};
TH1D* rathist8[ngenfiles][2] ={0};
TH1D* rathist9[nfiles][2] ={0};
TH1D* rathist10[nfiles][2] ={0};


const char* namex;
char namey[100];
const char* titlex;
char titley[100];

const char* namex2;
char namey2[100];
const char* titlex2;
char titley2[100];


int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return 1000;
}

void calc_deltax(TH1D* thin, double* param) {
  double avedev = 0.0;
  double errsq = 0.0;
  double chisq = 0.0;

  for (int ij=0; ij<thin->GetNbinsX(); ij++) {
    double val = thin->GetBinContent(ij+1);
    if (val >1) val=1./val;
    
    double err = max(min(1., thin->GetBinError(ij+1)), 1.e-6);
    if (err<=1.e-6 || err>=1) continue;
    avedev +=fabs(1 - val)/pow(err, 2.);
    errsq +=pow(1./err, 2.);
    chisq +=pow( (1 - val)/err, 2.);
  }
  avedev /=errsq; avedev = int(1000.*avedev+0.0005)/1000.;
  chisq = int(10.*chisq+0.05)/10.;
  param[0] = avedev;
  param[1] = chisq;

}

bool reject_for_rebin(int ix, int ivar) {
  // ix : bin, starting from underflow
  //ivar : Variable, thrust minor etc
  
  if (ivar==3 || ivar==6) { 
    if (ix<=0 || ix<=3 || ix==5 || ix==6 || ix==8 || ix==30) return true;
  } else if (ivar==9 || ivar==12 || ivar==18 || ivar==21) {
    if (ix==30) return true;
  }
  
  //  if (ix<=0) return true;
  //  if (ix<=4 || ix==6 || ix==7 || ix==8 || ix==10 || ix==11 || ix==13) return true;
  return false;

}


void hist_setbins(TH1D* th, int lm) {
  const int nbmx=100;
  int nbn=th->GetNbinsX();
  if (nbn>nbmx) nbn = nbmx;

  double dx=0.2*(nfiles/2.-lm)*th->GetBinWidth(nbn/2);
  double xval[nbmx];
  for (int ix=0; ix<=nbn; ix++) { //GMAA
    xval[ix] = th->GetBinLowEdge(ix+1)+dx;
  }
  th->SetBins(nbn, xval);
}
//Low PU run 193092

void convert_errorplot(TH1D* thin, int ifil, int ior) {
  int nbins = thin->GetNbinsX();
  for (int ix=1; ix<=nbins; ix++) {
    double val= thin->GetBinContent(ix);
    double err = thin->GetBinError(ix);
    if (err > val) cout <<"==       ================== "<<thin->GetName()<<" "<<ifil<<" "<<ior<<" "<<val <<" "<<err/max(1.e-6,val)<<endl;
    if (val>1.e-6) {
      thin->SetBinContent(ix, err/val);
    } else {
      thin->SetBinContent(ix, 0.0);
    }

    thin->SetBinError(ix, 0.0);
  }
  thin->SetLineColor(ifil);
  thin->SetLineStyle(2*ior+1);
  thin->SetLineWidth(2*ior+1);
}

//double arrayvar_first[5][3][8];
//double arrayvar_last[5][3][8];


int arrayvar_first[4][3][8] = {{{4,6,5,4,2,3,2,4},
                                {9,9,8,8,7,7,6,5},
                                {7,7,6,6,7,7,5,6}},  // thrustc
                               {{5,3,3,4,4,5,6,5},
                                {5,4,4,4,5,4,4,4},
                                {4,3,4,4,4,4,3,3}}, //t3mass
                                //{{11,9,9,8,7,7,6,5},  // For pt0 to pt4 remove one more bins in final plot == done
                               // {7,6,7,9,5,7,6,5},
                               // {7,6,7,9,5,7,6,5}},  // y3c
                               {{10,9,8,8,8,6,5,4},
                                {9,10,9,10,9,7,7,5},
                                {4,6,7,5,5,4,4,4}}, //broadt
                               {{8,7,7,7,7,4,5,5},
                                {10,7,6,5,5,5,4,4},
                                {6,5,6,6,5,5,4,4}}}; //ttmass

int arrayvar_last[4][3][8] =  {{{4,3,3,3,3,4,5,4},
                                {2,2,2,2,2,2,4,6},
                                {2,2,3,3,3,4,5,3}},  // thrustc
                               {{4,4,5,5,8,7,7,7},
                                {4,5,6,7,10,7,8,8},
                                {3,4,4,5,5,5,5,6}}, //t3mass
                               // {{3,6,7,7,5,7,6,5},
                               // {7,6,7,9,5,7,6,5},
                               // {7,6,7,9,5,7,6,5}},  // y3c
                               {{7,6,7,9,5,12,11,10},
                                {3,5,5,5,7,10,10,11},
                                {4,4,3,5,6,6,4,8}}, //broadt
                               {{7,6,7,5,7,7,9,8},
                                {6,6,7,7,9,7,9,9},
                                {5,5,6,5,6,6,7,7}}}; //ttmass



/*
  int arrayvar_first[5][3][8] ={{{1, 1, 1, 1, 1, 2, 1, 1},
  {14, 14, 1, 2, 12, 11, 11, 10},
  {7, 8, 7, 6, 7, 7, 2, 7}},
  {{1, 1, 1, 1, 1, 1, 2, 5},
  {1, 1, 1, 1, 1, 1, 2, 7},
  {2, 2, 3, 1, 1, 2, 1, 1}},
  {{12, 8, 9, 8, 7, 6, 6, 4},
  {5, 8, 4, 7, 4, 2, 2, 1},
  {3, 3, 2, 5, 2, 2, 4, 4}},
  {{11, 10, 9, 8, 7, 6, 6, 2},
  {12, 14, 13, 2, 12, 12, 11, 10},
  {5, 3, 3, 4, 4, 5, 5, 5}},
  {{4, 1, 2, 6, 5, 7, 6, 6},
  {9, 8, 8, 7, 7, 7, 7, 6},
  {5, 5, 4, 4, 5, 2, 4, 4}}};
  
  int arrayvar_last[5][3][8] = {{{31, 32, 32, 32, 31, 32, 32, 32},
  {29, 29, 29, 29, 29, 29, 29, 29},
  {20, 20, 20, 20, 20, 20, 20, 20}},
  {{36, 37, 37, 37, 37, 35, 36, 36},
  {36, 36, 37, 37, 37, 37, 37, 36},
  {18, 18, 17, 17, 17, 17, 18, 17}},
  {{18, 19, 20, 18, 19, 19, 19, 19},
  {16, 16, 17, 16, 17, 17, 16, 17},
  {9, 9, 8, 8, 8, 8, 8, 8}},
  {{31, 32, 32, 32, 32, 31, 31, 31},
  {32, 32, 32, 32, 32, 31, 31, 31},
  {17, 17, 16, 17, 16, 17, 16, 16}},
  {{30, 31, 32, 32, 32, 32, 32, 32},
  {29, 30, 31, 31, 31, 31, 31, 31},
  {16, 16, 18, 17, 17, 17, 17, 17}}};
  
*/

int arrayind[4] ={3,9,18,24};

TH1D* rebin_hist(TH1D* thin, int ijetpt, int ivar, int ityp) {

  //  return (TH1D*)thin->Clone();
  
  TH1D* thout;

  double yvl[100]={0.0};
  double erryvl[100]={0.0};
  int nbinx = thin->GetNbinsX();
  int indx=-1;
  for (int ij=0; ij<4; ij++) {
    //cout << " ivar " << ivar << " ; ij " << ij << " ; arrayind " << arrayind[ij] << endl;
    if (ivar==arrayind[ij]) { indx=ij; break;}
  }
  if (indx<0) return thout;
  
  int ifirst = 0;
  nxmod2= -1;
  if (nxmod2<0) { 
    for (int ij=0; ij<nbinx+2; ij++) {
      //cout<<"NBinX======="<<nbinx<< " ; ij= " << ij <<endl;
      if (ij <=arrayvar_first[indx][ityp][ijetpt]) {
	//cout <<" index= " << indx << " ij =" << ij <<"; ityp= "<< ityp << " ; ijetpt = " << ijetpt <<" iF = " << arrayvar_first[indx][ityp][ijetpt] << endl;
	//        xmod2[ifirst] = thin->GetBinLowEdge(ij);
	//        yvl[ifirst] +=thin->GetBinContent(ij);
	//        erryvl[ifirst] +=thin->GetBinError(ij)*thin->GetBinError(ij);	
      } else if ( ij > (nbinx+2 - arrayvar_last[indx][ityp][ijetpt])) {
	//cout << " ij = " << ij  << " iL = " << arrayvar_last[indx][ityp][ijetpt] << endl;
	//        yvl[ifirst+1] +=thin->GetBinContent(ij);
	//        erryvl[ifirst+1] +=thin->GetBinError(ij)*thin->GetBinError(ij);
      } else {
        xmod2[ifirst] = thin->GetBinLowEdge(ij);
	xmod2[ifirst+1] = thin->GetBinLowEdge(ij+1);
        yvl[ifirst+1] =thin->GetBinContent(ij);
        erryvl[ifirst+1] +=thin->GetBinError(ij)*thin->GetBinError(ij);
        ifirst++;
	  //cout << "ifirst = " << ifirst << "yval ==== " << yvl[ifirst]  << "Bin Content = " << thin->GetBinContent(ij) << " Error= " << erryvl[ifirst] <<endl;
      }
    }
    nxmod2 = ifirst;
    //cout << "New Bin No = " << nxmod2 << endl; 
  }

  namex = thin->GetName();
  sprintf(namey, "Roo_%s", namex);
  titlex = thin->GetTitle();
  sprintf(titley, "Roo_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod2);

  for (int ix=0; ix<nxmod2; ix++) {
    thout->SetBinContent(ix+1, yvl[ix+1]);
    //cout << "Rbin =" << thout->GetBinContent(ix+1);
    thout->SetBinError(ix+1, sqrt(erryvl[ix+1]));
  }

  thout->GetXaxis()->SetTitleFont(42);
  thout->GetXaxis()->SetLabelFont(42);
  thout->GetXaxis()->SetLabelSize(0.07);
  thout->GetXaxis()->SetLabelOffset(.01);

  thout->GetYaxis()->SetTitleFont(42);
  thout->GetYaxis()->SetLabelFont(42);
  thout->GetYaxis()->SetLabelSize(0.07);
  thout->GetYaxis()->SetLabelOffset(.01);

  thout->GetXaxis()->SetTitleOffset(0.9);
  thout->GetXaxis()->SetTitleSize(0.07);

  thout->SetLineWidth(1);
  thout->GetXaxis()->SetTitle(thout->GetTitle());

  return thout;

}

TH2D* rebin_hist2d(TH2D* thin, int ijetpt, int ivar) {
  //  return (TH2D*)thin->Clone();
  
  double zval[100][100]={0};
  double zerr[100][100]={0};
  int nbinx = thin->GetNbinsX();
  int nbiny = thin->GetNbinsY();

  for (int ij=0; ij<nbinx+2; ij++) {
    int ix = getbinid(thin->GetXaxis()->GetBinCenter(ij), nxmod2, xmod2);
    for (int jk=0; jk<nbiny+2; jk++) {
      if (ij==nbinx+1 && jk==nbiny+1) continue;
      int iy = getbinid(thin->GetYaxis()->GetBinCenter(jk), nxmod2, xmod2);
      if (ix<0 || iy<0) continue;
      if (ix>nxmod2 && iy>nxmod2) continue;
      double val=thin->GetBinContent(ij, jk);
      double err=thin->GetBinError(ij, jk);
      if (ix>nbinx+2 && iy<nbinx+2) {
	zval[nxmod2][iy] += val;
	zerr[nxmod2][iy] += err*err;
      } else if (ix<nbinx+2 && iy>nbinx+2) {
	zval[ix][nxmod2] += val;
	zerr[ix][nxmod2] += err*err;
      } else {
	zval[ix][iy] += val;
	zerr[ix][iy] += err*err;
      }
    }
  }
  namex = thin->GetName();
  sprintf(namey, "Roo_%s", namex);
  titlex = thin->GetTitle();
  sprintf(titley, "Roo_%s", titlex);
  TH2D* thout = new TH2D(namey, titley, nxmod2, xmod2, nxmod2, xmod2);
  
  for (int ix=0; ix<=nxmod2; ix++) {
    for (int iy=0; iy<=nxmod2; iy++) {
      thout->SetBinContent(ix+1, iy+1, zval[ix][iy]);
      thout->SetBinError(ix+1, iy+1, sqrt(zerr[ix][iy]));
    }
  }
  
  thout->GetXaxis()->SetTitleFont(42);
  thout->GetXaxis()->SetLabelFont(42);
  thout->GetXaxis()->SetLabelSize(0.07);
  thout->GetXaxis()->SetLabelOffset(.01);

  thout->GetYaxis()->SetTitleFont(42);
  thout->GetYaxis()->SetLabelFont(42);
  thout->GetYaxis()->SetLabelSize(0.07);
  thout->GetYaxis()->SetLabelOffset(.01);
  
  return thout;

}

void normlise_hist(TH1D* hist, double ntot) {
  if (ntot<1.) ntot = 1.;
  double nttox = max(1.,hist->Integral());
  double ratio = ntot/nttox;
  for (int ij=0; ij<hist->GetNbinsX(); ij++) {
    hist->SetBinContent(ij+1, ratio*hist->GetBinContent(ij+1));
    hist->SetBinError(ij+1, ratio*hist->GetBinError(ij+1));
  }
}

void reweight_hist(TH1D* th1, double* effi) {
  for (int ij=0; ij<th1->GetNbinsX(); ij++) {
    th1->SetBinContent(ij+1, th1->GetBinContent(ij+1)/max(0.1,effi[ij])); 
    //    th1->SetBinContent(ij+1, (th1->GetBinContent(ij+1))*(1+effi[ij]));
  }
}

void add_index(int ijk, TNamed* hist, int ndim) {

  namex = hist->GetName();
  sprintf(namey, "%i_%s", ijk, namex);
  hist->SetName(namey);
  
  titlex = hist->GetTitle();
  sprintf(titley, "%i_%s", ijk, titlex);
  hist->SetTitle(titley);
}

void add_nameindex(TNamed* hist, char ttl[]) {

  namex = hist->GetName();
  sprintf(namey, "%s_%s", ttl, namex);
  hist->SetName(namey);
  
  titlex = hist->GetTitle();
  sprintf(titley, "%s_%s", ttl, titlex);
  hist->SetTitle(titley);
}


void calc_chisq(double* chi2, TH1D* h_data, TH1D* h_mc) {
  const int nbinmx=100;
  double vdata[nbinmx]={0.};
  double vmc[nbinmx]={0.};
  double totdata=0.;
  double totmc = 0.;
  double errdata[nbinmx]={0.};
  double errmc[nbinmx]={0.};  
  double toterr[nbinmx]={0.}; 
  int nbin = h_mc->GetNbinsX();
  
  int nbndata = h_data->GetNbinsX();
  if (nbin <nbinmx) {
    for (int ij=0; ij<nbin; ij++) {
      
      vdata[ij] = h_data->GetBinContent(ij+1);
      totdata += vdata[ij]; 
      vmc[ij] = h_mc->GetBinContent(ij+1);
      totmc +=vmc[ij];

      errdata[ij] = h_data->GetBinError(ij+1);
      errmc[ij] = h_mc->GetBinError(ij+1);
    }

    if (totmc <1.) totmc = 1.;
    if (totdata <1.) totdata = 1.;
    for (int ij=0; ij<nbin; ij++) {
      vdata[ij] /=totdata;
      vmc[ij] /=totmc;
      errdata[ij] /=totdata;
      errmc[ij] /=totmc;
      toterr[ij] = sqrt(errdata[ij]*errdata[ij] + errmc[ij]*errmc[ij]);
    } 

    //Calculate chisquare
    chi2[0] =  chi2[1] = chi2[2] = 0;
    for (int ij=0; ij<nbin; ij++) {
      chi2[0] += pow((vdata[ij] - vmc[ij]), 2.0); // /max(1.e-6,errmc[ij]), 2.0); //17th Oct 2010
      chi2[1] += pow((vdata[ij] - vmc[ij])/max(1.e-6,errdata[ij]), 2.0);
      chi2[2] += pow((vdata[ij] - vmc[ij])/max(1.e-6,toterr[ij]), 2.0);
    }
  } else {
    cout <<"Increase nbinmx, which is equal to "<<nbinmx<<" less than histogrammes binning "<<nbin<<endl;

  }
}

void calc_chisq_svd(double* chi2, TH1D* h_unfold, TMatrixD cov, TMatrixD cov2,TH1D* h_mc) {

  double ratio = h_unfold->Integral()/max(1., h_mc->Integral());
  const int nbinmx=100;
  double vdata[nbinmx]={0.};
  double vmc[nbinmx]={0.};
  double errdata[nbinmx]={0.};
  double errmc[nbinmx]={0.};  
  double toterr[nbinmx]={0.}; 
  int nbin = h_mc->GetNbinsX();
  TVectorD diffdmc(nbin);

  int nbndata = h_unfold->GetNbinsX();
  if (nbin <nbinmx) {
    for (int ij=0; ij<nbin; ij++) {
      
      vdata[ij] = h_unfold->GetBinContent(ij+1);
      vmc[ij] = ratio*h_mc->GetBinContent(ij+1);

      errdata[ij] = h_unfold->GetBinError(ij+1);
      errmc[ij] = ratio*h_mc->GetBinError(ij+1);
      diffdmc(ij) = vdata[ij] - vmc[ij];
      toterr[ij] = sqrt(errdata[ij]*errdata[ij] + errmc[ij]*errmc[ij]);
    }
    //Calculate chisquare
    chi2[0] =  chi2[1] = chi2[2] = chi2[3] = chi2[4]=0;
    for (int ij=0; ij<nbin; ij++) {
      chi2[0] += pow((vdata[ij] - vmc[ij]), 2.0); // /max(1.e-6,errmc[ij]), 2.0); //17th Oct 2010
      chi2[1] += pow((vdata[ij] - vmc[ij])/max(1.e-6,errdata[ij]), 2.0);
      chi2[2] += pow((vdata[ij] - vmc[ij])/max(1.e-6,toterr[ij]), 2.0);
    } 

    int ncols=cov.GetNcols();

    for (int ix=0; ix<ncols; ix++) {
      cov(ix,ix) += errmc[ix]*errmc[ix];
      cov2(ix,ix) += errmc[ix]*errmc[ix];
    }

    cov.Invert();
    cov2.Invert();

    TVectorD covdx = cov*diffdmc;
    for (int ix=0; ix<cov.GetNcols(); ix++) { chi2[3] +=diffdmc(ix)*covdx(ix);}
    
    TVectorD covdx2 = cov2*diffdmc;
    for (int ix=0; ix< ncols; ix++) { chi2[4] +=diffdmc(ix)*covdx2(ix);}

  } else {
    cout <<"Increase nbinmx, which is equal to "<<nbinmx<<" less than histogrammes binning "<<nbin<<endl;

  }
}


double calc_chisq_sss(TH1D* h_unfold, TH1D* h_mc, TMatrixD cov) {

  double ratio = h_unfold->Integral()/max(1., h_mc->Integral());
  const int nbinmx=100;
  double vdata[nbinmx]={0.};
  double vmc[nbinmx]={0.};
  double errdata[nbinmx]={0.};
  double errmc[nbinmx]={0.};  
  double toterr[nbinmx]={0.}; 
  int nbin = h_mc->GetNbinsX();
  TVectorD diffdmc(nbin);

  int nbndata = h_unfold->GetNbinsX();
  double chi2 = 0.0;
  if (nbin <nbinmx) {
    for (int ij=0; ij<nbin; ij++) {
      
      vdata[ij] = h_unfold->GetBinContent(ij+1);
      vmc[ij] = ratio*h_mc->GetBinContent(ij+1);
      diffdmc(ij) = vdata[ij] - vmc[ij];
    }
    //Calculate chisquare
    
    TVectorD covdx = cov*diffdmc;
    for (int ix=0; ix<cov.GetNcols(); ix++) { chi2 +=diffdmc(ix)*covdx(ix);}
    
  } else {
    cout <<"Increase nbinmx, which is equal to "<<nbinmx<<" less than histogrammes binning "<<nbin<<endl;

  }
  return int(1000.*(chi2+0.0005))/1000.;
}



void calc_deltax_cov(double* delta, TH1D* h1, TH1D* h2, TMatrixD cov1, TMatrixD cov2) {

  //  int ncols=cov1.GetNcols();
  //  TVectorD diffdmc(ncols);
  //  TVectorD mcerror2(ncols);
  double totdata=0.;
  double totmc = 0.;

  const int nbinmx=100;
  double vdata[nbinmx]={0.};
  double vmc[nbinmx]={0.};

  double errdata[nbinmx]={0.};
  double errmc[nbinmx]={0.};  
  double errmc2[nbinmx]={0.};  
  double toterr[nbinmx]={0.}; 
  double totraterr[nbinmx]={0.}; 
  double rat[nbinmx]={0.};

  int nbin = h1->GetNbinsX();
  
  int nbndata = h1->GetNbinsX();
  TVectorD diffdmc(nbin);
  if (nbin <nbinmx) {
    for (int ij=0; ij<nbin; ij++) {
      
      vdata[ij] = h1->GetBinContent(ij+1);
      totdata += vdata[ij]; 
      vmc[ij] = h2->GetBinContent(ij+1);
      totmc +=vmc[ij];

      errdata[ij] = h1->GetBinError(ij+1);
      errmc[ij] = h2->GetBinError(ij+1);
    }

    if (totmc <1.) totmc = 1.;
    if (totdata <1.) totdata = 1.;

    double ratio = totdata/ max(1., totmc);
    for (int ij=0; ij<nbin; ij++) {
      diffdmc(ij) = vdata[ij] - ratio*vmc[ij];
      errmc2[ij] = ratio*errmc[ij]*ratio*errmc[ij];
      toterr[ij] = sqrt(errdata[ij]*errdata[ij] + errmc2[ij]);
    } 

    //Calculate chisquare
    //    double errsq = 0.;
    double chisq = 0.;

    for (int ij=0; ij<nbin; ij++) {
      if (errdata[ij] >1.e-6 || errmc[ij]>1.e-6) {
	
	chisq += pow( (vdata[ij] - ratio*vmc[ij])/max(1.e-6,toterr[ij]), 2.0);
      }
    }
    chisq = int(10.*(chisq+0.05))/10.;
    
    //    delta[0] = avedev;
    delta[1] = chisq;

    if (cov1.GetNcols() ==0) {
      int ncols = h1->GetNbinsX(); 
      cov1.ResizeTo(ncols, ncols);
      for (int ix=0; ix<ncols; ix++) {
	cov1(ix, ix) =errdata[ix]*errdata[ix];
      }
    }
    
    int ncols=cov2.GetNcols();
    
    if (ncols==0) {
      ncols = cov1.GetNcols();
      cov2.ResizeTo(ncols, ncols);
      for (int ix=0; ix<ncols; ix++) {
	cov2(ix, ix) =errmc2[ix];
      }
    } else { //If defined, just modified with ratio
      for (int ix=0; ix< ncols; ix++) {
	for (int iy=0; iy< ncols; iy++) {
	  cov2(ix, iy) *=ratio*ratio;
	}
      }
    }
    
    TMatrixD comcovariance = cov1 + cov2;

    comcovariance.Invert();
    
    TVectorD covdx = comcovariance*diffdmc;
    delta[2] = 0;
    for (int ix=0; ix< ncols; ix++) { delta[2] +=diffdmc(ix)*covdx(ix);}
    delta[2]= int(10.*delta[2]+0.05)/10.;

    for (int ij=0; ij<nbin; ij++) {
      vdata[ij] /=totdata;
      vmc[ij] /=totmc;
      errdata[ij] /=totdata;
      errmc[ij] /=totmc;
      
      rat[ij] = vdata[ij]/max(1.e-6,abs(vmc[ij]));	

      if (rat[ij] >1.0) rat[ij] = 1/rat[ij];
      if (rat[ij] <=0.0) rat[ij] = 1.e-6;

      totraterr[ij] = rat[ij] * sqrt (pow(errdata[ij]/max(1.e-6, abs(vdata[ij])),2.0) + pow(errmc[ij]/max(1.e-6, abs(vmc[ij])),2.0));
    } 

    //Calculate chisquare
    double avedev = 0.;
    double errsq = 0.;
    
    for (int ij=0; ij<nbin; ij++) {
      if ((vdata[ij] >1.e-6 && vmc[ij]>1.e-6) && (errdata[ij] >1.e-6 || errmc[ij]>1.e-6) && totraterr[ij] >1.e-6) {
	avedev +=abs( 1 - rat[ij])/pow(max(1.e-6,totraterr[ij]), 2.0);
	errsq +=pow(1./max(1.e-6,totraterr[ij]), 2.0);
      }
    }
    avedev /=errsq; avedev = int(1000.*avedev+0.0005)/1000.;
    
    delta[0] = avedev;
  } else {
    cout <<"Increase nbinmx, which is equal to "<<nbinmx<<" less than histogrammes binning "<<nbin<<endl;
  }
}

int subtract_background(TH2D* h2d_correl, TH1D* data, TH1D* mc1, TH1D* mc2, TH1D*mc3, TH1D* mc4, TH1D* mc5, double* fakerate, double* effi, int jk) {
  int nbinx = h2d_correl->GetNbinsX();
  int nbiny = h2d_correl->GetNbinsY();
  const int nbinmx=100;
  if (nbinx>nbinmx || nbiny > nbinmx) {
    cout <<"Increase nbinmx, which is "<< nbinmx<<" in subtract_background"<<endl;
    cout <<"Reconstructed objects are not pedestal subtracted"<<endl;
    return 0;
  }
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};
  if (jk==0) {
    for (int ix=0; ix<nbinx+1; ix++) {
      for (int iy=0; iy<nbiny+1; iy++) {
	if (ix==nbinx && iy==nbiny) continue;
	totalreco[ix] +=h2d_correl->GetBinContent(ix+1, iy+1);
	if (iy==nbiny) fakerate[ix] =h2d_correl->GetBinContent(ix+1, iy+1);
	
	totalgen[iy] +=h2d_correl->GetBinContent(ix+1, iy+1);
	if (ix==nbinx) effi[iy] =h2d_correl->GetBinContent(ix+1, iy+1);
	
	//      cout <<"ix "<< ix<<" "<<iy<<" "<<h2d_correl->GetBinContent(ix+1, iy+1)<<" "<<totalgen[iy]<<" "<<fakerate[ix]<<" "<<totalreco[ix]<<" "<<effi[iy]<<endl;
	//GMA 27th Jan 2013
	if (ix==nbinx || iy==nbiny) {
	  h2d_correl->SetBinContent(ix+1, iy+1, 0.0);
	  h2d_correl->SetBinError(ix+1, iy+1, 0.0);
	}
      }
    }
    
    for (int iy=0; iy<nbiny; iy++) {
      effi[iy] = (totalgen[iy] - effi[iy])/max(1.0, totalgen[iy]);
      if (effi[iy]<1.e-3) effi[iy]=1.e-3;
    }
    for (int ix=0; ix<nbinx; ix++){
      //    cout <<"fake "<< ix<<" "<< fakerate[ix]<<" "<<totalreco[ix]<<" "<<mc1->GetBinContent(ix+1)<<" "<<mc2->GetBinContent(ix+1)<<" "<<mc1->GetBinError(ix+1)<<" "<<mc2->GetBinError(ix+1)<<endl;
      
      fakerate[ix] /=max(1.0, totalreco[ix]);
    }
  }
  
  //121222
  if (pedsub) { 
    for (int ix=0; ix<nbinx; ix++) {
      if (jk<6) { 
	mc1->SetBinContent(ix+1, (1.- fakerate[ix])*mc1->GetBinContent(ix+1));
	mc1->SetBinError(ix+1, (1.- fakerate[ix])*mc1->GetBinError(ix+1));

	if (nfiles>=3) { 
	  mc2->SetBinContent(ix+1, (1.- fakerate[ix])*mc2->GetBinContent(ix+1));
	  mc2->SetBinError(ix+1, (1.- fakerate[ix])*mc2->GetBinError(ix+1));
	}	

	if (nfiles>=4) {
	  mc3->SetBinContent(ix+1, (1.- fakerate[ix])*mc3->GetBinContent(ix+1));
	  mc3->SetBinError(ix+1, (1.- fakerate[ix])*mc3->GetBinError(ix+1));
	}
	if (nfiles>=5) {
	  mc4->SetBinContent(ix+1, (1.- fakerate[ix])*mc4->GetBinContent(ix+1));
	  mc4->SetBinError(ix+1, (1.- fakerate[ix])*mc4->GetBinError(ix+1));
	}
	if (nfiles>=6) {
	  mc5->SetBinContent(ix+1, (1.- fakerate[ix])*mc5->GetBinContent(ix+1));
	  mc5->SetBinError(ix+1, (1.- fakerate[ix])*mc5->GetBinError(ix+1));
	}	
      } 
      if (jk==0 || jk==6 || jk==7) {
	data->SetBinContent(ix+1, (1.- fakerate[ix])*data->GetBinContent(ix+1));
	data->SetBinError(ix+1, (1.- fakerate[ix])*data->GetBinError(ix+1));
	//    cout <<"ix2   "<< ix<<" "<<fakerate[ix]<<" "<<mcreco->GetBinContent(ix+1)<<" "<<data->GetBinContent(ix+1)<<" "<<mcreco->GetBinError(ix+1)<<" "<<data->GetBinError(ix+1)<<endl;
      }
    }
  }
  return 1;
}



/*
void fill_histogramme(int lm, int jk, int ij, int mn, int kl, int nmcfill,int ii, int ijer) {
  // lm : input file number
  // jk : nextraUnfold
  // ij : jetptbn
  // mn : jetetamn
  // kl : Index of eventshape variables, for the time being, 3 or 6
  // ii : Typ 
  // nmcfill : fill Geninformation only once for two algorithms
  // ijer: jer up down index

  char histname[100];
  fileInput[lm]->cd();
  
  //  Only Reco input;
  // sprintf(histname, "%s_%s_c%i_e%i", varname[kl], jetsname[nAlgoorder[jk]], ij, mn);

  // sprintf(histname,"analyzeBasicPat/reco_typ_%i_pt%i_eta%i_%i", jk, ij, mn, kl);
  if(lm==0 || lm==3)   sprintf(histname,"analyzeBasicPat/reco_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl); 
  //else sprintf(histname,"reco_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
  else {
#ifdef JERUP     
    sprintf(histname,"recoreso_typ_%i_pt%i_eta%i_%i_%i", ii, ij, mn, kl,ijer);
#elif defined(JERDN)
    sprintf(histname,"recoreso_typ_%i_pt%i_eta%i_%i_%i", ii, ij, mn, kl,ijer);
#else
    sprintf(histname,"reco_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
#endif	
  }
  //if(lm==0 ) sprintf(histname,"analyzeBasicPat/reco_typ_%i_pt%i_eta%i_%i", jk, ij, mn, kl);

  cout << " OKKKK 1" << endl;
  // else sprintf(histname,"reco_typ_%i_pt%i_eta%i_%i", jk, ij, mn, kl); 
  //else  sprintf(histname,"reco_typ_%i_pt%i_eta%i_%i", norder[ii], ij, mn, kl);
  cout << "histname RECO ======="<<histname<<endl;
  tmphist2 = (TH1D*)fileInput[lm]->Get(histname);
  tmphist2x = (TH1D*)tmphist2->Clone();
  namex = tmphist2x->GetName();
  sprintf(namey, "%s_%s", modelname[lm], namex);
  tmphist2x->SetName(namey);
  
  titlex = tmphist2x->GetTitle();
  

  sprintf(titley, "%s_%s", modelname[lm], titlex);
  tmphist2x->SetTitle(titley);

  //  if (rebin>0) {tmphist2x->Rebin(rebin);}
  h_recoinput[lm] = rebin_hist(tmphist2x, ij, kl, ii);
  h_recoinputwobkgsub[lm] = (TH1D*)h_recoinput[lm]->Clone();
  //Generated input
 
  cout << " OKKKK 1" << endl;
 
  if (lm>0 && nmcfill==0) { //Only for MC and and only once for two jet types
    //  sprintf(histname, "%s_%s_c%i_e%i", varname[kl], jetsname[0], ij, mn);
    // sprintf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", norder[ii], ij, mn, kl); 
    //sprintf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", norder[ii], ij, mn, kl);
    //    sprintf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", jk, ij, mn, kl);
    if(lm==3) sprintf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
    else sprintf(histname,"gen_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
    cout << "histname GEN ========= "<<histname<<endl;
    tmphist1 = (TH1D*)fileInput[lm]->Get(histname);
    tmphist1x = (TH1D*)tmphist1->Clone();

    namex = tmphist1x->GetName();
    sprintf(namey, "%s_%s", modelname[lm], namex);
    tmphist1x->SetName(namey);
    titlex = tmphist1x->GetTitle();
    sprintf(titley, "%s_%s", modelname[lm], titlex);
    tmphist1x->SetTitle(titley);
    
    //    if (rebin>0) {tmphist1x->Rebin(rebin);}
    h_mcgeninput[lm] = (TH1D*)rebin_hist(tmphist1x, ij, kl, ii);
  }
  //Response matrix;

  if (lm==indexx && jk==0) { //2d only from PFjet
    
    //sprintf(histname, "%s_%s_Gen_c%i_e%i", varname[kl], jetsname[nAlgoorder[jk]], ij, mn);
    
    if(lm==3) sprintf(histname,"analyzeBasicPat/corr_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
    else sprintf(histname,"corr_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
    cout << "histname 2D"<<histname<< " ;ii = " << ii <<endl;
    //if (lm==1)   sprintf(histname,"corr_typ_%i_pt%i_eta%i_%i", jk, ij, mn, kl);

    //    sprintf(histname, "%s_%s_Gen_c%i_e%i", varname[kl], jetsname[igenres], ij, mn);
    //    if (rebin>0) { 
    //      h2d_mccorrel = rebin_hist2d ( (TH2D*)((TH2D*)fileInput[lm]->Get(histname))->Rebin2D(rebin, rebin), ij, kl);
    //    } else {
    h2d_mccorrel = rebin_hist2d((TH2D*)fileInput[lm]->Get(histname), ij, kl);
    // Study of revers smearing (Niki)
    //    sprintf(histname, "%s_%s_c%i_e%i", varname[kl], jetsname[igenres], ij, mn);
    //    h1d_mcreco = rebin_hist((TH1D*)fileInput[lm]->Get(histname), ij, kl);
    
    //    }
  }

  if (tmphist1) { delete tmphist1; tmphist1=0;}
  if (tmphist1x) { delete tmphist1x; tmphist1x=0;}
  
  if (tmphist2) { delete tmphist2; tmphist2=0;}
  if (tmphist2x) { delete tmphist2x; tmphist2x=0;}
}*/

void fill_genhistogramme(int lm, int ij, int mn, int kl, int ii) {
  // lm : input file number
  // ij : jetptbn
  // mn : jetetamn
  // kl : Index of eventshape variables, for the time being, 3 or 6
 // ii typ 
  char histname[100];
  fileInput[lm]->cd();
  
  //sprintf(histname, "%s_%s_c%i_e%i", varname[kl], jetsname[0], ij, mn);
   if(lm<30) sprintf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
   else if(lm==34) printf(histname,"analyzeBasicPat/gen_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
   else sprintf(histname,"gen_typ_%i_pt%i_eta%i_%i", ii, ij, mn, kl);
  cout << " HISTNAME==== " << histname << " ; lm = " <<lm << endl;
  cout << " Root FILE NAME==== " << fileInput[lm]->GetName() << endl;
  tmphist1 = (TH1D*)fileInput[lm]->Get(histname);
  tmphist1x = (TH1D*)tmphist1->Clone();
  
  namex = tmphist1x->GetName();
  sprintf(namey, "%s_%s", genmodelname[lm], namex);
  tmphist1x->SetName(namey);
  titlex = tmphist1x->GetTitle();
  sprintf(titley, "%s_%s", genmodelname[lm], titlex);
  tmphist1x->SetTitle(titley);
  
  //    if (rebin>0) {tmphist1x->Rebin(rebin);}
  h_mcgeninput[lm] = (TH1D*)rebin_hist(tmphist1x, ij, kl, ii);
  if (tmphist1) { delete tmphist1; tmphist1=0;}
  if (tmphist1x) { delete tmphist1x; tmphist1x=0;}
  
}


void hist_labels(TH1D* hist, int icol) {

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.08);
  hist->GetXaxis()->SetLabelOffset(.01);

  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.08);
  hist->GetYaxis()->SetLabelOffset(.01);

  hist->SetMarkerColor(icol);
  hist->SetMarkerStyle(22+icol);
  hist->SetMarkerSize(0.44);

}

int genhist() { //int main(argc, char**argv)
  
  gStyle->SetPaintTextFormat("5.2e"); //4.1f    
  //   minorc_JetJPT_Gen_c1_e1->Draw("text0")

  gStyle->SetOptLogy(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetStatColor(10);
  
  gStyle->SetCanvasColor(10);
  gStyle->SetOptStat(0); //1110);
  gStyle->SetOptTitle(1);

  gStyle->SetTitleW(.96);
  gStyle->SetTitleH(.06);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleX(.02);
  gStyle->SetTitleAlign(13);
  gStyle->SetTitleColor(10);
  //  gStyle->SetTitleOffset(-0.05);
  gStyle->SetTitleBorderSize(0); //1);
  gStyle->SetTitleFontSize(0.10);

  gStyle->SetPalette(1,0);
  gStyle->SetPadColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatBorderSize(1);

  gStyle->SetStatStyle(1001);
  gStyle->SetOptFit(101);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);

  gStyle->SetStatX(.99);
  gStyle->SetStatY(.99);
  gStyle->SetStatW(.45);
  gStyle->SetStatH(.16);
  gStyle->SetLabelSize(0.095,"XY");  
  gStyle->SetLabelOffset(0.21,"XYZ");
  gStyle->SetTitleSize(0.065,"XY");  
  gStyle->SetTitleOffset(0.06,"XYZ");
  gStyle->SetPadTopMargin(0.075); //2); //0.11 //0.09
  gStyle->SetPadBottomMargin(0.075);
  gStyle->SetPadLeftMargin(0.09);
  gStyle->SetPadRightMargin(0.02);
  //  gStyle->SetPadGridX(1);
  //  gStyle->SetPadGridY(1);
  //  gStyle->SetGridStyle(3);
  //  gStyle->SetNdivisions(606,"XY");

  gStyle->SetMarkerSize(0.44);
  gStyle->SetMarkerColor(2);
  gStyle->SetMarkerStyle(20);


  char mcrootfile[100];
  char datarootfile[100];
  char hist_mc[100];
  char hist_reco[100];
  char hist_unfol[100];


  //24th Jan 2013
  //Number of iterations/regularisation paramters in bayes and SVD 
  double dataerr[100]; //be sure that 
  double datval[100]; 
  double mxdaterr=0.0;
  double relsumerr=0.0;
  double relsum2err=0.0;

  int    mxdatbin=-1;
  TCanvas *c1;
  //  sprintf(titley, "_ab_%s_%s", modelname[indexx], varname[varindex[0]]);
  //  sprintf(titley, "_hhi3_%i_%i_%s_%i_c%i_e%i_130521_b", igenres, varindex[0], modelname[indexx], int(pedsub), ptstrt, etastrt);  
  sprintf(titley, "13Tev_%i_%i_%s_%i_c%i_e13_130529_kcov", igenres, varindex[0], modelname[indexx], int(pedsub), ptstrt);
  sprintf(namey, "outfile%s.txt", titley);
  sprintf(namey2, "outfile%s.inc", titley);  

  ofstream file_out(namey); // "outfilex.txt");
  ofstream file_out2(namey2);
  
  ofstream file_out3("chi.log");


  sprintf(namey, "outfile%s.ps", titley);
  TPostScript ps(namey, 111); // "outfilex.ps",111);  
  ps.Range(20,30); //ps.Range(10,20);
  ps.NewPage();  

  sprintf(namey, "outfile%s.root", titley);
  TFile* fileOut = new TFile(namey, "recreate");

  /*for (int ij=0; ij<nfiles; ij++) {
    fileInput[ij] = new TFile(rootname[ij], "read");
  }*/
  for (int ij=0; ij<ngenfiles; ij++) {
    fileInput[ij] = new TFile(genrootname[ij], "read");
  }
  
  //  TFile* fileInData = new TFile("all_run2010a_ch.root", "read"); //jetmetpromt4_upto144114_cf.root", "read");

  char histname[100];
  char histname2[100]; //for isdouble rat
  //  int rebin=3;


 

// for(int ii=0; ii<5 ; ii++){

// file_out3<<"                           " <<jetsname[ii] << endl;
 file_out3<< "                     #Delta            and                #chi^2"  << endl;

// c1 = new TCanvas("c1", "Statistics and efficiency", 600, 850);
// c1->cd(); 
 /*TPaveText pt(.1,.5,.9,.9);
 pt.AddText("New");
 pt.SetLabel("Born equation");
 pt.Draw(); */
 //TCanvas *c2 = new TCanvas("c2"); 
 //TLatex l;
 //l.DrawLatex(0.5,0.95," name ");
 //c1->Print("latex2.ps");
 // ps.NewPage();
 //c2->Print("latex2.ps");
 // c1->Update();


#ifdef JERUP
 int ijer=1;
#elif defined(JERDN)
 int ijer=2;
#else
 int ijer=0;
#endif 
 
 for (int  ij=ptstrt;  ij<njetptall; ij++) {//Pt bin
//      for (int  ij=4;  ij<5; ij++) {//Pt bin
   
   for (int mn=etastrt+1; mn<etastrt+2; mn++) { //eta bin
     //    for (int mn=etastrt; mn<etastrt+1; mn++) { //eta bin
     
     cout<< " Eta=======  " <<mn <<endl;
     
     for (int klx=0; klx<nvarx; klx++) { //variables
//     for (int klx=0; klx<1; klx++) { //variables
       int kl = varindex[klx];
       
       for(int ii=0; ii<3 ; ii++){//njtype Jets , Particles etc
	 
	 file_out3<<"                           " <<jetsname[ii] << endl;
	 file_out3<< "                     #Delta            and                #chi^2"  << endl;
	 
 
	 
	 
	 int natype=0; //fill Geninformation only once for four algorithm
	 
	 //Inintialization for lower and upper range
	 nxmod2 =-1;
	 
	     for (int lm=0; lm<ngenfiles; lm++) {
	       fill_genhistogramme(lm, ij, mn, kl,ii);
	     //  fill_genhistogramme(lm, ij, mn, kl, ii);
	     }
	   //}
		 
                for (int lm=0; lm<ngenfiles; lm++) {
                 cout << " OK got it 1" << endl; 
		  rathist8[lm][1] = (TH1D*)h_mcgeninput[lm]->Clone();
		  //		  for (int iy=0; iy<rathist8[lm][1]->GetNbinsX()+2; iy++) {rathist8[lm][1]->SetBinError(iy, 0.0);}
                 }
		 cout << " OK got it 2" << endl; 

	    //c1->Update();
	    //ps.NewPage();  


	    fileOut->cd();
             
            for (int lm=0; lm<ngenfiles; lm++) {
                if (rathist8[lm][1]) rathist8[lm][1]->Write();
            }
	    
	    for (int lm=0; lm<ngenfiles; lm++) {
	      for (int ix=0; ix<2; ix++) {
		if (rathist8[lm][ix]) {delete rathist8[lm][ix]; rathist8[lm][ix] = 0;}
	      }
	    }

	 for (int lm=0; lm<ngenfiles; lm++) {
	   if (h_mcgeninput[lm]) { delete h_mcgeninput[lm]; h_mcgeninput[lm]=0;}
	 }
       } //for (int klx=0; kl<nvarx; klx++) 
     } // for (int mn=0; mn<njetetamn; mn++) 
   } //  for (int ij=0; ij<njetptmn; ij++) 
   //file_out4<< jetsname[ii+1]  << " & " <<" &  "   <<"& " << "& \\\\"<<endl;
   //              file_out4<< "\\hline"<< endl;

 }


 
 ps.Close();
 //file_out3->Close();
}

/*

    vx = CompProd(vw, vxini);
      //TMatrixD X_inv1(_covmatrix.GetNrows(),_covmatrix.GetNrows());
    
      TMatrixD vrm1(_covmatrix.GetNrows(),_covmatrix.GetNrows());
      for(int i=0;i<vrm1.GetNcols();i++){
        for(int j=0;j<vrm1.GetNcols();j++){
          vrm1(i,j)=0.0;
          for(int k=0;k<vrm1.GetNcols();k++){
            vrm1(i,j)=vrm1(i,j)+Vort(k,i)*mC(k,j)/_xini->GetBinContent(j+1);
          }
        }

      _covmatrix_inverse.ResizeTo(_covmatrix.GetNrows(),_covmatrix.GetNrows());
      TMatrixD X_inv1(_covmatrix.GetNrows(),_covmatrix.GetNrows());
      for(int i=0;i<X_inv1.GetNcols();i++){
        for(int j=0;j<X_inv1.GetNcols();j++){
          X_inv1(i,j)=0.;
          _covmatrix_inverse(i,j)=0;
          for(int k=0;k<X_inv1.GetNcols();k++){
            X_inv1(i,j)+=vrm1(k,i)*vrm1(k,j)*ASV(k)*ASV(k);
            //X_inv2(i,j)+=1/(vxini->GetBinContent(i)*vxini->GetBinContent(j))*mA(k,i)*mA(k,j);
            _covmatrix_inverse(i,j)+=1/(vxini(i)*vxini(j))*mA(k,i)*mA(k,j);
          }
        }
      }
      // _covmatrix();


ij 15 -6.26667 -6.13333 0 0 0
ij 16 -6 -5.86667 5 2.23607 5
ij 17 -5.73333 -5.6 138 11.7473 138
ij 18 -5.46667 -5.33333 559 23.6432 559
ij 19 -5.2 -5.06667 1393 37.3229 1393
ij 20 -4.93333 -4.8 2449 49.4874 2449
ij 21 -4.66667 -4.53333 3372 58.0689 3372
ij 22 -4.4 -4.26667 3837 61.9435 3837
ij 23 -4.13333 -4 3777 61.4573 3837
ij 24 -3.86667 -3.73333 3503 59.1861 3837
ij 25 -3.6 -3.46667 2699 51.9519 3837
ij 26 -3.33333 -3.2 1775 42.1307 3837
ij 27 -3.06667 -2.93333 775 27.8388 3837
ij 28 -2.8 -2.66667 270 16.4317 3837
ij 29 -2.53333 -2.4 17 4.12311 3837
ij 30 -2.26667 -2.13333 0 0 3837
ij 31 -2 -1.86667 0 0 3837

nxmod2++++++++++++++++++++++++ 14
ij 0 -6.26667 -6.13333 0 0
ij 1 -6 -5.86667 5 2.23607
ij 2 -5.73333 -5.6 138 11.7473
ij 3 -5.46667 -5.33333 559 23.6432
ij 4 -5.2 -5.06667 1393 37.3229
ij 5 -4.93333 -4.8 2449 49.4874
ij 6 -4.66667 -4.53333 3372 58.0689
ij 7 -4.4 -4.26667 3837 61.9435
ij 8 -4.13333 -4 3777 61.4573
ij 9 -3.86667 -3.73333 3503 59.1861
ij 10 -3.6 -3.46667 2699 51.9519
ij 11 -3.33333 -3.2 1775 42.1307
ij 12 -3.06667 -2.93333 775 27.8388
ij 13 -2.8 -2.66667 270 16.4317
ij 14 -2.53333 -2.4 17 4.12311
xx2 
-2.26667
yy0 30 30 14 -6 -2.26667


anal_evt_unfold1

anal_evt_unfold

rm anal_evt_unfold1

make anal_evt_unfold1

anal_evt_unfold1

make clean

make

./anal_evt_unfold01

./anal_evt_unfold02

./anal_evt_unfold03

./anal_evt_unfold04

./anal_evt_unfold05

./anal_evt_unfold06

./anal_evt_unfold07

./anal_evt_unfold08

./anal_evt_unfold09

./anal_evt_unfold10

./anal_evt_unfold11

./anal_evt_unfold12


// to check the iteration number and bin number for each variable 
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_3_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_6_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_9_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_12_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_15_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_18_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_21_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_24_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_27_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_30_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_31_Pythia6_1_c0_e1.txt
grep "iter bayes_Data_PFjet_" outfile_hhi3_it10_33_Pythia6_1_c0_e1.txt



*/
