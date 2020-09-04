#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <TArrayC.h>
#include <string>
#include "TH1.h"
#include "TH2D.h"
#include <TH1D.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include "TUnfoldDensity.h"
#include "TUnfoldIterativeEM.h"

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldIterativeEM.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>



#define NOREG
//#define TIKH_SURE
//#define TIKH_LSCAN
//#define ITER_SURE
//#define ITER
//#define REBIN
//#define UOEx
#define CLOUSER

using namespace std;

static const auto feps = numeric_limits<float>::epsilon();


void testUnfold2c()
{
  // switch on histogram errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //Input Data and MC histogram
  //TFile *inputData=new TFile("PY8_Flat_UL17_Newbin_22Aug20.root");
 // TFile *inputData=new TFile("PY8_Flat_UL17_MissRM_22Aug20.root");
  //TFile *inputData=new TFile("Test_MC_QCD.root");
  TFile *inputData=new TFile("PY8_UL17_binned_22Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Binnned_MissRM_22Aug20.root");
  //TFile *inputData=new TFile("MG_UL17_MissRM_23Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Flat_25Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Binned_25Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Flat_31Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Flat_supoff_31Aug20.root");
  //TFile *inputData=new TFile("PY8_UL17_Binned_17Aug20_1.root");
  //TFile *inputData=new TFile("PY8_UL17_Binned_20Aug20.root");
  //TFile *inputData=new TFile("PY_UL17_Flat_header_31Aug20.root");
  //TFile *inputData=new TFile("Data_UL2017_JetHT_24Aug20.root");
  
  //TFile *inputMC=new TFile("Test_MC_QCD.root");
  //TFile *inputMC=new TFile("PY8_Flat_UL17_Newbin_22Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Binnned_MissRM_22Aug20.root");
  TFile *inputMC=new TFile("PY8_UL17_binned_22Aug20.root");
  //TFile *inputMC=new TFile("MG_UL17_MissRM_23Aug20.root");
  //TFile *inputMC=new TFile("PY8_CP5_UL_Reco2Gen_30July20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Binned_25Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Flat_31Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Flat_supoff_31Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Binned_27Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Binned_20Aug20.root");
  //TFile *inputMC=new TFile("PY8_UL17_Binned_17Aug20_2.root");
  //TFile *inputMC=new TFile("PY_UL17_Flat_header_31Aug20.root");
  //TFile *inputMC=new TFile("PY8_Binned_UL17_14Aug20_Fake0bins.root");
  
  //TFile *inputMC1=new TFile("MG_UL17_binned_6July20.root");
  TFile *inputMC1=new TFile("MG_UL17_MissRM_23Aug20.root");
  //TFile *inputMC1=new TFile("PY8_UL17_ULbinned_1July20.root");
  TFile *inputMC2=new TFile("HW7_UL17_Flat_23Aug20.root");
  //TFile *inputMC2=new TFile("HW7_Flat_Binned_8July20.root");
  
  //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  TFile *outputFile=new TFile("testunfold2c_unfolded.root","recreate");
  int irbin = 2;
  
  ofstream file;
  file.open("Tau_Value.txt");
  file <<"L Curve    "  <<"                   "<< "Scan Sure  "  <<"                "<<" Scan Tau "<<endl;
  
  int const type = 2;         // Jet & Charage particles
  int const itype[type]={0,1};   //{0}--->Jet ; {1}---> Charged Particles
  const  char* itypeN[type]={"Jets","Charged Particles"};
  char histname[100], name[100];
  const int nHLTmx=8; //HT2 Range
  const int njetetamn=1;  //eta value used 2.4
  //double etarange[njetetamn] ={2.4}; //
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass
  static const int nhist=10; //We need 8 But define 10
  const int njetptmn = nHLTmx;
  double leadingPtThreshold[njetptmn+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger
  const char* vartitle[nvar]={"Anti-Y_{23,C} ", "Anti-Y_{23,E} ", "Anti-Y_{23,R} ",
			      "#tau_{_{#perp}} ", "#tau_{_{#perp} _{   ,E}} ", "#tau_{_{#perp} _{   ,R}} ",
			      "T_{ m,C} ", "T_{ m,E} ", "T_{ m,R} ",
			      "#rho_{Tot} ", "#rho_{Tot,E} ", "#rho_{Tot,R} ",
                              "#rho_{H,C} ", "#rho_{H,E} ", "#rho_{H,R} ",
			      "Y_{23} ", "Y_{23,E} ", "Y_{23,R} ",
			      "B_{ T,C} ", "B_{ T,E} ", "B_{ T,R} ",
			      "B_{ W,C} ", "B_{ W,E} ", "B_{ W,R} ",
			      "#rho^{T}_{Tot} ", "#rho^{T}_{Tot,E} ", "#rho^{T}_{Tot,R} ",
			      "#rho^{T}_{H,C} ", "#rho^{T}_{H,E} ", "#rho^{T}_{H,R} ",
			      "S_{_{#perp} _{   ,C}}", "C-parameter_{C}"};
  
  //----------------------for reco and gen bin number
  int rnbinsx[type][nusedvar][nHLTmx];
  int gnbinsx[type][nusedvar][nHLTmx];
  
/*  // define the number of bin excluded form lower end
  int arrayvar_first[type][nusedvar][nHLTmx] = {{
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0}},
						{{0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0}}};
  
  // define the number of bin excluded from higher end
  int arrayvar_last[type][nusedvar][nHLTmx] =  {{
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0}},
						{{0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0}}};
  //gen bin 
  // define the number of bin excluded form lower end
int arrayvar_firstG[type][nusedvar][nHLTmx] = {{
                                               {0,0,0,0,0,0,0,0},
                                               {0,0,0,0,0,0,0,0},
                                               {0,0,0,0,0,0,0,0},
                                               {0,0,0,0,0,0,0,0},
                                               {0,0,0,0,0,0,0,0}},
					       {{0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0}}};
 
// define the number of bin excluded from higher end
 int arrayvar_lastG[type][nusedvar][nHLTmx] =  {{
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0},
                                                {0,0,0,0,0,0,0,0}},
						{{0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0},
						 {0,0,0,0,0,0,0,0}}};
 
 */
 TH1D *MC_Reco[type][nusedvar][njetptmn];  //Reconstructed MC
 TH1D *MC_fake[type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
 TH1D *MC_Gen[type][nusedvar][njetptmn];   //Generator MC
 TH1D *MC_miss[type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
 TH1D *Data_Reco[type][nusedvar][njetptmn];    //Reconstructed Data
 TH1D *PsudoData_Gen[type][nusedvar][njetptmn];    //Gen Level Psudo Data
 TH2D *h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco
 
 TH1D *MC_Reco1[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *MC_fake1[type][nusedvar][njetptmn];  //Fake :  Reco but No Gen Copy
 TH1D *MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy 
 TH1D *MC_miss1[type][nusedvar][njetptmn];   ///Miss:  No Reco but in Gen Copy
 TH1D *Data_Reco1[type][nusedvar][njetptmn];    //Reconstructed Data copy
 TH2D *h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy
 
 
 TH1D *MC_Reco2[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *MC_fake2[type][nusedvar][njetptmn];  //Fake :  Reco but No Gen Copy
 TH1D *MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
 TH1D *MC_miss2[type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen Copy
 TH1D *Data_Reco2[type][nusedvar][njetptmn];    //Reconstructed Data copy
 TH2D *h2dGenDetMC2[type][nusedvar][njetptmn];  //// MC generator Vs Reco copy
 
 TH1D *MC_Reco_re[type][nusedvar][njetptmn];  //Reco MC Rebins If needed 
 TH1D *MC_Gen_re[type][nusedvar][njetptmn];   //Generator Rebins If Needed
 TH1D *Data_Reco_re[type][nusedvar][njetptmn];    //Reconst Data with Rebin
 TH2D *h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco Rebin
 
 //---------------------MG5
 TH1D *MG5_MC_Reco[type][nusedvar][njetptmn];  //Reconstructed MC
 TH1D *MG5_MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
 TH2D *MG5_h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco
 
 TH1D *MG5_MC_Reco1[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *MG5_MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy
 TH2D *MG5_h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy
 
 
 TH1D *MG5_MC_Reco2[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *MG5_MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
 TH2D *MG5_h2dGenDetMC2[type][nusedvar][njetptmn];
 
 TH1D *MG5_MC_Reco_re[type][nusedvar][njetptmn];  //Reconstructed MC
 TH1D *MG5_MC_Gen_re[type][nusedvar][njetptmn];   //Generator level MC
 TH2D *MG5_h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco
 //----------------------HW
 TH1D *HW7_MC_Reco[type][nusedvar][njetptmn];  //Reconstructed MC
 TH1D *HW7_MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
 TH2D *HW7_h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco
 
 TH1D *HW7_MC_Reco1[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *HW7_MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy
 TH2D *HW7_h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy
 
 
 TH1D *HW7_MC_Reco2[type][nusedvar][njetptmn];  //Reconstructed MC copy
 TH1D *HW7_MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
 TH2D *HW7_h2dGenDetMC2[type][nusedvar][njetptmn];
 
 TH1D *HW7_MC_Reco_re[type][nusedvar][njetptmn];  //Reconstructed MC
 TH1D *HW7_MC_Gen_re[type][nusedvar][njetptmn];   //Generator level MC
 TH2D *HW7_h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco
 
 
 
 
 TH1D *unfold_NoReg[type][nusedvar][njetptmn];
 TH1D *unfold_Tau[type][nusedvar][njetptmn];
 TH1D *unfold_L[type][nusedvar][njetptmn];
 
 TH1D* hist_eff[type][nusedvar][njetptmn];
 TH1D* hist_fake[type][nusedvar][njetptmn];
 TH1D* hist_purity[type][nusedvar][njetptmn];
 TH1D* hist_stbl[type][nusedvar][njetptmn];
 
 TH1D* hist_eff1[type][nusedvar][njetptmn];
 TH1D* hist_fake1[type][nusedvar][njetptmn];
 TH1D* hist_purity1[type][nusedvar][njetptmn];
 TH1D* hist_stbl1[type][nusedvar][njetptmn];
 
 TH2D* COV_Mat_NoReg[type][nusedvar][njetptmn];
 TH2D* COV_Mat_Tau[type][nusedvar][njetptmn];
 TH2D* COV_Mat_L[type][nusedvar][njetptmn];
 
 TH2D* corr_mat_NoReg[type][nusedvar][njetptmn];
 
 
 /*TMatrixD covmatrix_NoReg[type][nusedvar][njetptmn];
   TMatrixD covmatrix_Tua[type][nusedvar][njetptmn];
   TMatrixD covmatrix_L[type][nusedvar][njetptmn];
   TMatrixD corrmatrix_NoReg[type][nusedvar][njetptmn];
   TMatrixD corrmatrix_Tua[type][nusedvar][njetptmn];
   TMatrixD corrmatrix_L[type][nusedvar][njetptmn];
 */
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);
 gStyle->SetPadRightMargin(0.06);
 gStyle->SetPadBottomMargin(0.18);
 gStyle->SetPadLeftMargin(0.19);
 gStyle->SetPadTopMargin(0.07);
 int font=42;
 gStyle->SetLabelFont(font,"XYZ");
 gStyle->SetLegendFont(font);
 gStyle->SetStatFont(font);
 gStyle->SetTextFont(font);
 gStyle->SetTitleFont(font,"XYZ");
 gStyle->SetTitleOffset(1.5,"y");
 gStyle->SetTitleOffset(1.2,"x");
 gStyle->SetTitleSize(0.8,"P");
 gStyle->SetTitleSize(0.06,"xy");
 gStyle->SetLabelSize(0.05,"xy");
 gStyle->SetLabelOffset(0.012,"xy");
 
 
 TDirectoryFile *inputDir0=new TDirectoryFile("Data","Inputs Data");
 TDirectoryFile *inputRebin0=new TDirectoryFile("RebinData","Inputs Data Rebin");
 
 TDirectoryFile *inputDir1=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
 TDirectoryFile *inputDir2=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
 TDirectoryFile *inputDir3=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
 
 TDirectoryFile *inputRebin1=new TDirectoryFile("RebinPythia8","Pythia8 MC Rebin");
 TDirectoryFile *inputRebin2=new TDirectoryFile("RebinMG8","Madgraph, MC Rebin");
 TDirectoryFile *inputRebin3=new TDirectoryFile("RebinHW7","Herwig7 MC Rebin");
 
 TDirectoryFile *foldpy8=new TDirectoryFile("fold","folded with Probablility matrix");
 TDirectoryFile *Unfold=new TDirectoryFile("Unfold","Unfolded, Refold, correlation");
 
 int subtract_background1(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
 int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
 void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
 void Condition (TH2 * RM, TH1* miss);
 TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
 TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco,TH1D* MC_gen, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8],int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8] );
 // TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar );
 TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
 
 
 //Read Input Data MC and Response matrix
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
       //Reco Data
       inputDir0->cd();
       sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       Data_Reco1[ity][ivar][ipt] = (TH1D*) inputData->Get(histname);
       //Data_Reco1[ity][ivar][ipt]->Write();
       //if(ity ==0){Data_Reco1[ity][ivar][ipt]->Rebin(2);}
       for (int ibin =1 ; ibin <  Data_Reco1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
       if(Data_Reco1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " Data Reco Bin is Zero for bin number : ******** "<<  ibin  << endl; }}
#ifdef CLOUSER
       sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       PsudoData_Gen[ity][ivar][ipt] = (TH1D*) inputData->Get(histname);
#else
       PsudoData_Gen[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
#endif

       inputDir1->cd();
       sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
       //if(ity ==0){MC_Reco1[ity][ivar][ipt] ->Rebin(2);}
       //MC_Reco1[ity][ivar][ipt]->Write();
       cout << histname << endl;
       int recobins = MC_Reco1[ity][ivar][ipt]->GetNbinsX();
       rnbinsx[ity][ivar][ipt]=MC_Reco1[ity][ivar][ipt]->GetNbinsX();
       if(recobins !=Data_Reco1[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
       
       //MC Fake
       sprintf(histname, "analyzeBasicPat/fake_reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MC_fake1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
       //MC_fake1[ity][ivar][ipt]->Write();
       //if(ity ==0){MC_fake1[ity][ivar][ipt]->Rebin(2);}
       cout << histname << endl;
       
       //Gen MC
       sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
       //if(ity ==0){MC_Gen1[ity][ivar][ipt]->Rebin(2);}
       //MC_Gen1[ity][ivar][ipt]->Write();
       int genbins = MC_Gen1[ity][ivar][ipt]->GetNbinsX();
       gnbinsx[ity][ivar][ipt]=MC_Gen1[ity][ivar][ipt]->GetNbinsX();
       for(int ibin = 1; ibin < MC_Gen1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
       if (MC_Gen1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  ibin  << endl; }}
       
       //MC miss 
       sprintf(histname, "analyzeBasicPat/miss_gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MC_miss1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
       // if(ity ==0){MC_miss1[ity][ivar][ipt]->Rebin(2);}
       
       
       //Response Matrix
       sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
       h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC->Get(histname);   //Xgen(coarse) , Yreco(fine)
       //h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)
       //if(recobins !=h2dGenDetMC1[ity][ivar][ipt]->GetNbinsX()){cout << "Reco Bin missmatch in response Matrix " << endl;}
       //if(genbins !=h2dGenDetMC1[ity][ivar][ipt]->GetNbinsY()){cout << "Gen Bin missmatch in response Matrix " << endl;}
       //cout << "Recobins = " <<  recobins  << "   Gen Bins = " << genbins << endl;

      //-------------------Replace Fake by Reco - RMx Projection------------------
      //-----------------------------------------------------------------------
       TH1D* RMxx= h2dGenDetMC1[ity][ivar][ipt]->ProjectionX();
       for (int i = 1; i <= MC_fake1[ity][ivar][ipt]->GetNbinsX(); ++i) {
         double content = MC_Reco1[ity][ivar][ipt]->GetBinContent(i);
         double factor = RMxx->GetBinContent(i);
         content -= factor;
         MC_fake1[ity][ivar][ipt]->SetBinContent(i, content);
       }


       Data_Reco2[ity][ivar][ipt] = (TH1D*)Data_Reco1[ity][ivar][ipt]->Clone();
       MC_Reco2[ity][ivar][ipt] = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
       MC_Gen2[ity][ivar][ipt] = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
       h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();
       
       MC_miss2[ity][ivar][ipt] = (TH1D*)MC_miss1[ity][ivar][ipt]->Clone();
       MC_fake2[ity][ivar][ipt] = (TH1D*)MC_fake1[ity][ivar][ipt]->Clone();
      




//A kind of way to remove underflow and overlow bins
#ifdef UOEx
       TH1D *NewData; TH1D *NewReco; TH1D *NewGen; TH1D *Newmiss; TH1D *Newfake;
       TH2D *NewRecoGen;
       NewData = (TH1D*)Data_Reco1[ity][ivar][ipt]->Clone();
       NewReco = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
       NewGen  = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
       NewRecoGen = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();
       
       Newmiss  = (TH1D*)MC_miss1[ity][ivar][ipt]->Clone();
       Newfake = (TH1D*)MC_fake1[ity][ivar][ipt]->Clone();
       
       int NbinxD = NewData->GetNbinsX();
       int NbinxR = NewReco->GetNbinsX();
       int NbinxG = NewGen->GetNbinsX();
       int NbinMx = NewRecoGen->GetNbinsX();
       int NbinMy = NewRecoGen->GetNbinsY();
       NewData->Reset(); NewReco->Reset(); NewGen->Reset(); NewRecoGen->Reset(); Newmiss->Reset(); Newfake->Reset();
       
       for(int ix=1; ix < NbinxD+1 ; ix++){
	 NewData->SetBinContent(ix,Data_Reco1[ity][ivar][ipt]->GetBinContent(ix));
	 NewData->SetBinError(ix, sqrt(Data_Reco1[ity][ivar][ipt]->GetBinError(ix)* Data_Reco1[ity][ivar][ipt]->GetBinError(ix)));
       }
       
       for(int ix=0; ix < NbinxR+2 ; ix++){
	 NewReco->SetBinContent(ix,MC_Reco1[ity][ivar][ipt]->GetBinContent(ix));
	 NewReco->SetBinError(ix, sqrt(MC_Reco1[ity][ivar][ipt]->GetBinError(ix)* MC_Reco1[ity][ivar][ipt]->GetBinError(ix)));
       }
       for(int ix=0; ix < NbinxG+2 ; ix++){
	 NewGen->SetBinContent(ix,MC_Gen1[ity][ivar][ipt]->GetBinContent(ix));
	 NewGen->SetBinError(ix, sqrt(MC_Gen1[ity][ivar][ipt]->GetBinError(ix)* MC_Gen1[ity][ivar][ipt]->GetBinError(ix)));
       }       
       for(int ix=0; ix < NbinMx+2 ; ix++){
	 for(int iy=0; iy < NbinMy+2 ; iy++){
	   NewRecoGen->SetBinContent(ix, iy, h2dGenDetMC1[ity][ivar][ipt]->GetBinContent(ix,iy));
	   NewRecoGen->SetBinError(ix, iy, sqrt(h2dGenDetMC1[ity][ivar][ipt]->GetBinError(ix,iy)* h2dGenDetMC1[ity][ivar][ipt]->GetBinError(ix,iy)));
	 }
       }
       //NewRecoGen->SetBinContent(0, 0,0.0); NewRecoGen->SetBinError(0,0,0.0);
       for(int ix=0; ix < NbinxR+2 ; ix++){
	 Newfake->SetBinContent(ix,MC_fake1[ity][ivar][ipt]->GetBinContent(ix));
	 Newfake->SetBinError(ix, sqrt(MC_fake1[ity][ivar][ipt]->GetBinError(ix)* MC_fake1[ity][ivar][ipt]->GetBinError(ix)));
       }
       for(int ix=0; ix < NbinxG+2 ; ix++){
	 Newmiss->SetBinContent(ix,MC_miss1[ity][ivar][ipt]->GetBinContent(ix));
	 Newmiss->SetBinError(ix, sqrt(MC_miss1[ity][ivar][ipt]->GetBinError(ix)* MC_miss1[ity][ivar][ipt]->GetBinError(ix)));
       }
       
       //Fill the MISS at reco underflow and gen underflow
       for(int iy=0; iy < NbinMy+2 ; iy++){
       }
       
       //if(ity==0){
       //Fill the fake in gen underflow
       for(int ix=0; ix < NbinMx+2 ; ix++){
       //NewRecoGen->SetBinError(ix, 0, MC_fake1[ity][ivar][ipt]->GetBinError(ix)* h2dGenDetMC1[ity][ivar][ipt]->GetBinError(ix,0));
  //     NewRecoGen->SetBinError(0, ix, MC_fake1[ity][ivar][ipt]->GetBinError(ix)* h2dGenDetMC1[ity][ivar][ipt]->GetBinError(0,ix));
       //NewRecoGen->SetBinError(ix, 0, 0.0);
       }
       //}
       
       Data_Reco1[ity][ivar][ipt] = (TH1D*)NewData->Clone();
       MC_Reco1[ity][ivar][ipt] = (TH1D*)NewReco->Clone();
       MC_Gen1[ity][ivar][ipt]  = (TH1D*)NewGen->Clone();
       h2dGenDetMC1[ity][ivar][ipt] = (TH2D*)NewRecoGen->Clone();
       
       MC_fake1[ity][ivar][ipt]  = (TH1D*)Newfake->Clone();
       MC_miss1[ity][ivar][ipt]  = (TH1D*)Newmiss->Clone();
       
#endif

//----------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20
//Saved in Pythia8 directory
       TH1* RMx = h2dGenDetMC1[ity][ivar][ipt]->ProjectionX();
       TH1* RMy = h2dGenDetMC1[ity][ivar][ipt]->ProjectionY();

       sprintf(name,"ProjectX_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       RMx->SetNameTitle(name,name);
       RMx->Write();

       sprintf(name,"Recominusfake_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       TH1* RecoFakeCorrect = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
       RecoFakeCorrect->SetNameTitle(name,name);
       

       sprintf(name,"ProjectY_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       RMy->SetNameTitle(name,name);
       RMy->Write();
       
       sprintf(name,"Genminusmiss_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       TH1* GenMissCorrect = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
       GenMissCorrect->SetNameTitle(name,name);

        for (int i = 1; i <= RecoFakeCorrect->GetNbinsX(); ++i) {
         double content = RecoFakeCorrect->GetBinContent(i);
         double factor = MC_fake1[ity][ivar][ipt]->GetBinContent(i);
         content -= factor;
	 cout << " fake Factor " << factor <<endl;
         RecoFakeCorrect->SetBinContent(i, content);
       }
         //RecoFakeCorrect->Divide(RMx);
        // RecoFakeCorrect->SetMinimum(0.5); MC_fake1[ity][ivar][ipt]->SetMaximum(1.6);

        RecoFakeCorrect->Write();
        for (int i = 1; i <= GenMissCorrect->GetNbinsX(); ++i) {
         double content = GenMissCorrect->GetBinContent(i);
         double factor = MC_miss1[ity][ivar][ipt]->GetBinContent(i);
         content -= factor;
	 cout << " Miss Factor " << factor << endl;
         GenMissCorrect->SetBinContent(i, content);
       }
         //GenMissCorrect->Divide(RMy);
        // GenMissCorrect->SetMinimum(0.85); MC_fake1[ity][ivar][ipt]->SetMaximum(1.15);
	 GenMissCorrect->Write();

//-----------------------------------------------------------------------------

//--------------------------------Fake and Miss rate : Thanks to Patrick  14Aug20

       MC_fake1[ity][ivar][ipt]->Divide(MC_fake1[ity][ivar][ipt], MC_Reco1[ity][ivar][ipt], 1, 1, "b");
       //MC_fake1[ity][ivar][ipt]->Divide(MC_Reco1[ity][ivar][ipt]);
       MC_miss1[ity][ivar][ipt]->Divide(MC_miss1[ity][ivar][ipt], MC_Gen1[ity][ivar][ipt], 1, 1, "b");
       //MC_miss1[ity][ivar][ipt]->Divide(MC_Gen1[ity][ivar][ipt]);
       

       MC_fake1[ity][ivar][ipt]->SetMinimum(-0.05); MC_fake1[ity][ivar][ipt]->SetMaximum(1.01);
       MC_miss1[ity][ivar][ipt]->SetMinimum(-0.05); MC_miss1[ity][ivar][ipt]->SetMaximum(1.01);

       Data_Reco[ity][ivar][ipt] =(TH1D*)Data_Reco1[ity][ivar][ipt]->Clone();
       MC_Reco[ity][ivar][ipt] = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
       MC_Gen[ity][ivar][ipt] = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
       MC_fake[ity][ivar][ipt] = (TH1D*)MC_fake1[ity][ivar][ipt]->Clone();
       MC_miss[ity][ivar][ipt] = (TH1D*)MC_miss1[ity][ivar][ipt]->Clone();
       h2dGenDetMC[ity][ivar][ipt] = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();
       
       inputDir0->cd();
       Data_Reco[ity][ivar][ipt]->Write();
       PsudoData_Gen[ity][ivar][ipt]->Write();
       inputDir1->cd();
       MC_Reco[ity][ivar][ipt]->Write();
       MC_Gen[ity][ivar][ipt]->Write();
       MC_fake[ity][ivar][ipt]->Write();
       MC_miss[ity][ivar][ipt]->Write();
       h2dGenDetMC[ity][ivar][ipt]->Write();

       
#ifdef REBIN
       const char *dataold; char datanew[100];
       dataold = Data_Reco1[ity][ivar][ipt]->GetName();
       sprintf(datanew, "Data_%s", dataold);  Data_Reco1[ity][ivar][ipt]->SetName(datanew);
       
       Data_Reco_re[ity][ivar][ipt] = rebin1d_hist(Data_Reco1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
       
       MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(MC_Reco1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
       //MC_Reco_re[ity][ivar][ipt]= (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
       //MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
       //MC_Gen_re[ity][ivar][ipt] = rebin1d_hist_gen(MC_Gen1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
       MC_Gen_re[ity][ivar][ipt] = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
       cout << " Rebin MC Gen bin : " <<  MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;
       
       
       h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(h2dGenDetMC1[ity][ivar][ipt], MC_Reco1[ity][ivar][ipt], MC_Gen1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last,arrayvar_firstG,arrayvar_lastG);
       cout << "ReBin REco Bin :"<< h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;
       
       //Define data MC to be used for unfolding
         Data_Reco[ity][ivar][ipt] = (TH1D*)Data_Reco_re[ity][ivar][ipt]->Clone();
         MC_Reco[ity][ivar][ipt] = (TH1D*)MC_Reco_re[ity][ivar][ipt]->Clone();
         MC_Gen[ity][ivar][ipt] = (TH1D*)MC_Gen_re[ity][ivar][ipt]->Clone();
         h2dGenDetMC[ity][ivar][ipt] = (TH2D*)h2dGenDetMC_re[ity][ivar][ipt]->Clone();
         cout << "Rebin for Pythia8 Done"
#else
	 
#endif
	 inputRebin0->cd();
	 Data_Reco[ity][ivar][ipt]->Write();
	 inputRebin1->cd();
	 MC_Reco[ity][ivar][ipt]->Write();
	 MC_fake[ity][ivar][ipt]->Write();
         MC_Gen[ity][ivar][ipt]->Write();
         MC_miss[ity][ivar][ipt]->Write();
         h2dGenDetMC[ity][ivar][ipt]->Write();
	 
	 
     }
   }
 }
 cout << "Histogram Read for Pythia8 Done " <<endl;
//------------------Read Madgraph
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
       inputDir2->cd();
       sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MG5_MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC1->Get(histname);
       //cout << histname << endl;
       MG5_MC_Reco1[ity][ivar][ipt]->Write();
       
       //Gen MC
       sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
       MG5_MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC1->Get(histname);
       //MC_Gen1[ity][ivar][ipt]->Rebin(2);
       MG5_MC_Gen1[ity][ivar][ipt]->Write();

       //Response Matrix
       sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
       MG5_h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC1->Get(histname);   //Xgen(coarse) , Yreco(fine)
       //MG5_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
       MG5_h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)
	 
       MG5_MC_Reco2[ity][ivar][ipt] = (TH1D*)MG5_MC_Reco1[ity][ivar][ipt]->Clone();
       MG5_MC_Gen2[ity][ivar][ipt] = (TH1D*)MG5_MC_Gen1[ity][ivar][ipt]->Clone();
       MG5_h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)MG5_h2dGenDetMC1[ity][ivar][ipt]->Clone();
	 
#ifdef REBIN
       MG5_MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(MG5_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
       //MG5_MC_Reco_re[ity][ivar][ipt]->Rebin(4);
	 
       //MG5_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(MG5_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
       MG5_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist_gen(MG5_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
       //MG5_MC_Gen_re[ity][ivar][ipt]->Rebin(2);
       //cout << " Rebin MC Gen bin : " <<  MG5_MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;
	 
       MG5_h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(MG5_h2dGenDetMC1[ity][ivar][ipt], MG5_MC_Reco1[ity][ivar][ipt],MG5_MC_Gen1[ity][ivar][ipt] ,ity, ipt, ivar,arrayvar_first,arrayvar_last,arrayvar_firstG,arrayvar_lastG);
       //MG5_h2dGenDetMC_re[ity][ivar][ipt] = h2dGenDetMC1[ity][ivar][ipt];
       //MG5_h2dGenDetMC_re[ity][ivar][ipt]->RebinY(2);
       //MG5_h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
       cout << "ReBin REco Bin :"<< MG5_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << MG5_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;
       //Define data MC to be used for unfolding
       MG5_MC_Reco[ity][ivar][ipt] = (TH1D*)MG5_MC_Reco_re[ity][ivar][ipt]->Clone();
       MG5_MC_Gen[ity][ivar][ipt] = (TH1D*)MG5_MC_Gen_re[ity][ivar][ipt]->Clone();
       MG5_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)MG5_h2dGenDetMC_re[ity][ivar][ipt]->Clone();
       cout << "Histogram MG Read and Rebin OK" <<endl;
#else
       MG5_MC_Reco[ity][ivar][ipt] = (TH1D*)MG5_MC_Reco1[ity][ivar][ipt]->Clone();
       MG5_MC_Gen[ity][ivar][ipt] = (TH1D*)MG5_MC_Gen1[ity][ivar][ipt]->Clone();
       MG5_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)MG5_h2dGenDetMC1[ity][ivar][ipt]->Clone();
#endif
       inputRebin2->cd();
       MG5_MC_Reco[ity][ivar][ipt]->Write();
       MG5_MC_Gen[ity][ivar][ipt]->Write();
       MG5_h2dGenDetMC[ity][ivar][ipt]->Write();
     }
   }
 }
 cout << "Histogram MG Read" <<endl;
 //-----------HW7
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
       inputDir3->cd();
       sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         HW7_MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC2->Get(histname);
	 //cout << histname << endl;
         HW7_MC_Reco1[ity][ivar][ipt]->Write();
         //Gen MC
         sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         HW7_MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC2->Get(histname);
	 //       HW7_MC_Gen1[ity][ivar][ipt]->Rebin(2);
         HW7_MC_Gen1[ity][ivar][ipt]->Write();
	 
         //Response Matrix
         sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
         HW7_h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC2->Get(histname);   //Xgen(coarse) , Yreco(fine)
	 //HW7_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
	 HW7_h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)

         HW7_MC_Reco2[ity][ivar][ipt] = (TH1D*)HW7_MC_Reco1[ity][ivar][ipt]->Clone();
         HW7_MC_Gen2[ity][ivar][ipt] = (TH1D*)HW7_MC_Gen1[ity][ivar][ipt]->Clone();
         HW7_h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)HW7_h2dGenDetMC1[ity][ivar][ipt]->Clone();
	 
	 
#ifdef REBIN
        HW7_MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(HW7_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
	// HW7_MC_Reco_re[ity][ivar][ipt]->Rebin(4);
	
	//  HW7_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(HW7_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
        HW7_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist_gen(HW7_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last);
	//HW7_MC_Gen_re[ity][ivar][ipt]->Rebin(2);
        //cout << " Rebin MC Gen bin : " <<  HW7_MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;
	
        HW7_h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(HW7_h2dGenDetMC1[ity][ivar][ipt], HW7_MC_Reco1[ity][ivar][ipt], HW7_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar,arrayvar_first,arrayvar_last,arrayvar_firstG,arrayvar_lastG);
	//HW7_h2dGenDetMC_re[ity][ivar][ipt] = h2dGenDetMC1[ity][ivar][ipt];
        //HW7_h2dGenDetMC_re[ity][ivar][ipt]->RebinY(2);
        //HW7_h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
	//cout << "ReBin REco Bin :"<< HW7_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << HW7_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;
	
	
        //Define data MC to be used for unfolding
	HW7_MC_Reco[ity][ivar][ipt] = (TH1D*)HW7_MC_Reco_re[ity][ivar][ipt]->Clone();
	HW7_MC_Gen[ity][ivar][ipt] = (TH1D*)HW7_MC_Gen_re[ity][ivar][ipt]->Clone();
	HW7_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)HW7_h2dGenDetMC_re[ity][ivar][ipt]->Clone();
	
	
        cout << "Histogram HW Read and Rebin OK" <<endl;
#else
        HW7_MC_Reco[ity][ivar][ipt] = (TH1D*)HW7_MC_Reco1[ity][ivar][ipt]->Clone();
        HW7_MC_Gen[ity][ivar][ipt] = (TH1D*)HW7_MC_Gen1[ity][ivar][ipt]->Clone();
        HW7_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)HW7_h2dGenDetMC1[ity][ivar][ipt]->Clone();
	 
#endif
        inputRebin3->cd();
        HW7_MC_Reco[ity][ivar][ipt]->Write();
        HW7_MC_Gen[ity][ivar][ipt]->Write();
	HW7_h2dGenDetMC[ity][ivar][ipt]->Write();
     }
   }
 }
 
 cout << "Histogram HW7 Read oK " <<endl;

 //------------------Fold check : Patrick 1 Sep 20
 //Get Probability Matrix  
 //(Multiply the gen level by the probability matrix. Of course, don't forget to account for miss and fake entries (if applicable).
foldpy8->cd();
  for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
      TH2D* RM  = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();
      TH1D* Reco  = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
      TH1D* Gen  = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
      TH1D* fake= (TH1D*)MC_fake2[ity][ivar][ipt]->Clone();
      TH1D* miss = (TH1D*)MC_miss2[ity][ivar][ipt]->Clone();
      
      RM->RebinY(2);Gen->Rebin(2);miss->Rebin(2);

      TH1D* Folded = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone(); Folded->Reset();
       
      Fold(RM, Reco, Gen, miss, fake, Folded);
      sprintf(name,"Fold_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
      Folded->SetNameTitle(name,name);

      Folded->Write();

     }
   }
  }
//Condition of Probability Matrix

for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
      TH2D* RM  = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();
      TH1D* miss = (TH1D*)MC_miss2[ity][ivar][ipt]->Clone();
      RM->RebinY(2); miss->Rebin(2);
      cout <<setw(2) << ity <<setw(5) <<ivar << setw(5) << ipt <<'\n';
      Condition(RM, miss);

     }
   }
}


 //--------------------------------------------------------------
 Unfold->cd();
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
       //Print the reco bins
       cout <<"Type = " << ity <<" var " << ivar << " HT2 =" << ipt << endl;
       //cout <<" Reco Bin : {"; 
       double rxbins[MC_Reco[ity][ivar][ipt]->GetNbinsX()+1]={0};
       for (int ix=0; ix<MC_Reco[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	 rxbins[ix] = MC_Reco[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
	// cout <<  rxbins[ix] <<", " ; 
       }
      // cout << endl; 
      // cout <<" Gen Bin : {"; 
       //Get Gen bins
       double gxbins[MC_Gen[ity][ivar][ipt]->GetNbinsX()+1]={0};
       for (int ix=0; ix<MC_Gen[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	 gxbins[ix] = MC_Gen[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
	// cout <<  gxbins[ix] <<", "; 
       }
      // cout << endl; 
      // cout <<" 2d Reco  Bin : {"; 
       //Get Gen bins
       for (int ix=0; ix<h2dGenDetMC[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
//	 cout <<  h2dGenDetMC[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1) <<" , ";
       }
     //  cout << endl; 
      // cout <<" 2D Gen  Bin : {"; 
       for (int ix=0; ix<h2dGenDetMC[ity][ivar][ipt]->GetNbinsY()+1; ix++) {
//	 cout <<  h2dGenDetMC[ity][ivar][ipt]->GetYaxis()->GetBinLowEdge(ix+1) <<" , ";
       }
       
       cout << endl; 
       sprintf(name,"Effi_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_eff[ity][ivar][ipt] = new TH1D(name,name,MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
       hist_eff[ity][ivar][ipt]->Sumw2();
       sprintf(name,"Fake_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_fake[ity][ivar][ipt] = new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_fake[ity][ivar][ipt]->Sumw2();
       sprintf(name,"Purity_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_purity[ity][ivar][ipt] =  new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_purity[ity][ivar][ipt]->Sumw2();
       sprintf(name,"stability_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_stbl[ity][ivar][ipt] =  new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_stbl[ity][ivar][ipt]->Sumw2();
       
       //-----------------------------------------------for cross check
       sprintf(name,"Effi1_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_eff1[ity][ivar][ipt] = new TH1D(name,name,MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
       hist_eff1[ity][ivar][ipt]->Sumw2();
       sprintf(name,"Fake1_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_fake1[ity][ivar][ipt] = new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_fake1[ity][ivar][ipt]->Sumw2();
       sprintf(name,"Purity1_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_purity1[ity][ivar][ipt] =  new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_purity1[ity][ivar][ipt]->Sumw2();
       sprintf(name,"stability1_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       hist_stbl1[ity][ivar][ipt] =  new TH1D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(),rxbins);
       hist_stbl1[ity][ivar][ipt]->Sumw2();
       //-----------------------------------------------for cross check
       
       sprintf(name,"No_Reg_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       unfold_NoReg[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
       unfold_NoReg[ity][ivar][ipt]->Sumw2();
       //unfold_NoReg[ity][ivar][ipt]->Rebin(2);
       
       sprintf(name,"TauScan_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       unfold_Tau[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
       unfold_Tau[ity][ivar][ipt]->Sumw2();
       //unfold_Tau[ity][ivar][ipt]->Rebin(2);
       
       sprintf(name,"Lcurve_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       unfold_L[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
       unfold_L[ity][ivar][ipt]->Sumw2();
       //unfold_L[ity][ivar][ipt]->Rebin(2);
       
       
       sprintf(name,"Cov_Matrix_NoReg_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       COV_Mat_NoReg[ity][ivar][ipt] = new TH2D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(), rxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
       COV_Mat_NoReg[ity][ivar][ipt]->Sumw2();
       //COV_Mat_NoReg[ity][ivar][ipt]->RebinY(2);
       
       sprintf(name,"Corr_Matrix_NoReg_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       corr_mat_NoReg[ity][ivar][ipt] = new TH2D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
       corr_mat_NoReg[ity][ivar][ipt]->Sumw2();
       //corr_mat_NoReg[ity][ivar][ipt]->RebinY(2);
       //corr_mat_NoReg[ity][ivar][ipt]->RebinX(2);
     }
   }
 }
 
 //cross check efficincy and purity calculation (Tunfold example code)
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       //----------------------------------------
       for(int binGen=0;binGen<= h2dGenDetMC[ity][ivar][ipt]->GetNbinsY()+1;binGen++) {
         double sum0=0.;
         double sum1=0.;
         for(int binRec=0;binRec<= h2dGenDetMC[ity][ivar][ipt]->GetNbinsX()+1;
             binRec++) {
	   //double c=  h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	   double c=  h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binRec,binGen);
	   sum0+=c;
	   if((binRec>0)&&(binRec<=h2dGenDetMC[ity][ivar][ipt]->GetNbinsX())) {
	     sum1+=c;
	   }
         }
         if(sum0>0.0) {
	   hist_eff1[ity][ivar][ipt]->SetBinContent(binGen,sum1/sum0);
         }
       }
       //---------------------------
       hist_eff1[ity][ivar][ipt]->SetMinimum(0.9); hist_eff1[ity][ivar][ipt]->SetMaximum(1.1);
       hist_eff1[ity][ivar][ipt]->Write();
       

       for(int binRec=0; binRec<= hist_purity1[ity][ivar][ipt]->GetNbinsX()+1; binRec++) {
         double sum=0.;
         for(int binGen=0; binGen<=hist_purity1[ity][ivar][ipt]->GetNbinsX()+1; binGen++) {
	   //sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	   sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binRec,binGen);
         }
         double p=0.;
         if(sum>0.0) {
	   p = h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binRec,binRec)/sum;
         }
         hist_purity1[ity][ivar][ipt]->SetBinContent(binRec,p);
       }

       hist_purity1[ity][ivar][ipt]->SetMinimum(-0.05); hist_purity1[ity][ivar][ipt]->SetMaximum(1.01);
       hist_purity1[ity][ivar][ipt]->Write();
     }
   }
  }
 
 //-------------------------- //cross check efficincy and purity calculation
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){

       double fakerate[100] = {0.};
       double eff[100] ;
       for(int jj=0; jj<100; jj++){
	 eff[jj] = 1.;
       }
       double purity[100] = {0.};
       double stbl[100] = {0.};
       
       //    lastbin[ij][jk] = LastBin_Counter(hist_data[ij][jk]);
       //    For Purity and Stability 
       subtract_background(h2dGenDetMC[ity][ivar][ipt],MC_Reco[ity][ivar][ipt],MC_Gen[ity][ivar][ipt],Data_Reco[ity][ivar][ipt],fakerate,eff,purity,stbl); //No correction to Reco or Gen
       
       for(int bn=0; bn<(MC_Gen[ity][ivar][ipt]->GetNbinsX()); bn++){
         //hist_eff[ity][ivar][ipt]->SetBinContent(bn+1,eff[bn+1]);
        // hist_fake[ity][ivar][ipt]->SetBinContent(bn+1,fakerate[bn+1]) ;
         hist_purity[ity][ivar][ipt]->SetBinContent(bn+1,purity[bn+1]);
         hist_stbl[ity][ivar][ipt]->SetBinContent(bn+1,stbl[bn+1]);
       }

       hist_fake[ity][ivar][ipt] = (TH1D*)MC_fake[ity][ivar][ipt]->Clone();
       hist_eff[ity][ivar][ipt]  = (TH1D*)MC_miss[ity][ivar][ipt]->Clone();
       
       hist_eff[ity][ivar][ipt]->SetMinimum(-0.01); hist_eff[ity][ivar][ipt]->SetMaximum(1.01);
       hist_fake[ity][ivar][ipt]->SetMinimum(-0.1); hist_fake[ity][ivar][ipt]->SetMaximum(1.01);
       hist_purity[ity][ivar][ipt]->SetMinimum(-0.01); hist_purity[ity][ivar][ipt]->SetMaximum(1.01);
       hist_stbl[ity][ivar][ipt]->SetMinimum(-0.01); hist_stbl[ity][ivar][ipt]->SetMaximum(1.01);
       
       hist_eff[ity][ivar][ipt]->Write();
       hist_fake[ity][ivar][ipt]->Write();
       hist_purity[ity][ivar][ipt]->Write();
       hist_stbl[ity][ivar][ipt]->Write();
     }
   }
 }
 
double bias[5]={0.0,0.0,0.0,0.0,0.0};
//double bias[5]={1.0,1.0,1.0,1.0,1.0};
 //int bias[5]={2.0,2.0,2.0,2.0,2.0};

 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){

       //     if (ivar==2) continue; 
       cout <<"type "<< ity << " : Variables " << ivar << " HT2 Bin : " << ipt << endl;
       file <<"["<< ity << "," << ivar << "," << ipt <<"] --->" << endl;
     
        
       //Get reco bins
       double rxbins[MC_Reco[ity][ivar][ipt]->GetNbinsX()+1]={0};
       for (int ix=0; ix<MC_Reco[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	 rxbins[ix] = MC_Reco[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1); }
       
       //Get Gen bins
       double gxbins[MC_Gen[ity][ivar][ipt]->GetNbinsX()+1]={0};
       for (int ix=0; ix<MC_Gen[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	 gxbins[ix] = MC_Gen[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1); }
       
       //Rebin for match the condition of reco vs gen bin 
          h2dGenDetMC[ity][ivar][ipt]->RebinY(irbin);
          MC_Gen[ity][ivar][ipt]->Rebin(irbin);
       //h2dGenDetMC[ity][ivar][ipt]->Rebin(1,2);
         



       TH2* hist_migrationCoarseFine_MC = (TH2D*)h2dGenDetMC[ity][ivar][ipt]->Clone();
       TH1* input = (TH1D*)Data_Reco[ity][ivar][ipt]->Clone();
       TH1* mcgen = (TH1D*)MC_Gen[ity][ivar][ipt]->Clone();
       
       
       //correction for Fake as background subtraction : Patrick
       TH1* mcbackground = (TH1D*)MC_fake[ity][ivar][ipt]->Clone();
       mcbackground->Multiply(input);

       double biasScale =bias[ivar] ;
       const char *REGULARISATION_DISTRIBUTION=0;
       const char *REGULARISATION_AXISSTEERING="*[UOB]";
       char Rhoname[100], Rhotitle[100], Lcursure[100], Lcurtitle[100], TauLsure[100], TuaLtitle[100], probMat[100], probmat_title[100],Rhoname2d[100],
	 Rhotitle2d[100],Ematrix[100],Ematrixtitle[100],foldback[100],foldback_title[100];
       
       // preserve the area
       //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
       TUnfoldDensity::EConstraint constraintMode= TUnfoldDensity::kEConstraintArea;
       //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
       //TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
       //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
       TUnfoldDensity::ERegMode regMode = TUnfoldDensity::kRegModeCurvature;
       //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
       TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
       
       //https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
       sprintf(name,"Data_covariance_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       TH2D* covarianceM = new  TH2D(name,name, input->GetNbinsX(), rxbins, input->GetNbinsX(), rxbins);
       covarianceM->Sumw2();
       for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
          double err = input->GetBinError(ix);
	  covarianceM->SetBinContent(ix,ix,err*err);
       }
       
       covarianceM->Write();
       
       double taumx = 0.; double taumi = 0.;
       double tau = 1e-4;
       
       //TUnfoldDensity class
       //No  regularisation --------------------------------
       TUnfoldDensity tunfoldNoRegularisation(hist_migrationCoarseFine_MC,
					      TUnfold::kHistMapOutputVert,
					      TUnfoldDensity::kRegModeNone, 
					      TUnfoldDensity::kEConstraintNone,
					      TUnfoldDensity::kDensityModeNone);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       
       // TUnfoldDensity::kDensityModeBinWidthAndUser);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       
       tunfoldNoRegularisation.SubtractBackground(mcbackground, "Background", 1.0, 0.0); // hist,name,scale, scale error 
       int status = tunfoldNoRegularisation.SetInput(input,biasScale,0,covarianceM);
       
       int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
       cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;
       //tunfoldNoRegularisation.SubtractBackground(mcbackground,"Background", 1.0,0.04); // hist,name,scale, scale error 
       
       //if(>=10000) { std::cout<<"Unfolding result may be wrong\n";  }
       
       
       //the initial bias vector is determined from the response matrix
       //but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
       //tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result
       
       //if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeBins(7,1,4,regMode);}  //Test for Regularization in few  unmatched bins
       //if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeCurvature(10,9,11,1.0,0.6);}
 
       //Choose a value of tau to unfold with tau,0 means no regularization
       tunfoldNoRegularisation.DoUnfold(0.0);//,input,biasScale);
       
       /* //Binmaps : Thanks to Suman(Tifr)
       int binnos = rnbinsx[ity][ivar][ipt];
       Int_t *binMap=new Int_t[binnos+2];
       for(Int_t i=1;i<=binnos;i++) binMap[i]=i;
       binMap[0]=-1;   binMap[binnos+1]=-1;
       */

       char unfoldhist[100], title[100], NoReg_RhoIJ[100], NoReg_RhoIJ_tit[100], NoReg_Ematrix[100], NoReg_Ematrix_tit[100];// NoReg_prob[100], NoReg_probtitle[100];
       sprintf(unfoldhist, "Tunfold_Noreg_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(title, "Tunfolded Noreg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(NoReg_RhoIJ, "Tunfold_Noreg_corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(NoReg_RhoIJ_tit, "2D correlation coefficients No Regularisation %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(NoReg_Ematrix, "Tunfold_Noreg_Emat_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(NoReg_Ematrix_tit, "EMatrix No Regularisation %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(probMat, "Tunfold_Noreg_probM_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       sprintf(probmat_title, "Probability matrix Noreg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(foldback, "Tunfold_NoReg_Refold_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(foldback_title, "TunfoldFolded back  NoReg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);


       // unfolding result, signal only
       TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
       //TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,"","*[UO]");//,"signal");
       //TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput("hist_PTunfolded_noRegularisation", "P_{T,unfolded} [GeV]","signal");

       TH2 *hist_RhoIJ_noRegularisation = tunfoldNoRegularisation.GetRhoIJtotal(NoReg_RhoIJ, NoReg_RhoIJ_tit);//,"signal");
       TH2 *hist_Rho2D_noRegularisation = tunfoldNoRegularisation.GetEmatrixTotal(NoReg_Ematrix, NoReg_Ematrix_tit);//,"signal");
       //tunfoldNoRegularisation.GetEmatrix(COV_Mat_NoReg[ity][ivar][ipt],binMap);//,"signal");
       TH1 *hist_foldedback_NoReg = tunfoldNoRegularisation.GetFoldedOutput(foldback, foldback_title);//,"signal");
       TH2 *hist_prob_Noreg = tunfoldNoRegularisation.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");

       // correction for Miss entries  : Partick 
       for (int i = 1; i <= hist_PTunfolded_noRegularisation->GetNbinsX(); ++i) {
	 double content = hist_PTunfolded_noRegularisation->GetBinContent(i);
	 double factor = 1;
	 factor += MC_miss1[ity][ivar][ipt]->GetBinContent(i);
         content *= factor;
	 hist_PTunfolded_noRegularisation->SetBinContent(i, content);
       }
      
       /*  //Quick check for closure Ratio Plots
       hist_PTunfolded_noRegularisation->Scale(1/(hist_PTunfolded_noRegularisation->Integral()));
       mcgen->Scale(1/mcgen->Integral());
       hist_PTunfolded_noRegularisation->Divide(mcgen);
       hist_PTunfolded_noRegularisation->SetMinimum(0.85); hist_PTunfolded_noRegularisation->SetMaximum(1.15);
       */

       hist_PTunfolded_noRegularisation->Write();
       hist_Rho2D_noRegularisation->Write();
       hist_RhoIJ_noRegularisation->Write();
       hist_foldedback_NoReg->Write();
       hist_prob_Noreg->Write();
       
      cout << "Without Regularization finish " << endl;
       //-----------------------------------------------------------------------// L curve Scan
       
       TUnfoldDensity tunfoldTikhonovLCurve
	 (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert,regMode, constraintMode,densityFlags,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       tunfoldTikhonovLCurve.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 
       tunfoldTikhonovLCurve.SetInput(input,biasScale);//,0.,covarianceM);
       // tunfoldTikhonovLCurve.SetBias(mcgen);
       int iBest_TikhonovLCurve=-1;
       double tauBest_TikhonovLCurve=-1.,DF_TikhonovLCurve=-1.;
       int iBest;
       int NPOINT_TikhonovLCurve=50;//200;
       TGraph *graph_LCurve_TikhonovLCurve;
       TSpline *spline_Curvature_TikhonovLCurve;
       taumx = 1e-5; taumi = 1e-9;
       //taumx = 0.0; taumi = 0.0;
       
       iBest=tunfoldTikhonovLCurve.ScanLcurve(NPOINT_TikhonovLCurve, taumi, taumx, &graph_LCurve_TikhonovLCurve ,0,0, &spline_Curvature_TikhonovLCurve );
	      
       sprintf(unfoldhist, "Tunfold_lscan_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       sprintf(title, "Tunfolded Tikhonov Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(Rhoname, "Tunfold_lscan_1dcorr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       sprintf(Rhotitle, "correlation coefficients Tikhonov Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(Rhoname2d, "Tunfold_lscan_corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(Rhotitle2d, "2D correlation coefficients Tikhonov Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(Ematrix, "Tunfold_lscan_Emat_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(Ematrixtitle, "covariance matrix all contributions Tikhonov Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(probMat, "Tunfold_lscan_probM_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       sprintf(probmat_title, "Probability matrix Tikhonov Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       sprintf(foldback, "Tunfold_lscan_Refold_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
       sprintf(foldback_title, "Folded back  Lcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
       
       iBest_TikhonovLCurve=iBest;
       tauBest_TikhonovLCurve = tunfoldTikhonovLCurve.GetTau();
       DF_TikhonovLCurve = tunfoldTikhonovLCurve.GetDF();
       
       cout << " Tau Value = " << tauBest_TikhonovLCurve << "  ibest ="<< iBest << endl;
       file << "ibest ="<< iBest << "/ Tau =" << tauBest_TikhonovLCurve <<"       ";
       
       TH1 *hist_PTunfolded_TikhonovLCurve = tunfoldTikhonovLCurve.GetOutput(unfoldhist, title,0,"*[UO]" ,true);
       //           TH1 *hist_PTunfolded_TikhonovLCurve = tunfoldTikhonovLCurve.GetOutput(unfoldhist, title);
       TH1 *hist_Rho_TikhonovLCurve = tunfoldTikhonovLCurve.GetRhoItotal(Rhoname, Rhotitle);//,"signal");
       TH2 *hist_RhoIJ_TikhonovLCurve = tunfoldTikhonovLCurve.GetRhoIJtotal(Rhoname2d, Rhotitle2d);//,"signal");
       TH2 *hist_Ematrix_TikhonovLCurve = tunfoldTikhonovLCurve.GetEmatrixTotal(Ematrix, Ematrixtitle);//"*[UO]");//,"signal");
       //  TH2 *hist_Ematrix_TikhonovLCurve = tunfoldTikhonovLCurve.GetEmatrix(EmatrixL, EmatrixtitleL);//,"signal");
       TH2 *hist_prob_TikhonovLCurve = tunfoldTikhonovLCurve.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");
       TH1 *hist_foldedback_Lcurve= tunfoldTikhonovLCurve.GetFoldedOutput(foldback, foldback_title);//,"signal");
       

       // correction for Miss entries
       for (int i = 1; i <= hist_PTunfolded_TikhonovLCurve->GetNbinsX(); ++i) {
         double content = hist_PTunfolded_TikhonovLCurve->GetBinContent(i);
         double factor = 1;
         factor += MC_miss[ity][ivar][ipt]->GetBinContent(i);
         content *= factor;
         hist_PTunfolded_TikhonovLCurve->SetBinContent(i, content);
       }

       hist_PTunfolded_TikhonovLCurve->Write();
       hist_Rho_TikhonovLCurve->Write();
       hist_RhoIJ_TikhonovLCurve->Write();
       hist_Ematrix_TikhonovLCurve->Write();
       hist_prob_TikhonovLCurve->Write();
       hist_foldedback_Lcurve->Write();
       
       sprintf(Lcursure, "Tunfold_lscan_Lcurve_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       sprintf(TauLsure, "Tau_TikhLscan_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       // save auxillary plots
       graph_LCurve_TikhonovLCurve->Write(Lcursure);
       spline_Curvature_TikhonovLCurve->Write(TauLsure);
       
       
       cout << "L curve Scan  ok" << endl;
       //----------------------------------------------//ScanSURE
       /*
	 {
	 TUnfoldDensity tunfoldTikhonovSURE
	 (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert, regMode, constraintMode, densityFlags, 0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
	 tunfoldTikhonovSURE.SetInput(input,biasScale);
	 tunfoldTikhonovSURE.SetBias(mcgen);
	 
	 int NPOINT_TikhonovSURE = 200;
         taumx = 1e-3; taumi = 1e-9;
	 //   taumx = 0.0; taumi = 0.0;
         TGraph *graph_logTauSURE_TikhonovSURE,*graph_dfChi2A_TikhonovSURE;
	 //const char *SCAN_DISTRIBUTION=0;
	 //const char *SCAN_AXISSTEERING=0;
	 //     int ibest_sure = tunfoldTikhonovSURE.ScanSURE (NPOINT_TikhonovSURE,0.,0., &graph_logTauSURE_TikhonovSURE, &graph_dfChi2A_TikhonovSURE, 0);
	 int ibest_sure = tunfoldTikhonovSURE.ScanSURE (NPOINT_TikhonovSURE,taumi,taumx, &graph_logTauSURE_TikhonovSURE, &graph_dfChi2A_TikhonovSURE, 0);
	 double tau_SURE = tunfoldTikhonovSURE.GetTau();
	 //     cout << " Ibest Sure Scan : "<<ibest_sure <<" SURE Tau Value : " << tau_SURE << endl;
	 //     file << " Ibest : "<<ibest_sure <<" //Tau = " << tau_SURE <<"      ";
	 sprintf(unfoldhist, "Tunfold_SURE_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(title, "Tunfolded Tikhonov SURE %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Rhoname, "Tunfold_SURE_1Dcorr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Rhotitle, "correlation coefficients Tikhonov: SURE %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Rhoname2d, "Tunfold_SURE_corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Rhotitle2d, "2D correlation coefficients Tikhonov: SURE %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Ematrix, "Tunfold_SURE_Emat_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Ematrixtitle, "covariance matrix all contributions Tikhonov: SURE %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(probMat, "Tunfold_SURE_probM_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(probmat_title, "Matrix of probabilities Tikhonov: SURE %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(foldback, "Tunfold_SURE_Refold_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(foldback_title, "Folded back  ScanSure %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 
	 
	 //unfolding result, signal only
	 TH1 *hist_PTunfolded_TikhonovSURE = tunfoldTikhonovSURE.GetOutput(unfoldhist, title);//,"signal");
	 TH1 *hist_Rho_TikhonovSURE = tunfoldTikhonovSURE.GetRhoItotal(Rhoname, Rhotitle);//,"signal");
	 TH2 *hist_Rhoij_TikhonovSURE = tunfoldTikhonovSURE.GetRhoIJtotal(Rhoname2d, Rhotitle2d);//,"signal");
	 TH2 *hist_Ematrix_TikhonovSURE = tunfoldTikhonovSURE.GetEmatrixTotal(Ematrix,Ematrixtitle);//,"signal");
	 TH2 *hist_Prob_TikhonovSURE = tunfoldTikhonovSURE.GetProbabilityMatrix(probMat, probmat_title);//,"signal");
	 TH1 *hist_foldedback_SURE = tunfoldTikhonovSURE.GetFoldedOutput(foldback, foldback_title);//,"signal");
	 
	 hist_PTunfolded_TikhonovSURE->Write();
	 hist_Rho_TikhonovSURE->Write();
	 hist_Rhoij_TikhonovSURE->Write();
	 hist_Ematrix_TikhonovSURE->Write();
	 hist_Prob_TikhonovSURE->Write();
	 hist_foldedback_SURE->Write();
	 
	 }
       */
//---------------------------------------------------Scan [Tau Minimization of a global correlation coefficient (Stefan Schmitt)]
 /*      {
	 TUnfoldDensity tunfoldTauSCAN
	   (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert,regMode, constraintMode,densityFlags);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
         tunfoldTauSCAN.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 

	 tunfoldTauSCAN.SetInput(input,biasScale);
	 //        tunfoldTauSCAN.SetBias(mcgen);
	 int NScan = 100;
	 //taumi= 1e-4; taumx = 1e-9;
	 taumi=0.; taumx =0.;
	 
	 TGraph *lCurve;
	 TSpline *Scanlcurve;
	 TSpline *logTauX,*logTauY;
	 const char *SCAN_DISTRIBUTION=0;
	 const char *SCAN_AXISSTEERING=0;//"*[UOB]";//"*[b]" ;
	 //int ibest_tau = tunfoldTauSCAN.ScanTau (NScan,taumi,taumx, &Scanlcurve,TUnfoldDensity::kEScanTauRhoAvgSys, SCAN_DISTRIBUTION, SCAN_AXISSTEERING,&lCurve);
	 //int ibest_tau = tunfoldTauSCAN.ScanTau (NScan,0.,0., &Scanlcurve,TUnfoldDensity::kEScanTauRhoAvg, SCAN_DISTRIBUTION, SCAN_AXISSTEERING,&lCurve);
	 //int ibest_tau = tunfoldTauSCAN.ScanTau (NScan,0.,0., &Scanlcurve,TUnfoldDensity::kEScanTauRhoMax, SCAN_DISTRIBUTION, SCAN_AXISSTEERING,&lCurve);
	 int ibest_tau = tunfoldTauSCAN.ScanTau (NScan,taumi,taumx, &Scanlcurve,TUnfoldDensity::kEScanTauRhoAvg, SCAN_DISTRIBUTION, SCAN_AXISSTEERING,&lCurve);
	 double tauScan = tunfoldTauSCAN.GetTau();
	 cout << " ibest : "<<ibest_tau <<" Scan Tau : " << tauScan << endl;
	 file << " ibest : "<<ibest_tau <<" / Tau : " << tauScan << endl;
	 sprintf(unfoldhist, "Tunfold_scantau_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(title, "Tunfolded Tikhonov ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Rhoname, "Tunfold_scantau_1dcorr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Rhotitle, "correlation coefficients ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Rhoname2d, "Tunfold_scantau_corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Rhotitle2d, "2D correlation coefficients ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(Ematrix, "Tunfold_scantau_Emat_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(Ematrixtitle, "covariance matrix all contributions ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(probMat, "Tunfold_scantau_probM_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(probmat_title, "Matrix of probabilities ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 sprintf(foldback, "Tunfold_scantau_Refold_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	 sprintf(foldback_title, "Folded back  ScanTau %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	 
	 //unfolding result
	 TH1 *hist_PTunfolded_TauScan = tunfoldTauSCAN.GetOutput(unfoldhist, title);//,"signal");
	 TH1 *hist_Rho_TauScan = tunfoldTauSCAN.GetRhoItotal(Rhoname, Rhotitle);//,"signal");
	 TH2 *hist_Rhoij_TauScan = tunfoldTauSCAN.GetRhoIJtotal(Rhoname2d, Rhotitle2d);//,"signal");
	 TH2 *hist_Ematrix_TauScan = tunfoldTauSCAN.GetEmatrixTotal(Ematrix,Ematrixtitle);//,"signal");
	 TH2 *hist_Prob_TauScan = tunfoldTauSCAN.GetProbabilityMatrix(probMat, probmat_title);//,"signal");
	 TH1 *hist_foldedback_TauScan = tunfoldTauSCAN.GetFoldedOutput(foldback, foldback_title);//,"signal");

         // correction for  Miss entries 
       for (int i = 1; i <= hist_PTunfolded_TauScan->GetNbinsX(); ++i) {
         double content = hist_PTunfolded_TauScan->GetBinContent(i);
         double factor = 1;
         factor += MC_miss[ity][ivar][ipt]->GetBinContent(i);
         content *= factor;
         hist_PTunfolded_TauScan->SetBinContent(i, content);
       }
   	 
	 //      Double_t t[1],x[1],y[1];
	 //      logTauX->GetKnot(ibest_tau,t[0],x[0]);
	 //      logTauY->GetKnot(ibest_tau,t[0],y[0]);
	 //        TGraph *bestLcurve=new TGraph(1,x,y);
	 //      TGraph *bestLogTauLogChi2=new TGraph(1,t,x);
	 
	 
	   
	   sprintf(title, "bestLcurve %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	   bestLcurve->SetTitle(title);
	   
	   sprintf(histname, "bestLcurve_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	   bestLcurve->SetName(histname);
	   
	   sprintf(histname, "Taulogchi2_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
	   bestLogTauLogChi2->SetName(histname);
	   
	   sprintf(title, "TaulogChi2 %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
	   bestLogTauLogChi2->SetTitle(title);
	   
	   
	   bestLogTauLogChi2->Write();
	   bestLcurve->Write();
	   //    Scanlcurve->SetTitle(title);
	   
	 Scanlcurve->Write();
	 
	 
	 hist_PTunfolded_TauScan->Write();
	 hist_Rho_TauScan->Write();
	 hist_Rhoij_TauScan->Write();
	 hist_Ematrix_TauScan->Write();
	 hist_Prob_TauScan->Write();
	 hist_foldedback_TauScan->Write();
       }
       */
       //-----------------------------------------End of Variables loop
     }
   }
 }
 
 
 delete outputFile;
 file.close();
}




//RooUnfold Subtract Background (Fake, miss ) method : Now replace by patrick code
int subtract_background1(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl) {  
  int nbinx = h2d_correl->GetNbinsX();
  int nbiny = h2d_correl->GetNbinsY();
  const int nbinmx = 100 ;
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};
  
  for (int ix=0; ix<nbinx+1; ix++) {
    for (int iy=0; iy<nbiny+1; iy++) {
      if(ix==0&&iy==0) continue ;
      totalreco[ix] +=h2d_correl->GetBinContent(ix, iy);
      if (iy==0) fakerate[ix-1] = h2d_correl->GetBinContent(ix,iy);      //fake from GEN Underflow bins
      totalgen[iy] +=h2d_correl->GetBinContent(ix, iy);
      if (ix==0) effi[iy-1] =h2d_correl->GetBinContent(ix, iy);         //miss from RECO Underflow bin 
    }//iy
  }//ix
  
  for (int iy=0; iy<nbiny; iy++) {
    effi[iy] = (totalgen[iy+1] - effi[iy])/max(1.e-10, totalgen[iy+1]);
    if(iy>10 && effi[iy]>0.1 && effi[iy]<0.9){
      //       effi[iy] = effi[iy-1]-0.001;
    }
    if (effi[iy]<1.e-6) effi[iy]=1.e-6;
  } //iy
  
  for (int ix=0; ix<nbinx; ix++){
    if(totalreco[ix+1]>1.e-6) {     fakerate[ix] /= totalreco[ix+1]; }
    else { fakerate[ix] = 0.; }
    
    // if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
  }//ix
  
  for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
    reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
    reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
  }
  
  for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
    data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
    data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
  }
  
  for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
    // gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]);
    // gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]);
  }
  
  double tot_reco_in[nbinmx] = {0.};
  double tot_gen_in[nbinmx] = {0.};
  
  for(int bn=0; bn<nbinx; bn++){
    for(int am=0; am<nbiny; am++){
      tot_gen_in[am]+=h2d_correl->GetBinContent(bn+1,am+1);
      tot_reco_in[bn]+=h2d_correl->GetBinContent(bn+1,am+1);
    }
  }
  
  for(int bn=0; bn<nbinx; bn++){
    purity[bn] = (tot_reco_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_reco_in[bn]:0.;
    stbl[bn] = (tot_gen_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_gen_in[bn]:0.;
  }
  return 0;
} //end substract func


//Copy of RooUnfold Subtract Background method with no correction to Reco or Gen --->Can be used  for purity and stability
int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl) {
  int nbinx = h2d_correl->GetNbinsX();
  int nbiny = h2d_correl->GetNbinsY(); 
  const int nbinmx = 100 ;
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};
  for (int ix=0; ix<nbinx+1; ix++) {
    for (int iy=0; iy<nbiny+1; iy++) {      
      if(ix==0&&iy==0) continue ;
      totalreco[ix] +=h2d_correl->GetBinContent(ix, iy);
      if (iy==0) fakerate[ix-1] = h2d_correl->GetBinContent(ix,iy);
      totalgen[iy] +=h2d_correl->GetBinContent(ix, iy);
      if (ix==0) effi[iy-1] =h2d_correl->GetBinContent(ix, iy);
      
    }//iy
  }//ix
  
  for (int iy=0; iy<nbiny; iy++) {
    effi[iy] = (totalgen[iy+1] - effi[iy])/max(1.e-10, totalgen[iy+1]);
    
   // if(iy>10 && effi[iy]>0.1 && effi[iy]<0.9){
   //       effi[iy] = effi[iy-1]-0.001;
   // }
    
  //  if (effi[iy]<1.e-6) effi[iy]=1.e-6;
   } //iy
  
  for (int ix=0; ix<nbinx; ix++){
    if(totalreco[ix+1]>1.e-6) {     fakerate[ix] /= totalreco[ix+1]; }
    else { fakerate[ix] = 0.; }
    
    // if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
  }//ix
  
 /// for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
    //        reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
    //        reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
 // }
  
 // for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
    //       data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
    //        data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
 // }
  
//  for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
    //  gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]) ;
    //  gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]) ;
//  }
  
  double tot_reco_in[nbinmx] = {0.};
  double tot_gen_in[nbinmx] = {0.};
  
  for(int bn=0; bn<nbinx; bn++){
    for(int am=0; am<nbiny; am++){
      tot_gen_in[am]+=h2d_correl->GetBinContent(bn+1,am+1);
      tot_reco_in[bn]+=h2d_correl->GetBinContent(bn+1,am+1);
    }
  }
  
  for(int bn=0; bn<nbinx; bn++){
    purity[bn] = (tot_reco_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_reco_in[bn]:0.;
    stbl[bn] = (tot_gen_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_gen_in[bn]:0.;
  }
  return 0;
} //end substract func



//1D Rebin function:Basically it cut off some bin you  don't needed
TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]){
  TH1D* thout;
  int nxmod2 = -1; double xmod2[300];
  const char *namex; const char *titlex;
  char namey[100], titley[100];
  
  double yvl[200]={0.0};
  double erryvl[200]={0.0};
  int nbinx = thin->GetNbinsX();
  int ifirst = 0;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
      } else {
        xmod2[ifirst] = thin->GetBinLowEdge(ij);
        xmod2[ifirst+1] = thin->GetBinLowEdge(ij+1);
        yvl[ifirst] =thin->GetBinContent(ij);
        erryvl[ifirst] +=thin->GetBinError(ij)*thin->GetBinError(ij);
        ifirst++;
      }
    }
    nxmod2 = ifirst;
  }
  
  double xmod3[nxmod2+1];
  for(int ii=0; ii < nxmod2+1 ; ii++){ xmod3[ii]=xmod2[ii]; }
  
  namex = thin->GetName();
  sprintf(namey, "rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "Rebin_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod3);
  
  for (int ix=0; ix<nxmod2+1; ix++) {
    thout->SetBinContent(ix+1, yvl[ix]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix]));     }
  //thout->SetBinContent(0,5);
  return thout;
}   //end of rebin_hist functiomC


//2D Rebin function:Basically it cut off some bin you  don't needed
TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco, TH1D* MC_gen, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8],int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8] ){
  TH2D* thout;
  cout << "Reco rebin:" <<endl;
  double xmod2[200]; double ymod2[200];
  const char* namex; const char* titlex;
  char namey[100],titley[100];
  
  double yvl[200][200]={0.0};
  double erryvl[200][200]={0.0};
  int nbinx = thin->GetNbinsX();
  int nbiny = thin->GetNbinsY();
  cout << "   X bin :=" << nbinx << "   Y bin  : = " <<nbiny <<endl; 
  int ifirst = 0; int jfirst =0;
  int  nreco= -1; int  ngen= -1;
  int xfirst =0;int yfirst =0;

 //Get Reco Level entries
  if (nreco<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
      } else {
        xmod2[ifirst] = MC_reco->GetBinLowEdge(ij);
        xmod2[ifirst+1] = MC_reco->GetBinLowEdge(ij+1);
        ifirst++;
      }
    }
  }
//Get Gen Level entries
  if (ngen<0) {
    for (int ij=0; ij<nbiny+2; ij++) {
      if (ij <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( ij > (nbiny - arrayvar_lastG[itype][ivar][ijetpt])) {
      } else {
        ymod2[jfirst] = MC_gen->GetBinLowEdge(ij);
        ymod2[jfirst+1] = MC_gen->GetBinLowEdge(ij+1);
        jfirst++;
      }
    }
  }
  if (xfirst==0 && yfirst==0) {
    for (int iy=0; iy<nbiny+2; iy++) {
      if (iy <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( iy > (nbiny - arrayvar_lastG[itype][ivar][ijetpt]) ) {
      } else {
	xfirst=0;
      for (int ix=0; ix<nbinx+2; ix++) {
        if (ix <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if (ix > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
      } else {
        yvl[xfirst][yfirst] = thin->GetBinContent(ix,iy);
        erryvl[xfirst][yfirst] += thin->GetBinError(ix,iy)*thin->GetBinError(ix,iy);
        xfirst++;
        }
      }
      yfirst++;
    }
  }
}

cout << "Reco Bins = " << ifirst <<"  " << xfirst << endl;
cout << "Gen Bins = " << jfirst <<"  " << yfirst << endl;

    nreco = ifirst;
    ngen = jfirst;
double xmod3[nreco+1];
for(int ii=0; ii < nreco+1 ; ii++){ xmod3[ii]=xmod2[ii]; }
double xmod4[ngen+1];
for(int ii=0; ii < ngen+1 ; ii++){  xmod4[ii]=ymod2[ii]; }

  namex = thin->GetName();
  sprintf(namey, "Rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "%s", titlex);
  thout = new TH2D(namey, titley, nreco, xmod3, ngen, xmod4);

  for (int ix=0; ix<nreco+1; ix++) {
  for (int iy=0; iy<ngen+1; iy++) {
    thout->SetBinContent(ix+1,iy+1, yvl[ix][iy]);
    thout->SetBinError(ix+1,iy+1, sqrt(erryvl[ix][iy]));
        }
     }
  return thout;
}   //end of rebin_hist functiomC


//1D Rebin function:Basically it cut off some bin you  don't needed
TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8]){
TH1D* thout;
int nxmod2 = -1; double xmod2[300];

const char *namex; const char *titlex;
char namey[100], titley[100];

  double yvl[200]={0.0};
  double erryvl[200]={0.0};
  int nbinx = thin->GetNbinsX();
  
  int ifirst = 0;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_lastG[itype][ivar][ijetpt])) {
      } else { 
        xmod2[ifirst] = thin->GetBinLowEdge(ij);
        xmod2[ifirst+1] = thin->GetBinLowEdge(ij+1);
        yvl[ifirst] =thin->GetBinContent(ij);
        erryvl[ifirst] +=thin->GetBinError(ij)*thin->GetBinError(ij);
        ifirst++;
      }
    }
    nxmod2 = ifirst;
  }
double xmod3[nxmod2+1];
for(int ii=0; ii < nxmod2+1 ; ii++){ xmod3[ii]=xmod2[ii]; }

  namex = thin->GetName();
  sprintf(namey, "rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "Rebin_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod3);

  for(int ix=0; ix<nxmod2+1; ix++) {
    thout->SetBinContent(ix+1, yvl[ix]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix]));     }
  //  thout->SetBinContent(0,5);
    return thout;
}   //end of rebin_hist functiomC


void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect){

//calculate Probability Matrix
for(int ij=1; ij<(HistoMatrix->GetNbinsY()+1); ij++){
double row_sum = 0.;
for(int jk=1; jk<(HistoMatrix->GetNbinsX()+1); jk++){
  row_sum+=HistoMatrix->GetBinContent(jk,ij);
}//jk
if(row_sum>1.e-10){
 for(int jk=1; jk<(HistoMatrix->GetNbinsX()+1); jk++){
   HistoMatrix->SetBinContent(jk,ij,HistoMatrix->GetBinContent(jk,ij)*1./row_sum) ; //Probability
  }//jk
}
}//ij

//folding gen level to Reco
for(int i=1;i<HistReco->GetNbinsX()+1;i++){
     double sum=0.; double Err =0.;
       for(int j=1;j<HistGen->GetNbinsX()+1;j++){
       double misscorr = (HistGen->GetBinContent(j))-(miss->GetBinContent(j)); //Miss correction
       sum += HistoMatrix->GetBinContent(i,j)*misscorr;
       Err += (HistoMatrix->GetBinContent(i,j)*HistGen->GetBinError(j))*(HistoMatrix->GetBinContent(i,j)*(HistGen->GetBinError(j)));
           }
      sum = sum +(fake->GetBinContent(i)); //fake correction
      HistoCorrect->SetBinContent(i,sum);
      HistoCorrect->SetBinError(i,sqrt(Err));
    	}
}//end Fold

//Condition number calculation
void Condition (TH2* RM, TH1* miss){
    const int Nx = RM->GetNbinsX(),
              Ny = RM->GetNbinsY();
    cout << Nx <<"  " << Ny << endl;
    if (Ny*2 != Nx) { cout << Nx << ' ' << Ny << endl;  return; }

    TH1D* RMy = RM->ProjectionY("RMx", 0, -1); //Gen Projection

    // normalisation & condition
    TMatrixD m(Ny,Ny);
    for (int i = 1; i <= Ny; ++i) {
        double normalisation = RMy->GetBinContent(i);
               normalisation += miss->GetBinContent(i);
        if (normalisation > 0)
        for (int j = 1; j <= Nx; ++j) {
            double content = RM->GetBinContent(j,i);
            content /= normalisation;
            m((j-1)/2,i-1) += content;
        }
    }
    TDecompSVD svd(m);
    TVectorD v = svd.GetSig();

    double Max = v[0];
    for(int k = 0; k < Ny; ++k) {
        if (abs(v[k]) < feps) break;
        cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/v[k] << '\n';
    }
}



