// Authors: Suman Kumar Kundu
//#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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
#include <TPad.h>
#include <TLine.h>
#include <TRandom.h>

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldResponse.h"
#include "TLegend.h"


#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"

#define REBIN
//#define Roo

using namespace std;

//#endif


void RooUnfoldExample(){
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

const int svd_par[type][nusedvar][njetptmn]= {{{4,5,5,6,6,7,6,6},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7}},
                                              {{5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7}}};

const int bayes_par[type][nusedvar][njetptmn]= {{{4,4,4,4,4,4,4,4},
                                              {5,5,5,5,5,7,7,7},
                                              {3,3,3,4,4,4,4,4},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7}},
                                              {{5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7},
                                              {5,5,5,5,5,7,7,7}}};
//Define all required function

//int getbinid(double val, int nbmx, double* array);
TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar);
TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco, int itype, int ijetpt, int ivar);

void subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
void Fold(TH2D* HistoMatrix, TH1D* HistoGen, TH1D* HistoUnfold, TH1D* HistoCorrect);
void ReFold(TH2D* HistoMatrix,  TH1D* HistoGen, TH1D* HistoUnfold, TH1D* HistoCorrect);
void BackFold_Gen(RooUnfoldResponse response_b, TH1D* hist_gen, TH1D* hist_back,double *eff,double* fakerate);
void BackFold_Bayesian(RooUnfoldResponse response_b, RooUnfoldBayes unfold, TH1D* hist_back,double *eff,double* fakerate);
void BackFold_Svd(RooUnfoldResponse response_b, RooUnfoldSvd unfold, TH1D* hist_back,double *eff,double* fakerate);
void BackFold_BinbyBin(RooUnfoldResponse response_b, RooUnfoldBinByBin unfold, TH1D* hist_back,double *eff,double* fakerate);


//Input Data and MC histogram 
 //TFile *inputData=new TFile("PY8_UL2017_120Bin_14june20.root");
//  TFile *inputData=new TFile("PY8_UL17_ULbinned_1July20.root");
     TFile *inputData=new TFile("DATA_JetHT_UL17_ULbin_4July.root");
   // TFile *inputData=new TFile("PY8_UL17_binnedUL_1July20_NOweight.root");

 //TFile *inputMC=new TFile("PY8_UL2017_120Bin_14june20.root");
   TFile *inputMC=new TFile("PY8_UL17_ULbinned_1July20.root");
  // TFile *inputMC=new TFile("PY8_UL17_binnedUL_1July20_NOweight.root");
  //TFile *inputMC=new TFile("PY8_C5_1523_14PU_GT102X_3sig_16April.root");

  //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
    TFile *outputFile=new TFile("Roounfold_closure.root","recreate");


         //----------------------for reco and gen bin number
  int rnbinsx[type][nusedvar][nHLTmx];
  int gnbinsx[type][nusedvar][nHLTmx];

  ///////////////////Modify /////////////////////////////

   TH1D *MC_Reco[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco

   TH1D *MC_Reco1[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco1[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco


   TH1D *MC_Reco2[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco2[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC2[type][nusedvar][njetptmn];

   TH1D *MC_Reco_re[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MC_Gen_re[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco_re[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco


   TH1D *unfold_svd[type][nusedvar][njetptmn];
   TH1D *unfold_bayes[type][nusedvar][njetptmn];
   TH1D *unfold_BbB[type][nusedvar][njetptmn];

   TH1D *unfold_bayes_back[type][nusedvar][njetptmn];
   TH1D *unfold_svd_back[type][nusedvar][njetptmn];
   TH1D *unfold_BbB_back[type][nusedvar][njetptmn];


   TH1D* hist_eff[type][nusedvar][njetptmn];
   TH1D* hist_fake[type][nusedvar][njetptmn];
   TH1D* hist_purity[type][nusedvar][njetptmn];
   TH1D* hist_stbl[type][nusedvar][njetptmn];

   TH1D* hist_eff1[type][nusedvar][njetptmn];
   TH1D* hist_fake1[type][nusedvar][njetptmn];
   TH1D* hist_purity1[type][nusedvar][njetptmn];
   TH1D* hist_stbl1[type][nusedvar][njetptmn];

   TH2D* COV_Mat_svd[type][nusedvar][njetptmn];
   TH2D* COV_Mat_bayes[type][nusedvar][njetptmn];
   TH2D* COV_Mat_BbB[type][nusedvar][njetptmn];

   TMatrixD covmatrix_bayes[type][nusedvar][njetptmn];
   TMatrixD covmatrix_svd[type][nusedvar][njetptmn];
   TMatrixD covmatrix_BbB[type][nusedvar][njetptmn];

   TMatrixD corrmatrix_bayes[type][nusedvar][njetptmn];
   TMatrixD corrmatrix_svd[type][nusedvar][njetptmn];
   TMatrixD corrmatrix_BbB[type][nusedvar][njetptmn];


   TDirectoryFile *inputDir1=new TDirectoryFile("Data","Inputs Data");
   TDirectoryFile *inputDir2=new TDirectoryFile("MC","Inputs MC and Probability Matrix");
   TDirectoryFile *inputRebin1=new TDirectoryFile("RebinData","Inputs Data Rebin");
   TDirectoryFile *inputRebin2=new TDirectoryFile("RebinMC","Inputs MC Rebin");
   TDirectoryFile *Unfold=new TDirectoryFile("Unfold","Unfolded, Refold, correlation,effi,purity,COV");



//Read Input Data MC and Response matrix 
   for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){

         //Reco Data
         inputDir1->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         Data_Reco1[ity][ivar][ipt] = (TH1D*) inputData->Get(histname);
 //      Data_Reco1[ity][ivar][ipt]->Rebin(2);
         Data_Reco1[ity][ivar][ipt]->Write();
       //for (int ibin =1 ; ibin <  Data_Reco1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
       // if (Data_Reco1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " Data Reco Bin is Zero for bin number : ******** "<<  ibin  << endl; }}

         inputDir2->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
         cout << histname << endl;
         int recobins = MC_Reco1[ity][ivar][ipt]->GetNbinsX();
         rnbinsx[ity][ivar][ipt]=MC_Reco1[ity][ivar][ipt]->GetNbinsX();
         if(recobins !=Data_Reco1[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
         MC_Reco1[ity][ivar][ipt]->Write();
         //for(int ibin = 1; ibin < MC_Reco1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
        // if (MC_Reco1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " MC reco Bin is Zero for bin number :** "<<  ibin  << "  Bin lowedge "<< MC_Reco1[ity][ivar][ipt]->GetBinLowEdge(ibin)<< endl; }}

         //Gen MC
         sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
  //       MC_Gen1[ity][ivar][ipt]->Rebin(2);
         MC_Gen1[ity][ivar][ipt]->Write();
         int genbins = MC_Gen1[ity][ivar][ipt]->GetNbinsX();
         gnbinsx[ity][ivar][ipt]=MC_Gen1[ity][ivar][ipt]->GetNbinsX();
         //for(int ibin = 1; ibin < MC_Gen1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
        // if (MC_Gen1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  ibin  << endl; }}


         //Response Matrix
         sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
         h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC->Get(histname);   //Xgen(coarse) , Yreco(fine)
 //        h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
         h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)
         if(recobins !=h2dGenDetMC1[ity][ivar][ipt]->GetNbinsX()){cout << "Reco Bin missmatch in response Matrix " << endl;}
         if(genbins !=h2dGenDetMC1[ity][ivar][ipt]->GetNbinsY()){cout << "Gen Bin missmatch in response Matrix " << endl;}
         cout << "Recobins = " <<  recobins  << "   Gen Bins = " << genbins << endl;

         Data_Reco2[ity][ivar][ipt] = (TH1D*)Data_Reco1[ity][ivar][ipt]->Clone();
         MC_Reco2[ity][ivar][ipt] = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
         MC_Gen2[ity][ivar][ipt] = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
         h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();


#ifdef REBIN

        Data_Reco_re[ity][ivar][ipt] = rebin1d_hist(Data_Reco1[ity][ivar][ipt], ity, ipt, ivar);
       // Data_Reco_re[ity][ivar][ipt]->Rebin(4);
        cout << "Rebin Data Reco Bin : " << Data_Reco_re[ity][ivar][ipt]->GetNbinsX() << endl;

        inputRebin2->cd();
        MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
       // MC_Reco_re[ity][ivar][ipt]->Rebin(4);

        MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
       // MC_Gen_re[ity][ivar][ipt]->Rebin(8);
        cout << " Rebin MC Gen bin : " <<  MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;

        h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(h2dGenDetMC1[ity][ivar][ipt], MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
       // h2dGenDetMC_re[ity][ivar][ipt]->RebinY(8);
       // h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
        cout << "ReBin REco Bin :"<< h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;



//Define data MC to be used for unfolding
         Data_Reco[ity][ivar][ipt] = (TH1D*)Data_Reco_re[ity][ivar][ipt]->Clone();
         MC_Reco[ity][ivar][ipt] = (TH1D*)MC_Reco_re[ity][ivar][ipt]->Clone();
         MC_Gen[ity][ivar][ipt] = (TH1D*)MC_Gen_re[ity][ivar][ipt]->Clone();
         h2dGenDetMC[ity][ivar][ipt] = (TH2D*)h2dGenDetMC_re[ity][ivar][ipt]->Clone();


#else
         Data_Reco[ity][ivar][ipt] = (TH1D*)Data_Reco1[ity][ivar][ipt]->Clone();
         MC_Reco[ity][ivar][ipt] = (TH1D*)MC_Reco1[ity][ivar][ipt]->Clone();
         MC_Gen[ity][ivar][ipt] = (TH1D*)MC_Gen1[ity][ivar][ipt]->Clone();
         h2dGenDetMC[ity][ivar][ipt] = (TH2D*)h2dGenDetMC1[ity][ivar][ipt]->Clone();

#endif
         inputRebin1->cd();
         Data_Reco[ity][ivar][ipt]->Write();
         inputRebin2->cd();
         MC_Reco[ity][ivar][ipt]->Write();
         MC_Gen[ity][ivar][ipt]->Write();
         h2dGenDetMC[ity][ivar][ipt]->Write();
       }
     }
   }
   cout << "Histogram Read and Rebin oK " <<endl;


 Unfold->cd();
 for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){

 //Get reco bins
 double rxbins[MC_Reco[ity][ivar][ipt]->GetNbinsX()+1]={0};
 for (int ix=0; ix<MC_Reco[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
          rxbins[ix] = MC_Reco[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
        }

 //Get Gen bins
 double gxbins[MC_Gen[ity][ivar][ipt]->GetNbinsX()+1]={0};
   for (int ix=0; ix<MC_Gen[ity][ivar][ipt]->GetNbinsX()+1; ix++) {
     gxbins[ix] = MC_Gen[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
   }

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

   sprintf(name,"unfold_SVD_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_svd[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_svd[ity][ivar][ipt]->Sumw2();
   //unfold_svd[ity][ivar][ipt]->Rebin(2);

   sprintf(name,"unfold_bayes_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_bayes[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_bayes[ity][ivar][ipt]->Sumw2();
   //unfold_bayes[ity][ivar][ipt]->Rebin(2);

   sprintf(name,"unfold_BbB_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_BbB[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_BbB[ity][ivar][ipt]->Sumw2();
   //unfold_BbB[ity][ivar][ipt]->Rebin(2);

   sprintf(name,"unfold_back_SVD_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_svd_back[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_svd_back[ity][ivar][ipt]->Sumw2();
   //unfold_svd_back[ity][ivar][ipt]->Rebin(2);

   sprintf(name,"unfold_back_bayes_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_bayes_back[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_bayes_back[ity][ivar][ipt]->Sumw2();
   //unfold_bayes_back[ity][ivar][ipt]->Rebin(2);

   sprintf(name,"unfold_back_BbB_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   unfold_BbB_back[ity][ivar][ipt] = new TH1D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(),gxbins);
   unfold_BbB_back[ity][ivar][ipt]->Sumw2();


   sprintf(name,"Cov_Matrix_svd_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   COV_Mat_svd[ity][ivar][ipt] = new TH2D(name,name, MC_Reco[ity][ivar][ipt]->GetNbinsX(), rxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
   COV_Mat_svd[ity][ivar][ipt]->Sumw2();
 //COV_Mat_svd[ity][ivar][ipt]->RebinY(2);

   sprintf(name,"Corr_Matrix_bayes_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   COV_Mat_bayes[ity][ivar][ipt] = new TH2D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
   COV_Mat_bayes[ity][ivar][ipt]->Sumw2();
 //COV_Ma[ity][ivar][ipt]->RebinY(2);
 //corr_mat_NoReg[ity][ivar][ipt]->RebinX(2);
 
   

   sprintf(name,"Corr_Matrix_BbB_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   COV_Mat_BbB[ity][ivar][ipt] = new TH2D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
   COV_Mat_BbB[ity][ivar][ipt]->Sumw2();
       }
     }
 }
//Recheck the efficiency and purity same as Tufold code
 for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){
 //----------------------------------------
        for(int binGen=0;binGen<= h2dGenDetMC[ity][ivar][ipt]->GetNbinsY()+1;binGen++) {
         double sum0=0.;
         double sum1=0.;
         for(int binRec=0;binRec<= h2dGenDetMC[ity][ivar][ipt]->GetNbinsX()+1;binRec++) {
            double c=  h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
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
 hist_eff1[ity][ivar][ipt]->SetMinimum(-0.01); hist_eff1[ity][ivar][ipt]->SetMaximum(1.01);
hist_eff1[ity][ivar][ipt]->Write();


for(int binRec=0; binRec<= hist_purity1[ity][ivar][ipt]->GetNbinsX()+1; binRec++) {
         double sum=0.;
         for(int binGen=0; binGen<=hist_purity1[ity][ivar][ipt]->GetNbinsX()+1; binGen++) {
            sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
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
//--------------------------------------------------//Recheck the efficiency and purity same as Tufold code

//Start Unfolding
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
       subtract_background(h2dGenDetMC[ity][ivar][ipt],MC_Reco[ity][ivar][ipt],MC_Gen[ity][ivar][ipt],Data_Reco[ity][ivar][ipt],fakerate,eff,purity,stbl);

       for(int bn=0; bn<(MC_Gen[ity][ivar][ipt]->GetNbinsX()); bn++){
         hist_eff[ity][ivar][ipt]->SetBinContent(bn+1,eff[bn+1]);
         hist_fake[ity][ivar][ipt]->SetBinContent(bn+1,fakerate[bn+1]);
         hist_purity[ity][ivar][ipt]->SetBinContent(bn+1,purity[bn+1]);
         hist_stbl[ity][ivar][ipt]->SetBinContent(bn+1,stbl[bn+1]);
       }

       hist_eff[ity][ivar][ipt]->SetMinimum(-0.01); hist_eff[ity][ivar][ipt]->SetMaximum(1.01);
       hist_fake[ity][ivar][ipt]->SetMinimum(-0.01); hist_fake[ity][ivar][ipt]->SetMaximum(1.01);
       hist_purity[ity][ivar][ipt]->SetMinimum(-0.01); hist_purity[ity][ivar][ipt]->SetMaximum(1.01);
       hist_stbl[ity][ivar][ipt]->SetMinimum(-0.01); hist_stbl[ity][ivar][ipt]->SetMaximum(1.01);

       hist_eff[ity][ivar][ipt]->Write();
       hist_fake[ity][ivar][ipt]->Write();
       hist_purity[ity][ivar][ipt]->Write();
       hist_stbl[ity][ivar][ipt]->Write();


//Define Roounfold Response matrix
RooUnfoldResponse response_b;
sprintf(name,"unfold_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
response_b = RooUnfoldResponse(MC_Reco[ity][ivar][ipt], MC_Gen[ity][ivar][ipt] , h2dGenDetMC[ity][ivar][ipt], name, name);

TH1D *hist_input;

hist_input = (TH1D*)Data_Reco[ity][ivar][ipt]->Clone();

int iter = bayes_par[ity][ivar][ipt];
RooUnfoldBayes unfoldBayes(&response_b,hist_input,iter,false);

RooUnfoldSvd unfoldSvd(&response_b,hist_input,svd_par[ity][ivar][ipt]);

RooUnfoldBinByBin unfoldBbB(&response_b,hist_input);


covmatrix_bayes[ity][ivar][ipt].ResizeTo(h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(), h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());
covmatrix_bayes[ity][ivar][ipt] = unfoldBayes.Ereco(RooUnfold::kCovariance); //Get Covarince matrix
covmatrix_svd[ity][ivar][ipt].ResizeTo(h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(), h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());
covmatrix_svd[ity][ivar][ipt] = unfoldSvd.Ereco(RooUnfold::kCovariance);
covmatrix_BbB[ity][ivar][ipt].ResizeTo( h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(), h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());
covmatrix_BbB[ity][ivar][ipt] = unfoldBbB.Ereco(RooUnfold::kCovariance);

//cout<<hist_input->GetName()<<" Chisquare "<<unfoldBayes.Chi2(hist_gen[ij][jk],RooUnfold::kCovariance)<<endl;

corrmatrix_bayes[ity][ivar][ipt].ResizeTo(  h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(),  h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());
corrmatrix_svd[ity][ivar][ipt].ResizeTo( h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(),  h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());
corrmatrix_BbB[ity][ivar][ipt].ResizeTo( h2dGenDetMC[ity][ivar][ipt]->GetNbinsX(),  h2dGenDetMC[ity][ivar][ipt]->GetNbinsY());

for(int ix=0; ix<(h2dGenDetMC[ity][ivar][ipt]->GetNbinsX()); ix++){
    for(int iy=0; iy<(h2dGenDetMC[ity][ivar][ipt]->GetNbinsY()); iy++){
                corrmatrix_bayes[ity][ivar][ipt](ix,iy) = covmatrix_bayes[ity][ivar][ipt](ix,iy) * 1./sqrt(covmatrix_bayes[ity][ivar][ipt](ix,ix)*covmatrix_bayes[ity][ivar][ipt](iy,iy));
                corrmatrix_svd[ity][ivar][ipt](ix,iy) = covmatrix_svd[ity][ivar][ipt](ix,iy) * 1./sqrt(covmatrix_svd[ity][ivar][ipt](ix,ix)*covmatrix_svd[ity][ivar][ipt](iy,iy));
                corrmatrix_BbB[ity][ivar][ipt](ix,iy) = covmatrix_BbB[ity][ivar][ipt](ix,iy) * 1./sqrt(covmatrix_BbB[ity][ivar][ipt](ix,ix)*covmatrix_BbB[ity][ivar][ipt](iy,iy));
                }
        }

 corrmatrix_bayes[ity][ivar][ipt].Write();
 corrmatrix_bayes[ity][ivar][ipt].Write();
 corrmatrix_BbB[ity][ivar][ipt].Write();

unfold_bayes[ity][ivar][ipt] = (TH1D*) unfoldBayes.Hreco(RooUnfold::kCovariance)->Clone();
sprintf(name,"bayes_%s", unfold_bayes[ity][ivar][ipt]->GetName());
unfold_bayes[ity][ivar][ipt]->SetName(name);
for(int bn=0; bn<(unfold_bayes[ity][ivar][ipt]->GetNbinsX()); bn++){
//eff[bn] = 1.;
/*
hist_bayes[ij][jk]->SetBinContent(bn+1,(hist_bayes[ij][jk]->GetBinContent(bn+1))*1./eff[bn]);
hist_bayes[ij][jk]->SetBinError(bn+1,(hist_bayes[ij][jk]->GetBinError(bn+1))*1./sqrt(eff[bn]));
*/
}
unfold_bayes[ity][ivar][ipt]->Write();

unfold_svd[ity][ivar][ipt]  =  (TH1D*) unfoldSvd.Hreco(RooUnfold::kCovariance)->Clone();
sprintf(name,"svd_%s", unfold_svd[ity][ivar][ipt]->GetName());
unfold_svd[ity][ivar][ipt]->SetName(name);
/* for(int bn=0; bn<(hist_svd[ij][jk]->GetNbinsX()); bn++){
//eff[bn] = 1.;
hist_svd[ij][jk]->SetBinContent(bn+1,(hist_svd[ij][jk]->GetBinContent(bn+1))*1./eff[bn]);
hist_svd[ij][jk]->SetBinError(bn+1,(hist_svd[ij][jk]->GetBinError(bn+1))*1./(eff[bn]));
}*/
unfold_svd[ity][ivar][ipt]->Write();

unfold_BbB[ity][ivar][ipt]  =  (TH1D*) unfoldBbB.Hreco(RooUnfold::kCovariance)->Clone();
sprintf(name,"BbB_%s",unfold_BbB[ity][ivar][ipt]->GetName());
unfold_BbB[ity][ivar][ipt]->SetName(name);

unfold_BbB[ity][ivar][ipt]->Write();
/*
Fold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt], unfold_bayes[ity][ivar][ipt], unfold_bayes_back[ity][ivar][ipt]);
Fold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt], unfold_svd[ity][ivar][ipt], unfold_svd_back[ity][ivar][ipt]);
Fold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt], unfold_BbB[ity][ivar][ipt], unfold_BbB_back[ity][ivar][ipt]);
*/
/*
ReFold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt], unfold_bayes[ity][ivar][ipt], unfold_bayes_back[ity][ivar][ipt]);
ReFold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt], unfold_svd[ity][ivar][ipt], unfold_svd_back[ity][ivar][ipt]);
ReFold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt],  unfold_BbB[ity][ivar][ipt], unfold_BbB_back[ity][ivar][ipt]);
*/

BackFold_Bayesian(response_b,unfoldBayes, unfold_bayes_back[ity][ivar][ipt],eff,fakerate);
BackFold_Svd(response_b,unfoldSvd,unfold_svd_back[ity][ivar][ipt],eff,fakerate);
BackFold_BinbyBin(response_b,unfoldBbB, unfold_BbB_back[ity][ivar][ipt],eff,fakerate);

//Fold(h2dGenDetMC[ity][ivar][ipt], MC_Gen[ity][ivar][ipt],MC_Gen[ity][ivar][ipt],hist_gen_fold[ij][jk]);
//BackFold_Gen(response_b,hist_gen[ij][jk],hist_gen_fold[ij][jk],eff,fakerate);
//hist_gen_fold[ij][jk]->Divide(hist_reco2[ij][jk]) ;
unfold_bayes_back[ity][ivar][ipt]->Write();
unfold_svd_back[ity][ivar][ipt]->Write();
unfold_BbB_back[ity][ivar][ipt]->Write();


}
   }
}



   
   
  delete outputFile;
}//end of RooUnfoldExample


void subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl) {

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

        if (ix==0 || iy==0) {
          h2d_correl->SetBinContent(ix, iy, 0.0);
          h2d_correl->SetBinError(ix, iy, 0.0);
        }
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

 if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
 }//ix

 for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
        reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
        reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
  }

 for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
        data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
        data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
  }
/*
 for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
       gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]) ;
       gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]) ;
  }
*/
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
} //end substract func



//1D Rebin function:Basically it cut off some bin you  don't needed
TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar ){
TH1D* thout;
cout << "Eventshape histogram rebined" <<endl;
int nxmod2 = 0;
double xmod2[300];
int numHT =8;
int numvar =5;
int numtype =2;

const char *namex;
const char *titlex;
char namey[100], titley[100];

/*
// define the number of bin excluded form lower end
int arrayvar_first[numtype][numvar][numHT] = {{
	                    {24,24,24,24,24,24,24,24},
	                    {24,24,24,24,24,24,24,24},  // thrustc
                            {52,42,28,24,24,24,24,24}, //t3mass
                            {64,62,60,50,50,50,44,24}, //broadt
                            {24,24,24,24,24,24,24,24}},
                            {{50,50,46,44,44,42,46,46},
		            {22,22,22,22,22,22,22,22},  // thrustc
                            {22,22,22,22,22,22,22,22}, //t3mass
                            {26,24,24,22,22,22,46,52}, //broadt
                            {22,22,22,22,22,22,22,22}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
	                       {2,2,2,2,2,2,2,2},  // thrustc
	                       {2,2,2,2,2,2,2,2},  // thrustc
                               {10,10,6,6,2,2,2,2}, //t3mass
                               {0,0,0,0,0,0,0,0}, //broadt
                               {0,0,0,0,0,0,0,0}},
                               {{2,2,2,2,2,2,2,2},  // thrustc
                               {0,0,0,0,0,0,0,0},  // thrustc
                               {16,8,14,14,16,14,12,12}, //t3mass
                               {0,0,0,0,0,0,0,0}, //broadt
                               {0,0,0,0,0,0,0,0}}}; //ttmass


*/
/*
// define the number of bin excluded form lower end
int arrayvar_first[numtype][numvar][numHT] = {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}},
                            {{0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}},
                            {{0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass
*/


// define the number of bin excluded form lower end
int arrayvar_first[numtype][numvar][numHT] = {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {4,2,1,0,0,0,0,0},
                            {8,8,7,5,5,5,4,0},
                            {0,0,0,0,0,0,0,0}},
                            {{2,1,1,0,0,0,0,1},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {1,1,1,0,0,0,4,6},
                            {0,0,0,0,0,0,0,0}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,2,2,4},
                            {2,3,1,0,0,0,0,0},
                            {0,0,1,1,1,1,0,2},
                            {0,0,0,0,0,0,0,0}},
                            {{2,1,1,0,0,0,0,1},
                            {2,2,2,0,0,0,2,2},
                            {0,0,0,0,0,0,0,0},
                            {1,1,1,2,0,2,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass



int arrayind[numvar] ={3,9,15,18,24};

  double yvl[200]={0.0};
  double erryvl[200]={0.0};
  int nbinx = thin->GetNbinsX();
  int indx=-1;
  /*for (int ij=0; ij<numvar; ij++) {
    if (ivar==arrayind[ij]) { indx=ij; break;}
  }
  if (indx<0)  cout << " No Matching Variables" << endl;*/
  indx=ivar;
  int ifirst = 0;
  nxmod2= -1;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][indx][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][indx][ijetpt])) {
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
for(int ii=0; ii < nxmod2+1 ; ii++){
      xmod3[ii]=xmod2[ii];        }

  namex = thin->GetName();
  sprintf(namey, "rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "Rebin_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod3);

  for (int ix=0; ix<nxmod2+1; ix++) {
    thout->SetBinContent(ix+1, yvl[ix]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix]));
     }

  cout << "final Bin :"<< thout->GetNbinsX()<< endl;
  return thout;

}   //end of rebin_hist functiomC

//2D Rebin function:Basically it cut off some bin you  don't needed
TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco, int itype, int ijetpt, int ivar){
TH2D* thout;
cout << "Eventshape histogram rebined 2d" <<endl;
int nxmod2 = 0;
double xmod2[200];

int numHT =8;
int numvar =5;
int numtype =2;

const char* namex;
char namey[100],titley[100];
const char* titlex;

// define the number of bin excluded form lower end
// define the number of bin excluded form lower end

/*
int arrayvar_first[numtype][numvar][numHT] = {{
                            {24,24,24,24,24,24,24,24},
                            {24,24,24,24,24,24,24,24},  // thrustc
                            {52,42,28,24,24,24,24,24}, //t3mass
                            {64,62,60,50,50,50,44,24}, //broadt
                            {24,24,24,24,24,24,24,24}},
                            {{50,50,46,44,44,42,46,46},
                            {22,22,22,22,22,22,22,22},  // thrustc
                            {22,22,22,22,22,22,22,22}, //t3mass
                            {26,24,24,22,22,22,46,52}, //broadt
                            {22,22,22,22,22,22,22,22}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
                               {2,2,2,2,2,2,2,2},  // thrustc
                               {2,2,2,2,2,2,2,2},  // thrustc
                               {10,10,6,6,2,2,2,2}, //t3mass
                               {0,0,0,0,0,0,0,0}, //broadt
                               {0,0,0,0,0,0,0,0}},
                               {{2,2,2,2,2,2,2,2},  // thrustc
                               {0,0,0,0,0,0,0,0},  // thrustc
                               {16,8,14,14,16,14,12,12}, //t3mass
                               {0,0,0,0,0,0,0,0}, //broadt
                               {0,0,0,0,0,0,0,0}}}; //ttmass

*/
/*
// define the number of bin excluded form lower end
int arrayvar_first[numtype][numvar][numHT] = {{
                            {2,2,2,2,2,2,2,2},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}},
                            {{0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
                            {1,1,1,1,1,1,1,1},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}},
                            {{0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass
*/
// define the number of bin excluded form lower end
int arrayvar_first[numtype][numvar][numHT] = {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {4,2,1,0,0,0,0,0},
                            {8,8,7,5,5,5,4,0},
                            {0,0,0,0,0,0,0,0}},
                            {{2,1,1,0,0,0,0,1},
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0},
                            {1,1,1,0,0,0,4,6},
                            {0,0,0,0,0,0,0,0}}}; //ttmass

// define the number of bin excluded from higher end
int arrayvar_last[numtype][numvar][numHT] =  {{
                            {0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,2,2,4},
                            {2,3,1,0,0,0,0,0},
                            {0,0,1,1,1,1,0,2},
                            {0,0,0,0,0,0,0,0}},
                            {{2,1,1,0,0,0,0,1},
                            {2,2,2,0,0,0,2,2},
                            {0,0,0,0,0,0,0,0},
                            {1,1,1,2,0,2,0,0},
                            {0,0,0,0,0,0,0,0}}}; //ttmass


int arrayind[numvar] ={3,9,15,18,24};

  double yvl[200][200]={0.0};
  double erryvl[200][200]={0.0};
  int nbinx = thin->GetNbinsX();
  int indx=-1;
 /*for (int ij=0; ij<numvar; ij++) {
    if (ivar==arrayind[ij]) { indx=ij; break;}
  }
  if (indx<0)  cout << " No Matching Variables" << endl;*/
  indx =ivar;
  int ifirst = 0;
  int jfirst =0;
  nxmod2= -1;

  int kfirst =0;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][indx][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][indx][ijetpt])) {
      } else {
        xmod2[kfirst] = MC_reco->GetBinLowEdge(ij);
        xmod2[kfirst+1] = MC_reco->GetBinLowEdge(ij+1);
        kfirst++;
//        xbin[][]

      }
    }
  }

  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][indx][ijetpt] ) {
      } else if ( ij > (nbinx - arrayvar_last[itype][indx][ijetpt]) ) {
      } else {
	jfirst=0;
    for (int jk=0; jk<nbinx+2; jk++) {
      if (jk <=arrayvar_first[itype][indx][ijetpt]) {
      } else if (jk > (nbinx - arrayvar_last[itype][indx][ijetpt])) {
      } else {
  //      xmod2[ifirst] = MC_reco->GetBinLowEdge(ij);
   //     xmod2[ifirst+1] = MC_reco->GetBinLowEdge(ij);
        yvl[ifirst][jfirst] = thin->GetBinContent(ij,jk);
        erryvl[ifirst][jfirst] += thin->GetBinError(ij,jk)*thin->GetBinError(ij,jk);
        jfirst++;
        }
      }
      ifirst++;
    }
  }
    nxmod2 = ifirst;
  }
double xmod3[nxmod2+1];
for(int ii=0; ii < nxmod2+1 ; ii++){
      xmod3[ii]=xmod2[ii];        }


  namex = thin->GetName();
  sprintf(namey, "Rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "%s", titlex);
  thout = new TH2D(namey, titley, nxmod2, xmod3,nxmod2, xmod3);

  for (int ix=0; ix<nxmod2+1; ix++) {
  for (int iy=0; iy<nxmod2+1; iy++) {
    thout->SetBinContent(ix+1,iy+1, yvl[ix][iy]);
    thout->SetBinError(ix+1,iy+1, sqrt(erryvl[ix][iy]));
        }
     }

//  cout << "final Bin :"<< thout->GetNbinsX()<< endl;
  return thout;

}   //end of 2d rebin_hist functiomC

void Fold(TH2D* HistoMatrix, TH1D* HistoGen, TH1D* HistoUnfold, TH1D* HistoCorrect){

  double probFake=0.;
  double AwayMatrixTot=0.;

   for(int i=1;i<HistoGen->GetXaxis()->GetNbins()+1;i++){
    double Inside=0.;
    double ErrorP=0.;
    double sum = 0.;
    for(int j=1;j<HistoGen->GetXaxis()->GetNbins()+1;j++){
      probFake=0;
      AwayMatrixTot=HistoMatrix->GetBinContent(i,j)/HistoGen->GetBinContent(j);
      if(HistoGen->GetBinContent(j)==0) AwayMatrixTot=1;
      Inside=Inside+HistoUnfold->GetBinContent(j)*AwayMatrixTot;//*eff[j];///(1-probFake);
      ErrorP=ErrorP+HistoUnfold->GetBinError(j)*AwayMatrixTot/(1-probFake)*HistoUnfold->GetBinError(j)*AwayMatrixTot;
      sum+= +AwayMatrixTot;
    }
    HistoCorrect->SetBinContent(i,Inside);
    double Error=sqrt(ErrorP);
    //double Error=ErrorP;
    HistoCorrect->SetBinError(i,Error);
  }
}

void ReFold(TH2D* HistoMatrix,  TH1D* HistoGen, TH1D* HistoUnfold, TH1D* HistoCorrect){

for(int ij=0; ij<(HistoMatrix->GetNbinsX()); ij++){
double row_sum = 0.;
 for(int jk=0; jk<(HistoMatrix->GetNbinsX()); jk++){
  row_sum+=HistoMatrix->GetBinContent(jk+1,ij+1);
}//jk
if(row_sum>1.e-10){
 for(int jk=0; jk<(HistoMatrix->GetNbinsX()); jk++){
   HistoMatrix->SetBinContent(jk+1,ij+1,HistoMatrix->GetBinContent(jk+1,ij+1)*1./row_sum) ;
  }//jk
 }
}//ij

for(int ij=0; ij<(HistoMatrix->GetNbinsX()); ij++){
 double sum =0.;
 double error = 0.;
 for(int jk=0; jk<(HistoMatrix->GetNbinsX()); jk++){
  sum+= HistoMatrix->GetBinContent(ij+1,jk+1) * HistoUnfold->GetBinContent(jk+1) ;
  error+=HistoMatrix->GetBinContent(ij+1,jk+1) * HistoUnfold->GetBinError(jk+1);
 }//jk
   HistoCorrect->SetBinContent(ij+1,sum);
    HistoCorrect->SetBinError(ij+1,error);
}//ij

}

void BackFold_Gen(RooUnfoldResponse response_b, TH1D* hist_gen, TH1D* hist_back,double *eff,double* fakerate){
for(int am=0; am<(hist_gen->GetNbinsX()); am++){
hist_gen->SetBinContent(am+1,/*eff[am]*/hist_gen->GetBinContent(am+1)) ;
}
TMatrixD rm_mat_bfold = response_b.Mresponse();
for(int bn=0; bn<(hist_gen->GetNbinsX()); bn++){
  double sum = 0.;
  double err_sq=0.;
  for(int am=0; am<(hist_gen->GetNbinsX()); am++){
   sum+=rm_mat_bfold[bn][am]*hist_gen->GetBinContent(am+1);
   err_sq+=pow(rm_mat_bfold[bn][am]*hist_gen->GetBinError(am+1),2.);
  }//am
    err_sq = sqrt(err_sq);
    hist_back->SetBinContent(bn+1,sum*1./(1-fakerate[bn]));
    hist_back->SetBinError(bn+1,err_sq*1./sqrt(1-fakerate[bn]));
 }//bn
}

void BackFold_Bayesian(RooUnfoldResponse response_b, RooUnfoldBayes unfold, TH1D* hist_back,double *eff,double* fakerate){

TMatrixD rm_mat_bfold = response_b.Mresponse();
TVectorD bfold_proxy = unfold.Vreco();
TVectorD bfold_err_proxy = unfold.ErecoV(RooUnfold::kCovariance);
/*
for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
bfold_proxy[bn]*=eff[bn];
bfold_err_proxy[bn]*=eff[bn];
}
*/
TVectorD bfold = rm_mat_bfold * bfold_proxy;
TVectorD bfold_err = rm_mat_bfold * bfold_err_proxy;

for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
 hist_back->SetBinContent(bn+1,bfold[bn]*1./(1-fakerate[bn])) ;
 hist_back->SetBinError(bn+1,bfold_err[bn]*1./sqrt(1-fakerate[bn])) ;
}

double sum=0;
for(int ij=0; ij < (hist_back->GetNbinsX()); ij++){
sum+=rm_mat_bfold[5][ij];
}
}

void BackFold_Svd(RooUnfoldResponse response_b, RooUnfoldSvd unfold, TH1D* hist_back,double *eff,double* fakerate){

TMatrixD rm_mat_bfold = response_b.Mresponse();
TVectorD bfold_proxy = unfold.Vreco();
TVectorD bfold_err_proxy = unfold.ErecoV();
/*
for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
bfold_proxy[bn]*=eff[bn];
bfold_err_proxy[bn]*=eff[bn];
}
*/
TVectorD bfold = response_b.Mresponse() * bfold_proxy;//unfold.Vreco();
TVectorD bfold_err = rm_mat_bfold * bfold_err_proxy;

for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
  hist_back->SetBinContent(bn+1,bfold[bn]*1./(1-fakerate[bn])) ;
  hist_back->SetBinError(bn+1,bfold_err[bn]*1./sqrt(1-fakerate[bn])) ;
 }
}

void BackFold_BinbyBin(RooUnfoldResponse response_b, RooUnfoldBinByBin unfold, TH1D* hist_back,double *eff,double* fakerate){

TMatrixD rm_mat_bfold = response_b.Mresponse();
TVectorD bfold_proxy = unfold.Vreco();
TVectorD bfold_err_proxy = unfold.ErecoV();
/*
for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
bfold_proxy[bn]*=eff[bn];
bfold_err_proxy[bn]*=eff[bn];
}
*/
TVectorD bfold = response_b.Mresponse() * bfold_proxy; //unfold.Vreco();
TVectorD bfold_err = rm_mat_bfold * bfold_err_proxy;

for(int bn=0; bn < (hist_back->GetNbinsX()); bn++){
  hist_back->SetBinContent(bn+1,bfold[bn]*1./(1-fakerate[bn])) ;
  hist_back->SetBinError(bn+1,bfold_err[bn]*1./sqrt(1-fakerate[bn])) ;
 }

}

int LastBin_Counter(TH1D *h1){
int lastbin = 0;
for(int bn=0; bn<(h1->GetNbinsX()); bn++){
if((h1->GetBinContent(bn+1)) > 1.e-12){
lastbin = bn+1;
}
else { break;}
}
return lastbin;
}


#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
