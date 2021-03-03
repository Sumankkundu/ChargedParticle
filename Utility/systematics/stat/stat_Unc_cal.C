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
#include <TProfile.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>

using namespace std;

double ComputeStdDev (int n, double *x, double cv){
double sum=0;
for(int ij=2; ij<n; ij++){
sum +=pow(x[ij]-cv, 2.0); }
sum /=n-1;
return sqrt(sum);    }

void stat_Unc_cal(){
  int const nPDF = 100;
  int const nmc=4; //Number of MC -> 0,1,2 PY8, MG, HW7
  int const umc=0; //Which MC will used for Unfold

  int irbin = 1; //Rebin
  int const untype = 1;  // 0 for 2D , 1 for 2D 
  
  char histname[100], name[100], title[100], Axisname[100];
  const int nHLTmx=8; //HT2 Range
  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =2;
  Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass
  static const int nhist=10; //We need 8 But define 10
  const int njetptmn = nHLTmx;
  string NominalDir = "Unfold2D";
  string PDFdir = "MG8";

  ofstream file_out("jerunc_eta0_test.inc");

 // double mean[100]; double poserr[100];  double negerr[100];
  TH1F* fdfz;  TH1F* fdfx;  TH1F* fdfy;

//--------------------------Functions
  void setgstyle();
  void Integralhist(TH1 *hist);
  void Normalise(TH1* h, TH2* covmat);
  void Extract(TH1* global2d, TUnfoldBinning* Bin , char* Axisname, bool iw=0);
  TH1D* ReadHist1D(string name,TFile* root);
  TH2D* ReadHist2D(string name,TFile* root);
  void HT2_NormalV3(TH1* global2d, TUnfoldBinning* Bin, char* Axisname, int nht,  bool iw=0);
//-----------------------------------------

  TH1D *Nominal[ntype][nusedvar];                    //Nominal Data Unfolded 
  TH1D *Nominal_HT[ntype][nusedvar][nHLTmx];                    //Nominal Data Unfolded 
  TH1D *stat_err[ntype][nusedvar][nHLTmx];                    //Nominal Data Unfolded 
  TH1D *PDF_gen[nPDF][ntype][nusedvar];        //JEC Data Unfolded
  TH1D *PDF_gen_HT[nPDF][ntype][nusedvar][nHLTmx];        //JEC Data Unfolded

  TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/JEC_Tunfold_19Jan21/Unfolded_Result.root");  // Unfolded data 
  
  TFile *stat_Result =new TFile("stat_Result.root","recreate");   

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
      Nominal[ity][ivar] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unfoldroot)->Clone();
      //HT2_NormalV3(Data_Reco[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      Nominal_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
      stat_err[ity][ivar][ipt]= (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();       
      stat_err[ity][ivar][ipt]->Reset();
      stat_err[ity][ivar][ipt]->SetNameTitle(Form("stat_erro_%i_eta0_%i_pt%i" ,ity, var[ivar], ipt),Form("PDF Relative %i_eta0_%i_pt%i", ity, var[ivar], ipt));      
      }
       
       }
    }

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       double meanerr[100]={0};
       for(int ipt=0; ipt < nHLTmx; ipt++){
  
       TH1D* default_hist = (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();
     
       cout <<" ipt: " << ipt << " ivar: " << ivar << " ity: " << ity << endl;
       double mean[100]={0};


       int nbins = default_hist->GetNbinsX();
       double histmean = default_hist->GetMean(); 
       for (int ix=1; ix<nbins+1; ix++) { mean[ix] = default_hist->GetBinContent(ix);
       }
       int xx =0; 
       for (int ix=1; ix<nbins+1; ix++){
        xx =0;
       TH1D* hist = (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();
       double val1 = hist->GetBinContent(ix);
       double error = hist->GetBinError(ix);
 
       stat_err[ity][ivar][ipt]->SetBinContent( ix, abs(error)/val1); 
        cout << "Relative error :" << abs(error)/val1 ;
        } //nbinsx
       stat_err[ity][ivar][ipt]->Write();
     }//ipt
  }//ivar
}//ity
//------------------------
delete stat_Result;
}

//Functions :
TH1D* ReadHist1D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH1D* hist=(TH1D*)root->Get(histname);
//hist->Write();
return hist;
}
//-----------------------------------------------------------
TH2D* ReadHist2D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH2D* hist=(TH2D*)root->Get(histname);
//hist->Write();
return hist;
}


