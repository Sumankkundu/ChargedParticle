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

void PDF_Unc_cal(){
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
  TH1D *PDFerr_HT[ntype][nusedvar][nHLTmx];                    //Nominal Data Unfolded 
  TH1D *PDF_gen[nPDF][ntype][nusedvar];        //JEC Data Unfolded
  TH1D *PDF_gen_HT[nPDF][ntype][nusedvar][nHLTmx];        //JEC Data Unfolded

  TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/JEC_Tunfold_19Jan21/Unfolded_Result.root");  // Unfolded data 
  TFile *PFD_file = TFile::Open("Unfolded_Result.root");  // Unfolded data 
  
  TFile *PDF_Result =new TFile("PDF_Result.root","recreate");   

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
      Nominal[ity][ivar] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unfoldroot)->Clone();
      //HT2_NormalV3(Data_Reco[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      Nominal_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
      PDFerr_HT[ity][ivar][ipt]= (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();       
      PDFerr_HT[ity][ivar][ipt]->Reset();
      PDFerr_HT[ity][ivar][ipt]->SetNameTitle(Form("PDF_erro_%i_eta0_%i_pt%i" ,ity, var[ivar], ipt),Form("PDF Relative %i_eta0_%i_pt%i", ity, var[ivar], ipt));      
      }
      for(int ipdf=0; ipdf <nPDF; ipdf++){
      cout << ipdf<< "  :  " ;
//    PDF_gen[ipdf][ity][ivar]= (TH1D*)ReadHist1D(PDFdir+"/Edd_genpdf_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", PFD_file)->Clone();
        
      for(int ipt=0; ipt < nHLTmx; ipt++){
      PDF_gen_HT[ipdf][ity][ivar][ipt]= (TH1D*)ReadHist1D(PDFdir+"/Edd_genpdf_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_"+to_string(ipdf+1)+"_pt"+to_string(ipt), PFD_file)->Clone();
             }
      //HT2_NormalV3(Data_Reco_JEC[ity][ivar][ijec], binsRec[ity][ivar], Axisname,nHLTmx);   
	  }
       }
    }

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       double meanerr[100]={0};
       for(int ipt=0; ipt < nHLTmx; ipt++){
  
       TH1D* default_hist = (TH1D*)PDF_gen_HT[0][ity][ivar][ipt]->Clone();
     
       cout <<" ipt: " << ipt << " ivar: " << ivar << " ity: " << ity << endl;
       double mean[100]={0};


       int nbins = default_hist->GetNbinsX();
       double histmean = default_hist->GetMean(); 
       for (int ix=1; ix<nbins+1; ix++) { mean[ix] = default_hist->GetBinContent(ix);
       cout <<" bin conternt " <<  mean[ix];
       }
       int xx =0; 
       for (int ix=1; ix<nbins+1; ix++){
        xx =0;
       double pdfbincontent[nPDF];
       for(int ipdf=1; ipdf < nPDF; ipdf++){
       //cout << " ipdf : " << ipdf ;
       TH1D* pdf_hist = (TH1D*)PDF_gen_HT[ipdf][ity][ivar][ipt]->Clone();
       double val1 = pdf_hist->GetBinContent(ix);
//       cout << " val1 = " << val1 ;
       if(mean[ix]!= 0) pdfbincontent[ipdf] = val1/mean[ix];
       //if(mean[ix]!= 0) pdfbincontent[ipdf] = val1;
       
             xx += pow(abs(val1-mean[ix]),2.0);
              //cout << "xx = " << sqrt(xx/100) << endl;
            }//nPDf
             xx /=100;
              cout << "xx = " << xx << "  sqrt :  "<< sqrt(abs(xx))<<" : Relative:  " <<sqrt(abs(xx))/mean[ix] << endl;

       PDFerr_HT[ity][ivar][ipt]->SetBinContent( ix, sqrt(abs(xx))/mean[ix]); 

      
     double xgError = ComputeStdDev (nPDF, pdfbincontent, 1);
     cout << " xgError = " << xgError << endl;
//     cout << " Rel Error = " << mean[ix] << "  " << xgError/mean[ix] << endl;
        } //nbinsx
       PDFerr_HT[ity][ivar][ipt]->Write();
     }//ipt
  }//ivar
}//ity
//------------------------
delete PDF_Result;
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


