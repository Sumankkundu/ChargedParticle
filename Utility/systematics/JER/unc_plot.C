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

#define JER
using namespace std;

void unc_plot(){

	int const njer = 2;
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

  bool isstat =1;  int irbin=1;
  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={2,4,6,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[nusedvar]={3,9,15,18,24};
  Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000};

  const  char* itypeN[ntype]={"Jets","Charged Particles"};
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[8]={"83 < H_{T,2} < 109", "109 < H_{T,2} < 172", "172 < H_{T,2} < 241", "241 < H_{T,2} < 309","309 < H_{T,2} < 377", "377 < H_{T,2} < 462", "462 < H_{T,2} <570","H_{T,2} > 570"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};

  TH1D* ReadHist1D(string name,TFile* root, int irbin=1);
  TH2D* ReadHist2D(string name,TFile* root, int irbin=1);
  TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4="");
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl);


  TH1D *JEC_Unc[ntype][nusedvar][nHLTmx];
  TH1D *JER_Unc[ntype][nusedvar][nHLTmx];
  TH1D *PDF_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Unfold_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Unfold_Unc[ntype][nusedvar][nHLTmx];

  TFile *JEC_root = TFile::Open("JEC_Result.root");
  TFile *JER_root = TFile::Open("JER_Result.root");
  TFile *PDF_root = TFile::Open("PDF_Result.root");


  for(int ity=0; ity <ntype; ity++){
      for(int ivar =0;ivar < nusedvar; ivar++){
          for(int ipt=0; ipt < nHLTmx; ipt++){
          JEC_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("jec_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),JEC_root);
          JER_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("jer_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),JER_root);
          PDF_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("jer_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),JER_root);
      
 	  }
      
      }
    }





  //----------------------------------------------------------
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "JEC Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *jec_unc = (TH1D*)JEC_Unc[ity][ivar][ipt]->Clone();
     jec_unc->Draw()

     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "Unc_plot.pdf("); sprintf(pdfname1, "Unc_plot.pdf");sprintf(pdfname2, "Unc_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
     cptr->Clear();
        }
      }
   }






}

TH1D* ReadHist1D(string name, TFile* root, int irbin=1){
TString histname = name; cout << histname<< endl;
TH1D* hist=(TH1D*)root->Get(histname);
hist->Rebin(irbin);
//hist->Write();
return hist;
}

TH2D* ReadHist2D(string name, TFile* root, int irbin=1){
TString histname = name; cout << histname<<endl;
TH2D* hist=(TH2D*)root->Get(histname);
hist->RebinY(irbin);
//hist->Write();
return hist;
}

