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


#define NOREG
//#define TIKH_SURE
//#define TIKH_LSCAN
//#define ITER_SURE
//#define ITER
#define REBIN

using namespace std;
void testUnfold2c()
{
// switch on histogram errors
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();

  //Input Data and MC histogram 
  //TFile *inputData=new TFile("PY8_UL2017_120Bin_14june20.root");
 // TFile *inputData=new TFile("PY8_UL17_ULbinned_1July20.root");
  //  TFile *inputData=new TFile("MG_UL17_binned_6July20.root");
    TFile *inputData=new TFile("DATA_JetHT_UL17_ULbin_4July.root");
  //TFile *inputData=new TFile("PY8_UL17_binnedUL_1July20_NOweight.root");
  //TFile *inputData=new TFile("HW7_Flat_Binned_8July20.root");

// TFile *inputMC=new TFile("PY8_CP5_UL_Binned_20May.root");
// TFile *inputMC=new TFile("PY8_UL2017_120Bin_14june20.root");
    TFile *inputMC=new TFile("PY8_UL17_ULbinned_1July20.root");
 // TFile *inputMC=new TFile("MG_UL17_binned_6July20.root");
 //TFile *inputMC=new TFile("PY8_UL17_binnedUL_1July20_NOweight.root");
  //TFile *inputMC=new TFile("HW7_Flat_Binned_8July20.root");
 

  TFile *inputMC1=new TFile("MG_UL17_binned_6July20.root");
 // TFile *inputMC1=new TFile("PY8_UL17_ULbinned_1July20.root");
  TFile *inputMC2=new TFile("HW7_Flat_Binned_8July20.root");

  //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  TFile *outputFile=new TFile("testunfold2c_unfolded.root","recreate");
  
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
  
  ///////////////////Modify /////////////////////////////
   
   TH1D *MC_Reco[type][nusedvar][njetptmn];  //Reconstracted MC 
   TH1D *MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco
   
   TH1D *MC_Reco1[type][nusedvar][njetptmn];  //Reconstracted MC copy
   TH1D *MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy
   TH1D *Data_Reco1[type][nusedvar][njetptmn];    //Reconstracted Data copy
   TH2D *h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy


   TH1D *MC_Reco2[type][nusedvar][njetptmn];  //Reconstracted MC copy 
   TH1D *MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
   TH1D *Data_Reco2[type][nusedvar][njetptmn];    //Reconstracted Data copy
   TH2D *h2dGenDetMC2[type][nusedvar][njetptmn];
   
   TH1D *MC_Reco_re[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MC_Gen_re[type][nusedvar][njetptmn];   //Generator level MC
   TH1D *Data_Reco_re[type][nusedvar][njetptmn];    //Reconstracted Data
   TH2D *h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco

  //---------------------MG5
  TH1D *MG5_MC_Reco[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MG5_MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
   TH2D *MG5_h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco

   TH1D *MG5_MC_Reco1[type][nusedvar][njetptmn];  //Reconstracted MC copy
   TH1D *MG5_MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy
   TH2D *MG5_h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy


   TH1D *MG5_MC_Reco2[type][nusedvar][njetptmn];  //Reconstracted MC copy
   TH1D *MG5_MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
   TH2D *MG5_h2dGenDetMC2[type][nusedvar][njetptmn];

   TH1D *MG5_MC_Reco_re[type][nusedvar][njetptmn];  //Reconstracted MC
   TH1D *MG5_MC_Gen_re[type][nusedvar][njetptmn];   //Generator level MC
   TH2D *MG5_h2dGenDetMC_re[type][nusedvar][njetptmn];   // MC generator Vs Reco
//----------------------HW
   TH1D *HW7_MC_Reco[type][nusedvar][njetptmn];  //Reconstracted MC 
   TH1D *HW7_MC_Gen[type][nusedvar][njetptmn];   //Generator level MC
   TH2D *HW7_h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco

   TH1D *HW7_MC_Reco1[type][nusedvar][njetptmn];  //Reconstracted MC copy
   TH1D *HW7_MC_Gen1[type][nusedvar][njetptmn];   //Generator level MC copy
   TH2D *HW7_h2dGenDetMC1[type][nusedvar][njetptmn];   // MC generator Vs Reco copy


   TH1D *HW7_MC_Reco2[type][nusedvar][njetptmn];  //Reconstracted MC copy 
   TH1D *HW7_MC_Gen2[type][nusedvar][njetptmn];   //Generator level MC copy
   TH2D *HW7_h2dGenDetMC2[type][nusedvar][njetptmn];

   TH1D *HW7_MC_Reco_re[type][nusedvar][njetptmn];  //Reconstracted MC
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
   
   TDirectoryFile *inputDir1=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
   TDirectoryFile *inputDir2=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
   TDirectoryFile *inputDir3=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
   
   TDirectoryFile *inputRebin0=new TDirectoryFile("RebinData","Inputs Data Rebin");
   TDirectoryFile *inputRebin1=new TDirectoryFile("RebinPythia8","Pythia8 MC Rebin");
   TDirectoryFile *inputRebin2=new TDirectoryFile("RebinMG8","Madgraph, MC Rebin");
   TDirectoryFile *inputRebin3=new TDirectoryFile("RebinHW7","Herwig7 MC Rebin");
   
   TDirectoryFile *Unfold=new TDirectoryFile("Unfold","Unfolded, Refold, correlation");

   int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
   TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar);
   TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco, int itype, int ijetpt, int ivar);


//Read Input Data MC and Response matrix 
   for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){

         //Reco Data
         inputDir0->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         Data_Reco1[ity][ivar][ipt] = (TH1D*) inputData->Get(histname);
 //      Data_Reco1[ity][ivar][ipt]->Rebin(2);
         Data_Reco1[ity][ivar][ipt]->Write();
       for (int ibin =1 ; ibin <  Data_Reco1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
         if (Data_Reco1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " Data Reco Bin is Zero for bin number : ******** "<<  ibin  << endl; }}


	 inputDir1->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
         cout << histname << endl;
         int recobins = MC_Reco1[ity][ivar][ipt]->GetNbinsX();
         rnbinsx[ity][ivar][ipt]=MC_Reco1[ity][ivar][ipt]->GetNbinsX(); 
         if(recobins !=Data_Reco1[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
	 

	 MC_Reco1[ity][ivar][ipt]->Write();
         for(int ibin = 1; ibin < MC_Reco1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
         if (MC_Reco1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " MC reco Bin is Zero for bin number :** "<<  ibin  << "  Bin lowedge "<< MC_Reco1[ity][ivar][ipt]->GetBinLowEdge(ibin)<< endl; }}

         //Gen MC
         sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC->Get(histname);
//         MC_Gen1[ity][ivar][ipt]->Rebin(2);
         MC_Gen1[ity][ivar][ipt]->Write();
         int genbins = MC_Gen1[ity][ivar][ipt]->GetNbinsX();
         gnbinsx[ity][ivar][ipt]=MC_Gen1[ity][ivar][ipt]->GetNbinsX();
	 for(int ibin = 1; ibin < MC_Gen1[ity][ivar][ipt]->GetNbinsX()+1; ibin++ ){
         if (MC_Gen1[ity][ivar][ipt]->GetBinContent(ibin) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  ibin  << endl; }}


         //Response Matrix
         sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
         h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC->Get(histname);   //Xgen(coarse) , Yreco(fine)
//         h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
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
   // 	Data_Reco_re[ity][ivar][ipt]->Rebin(4);
    	cout << "Rebin Data Reco Bin : " << Data_Reco_re[ity][ivar][ipt]->GetNbinsX() << endl; 
 	
	MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
 //	MC_Reco_re[ity][ivar][ipt]->Rebin(4);
  	
	MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
 //	MC_Gen_re[ity][ivar][ipt]->Rebin(2);
  	cout << " Rebin MC Gen bin : " <<  MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;
  	
	h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(h2dGenDetMC1[ity][ivar][ipt], MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
//	h2dGenDetMC_re[ity][ivar][ipt] = h2dGenDetMC1[ity][ivar][ipt];
  //	h2dGenDetMC_re[ity][ivar][ipt]->RebinY(2);
//  	h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
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
	 inputRebin0->cd();
	 Data_Reco[ity][ivar][ipt]->Write();
         
	 inputRebin1->cd();
	 MC_Reco[ity][ivar][ipt]->Write();
         MC_Gen[ity][ivar][ipt]->Write();
         h2dGenDetMC[ity][ivar][ipt]->Write();


       }
     }
   }
   cout << "Histogram Read and Rebin oK " <<endl;      
//------------------Read MG
   for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){

         inputDir2->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MG5_MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC1->Get(histname);
         cout << histname << endl;
         MG5_MC_Reco1[ity][ivar][ipt]->Write();

         //Gen MC
         sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         MG5_MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC1->Get(histname);
//       MC_Gen1[ity][ivar][ipt]->Rebin(2);
         MG5_MC_Gen1[ity][ivar][ipt]->Write();

         //Response Matrix
         sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
         MG5_h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC1->Get(histname);   //Xgen(coarse) , Yreco(fine)
//       MG5_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
         MG5_h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)

         MG5_MC_Reco2[ity][ivar][ipt] = (TH1D*)MG5_MC_Reco1[ity][ivar][ipt]->Clone();
         MG5_MC_Gen2[ity][ivar][ipt] = (TH1D*)MG5_MC_Gen1[ity][ivar][ipt]->Clone();
         MG5_h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)MG5_h2dGenDetMC1[ity][ivar][ipt]->Clone();


#ifdef REBIN
        MG5_MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(MG5_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
 //     MG5_MC_Reco_re[ity][ivar][ipt]->Rebin(4);

        MG5_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(MG5_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
 //     MG5_MC_Gen_re[ity][ivar][ipt]->Rebin(2);
        cout << " Rebin MC Gen bin : " <<  MG5_MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;

        MG5_h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(MG5_h2dGenDetMC1[ity][ivar][ipt], MG5_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
//      MG5_h2dGenDetMC_re[ity][ivar][ipt] = h2dGenDetMC1[ity][ivar][ipt];
  //    MG5_h2dGenDetMC_re[ity][ivar][ipt]->RebinY(2);
//      MG5_h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
        cout << "ReBin REco Bin :"<< MG5_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << MG5_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;
//Define data MC to be used for unfolding
         MG5_MC_Reco[ity][ivar][ipt] = (TH1D*)MG5_MC_Reco_re[ity][ivar][ipt]->Clone();
         MG5_MC_Gen[ity][ivar][ipt] = (TH1D*)MG5_MC_Gen_re[ity][ivar][ipt]->Clone();
         MG5_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)MG5_h2dGenDetMC_re[ity][ivar][ipt]->Clone();


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
  
//-----------HW7
   for(int ity=0; ity <type; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){

         inputDir3->cd();
         sprintf(histname, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         HW7_MC_Reco1[ity][ivar][ipt] = (TH1D*) inputMC2->Get(histname);
         cout << histname << endl;
         HW7_MC_Reco1[ity][ivar][ipt]->Write();

         //Gen MC
         sprintf(histname, "analyzeBasicPat/gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
         HW7_MC_Gen1[ity][ivar][ipt] = (TH1D*) inputMC2->Get(histname);
//       HW7_MC_Gen1[ity][ivar][ipt]->Rebin(2);
         HW7_MC_Gen1[ity][ivar][ipt]->Write();

         //Response Matrix
         sprintf(histname, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
         HW7_h2dGenDetMC1[ity][ivar][ipt] = (TH2D*) inputMC2->Get(histname);   //Xgen(coarse) , Yreco(fine)
	 //       HW7_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
	 HW7_h2dGenDetMC1[ity][ivar][ipt]->Write();   //Xgen(coarse) , Yreco(fine)

         HW7_MC_Reco2[ity][ivar][ipt] = (TH1D*)HW7_MC_Reco1[ity][ivar][ipt]->Clone();
         HW7_MC_Gen2[ity][ivar][ipt] = (TH1D*)HW7_MC_Gen1[ity][ivar][ipt]->Clone();
         HW7_h2dGenDetMC2[ity][ivar][ipt] = (TH2D*)HW7_h2dGenDetMC1[ity][ivar][ipt]->Clone();


#ifdef REBIN
        HW7_MC_Reco_re[ity][ivar][ipt] = rebin1d_hist(HW7_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
 //     HW7_MC_Reco_re[ity][ivar][ipt]->Rebin(4);

        HW7_MC_Gen_re[ity][ivar][ipt] = rebin1d_hist(HW7_MC_Gen1[ity][ivar][ipt], ity, ipt, ivar);
 //     HW7_MC_Gen_re[ity][ivar][ipt]->Rebin(2);
        cout << " Rebin MC Gen bin : " <<  HW7_MC_Gen_re[ity][ivar][ipt]->GetNbinsX() << endl;

        HW7_h2dGenDetMC_re[ity][ivar][ipt] = rebin2d_hist(HW7_h2dGenDetMC1[ity][ivar][ipt], HW7_MC_Reco1[ity][ivar][ipt], ity, ipt, ivar);
//      HW7_h2dGenDetMC_re[ity][ivar][ipt] = h2dGenDetMC1[ity][ivar][ipt];
  //    HW7_h2dGenDetMC_re[ity][ivar][ipt]->RebinY(2);
//      HW7_h2dGenDetMC_re[ity][ivar][ipt]->RebinX(4);
        cout << "ReBin REco Bin :"<< HW7_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsX() <<" Rebin Gen Bin : " << HW7_h2dGenDetMC_re[ity][ivar][ipt]->GetNbinsY()<< endl;



//Define data MC to be used for unfolding
         HW7_MC_Reco[ity][ivar][ipt] = (TH1D*)HW7_MC_Reco_re[ity][ivar][ipt]->Clone();
         HW7_MC_Gen[ity][ivar][ipt] = (TH1D*)HW7_MC_Gen_re[ity][ivar][ipt]->Clone();
         HW7_h2dGenDetMC[ity][ivar][ipt] = (TH2D*)HW7_h2dGenDetMC_re[ity][ivar][ipt]->Clone();


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


//-----------------------------------------------  
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
   //  COV_Mat_NoReg[ity][ivar][ipt]->RebinY(2);
   
   sprintf(name,"Corr_Matrix_NoReg_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
   corr_mat_NoReg[ity][ivar][ipt] = new TH2D(name,name, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins, MC_Gen[ity][ivar][ipt]->GetNbinsX(), gxbins);
   corr_mat_NoReg[ity][ivar][ipt]->Sumw2();
   //corr_mat_NoReg[ity][ivar][ipt]->RebinY(2);
   // corr_mat_NoReg[ity][ivar][ipt]->RebinX(2);
   
   
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
       subtract_background(h2dGenDetMC[ity][ivar][ipt],MC_Reco[ity][ivar][ipt],MC_Gen[ity][ivar][ipt],Data_Reco[ity][ivar][ipt],fakerate,eff,purity,stbl) ;
       
       for(int bn=0; bn<(MC_Gen[ity][ivar][ipt]->GetNbinsX()); bn++){
         hist_eff[ity][ivar][ipt]->SetBinContent(bn+1,eff[bn+1]);
         hist_fake[ity][ivar][ipt]->SetBinContent(bn+1,fakerate[bn+1]) ;
         hist_purity[ity][ivar][ipt]->SetBinContent(bn+1,purity[bn+1]);
         hist_stbl[ity][ivar][ipt]->SetBinContent(bn+1,stbl[bn+1]);
       }
       
       hist_eff[ity][ivar][ipt]->SetMinimum(-0.01); hist_eff[ity][ivar][ipt]->SetMaximum(1.2);
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
 
int bias[5]={1,1,1,1,1};
//int bias[5]={2,2,2,2,2};
 
 for(int ity=0; ity <type; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
     for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       cout <<"type "<< ity << " : Variables " << ivar << " HT2 Bin : " << ipt << endl;
       file <<"["<< ity << "," << ivar << "," << ipt <<"] --->" << endl;

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

     
       h2dGenDetMC[ity][ivar][ipt]->RebinY(2);
       MC_Gen[ity][ivar][ipt]->Rebin(2);
      // h2dGenDetMC[ity][ivar][ipt]->Rebin(1,2);

       TH2* hist_migrationCoarseFine_MC = (TH2D*)h2dGenDetMC[ity][ivar][ipt]->Clone();
       TH1* input = (TH1D*)Data_Reco[ity][ivar][ipt]->Clone();
       TH1* mcgen = (TH1D*)MC_Gen[ity][ivar][ipt]->Clone();
       
      int biasScale =bias[ivar] ;
      const char *REGULARISATION_DISTRIBUTION=0;
      const char *REGULARISATION_AXISSTEERING="*[UOB]";
      char Rhoname[100], Rhotitle[100], Lcursure[100], Lcurtitle[100], TauLsure[100], TuaLtitle[100], probMat[100], probmat_title[100],Rhoname2d[100],
                Rhotitle2d[100],Ematrix[100],Ematrixtitle[100],foldback[100],foldback_title[100];

       // preserve the area
  //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
  TUnfoldDensity::EConstraint constraintMode= TUnfoldDensity::kEConstraintArea;
//  TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
//  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfoldDensity::ERegMode regMode = TUnfoldDensity::kRegModeCurvature;
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  
  //this defines the data covariance matrix
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
       //TUnfold density class	   
      //No  regularisation -------------------------------- 
       TUnfoldDensity tunfoldNoRegularisation(hist_migrationCoarseFine_MC,
					      TUnfoldDensity::kHistMapOutputVert, TUnfoldDensity::kRegModeNone,
					      TUnfoldDensity::kEConstraintNone,TUnfoldDensity::kDensityModeNone);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       tunfoldNoRegularisation.SetInput(input,biasScale,0.,covarianceM);

       //the initial bias vector is determined from the response matrix 
       //but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
       tunfoldNoRegularisation.SetBias(mcgen);   //not much affect result

       //Choose a value of tau to unfold with tau,0 means no regularization
       tunfoldNoRegularisation.DoUnfold(0.0,input,biasScale);
       
       int binnos = rnbinsx[ity][ivar][ipt];
       Int_t *binMap=new Int_t[binnos+2];
       for(Int_t i=1;i<=binnos;i++) binMap[i]=i;
       binMap[0]=-1;   binMap[binnos+1]=-1;
       
       
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
       TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
       //TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,"","*[UO]");//,"signal");
    // TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput("hist_PTunfolded_noRegularisation", "P_{T,unfolded} [GeV]","signal");
       hist_PTunfolded_noRegularisation->Write();
       
       TH2 *hist_RhoIJ_noRegularisation = tunfoldNoRegularisation.GetRhoIJtotal(NoReg_RhoIJ, NoReg_RhoIJ_tit);//,"signal");
       TH2 *hist_Rho2D_noRegularisation = tunfoldNoRegularisation.GetEmatrixTotal(NoReg_Ematrix, NoReg_Ematrix_tit);//,"signal");
       tunfoldNoRegularisation.GetEmatrix(COV_Mat_NoReg[ity][ivar][ipt],binMap);//,"signal");
       TH1 *hist_foldedback_NoReg = tunfoldNoRegularisation.GetFoldedOutput(foldback, foldback_title);//,"signal");
       TH2 *hist_prob_Noreg = tunfoldNoRegularisation.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");

       // hist_Rho_noRegularisation->Write();
       hist_Rho2D_noRegularisation->Write();
       hist_RhoIJ_noRegularisation->Write();
       hist_foldedback_NoReg->Write();
       hist_prob_Noreg->Write();
       for(int ix=0; ix<(h2dGenDetMC[ity][ivar][ipt]->GetNbinsX()); ix++){
	 for(int iy=0; iy<(h2dGenDetMC[ity][ivar][ipt]->GetNbinsY()); iy++){
	   corr_mat_NoReg[ity][ivar][ipt]->SetBinContent(ix,iy,COV_Mat_NoReg[ity][ivar][ipt]->GetBinContent(ix,iy) * 1./sqrt(COV_Mat_NoReg[ity][ivar][ipt]->GetBinContent(ix,ix)*COV_Mat_NoReg[ity][ivar][ipt]->GetBinContent(iy,iy)));
	 }
       }
       
       corr_mat_NoReg[ity][ivar][ipt]->Write();
       COV_Mat_NoReg[ity][ivar][ipt]->Write();
       
//-----------------------------------------------------------------------// L curve Scan
      
           TUnfoldDensity tunfoldTikhonovLCurve
             (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert,regMode, constraintMode,densityFlags,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
              tunfoldTikhonovLCurve.SetInput(input,biasScale,0.,covarianceM); 
              tunfoldTikhonovLCurve.SetBias(mcgen);
	      int iBest_TikhonovLCurve=-1;
              double tauBest_TikhonovLCurve=-1.,DF_TikhonovLCurve=-1.;
              int iBest;
	      int NPOINT_TikhonovLCurve=200;
              TGraph *graph_LCurve_TikhonovLCurve;
              TSpline *spline_Curvature_TikhonovLCurve;
   //          taumx = 1e-5; taumi = 1e-9; 
            taumx = 0.0; taumi = 0.0; 
	  
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
           cout << " Ibest Sure Scan : "<<ibest_sure <<" SURE Tau Value : " << tau_SURE << endl;
           file << " Ibest : "<<ibest_sure <<" //Tau = " << tau_SURE <<"      ";
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
//---------------------------------------------------Scan [Tau Minimization of a global correlation coefficient (Stefan Schmitt)]
{
          TUnfoldDensity tunfoldTauSCAN
              (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert,regMode, constraintMode,densityFlags);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
        tunfoldTauSCAN.SetInput(input,biasScale);
        tunfoldTauSCAN.SetBias(mcgen);
       	int NScan = 200;
      taumi= 1e-4; taumx = 1e-9;
      //  taumi=0.; taumx =0.;

       TGraph *lCurve;
       TSpline *Scanlcurve;
       TSpline *logTauX,*logTauY;
       const char *SCAN_DISTRIBUTION=0;
       const char *SCAN_AXISSTEERING="*[UOB]";/*"*[b]"*/ ;
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

  //      Double_t t[1],x[1],y[1];
  //      logTauX->GetKnot(ibest_tau,t[0],x[0]);
  //      logTauY->GetKnot(ibest_tau,t[0],y[0]);
//        TGraph *bestLcurve=new TGraph(1,x,y);
  //      TGraph *bestLogTauLogChi2=new TGraph(1,t,x);
        
/*

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
*/
    Scanlcurve->Write();


           hist_PTunfolded_TauScan->Write();
           hist_Rho_TauScan->Write();
           hist_Rhoij_TauScan->Write();
           hist_Ematrix_TauScan->Write();
           hist_Prob_TauScan->Write();
           hist_foldedback_TauScan->Write();
}

 //-----------------------------------------End of Variables loop
             }
          }
      }
  









delete outputFile;
  file.close();
} 


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
                            {2,4,1,0,0,0,0,0},
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
                            {2,4,1,0,0,0,0,0},
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

}   //end of rebin_hist functiomC




