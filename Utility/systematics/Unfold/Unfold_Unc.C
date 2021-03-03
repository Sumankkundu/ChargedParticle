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

void Unfold_Unc(){
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
  string NominalDir = "Unfold2D";
  string PyDir = "Pythia8";
  string MGDir = "MG8";
  string HWDir = "HW7";
  string UnfDir = "Unfold2D";


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
  float Diff_HT[6][ntype][nusedvar][nHLTmx];

  TH1D *PY_HT[ntype][nusedvar][nHLTmx];                //Pythia8 Gen Level
  TH1D *MG_HT[ntype][nusedvar][nHLTmx];                //MG Gen Level
  TH1D *HW_HT[ntype][nusedvar][nHLTmx];                //HW Gen Level 
  
  TH1D *Unfolded_HT[ntype][nusedvar][nHLTmx];                //Unfolded Result unfolded data by pythia
  
  TH1D *JERRelerr_HT[ntype][nusedvar][nHLTmx];              //Nominal Data Unfolded 
  TH1D *Unc_err_HT[njer][ntype][nusedvar];                  //JEC Data Unfolded

  TH1D *Unfold_PY_MG[ntype][nusedvar][nHLTmx];        //PY unfold by MG
  TH1D *Unfold_PY_HW[ntype][nusedvar][nHLTmx];        //PY unfold by HW
  TH1D *Unfold_MG_PY[ntype][nusedvar][nHLTmx];        //MG unfold by PY
  TH1D *Unfold_MG_HW[ntype][nusedvar][nHLTmx];        //MG unfold by HW
  TH1D *Unfold_HW_PY[ntype][nusedvar][nHLTmx];        //HW unfold by PY
  TH1D *Unfold_HW_MG[ntype][nusedvar][nHLTmx];        //HW unfold by MG

  //TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/JEC_Tunfold_19Jan21/Unfolded_Result.root");  // Unfolded data 
  TFile *MC_file = TFile::Open("Unfolded_Data_PY8_Nominal.root");  // Unfolded data 
  TFile *unfold_file = TFile::Open("Unfolded_Data_PY8_Nominal.root");  // Unfolded data 
  
  TFile *unfold_py_mg_file = TFile::Open("Unfolded_PY8_MG_17Jan21.root"); 
  TFile *unfold_py_hw_file = TFile::Open("Unfolded_PY8_HW_14Jan21.root");
  
  TFile *unfold_mg_py_file = TFile::Open("Unfolded_MG_Py8_14Jan21.root");
  TFile *unfold_mg_hw_file = TFile::Open("Unfolded_PY8_HW_14Jan21.root");
  
  TFile *unfold_hw_py_file = TFile::Open("Unfolded_HW_Py8_14Jan21.root");
  TFile *unfold_hw_mg_file = TFile::Open("Unfolded_HW_MG_17Jan21.root");
 
  TFile *Unf_Result =new TFile("Unf_Result.root","recreate");   

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt=0; ipt < nHLTmx; ipt++){
      //Data Unfolded 

      Unfolded_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(UnfDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), MC_file);

      //MC Gen file
      PY_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(PyDir+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), MC_file);
      MG_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(MGDir+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), MC_file);
      HW_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(HWDir+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), MC_file);

      //Unfolded file
      Unfold_PY_MG[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_py_mg_file);
      Unfold_PY_HW[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_py_hw_file);      
      Unfold_MG_PY[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_mg_py_file);
      Unfold_MG_HW[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_mg_hw_file);
      Unfold_HW_PY[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_hw_py_file);
      Unfold_HW_MG[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfold_hw_mg_file);
      //differ
      for(int iun=0; iun <6 ; iun++){
      Diff_HT[iun][ity][ivar][ipt]= 0.0;
      }
      //Uncertainty file
      
      Unc_err_HT[ity][ivar][ipt]= (TH1D*)PY_HT[ity][ivar][ipt]->Clone();       
      Unc_err_HT[ity][ivar][ipt]->Reset();
      Unc_err_HT[ity][ivar][ipt]->SetNameTitle(Form("unfold_erro_%i_eta0_%i_pt%i", ity, var[ivar], ipt),Form("Unfold Relative %i_eta0_%i_pt%i", ity, var[ivar], ipt));    
           }       
         }
      }

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt=0; ipt < nHLTmx; ipt++){
       
       Integralhist(Unfolded_HT[ity][ivar][ipt]);

       Integralhist(PY_HT[ity][ivar][ipt]);
       Integralhist(MG_HT[ity][ivar][ipt]);
       Integralhist(HW_HT[ity][ivar][ipt]);


       Integralhist(Unfold_PY_MG[ity][ivar][ipt]);
       Integralhist(Unfold_PY_HW[ity][ivar][ipt]);
       Integralhist(Unfold_MG_PY[ity][ivar][ipt]);
       Integralhist(Unfold_MG_HW[ity][ivar][ipt]);
       Integralhist(Unfold_HW_PY[ity][ivar][ipt]);
       Integralhist(Unfold_HW_MG[ity][ivar][ipt]);

       }
   }
}




//error for PY8  
for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt=0; ipt < nHLTmx; ipt++){
        
        int nbins = PY_HT[ity][ivar][ipt]->GetNbinsX();
         
        for (int ix=1; ix<nbins+1; ix++) {
        double nom_unf =  Unfolded_HT[ity][ivar][ipt]->GetBinContent(ix);
	
        double nom_py =  PY_HT[ity][ivar][ipt]->GetBinContent(ix);
        double nom_mg =  MG_HT[ity][ivar][ipt]->GetBinContent(ix);
        double nom_hw =  HW_HT[ity][ivar][ipt]->GetBinContent(ix);
        	
	double py1 =  Unfold_PY_MG[ity][ivar][ipt]->GetBinContent(ix);
	double py2 =  Unfold_PY_HW[ity][ivar][ipt]->GetBinContent(ix);

	double MG1 =  Unfold_MG_PY[ity][ivar][ipt]->GetBinContent(ix);
        double MG2 =  Unfold_MG_HW[ity][ivar][ipt]->GetBinContent(ix);

	double hw1 =  Unfold_HW_PY[ity][ivar][ipt]->GetBinContent(ix);
        double hw2 =  Unfold_HW_MG[ity][ivar][ipt]->GetBinContent(ix);

        Diff_HT[0][ity][ivar][ipt]= nom_py - py1;
        Diff_HT[1][ity][ivar][ipt]= nom_py - py2;
        Diff_HT[2][ity][ivar][ipt]= nom_mg - MG1;
        Diff_HT[3][ity][ivar][ipt]= nom_mg - MG2;
        Diff_HT[4][ity][ivar][ipt]= nom_hw - hw1;
        Diff_HT[5][ity][ivar][ipt]= nom_hw - hw2;

	double errmax =0.0;
        errmax = Diff_HT[0][ity][ivar][ipt];
	//find Max error
	for(int ixx=0; ixx < 6 ; ixx++){	
	if(errmax < Diff_HT[ixx][ity][ivar][ipt]) errmax = Diff_HT[ixx][ity][ivar][ipt];  }                                     	
        //double relerr = errmax/nom_py;
        double relerr = errmax/nom_unf;
	cout << "relerr = " <<   relerr << endl ;
        Unc_err_HT[ity][ivar][ipt]->SetBinContent(ix, relerr);	 
	  }
       cout << endl;
       Unc_err_HT[ity][ivar][ipt]->Write();
//---------------------------------------------
      } 

    }
  }






/*	

       cout << "relative error : ";
       for (int ix=1; ix<nbins+1; ix++) {
       xpos[ix] = default_hist->GetBinCenter(ix);
       poserr[ix] = sqrt(poserr[ix]);
       negerr[ix] = sqrt(negerr[ix]);
       relerr[ix] = max(poserr[ix]/max(1.e-6,mean[ix]), negerr[ix]/max(1.e-6,mean[ix]));
       // cout << " error : " << poserr[ix] << "  " << negerr[ix];
       cout  << poserr[ix]/max(1.e-6,mean[ix]) << "   " << negerr[ix]/max(1.e-6,mean[ix])<<" : Max " << relerr[ix]  << "  ";
       JERRelerr_HT[ity][ivar][ipt]->SetBinContent( ix, relerr[ix]); }
       JERRelerr_HT[ity][ivar][ipt]->Write();
*/

delete Unf_Result;

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

void Integralhist(TH1 *hist){ hist->Scale(1/(hist->Integral()));}

