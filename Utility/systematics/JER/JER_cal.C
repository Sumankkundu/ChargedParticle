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

void JER_cal(){
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
  string JERDir = "Unfold2D";

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
  TH1D *JERRelerr_HT[ntype][nusedvar][nHLTmx];                    //Nominal Data Unfolded 
  TH1D *Unfold_JER[njer][ntype][nusedvar];        //JEC Data Unfolded
  TH1D *Unfold_JER_HT[njer][ntype][nusedvar][nHLTmx];        //JEC Data Unfolded

//  TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/JEC_Tunfold_19Jan21/Unfolded_Result.root");  // Unfolded data 
  TFile *Unfoldroot = TFile::Open("Unfolded_Data_PY8_Nominal.root");  // Unfolded data 
  //TFile *Unfoldroot = TFile::Open("Unfolded_Result_MG_JER_0_6Feb21.root");  // Unfolded data 
  //TFile *Unfoldroot = TFile::Open("Unfolded_Result_MG_JERV3_0_6Feb21.root");  // Unfolded data 
  TFile *JER_file[2];
 
  //JER_file[0] = TFile::Open("Unfolded_Result_MG_Reso_1_RM_3Feb21.root");  // Unfolded data 
 //JER_file[1] = TFile::Open("Unfolded_Result_MG_Reso_2_RM_3Feb21.root");  // Unfolded data 
   JER_file[0] = TFile::Open("Unfolded_Result_JER_RM_1_3Feb21.root");  // Unfolded data 
   JER_file[1] = TFile::Open("Unfolded_Result_PY8_JER2_19Feb21.root");  // Unfolded data 
 
   //JER_file[0] = TFile::Open("Unfolded_Result_MG_JER_1_6Feb21.root");  // Unfolded data 
   //JER_file[1] = TFile::Open("Unfolded_Result_MG_JER_2_6Feb21.root");  // Unfolded data 
  TFile *JER_Result =new TFile("JER_Result.root","recreate");   

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
      Nominal[ity][ivar]= (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unfoldroot)->Clone();
      //HT2_NormalV3(Data_Reco[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      Nominal_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(NominalDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
      JERRelerr_HT[ity][ivar][ipt]= (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();       
      JERRelerr_HT[ity][ivar][ipt]->Reset();
      JERRelerr_HT[ity][ivar][ipt]->SetNameTitle(Form("jer_erro_%i_eta0_%i_pt%i" ,ity, var[ivar], ipt),Form("JER Relative %i_eta0_%i_pt%i", ity, var[ivar], ipt));      }

      for(int ijer=0; ijer <njer; ijer++){
      cout << ijer<< "  :  " ;
      Unfold_JER[ijer][ity][ivar]= (TH1D*)ReadHist1D(JERDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", JER_file[ijer])->Clone();
        
      for(int ipt=0; ipt < nHLTmx; ipt++){
      Unfold_JER_HT[ijer][ity][ivar][ipt]= (TH1D*)ReadHist1D(JERDir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), JER_file[ijer])->Clone();
           } 
      //HT2_NormalV3(Data_Reco_JEC[ity][ivar][ijec], binsRec[ity][ivar], Axisname,nHLTmx);   
        }
      }
    }
for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       double meanerr[100]={0};
       for(int ipt=0; ipt < nHLTmx; ipt++){
  
       TH1D* default_hist = (TH1D*)Nominal_HT[ity][ivar][ipt]->Clone();
     
       cout <<" ipt: " << ipt << " ivar: " << ivar << " ity: " << ity << endl;
       double xpos[100]={0};
       double xerr[100]={0};
       double mean[100]={0};

       double nnmean[100]={0};
       double nnrms[100]={0};
 
       double poserr[100]={0};
       double negerr[100]={0};

       double relerr[100]={0};

       int nbins = default_hist->GetNbinsX();
       double histmean = default_hist->GetMean(); 
       for (int ix=1; ix<nbins+1; ix++) { mean[ix] = default_hist->GetBinContent(ix);
cout << mean[ix] << " ";
       }

     
       //cout << " JEC : " << ijec ;

     TH1D* up_hist = (TH1D*)Unfold_JER_HT[0][ity][ivar][ipt]->Clone();
     TH1D* down_hist = (TH1D*)Unfold_JER_HT[1][ity][ivar][ipt]->Clone();
     
     double binwid = default_hist->GetBinWidth(nbins/2);
     
     if (binwid<0) continue;
     if (nbins>99) nbins=99;
     double histmean1 = up_hist->GetMean();
     double histmean2 = down_hist->GetMean();

     for (int ix=1; ix<nbins+1; ix++) {
     //double thewid = binwid/default_hist->GetBinWidth(ix);
     
     double val1 =  up_hist->GetBinContent(ix);
     double val2 =  down_hist->GetBinContent(ix);
//     cout << endl<<mean[ix] <<"  " << val1  << "  " << val2 << " " ;
     double dif1 = val1 - mean[ix];
     double dif2 = val2 - mean[ix];
    
     cout << "Relative error " << abs(dif1)/mean[ix] << "  " <<abs(dif2)/mean[ix]<<endl;
if (dif1*dif2<0) {
      if (dif1>0) {
            poserr[ix] += dif1*dif1; negerr[ix] += dif2*dif2;
        } else {
            poserr[ix] += dif2*dif2; negerr[ix] += dif1*dif1;
        }
     } else {
        if (dif1>0) {
           poserr[ix] += pow(TMath::Max(dif1, dif2),2);
         } else {
           negerr[ix] += pow(TMath::Min(dif1, dif2),2);
         }
       }
     }
//-------------------
 /*    double meandiff_up=abs(abs(histmean)-abs(histmean1));
       double meandiff_down=abs(abs(histmean)-abs(histmean2));
       if(meandiff_up>meandiff_down) meanerr[ipt]=meandiff_up;
       else meanerr[ipt]=meandiff_down;*/
	

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

cout << endl  ;
cout << endl  ;
     }//ipt
  }//ivar
}//ity
//------------------------

delete JER_Result;

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


