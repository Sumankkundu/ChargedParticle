//author : S.k.kundu
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraph.h"
#include "TObject.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TH2D.h>
#include <TTree.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include "CMS_lumi.C"

#define PYTHIA
#define MADGRAPH
#define Herwig
//#define REBIN
#define CLOUSER
//#define UFAxis
#define TRUEAXIS
void Unfoldplot(){
  
  static const int unfold_ty =1; //Unfold method to be plots 
  static const int nHLTmx=8; //HT2 Range
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =2;
  static const int nmc = 3;
  static const int ndd = 3; //0 : root Hist, 1: 1D , 2: 2D
  static const int umc = 0; //0 for Py8, 1 for MG , 2 for Herwig : which MC have used in Unfolding

  bool isstat =1;  int irbin=1;

  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={1,2,4,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[nusedvar]={3,9,15,18,24};
  Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000};
  
  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  
  const  char* ddtag[3]={"","d_","dd_"};
  const  char* trueAxis[3]={"","E","E"};
  const  char* histtag[3]={"","d_","dd_"};
  const  char* Datadirtag[3]={"Data","Data1D","Data2D"};
  const  char* foldtag[3]={"Folded","Folded1D","Folded2D"};
  const  char* dirtag[3]={"Unfold","Unfold1D","Unfold2D"};
  
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[8]={"73 < H_{T,2} < 93", "93 < H_{T,2} < 165", "165 < H_{T,2} < 225", "225 < H_{T,2} < 298", "298 < H_{T,2} < 365", "365 < H_{T,2} <452", "452 < H_{T2} <557","H_{T,2} >557"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};
  
  const  char* regN[4]={"TUnfold_NoReg_typ_","Tunfold_lscan_typ_","Tunfold_scantau_typ_","Tunfold_SURE_typ_"};
  const  char* reg_refold[4]={"Refold_NoReg_typ_","Tunfold_lscan_Refold_typ_","Tunfold_scantau_Refold_typ_","Tunfold_SURE_Refold_typ_"};
  const  char* CORR[4]={"corr_NoReg_typ_","Tunfold_lscan_corr_typ_","Tunfold_scantau_corr_typ_","Tunfold_SURE_corr_typ_"};
  const  char* COVN[4]={"Emat_NoReg_typ_","Tunfold_lscan_Emat_typ_","Tunfold_scantau_Emat_typ_","Tunfold_SURE_Emat_typ_"};
  const  char* ProbN[4]={"Prob_NoReg_typ_","Tunfold_lscan_probM_typ_","Tunfold_scantau_probM_typ_","Tunfold_SURE_probM_typ_"};
  const  char* dirname[3]={"Pythia8","MG8","HW7"};
  const  char* rebindirname[3]={"RebinPythia8","RebinMG8","RebinHW7"};

  const  char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const  char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};
  const  char* mcname[3]={"Pythia8 CP5 Tune","Madgraph","Herwig++"};
  const  char* mcnamerco[3]={"Pythia8 RECO","Madgraph RECO","Herwig++ RECO"};
  const  char* mcnamegen[3]={"Pythia8 GEN","Madgraph GEN","Herwig++ GEN"};
  const  char* itypeN[ntype]={"Jets","Charged Particles"}; 
  const  char* DataEra[3]={"Data RECO","Data RECO","Data RECO"};
  const  char* UndoldEra[3]={"Unfold 2016","Unfold 2017","Unfold 2018"};
  const  char* RefoldEra[5]={"Refold Pythia8(No Regularisation) ","Refold Pythia8(L-Curve scan)","Refold Pythia8(Scan Tau)","Refold Pythia8(ScanSURE)","Refold Pythia8(Iterative method)"};
//const  char* RefoldEra[5]={"Refold Madgraph(No Regularisation) ","Refold Madgraph(L-Curve scan)","Refold Madgraph(Scan Tau)","Refold Madgraph(ScanSURE)","Refold Pythia8(Iterative method)"};
  const  char* Unfoldtype[5]={"TUnfold(No Regularisation) ","TUnfold(L-Curve scan)","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
  const  char* closuretype[5]={"Unfolded Pythia8(No Regularisation) ","Unfolded Pythia8 (L-curve scan)","Unfolded Pythia8(Scan Tau)","Unfolded Pythia8 (ScanSURE)","Unfolded(Iterative method)"};
 // const  char* closuretype[5]={"Unfolded Madgraph(No Regularisation) ","Unfolded Madgraph (L-curve scan)","Unfolded Madgraph(Scan Tau)","Unfolded Madgraph (ScanSURE)","Unfolded(Iterative method)"};
  const  char* Modelnm[3]={"Pythia8","Madgraph","Herwig"};
  const  char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const  char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  static const int iera = 1;
  int iPeriod = 0;  int iPos=10 ;
  
  //******************************************************************************
//TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_1D/TUnfolding/testunfold2c_unfolded.root");  // Unfolded data 
TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D/TUnfolding/Unfolded_Result.root");  // Unfolded data 
//TFile *Unfoldroot = TFile::Open("testunfold2c_unfolded.root");  // Unfolded data 

  //--------------------------------------Function declearation
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_can1(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);


for(int idd=0; idd <ndd; idd++){
  TH1D *MC_gen[nmc][ntype][nusedvar][nHLTmx];   //MC gen
  TH1D *MC_gen_miss[nmc][ntype][nusedvar][nHLTmx]; //MC Gen-miss
  TH1D *Ex_MC_gen[nmc][ntype][nusedvar][nHLTmx]; //Truth Dist Gen
  TH1D *Ex_MC_gen_miss[nmc][ntype][nusedvar][nHLTmx];  //Truth Gen-miss
  
  TH1D *MC_reco[nmc][ntype][nusedvar][nHLTmx];   //MC Reco 
  TH1D *MC_reco_fake[nmc][ntype][nusedvar][nHLTmx]; //MC Reco -Fake
  TH1D *Ex_MC_reco[nmc][ntype][nusedvar][nHLTmx];   //Truth Dist Reco
  TH1D *Ex_MC_reco_fake[nmc][ntype][nusedvar][nHLTmx];  //Truth Reco -Fake
  
  TH2D *MC_Res[nmc][ntype][nusedvar][nHLTmx];   //RM
#ifdef CLOUSER  
  TH1D *Psudo_Data_gen[ntype][nusedvar][nHLTmx];  
  TH1D *Psudo_Data_reco_fake[ntype][nusedvar][nHLTmx];  
  TH1D *Ex_Psudo_Data_gen[ntype][nusedvar][nHLTmx];
  TH1D *Ex_Psudo_Data_gen_miss[ntype][nusedvar][nHLTmx];
  TH1D *Ex_Psudo_Data_reco_fake[ntype][nusedvar][nHLTmx];
#endif

  TH1D *Data_reco[ntype][nusedvar][nHLTmx];
  TH1D *Ex_Data_reco[ntype][nusedvar][nHLTmx];
  
  TH1D *hist_eff[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_fake[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_purity[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_stbl[nmc][ntype][nusedvar][nHLTmx];

  TH1D *Unfold[unfold_ty][ntype][nusedvar][nHLTmx];
  TH1D *Refold[unfold_ty][ntype][nusedvar][nHLTmx];

  TH1D *Ex_Unfold[unfold_ty][ntype][nusedvar][nHLTmx];
  TH1D *Ex_Refold[unfold_ty][ntype][nusedvar][nHLTmx];

  TH2D *Corr[unfold_ty][ntype][nusedvar][nHLTmx];
  TH2D *Prob[unfold_ty][ntype][nusedvar][nHLTmx];
  TH2D *Ematrix[unfold_ty][ntype][nusedvar][nHLTmx];


for(int ity=0; ity <ntype; ity++){
      for(int ivar =0;ivar < nusedvar; ivar++){
           for(Int_t ipt =0; ipt < nHLTmx ; ipt++){     
#ifdef CLOUSER
          sprintf(histname, "%s/%sgen_typ_%i_pt%i_eta0_%i",Datadirtag[idd],ddtag[idd],ity, ipt, var[ivar]);
          Psudo_Data_gen[ity][ivar][ipt] =(TH1D*)Unfoldroot->Get(histname);    Psudo_Data_gen[ity][ivar][ipt]->Rebin(irbin);
          cout << histname << endl;
#ifdef TRUEAXIS
if(idd >=1){
	  sprintf(histname, "%s/%s%sgen_typ_%i_pt%i_eta0_%i",Datadirtag[idd],trueAxis[idd],ddtag[idd],ity, ipt, var[ivar]);
	  cout << histname << endl;
          if(idd==1){Ex_Psudo_Data_gen[ity][ivar][ipt] =(TH1D*)Unfoldroot->Get(histname); }
          if(idd==2){TH2D* temp_Psudo_Data_gen =(TH2D*)Unfoldroot->Get(histname); Ex_Psudo_Data_gen[ity][ivar][ipt] = temp_Psudo_Data_gen->ProjectionX();}
	  Ex_Psudo_Data_gen[ity][ivar][ipt]->Rebin(irbin);          cout << histname << endl;

          sprintf(histname, "%s/%s%sRecominusfake_%i_pt%i_eta0_%i",Datadirtag[idd],trueAxis[idd],ddtag[idd],ity, ipt, var[ivar]);
          cout << histname << endl;
          if(idd==1){Ex_Psudo_Data_reco_fake[ity][ivar][ipt] =(TH1D*)Unfoldroot->Get(histname); }
          if(idd==2){TH2D* temp_Psudo_reco_fake =(TH2D*)Unfoldroot->Get(histname); Ex_Psudo_Data_reco_fake[ity][ivar][ipt] = temp_Psudo_reco_fake->ProjectionX();}
          Ex_Psudo_Data_reco_fake[ity][ivar][ipt]->Rebin(irbin);          cout << histname << endl;

          sprintf(histname, "%s/%s%sGenminusmiss_%i_pt%i_eta0_%i",Datadirtag[idd],trueAxis[idd],ddtag[idd],ity, ipt, var[ivar]);
          cout << histname << endl;
          if(idd==1){Ex_Psudo_Data_gen_miss[ity][ivar][ipt] =(TH1D*)Unfoldroot->Get(histname); }
          if(idd==2){TH2D* temp_Psudo_gen_miss =(TH2D*)Unfoldroot->Get(histname); Ex_Psudo_Data_gen_miss[ity][ivar][ipt] = temp_Psudo_gen_miss->ProjectionX();}
          Ex_Psudo_Data_gen_miss[ity][ivar][ipt]->Rebin(irbin);          cout << histname << endl;
}
#endif
#endif

//Exclude Underflow overlow
#ifdef UFAxis
          TH1D *NewData;          NewData = (TH1D*)Psudo_Data_gen[ity][ivar][ipt]->Clone();
	  int NbinxD = NewData->GetNbinsX();          NewData->Reset();
	  for(int ix=1; ix < NbinxD+1 ; ix++){
	  NewData->SetBinContent(ix,Psudo_Data_gen[ity][ivar][ipt]->GetBinContent(ix));
	  NewData->SetBinError(ix, sqrt(Psudo_Data_gen[ity][ivar][ipt]->GetBinError(ix)*Psudo_Data_gen[ity][ivar][ipt]->GetBinError(ix)));}
#endif
	}
      }
    }


  for(int  imc =0; imc < nmc ; imc++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	for(Int_t ipt =0; ipt < nHLTmx ; ipt++){     
	  sprintf(histname, "%s/%sgen_typ_%i_pt%i_eta0_%i",dirname[imc],ddtag[idd], ity, ipt, var[ivar]); 
	  MC_gen[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  MC_gen[imc][ity][ivar][ipt]->Rebin(irbin);
	  
	  sprintf(histname, "%s/%sreco_typ_%i_pt%i_eta0_%i", dirname[imc],ddtag[idd], ity,  ipt, var[ivar]); 
	  MC_reco[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);	  cout << histname << endl;

	  //Fake and Miss corrected Reco and Gen
          sprintf(histname, "%s/%sGenminusmiss_%i_pt%i_eta0_%i",dirname[imc],ddtag[idd], ity, ipt, var[ivar]);
          MC_gen_miss[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);     MC_gen_miss[imc][ity][ivar][ipt]->Rebin(irbin);

          sprintf(histname, "%s/%sRecominusfake_%i_pt%i_eta0_%i", dirname[imc],ddtag[idd], ity,  ipt, var[ivar]);
          MC_reco_fake[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);       cout << histname << endl;
#ifdef TRUEAXIS
if(idd>=1){
          sprintf(histname, "%s/%s%sgen_typ_%i_pt%i_eta0_%i",dirname[imc],trueAxis[idd],ddtag[idd], ity, ipt, var[ivar]); 
          if(idd==1){Ex_MC_gen[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);}     
          if(idd==2){TH2D* temp_MC_gen = (TH2D*) Unfoldroot->Get(histname); Ex_MC_gen[imc][ity][ivar][ipt] = temp_MC_gen->ProjectionX();    } 
	  Ex_MC_gen[imc][ity][ivar][ipt]->Rebin(irbin);

          sprintf(histname, "%s/%s%sreco_typ_%i_pt%i_eta0_%i", dirname[imc], trueAxis[idd],ddtag[idd], ity,  ipt, var[ivar]);
          if(idd==1){Ex_MC_reco[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname); cout << histname << endl;}
          if(idd==2){TH2D* temp_MC_reco = (TH2D*) Unfoldroot->Get(histname);   Ex_MC_reco[imc][ity][ivar][ipt] = temp_MC_reco->ProjectionX();  cout << histname << endl;}
      

          sprintf(histname, "%s/%s%sGenminusmiss_%i_pt%i_eta0_%i",dirname[imc],trueAxis[idd],ddtag[idd], ity, ipt, var[ivar]);
          if(idd==1){Ex_MC_gen_miss[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);}
          if(idd==2){TH2D* temp_MC_gen_miss = (TH2D*) Unfoldroot->Get(histname); Ex_MC_gen_miss[imc][ity][ivar][ipt] = temp_MC_gen_miss->ProjectionX();    }
          Ex_MC_gen_miss[imc][ity][ivar][ipt]->Rebin(irbin);

          sprintf(histname, "%s/%s%sRecominusfake_%i_pt%i_eta0_%i", dirname[imc], trueAxis[idd],ddtag[idd], ity,  ipt, var[ivar]);
          if(idd==1){Ex_MC_reco_fake[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname); cout << histname << endl;}
          if(idd==2){TH2D* temp_MC_reco_fake = (TH2D*) Unfoldroot->Get(histname);   Ex_MC_reco_fake[imc][ity][ivar][ipt] = temp_MC_reco_fake->ProjectionX();  cout << histname << endl;}
          }
#endif
	  //-------------------------------------Read Reso--------------------------------------
	  sprintf(histname, "%s/%scorr_typ_%i_pt%i_eta0_%i",dirname[imc], ddtag[idd],ity,  ipt, var[ivar]); 
	  MC_Res[imc][ity][ivar][ipt] = (TH2D*) Unfoldroot->Get(histname);	  cout << histname << endl;
	  
	  //------------------------------------------------------
	  sprintf(histname,"%s/%smiss_rate_%i_pt%i_eta0_%i",dirname[imc],ddtag[idd], ity, ipt, var[ivar]);
	  hist_eff[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);  cout << histname << endl;
	  sprintf(histname,"%s/%sfake_rate_%i_pt%i_eta0_%i",dirname[imc],ddtag[idd], ity, ipt, var[ivar]);
	  hist_fake[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname); cout << histname << endl;
	  sprintf(histname,"%s/%sPurity_%i_pt%i_eta0_%i", dirname[imc],ddtag[idd], ity, ipt, var[ivar]);
	  hist_purity[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname); cout << histname << endl;
	  sprintf(histname,"%s/%sstability_%i_pt%i_eta0_%i", dirname[imc],ddtag[idd], ity, ipt, var[ivar]);
	  hist_stbl[imc][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname); cout << histname << endl;
	  
	}
      }
    }
  }  //end of one MCinput root file  reading
  
  //----------------------------------------------------------------------------//Read Data
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(Int_t ipt =0; ipt < nHLTmx ; ipt++){
	sprintf(histname, "%s/%sreco_typ_%i_pt%i_eta0_%i",Datadirtag[idd], ddtag[idd], ity, ipt, var[ivar]); //reco_typ_1_pt6_eta0_15
	Data_reco[ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);	cout << histname << endl;     

#ifdef TRUEAXIS
if(idd>=1){
        sprintf(histname, "%s/%s%sreco_typ_%i_pt%i_eta0_%i",Datadirtag[idd], trueAxis[idd],ddtag[idd], ity, ipt, var[ivar]);
        if(idd==1){Ex_Data_reco[ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);  cout << histname << endl;}
        if(idd==2){TH2D* temp_Data_reco = (TH2D*) Unfoldroot->Get(histname); Ex_Data_reco[ity][ivar][ipt] = temp_Data_reco->ProjectionX();  cout << histname << endl;}
}
#endif
      }
    }
  }
  
  //---------------------------------------------------------//Read Unfolded
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
	  sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirtag[idd] ,ddtag[idd],regN[iun],ity,ipt, var[ivar]); 
	  Unfold[iun][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);	  cout << histname << endl;
#ifdef TRUEAXIS
if(idd>=1){
          sprintf(histname, "%s/%s%s%s%i_pt%i_eta0_%i",dirtag[idd], trueAxis[idd],ddtag[idd],regN[iun], ity, ipt, var[ivar]);
          if(idd==1){Ex_Unfold[iun][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);   cout << histname << endl;}
          if(idd==2){TH2D* temp_Unfold = (TH2D*) Unfoldroot->Get(histname); Ex_Unfold[iun][ity][ivar][ipt] = temp_Unfold->ProjectionX();  cout << histname << endl;}
	}
#endif
//Exclude Underflow and overflow bins
#ifdef UFAxis
	  TH1D *NewData;   NewData = (TH1D*) Unfold[iun][ity][ivar][ipt]->Clone(); int NbinxD = NewData->GetNbinsX();  NewData->Reset();
          for(int ix=1; ix < NbinxD+1 ; ix++){
          NewData->SetBinContent(ix, Unfold[iun][ity][ivar][ipt]->GetBinContent(ix));
          NewData->SetBinError(ix, sqrt( Unfold[iun][ity][ivar][ipt]->GetBinError(ix)* Unfold[iun][ity][ivar][ipt]->GetBinError(ix)));       }
#endif
	  sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirtag[idd] ,ddtag[idd],CORR[iun],ity,ipt, var[ivar]); 
          Corr[iun][ity][ivar][ipt]=(TH2D*) Unfoldroot->Get(histname);
	  sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirtag[idd] ,ddtag[idd],COVN[iun],ity,ipt, var[ivar]); 
           Ematrix[iun][ity][ivar][ipt]= (TH2D*) Unfoldroot->Get(histname);
	  sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirtag[idd] ,ddtag[idd],ProbN[iun],ity,ipt, var[ivar]); 
	   Prob[iun][ity][ivar][ipt] = (TH2D*) Unfoldroot->Get(histname);
	}
      }
    } //unfolded histUnfold
  }
//------------------------------------------------------------//Read Folded Back
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
	  sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirtag[idd] ,ddtag[idd],reg_refold[iun],ity,ipt, var[ivar]); //Tunfolded_typ_0_pt0_eta0_3
	  Refold[iun][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);	  cout << histname << endl;
#ifdef TRUEAXIS
if(idd>=1){
          sprintf(histname, "%s/%s%s%s%i_pt%i_eta0_%i",dirtag[idd], trueAxis[idd],ddtag[idd],reg_refold[iun], ity, ipt, var[ivar]);
           if(idd==1){Ex_Refold[iun][ity][ivar][ipt] = (TH1D*) Unfoldroot->Get(histname);  cout << histname << endl;}
           if(idd==2){TH2D* temp_Refold = (TH2D*) Unfoldroot->Get(histname); Ex_Refold[iun][ity][ivar][ipt] = temp_Refold->ProjectionX();  cout << histname << endl;}
      }
#endif
	}
      }
    } //unfolded histUnfold
  }
  
  cout <<" Read histogram OK " << endl;
  //PLot canvas declear
  TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 700,600 );  //for Reco
  TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,575 );  //for Corr
  TCanvas *cpt6 = new TCanvas("cpt6", "canvas6", 600,575 );  //for Prob
  TCanvas *cpt7 = new TCanvas("cpt7", "canvas7", 600,575 );  //for COV
  TCanvas *cpt8 = new TCanvas("cpt8", "canvas8", 600,575 );  //for Response
  TCanvas *cpt9 = new TCanvas("cpt9", "canvas9", 600,575 );  //for Projection
  
  //Reco comparison 
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){
	TH1D *MyHist  = (TH1D*) Data_reco[ity][ivar][ipt]->Clone();
	Integralhist(MyHist);
	//divBybinWidth(MyHist);
	Myplotset(MyHist,0,0);
	
	if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
	else if(ipt==7){ sprintf(Title,"%s:     H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
	sprintf(Yaxis," %s" ,Esvlogy[ivar]);
	MyHist->SetTitle(Title);
	MyHist->GetXaxis()->SetTitle("");
	MyHist->GetYaxis()->SetTitle(Yaxis);
	
	TH1D *MC_input[nmc];
	const char *MCinput_index[nmc+1];
	const char *data_index[1];
	for(int iout = 0 ; iout < nmc ; iout++){
	  MC_input[iout] = (TH1D*) MC_reco[iout][ity][ivar][ipt]->Clone(); 
	  Integralhist(MC_input[iout]);
	 //divBybinWidth(MC_input[iout]);
	  MCinput_index[iout]= mcnamerco[iout]; }
	  data_index[0]= DataEra[1]; 
	
	char lplot_xtitle[100];
	sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
	//float ratio_range1[2]={1.2,0.9};
	int num1[2]={nmc,1} ;
	float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	
	cpt0->cd();
	SetMycanvas(cpt0,0,0,0,0,0);

	cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
	CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	
	sprintf(pdfname, "%sRecoEVS_Plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sRecoEVS_Plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sRecoEVS_Plot.pdf)",ddtag[idd]); 
	if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
	}else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
	}else{  cpt0->Print(pdfname1,"pdf");};
      }
    }  
  }//-------------------end of RECO PLOT    
//------------------------- Reco Projection comparison-----------------------
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0;ivar < nusedvar; ivar++){
         for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirname[umc],ddtag[idd],"ProjectX_",ity,ipt, var[ivar]);    TH1* RMx=(TH1D*) Unfoldroot->Get(histname);
      sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirname[umc],ddtag[idd],"Recominusfake_",ity,ipt, var[ivar]);    TH1* RecoFake=(TH1D*) Unfoldroot->Get(histname);
      cpt9->cd();      SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);
      
      TH1D *MyHist  = (TH1D*)RMx->Clone();
      Integralhist(MyHist);
        //divBybinWidth(MyHist);
        Myplotset(MyHist,0,0);

        if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title,"%s:     H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
        sprintf(Yaxis," %s" ,Esvlogy[ivar]);
        MyHist->SetTitle(Title);
        MyHist->GetXaxis()->SetTitle("");
        MyHist->GetYaxis()->SetTitle(Yaxis);
        int imc =1;
        TH1D *MC_input[imc];
        const char *MCinput_index[imc+1];
        const char *data_index[1];
        for(int iout = 0 ; iout < imc ; iout++){
          MC_input[iout] = (TH1D*) RecoFake->Clone();
          Integralhist(MC_input[iout]);
         // divBybinWidth(MC_input[iout]);
          MCinput_index[iout]= "Reco #minus Fake"; }
          data_index[0]= "RM ProjectionX";

        char lplot_xtitle[100];
        sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
        //float ratio_range1[2]={1.2,0.9};
        int num1[2]={imc,1} ;
        float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};

        cpt9->Clear();  cpt9->cd();
        SetMycanvas(cpt9,0,0,0,0,0);

        cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();
        sprintf(pdfname, "%sRM_ProjectX_Plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sRM_ProjectX_Plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sRM_ProjectX_Plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt9->Print(pdfname2,"pdf");
        }else{  cpt9->Print(pdfname1,"pdf");};
      }
    }
  }
//--------------------------Gen Projection comparison------------------
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirname[umc],ddtag[idd],"ProjectY_",ity,ipt, var[ivar]);  TH1* RMx=(TH1D*) Unfoldroot->Get(histname);
      sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",dirname[umc],ddtag[idd],"Genminusmiss_",ity,ipt, var[ivar]); TH1* GenMiss=(TH1D*) Unfoldroot->Get(histname);
      cpt9->cd();
      SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);

      TH1D *MyHist  = (TH1D*)RMx->Clone();
      Integralhist(MyHist);
        //divBybinWidth(MyHist);
        Myplotset(MyHist,0,0);

        if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title,"%s:     H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
        sprintf(Yaxis," %s" ,Esvlogy[ivar]);
        MyHist->SetTitle(Title);
        MyHist->GetXaxis()->SetTitle("");
        MyHist->GetYaxis()->SetTitle(Yaxis);
        int imc =1;
        TH1D *MC_input[imc];
        const char *MCinput_index[imc+1];
        const char *data_index[1];
        for(int iout = 0 ; iout < imc ; iout++){
          MC_input[iout] = (TH1D*) GenMiss->Clone();
          Integralhist(MC_input[iout]);
         // divBybinWidth(MC_input[iout]);
          MCinput_index[iout]= "Gen #minus Miss"; }
          data_index[0]= "RM ProjectionY";

        char lplot_xtitle[100];
        sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
        //float ratio_range1[2]={1.2,0.9};
        int num1[2]={imc,1} ;
        float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.15,0.85};

        cpt9->Clear();        
	cpt9->cd();
        SetMycanvas(cpt9,0,0,0,0,0);

        cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

        sprintf(pdfname, "%sRM_ProjectY_Plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sRM_ProjectY_Plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sRM_ProjectY_Plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt9->Print(pdfname2,"pdf");
        }else{  cpt9->Print(pdfname1,"pdf");};
      }
    }
  }
//--------------------------------Fold comparison---------------
//--------------------------------------------------------------
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){

      sprintf(histname, "%s/%s%s%i_pt%i_eta0_%i",foldtag[idd],ddtag[idd],"Fold_",ity,ipt, var[ivar]);  TH1D* GenFold=(TH1D*)Unfoldroot->Get(histname);
        
        TH1D *MyHist  = (TH1D*)GenFold->Clone();
        Integralhist(MyHist);
        //divBybinWidth(MyHist);
        Myplotset(MyHist,0,0);

        if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title,"%s:     H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
        sprintf(Yaxis," %s" ,Esvlogy[ivar]);
        MyHist->SetTitle(Title);
        MyHist->GetXaxis()->SetTitle("");
        MyHist->GetYaxis()->SetTitle(Yaxis);
        int imc =1;
        TH1D *MC_input[imc];
        const char *MCinput_index[imc+1];
        const char *data_index[1];
        for(int iout = 0 ; iout < imc ; iout++){
        //sprintf(histname, "Pythia8/%s%s%i_pt%i_eta0_%i",ddtag[idd],"Recominusfake_",ity,ipt, var[ivar]); //Tunfolded_typ_0_pt0_eta0_3
        //MC_input[iout]=(TH1D*) Unfoldroot->Get(histname);   //Genfold  are already corrected for fake 
          MC_input[iout] = (TH1D*) MC_reco[umc][ity][ivar][ipt]->Clone(); 
          Integralhist(MC_input[iout]);
         // divBybinWidth(MC_input[iout]);
          MCinput_index[iout]= "Pythia RECO"; }
          data_index[0]= "Folded";

        char lplot_xtitle[100];
        sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
        //float ratio_range1[2]={1.2,0.9};
        int num1[2]={imc,1} ;
        float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
        
	cpt9->Clear();
        cpt9->cd();
        SetMycanvas(cpt9,0,0,0,0,0);

        cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

        sprintf(pdfname,"%sGenfold_Plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sGenfold_Plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sGenfold_Plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt9->Print(pdfname2,"pdf");
        }else{  cpt9->Print(pdfname1,"pdf");};

      }
    }  //end of phase space cut and variable loop
  }
#ifdef TRUEAXIS
if(idd>=1){
//--------------------------------------------------------------
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){

      TH1D *GenFold;
      sprintf(histname, "%s/E%s%s%i_pt%i_eta0_%i",foldtag[idd],ddtag[idd],"Fold_",ity,ipt, var[ivar]); cout <<histname<< endl; 
      if(idd==1){GenFold = (TH1D*)Unfoldroot->Get(histname); }
      if(idd==2){TH2D* GenFold2D=(TH2D*)Unfoldroot->Get(histname); GenFold = GenFold2D->ProjectionX();}
      
        TH1D *MyHist  = (TH1D*)GenFold->Clone();
        Integralhist(MyHist);
        //divBybinWidth(MyHist);
        Myplotset(MyHist,0,0);

        if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title,"%s:     H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
        sprintf(Yaxis," %s" ,Esvlogy[ivar]);
        MyHist->SetTitle(Title);
        MyHist->GetXaxis()->SetTitle("");
        MyHist->GetYaxis()->SetTitle(Yaxis);
        int imc =1;
        TH1D *MC_input[imc];
        const char *MCinput_index[imc+1];
        const char *data_index[1];
        for(int iout = 0 ; iout < imc ; iout++){
        sprintf(histname, "Pythia8/%s%s%i_pt%i_eta0_%i",ddtag[idd],"Recominusfake_",ity,ipt, var[ivar]);
        //MC_input[iout]=(TH1D*) Unfoldroot->Get(histname);
        MC_input[iout] = (TH1D*) Ex_MC_reco[iout][ity][ivar][ipt]->Clone();
        Integralhist(MC_input[iout]);
        //divBybinWidth(MC_input[iout]);
        MCinput_index[iout]= "Pythia RECO"; }
        data_index[0]= "Folded";

        char lplot_xtitle[100];
        sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
        //float ratio_range1[2]={1.2,0.9};
        int num1[2]={imc,1} ;
        float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};

        cpt9->Clear();
        cpt9->cd();
        SetMycanvas(cpt9,0,0,0,0,0);

        cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();
        sprintf(pdfname, "%s_TrueGenfold_Plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%s_TrueGenfold_Plot.pdf",ddtag[idd]); sprintf(pdfname2, "%s_TrueGenfold_Plot.pdf)",ddtag[idd]);
        if(ity==0 && ivar==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt9->Print(pdfname2,"pdf");
        }else{  cpt9->Print(pdfname1,"pdf");};
      }
    }  
  }
}//if(idd>=0)
#endif
 
 
//Response Matrix------------------------------------------------------
for(int  imc =0; imc < nmc ; imc++){
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){
	cpt8->cd();
//	MC_Res[imc][ity][ivar][ipt]->RebinY(irbin);
        char lplot_xtitle[100]; char lplot_ytitle[100];
        sprintf(lplot_xtitle, "RECO     %s",Esvlogx[ivar]); sprintf(lplot_ytitle, "GEN      %s",Esvlogx[ivar]);	
	double titoff1[3]={1.2,1.3,1.0};
        double titsize1[3] ={0.035,0.035,0.035};
	
	SetMycanvas(cpt8,0,0.1,0.15,0.05,0.1);
        Set2dHist( MC_Res[imc][ity][ivar][ipt],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
        
	MC_Res[imc][ity][ivar][ipt]->Draw("colz");
        
	if(ipt<=6){sprintf(Title,"%i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title," H_{T,2} > %i %s", HT2range[ipt] ,"GeV/c" );}
         
	TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
	CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,itypeN[ity] , "");leg1->SetTextColor(-8);leg1->Draw();
        CMS_lumi( cpt8, iPeriod, iPos ); cpt8->Update();
	sprintf(pdfname, "%sResponse_Mat%i.pdf(",ddtag[idd],imc); sprintf(pdfname1, "%sResponse_Mat%i.pdf",ddtag[idd],imc); sprintf(pdfname2, "%sResponse_Mat%i.pdf)",ddtag[idd],imc);
        if(ity==0 && ivar==0 && ipt ==0){cpt8->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt8->Print(pdfname2,"pdf");
        }else{  cpt8->Print(pdfname1,"pdf");};
        }
    }
  }
}
  cout <<" RECO PLOT OK " << endl;
//------------------------------------------------------Efficiency, Purity,Fake rate, stability
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
  TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 600,600 );  //for 
  TCanvas *cpt2 = new TCanvas("cpt2", "canvas2", 600,600 );  //for ESVs
  TCanvas *cpt3 = new TCanvas("cpt3", "canvas3", 600,600 );  //for ESVs
  TCanvas *cpt4 = new TCanvas("cpt4", "canvas4", 800,800 );  //for 

   	TLegend *leg2 = new TLegend(0.4,0.5,0.7,0.8);
	CTLegend(leg2," ","");
	TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        CTLegend(leg1,"", itypeN[ity]); 

	for(int ipt = 0 ; ipt <nHLTmx; ipt++){
        cpt1->cd();
	
	if(ipt<=6){sprintf(Title,"%i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title," H_{T,2} > %i %s", HT2range[ipt] ,"GeV/c" );}

        for (int i = 1; i <= hist_eff[umc][ity][ivar][ipt]->GetNbinsX(); ++i) {
         double content = 1 - hist_eff[umc][ity][ivar][ipt]->GetBinContent(i);
         hist_eff[umc][ity][ivar][ipt]->SetBinContent(i, content);
       }
         hist_eff[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_eff[umc][ity][ivar][ipt]->SetMaximum(1.01);
         hist_fake[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_fake[umc][ity][ivar][ipt]->SetMaximum(1.01);

        SetMycanvas(cpt1,0,0.1,0.15,0.05,0.12);
        Myplotset(hist_eff[umc][ity][ivar][ipt],Esvlogx[ivar],"Efficiency");
        hist_eff[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
	hist_eff[umc][ity][ivar][ipt]->Draw("same e1"); leg1->Draw();
   
        leg2->AddEntry(hist_eff[0][ity][ivar][ipt], Title ,"lp");
	if(ipt==7){leg2->Draw();}

        cpt2->cd();
//	SetMycanvas(cpt2);
        SetMycanvas(cpt2,0,0.1,0.15,0.05,0.12);
 	Myplotset(hist_purity[umc][ity][ivar][ipt],Esvlogx[ivar],"Purity");
 	hist_purity[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
        hist_purity[umc][ity][ivar][ipt]->Draw("same"); leg1->Draw();
	 if(ipt==7){leg2->Draw();}
	cpt2->Update();

        cpt3->cd();
        SetMycanvas(cpt3,0,0.1,0.1,0.05,0.12);
        Myplotset(hist_fake[umc][ity][ivar][ipt],Esvlogx[ivar],"Fake rate");
        hist_fake[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
        hist_fake[umc][ity][ivar][ipt]->Draw("same e1");leg1->Draw();
        if(ipt==7){leg2->Draw();}
         cpt3->Update();
	cpt4->cd();
        SetMycanvas(cpt4,0,0.1,0.1,0.05,0.12);
        Myplotset(hist_stbl[umc][ity][ivar][ipt],Esvlogx[ivar],"Stability");
        hist_stbl[umc][ity][ivar][ipt]->SetLineColor(color[ipt]); 
        hist_stbl[umc][ity][ivar][ipt]->Draw("same"); leg1->Draw();
        if(ipt==7){leg2->Draw();}
   	 cpt4->Update();
    }
    sprintf(pdfname, "%seffi_plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%seffi_plot.pdf",ddtag[idd]); sprintf(pdfname2, "%seffi_plot.pdf)",ddtag[idd]);
        if(ity==0 && ivar==0 ){cpt1->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt1->Print(pdfname2,"pdf");
        }else{  cpt1->Print(pdfname1,"pdf");};
	cpt1->Clear();
    sprintf(pdfname, "%spuri_plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%spuri_plot.pdf",ddtag[idd]); sprintf(pdfname2, "%spuri_plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 ){cpt2->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt2->Print(pdfname2,"pdf");
        }else{  cpt2->Print(pdfname1,"pdf");};
	cpt2->Clear();
     sprintf(pdfname, "%sfake_plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sfake_plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sfake_plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 ){cpt3->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4) {cpt3->Print(pdfname2,"pdf");
        }else{  cpt3->Print(pdfname1,"pdf");};
	cpt3->Clear();
     sprintf(pdfname, "%sstab_plot.pdf(",ddtag[idd]); sprintf(pdfname1, "%sstab_plot.pdf",ddtag[idd]); sprintf(pdfname2, "%sstab_plot.pdf)",ddtag[idd]); 
        if(ity==0 && ivar==0 ){cpt4->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt4->Print(pdfname2,"pdf");
        }else{  cpt4->Print(pdfname1,"pdf");};
	cpt4->Clear();
    
     }
  }   

cout << "Effi, Stabilty, fake, purity " <<endl;
  
  //Unfold comparisoin
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	for(int ipt = 0 ; ipt <nHLTmx; ipt++){
	  
	  TH1D *MyHist  = (TH1D*)Unfold[iun][ity][ivar][ipt]->Clone();
	  Integralhist(MyHist);
	  //divBybinWidth(MyHist);
	  Myplotset(MyHist,0,0);
	  if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
	  else if(ipt==7){ sprintf(Title,"%s:       H_{T,2} > %i %s",itypeN[ity] ,  HT2range[ipt] ,"GeV/c" );}
	  sprintf(Yaxis," %s" ,Esvlogy[ivar]);
	  MyHist->SetTitle(Title);
	  MyHist->GetXaxis()->SetTitle("");
	  MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  TH1D *MC_input[nmc];
	  const char *MCinput_index[nmc], *data_index[1];
	  for(int iout = 0 ; iout < nmc ; iout++){
	    MC_input[iout] = (TH1D*) MC_gen[iout][ity][ivar][ipt]->Clone();
	    Integralhist(MC_input[iout]);
	    //divBybinWidth(MC_input[iout]);
	    
	    MCinput_index[iout]= mcnamegen[iout]; }
	 
	  data_index[0]= Unfoldtype[iun]; 
	  //data_index[0]= UndoldEra[1]; 
	  char lplot_xtitle[100];
	  sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
	  //  float ratio_range1[2]={1.2,0.9};
	  int num1[2]={nmc,1} ;
	  float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  cpt0->cd();
          SetMycanvas(cpt0,0,0,0,0,0);
	  cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  sprintf(pdfname, "%sTUnfold_plot_%i.pdf(" ,ddtag[idd],iun); sprintf(pdfname1, "%sTUnfold_plot_%i.pdf" ,ddtag[idd],iun);sprintf(pdfname2, "%sTUnfold_plot_%i.pdf)" ,ddtag[idd],iun);
	  if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
	  }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
	  }else{cpt0->Print(pdfname,"pdf");};
	 
       cpt5->cd();
       double titoff1[3]={1.2,1.3,1.0};
       double titsize1[3] ={0.035,0.035,0.035};
       SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       gStyle->SetPaintTextFormat( "4.2f");
       Set2dHist(Corr[iun][ity][ivar][ipt],lplot_xtitle, lplot_xtitle,"correlation coefficients", titoff1, titsize1);
       Corr[iun][ity][ivar][ipt]->Draw("colz ");
       //Corr[iun][ity][ivar][ipt]->Draw("colz text");
       if(ipt<=6){sprintf(Title,"%i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        else if(ipt==7){ sprintf(Title," H_{T,2} > %i %s", HT2range[ipt] ,"GeV/c" );}

        TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
        CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,itypeN[ity] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       sprintf(pdfname, "%sTUnfold_corr_%i.pdf(" ,ddtag[idd],iun); sprintf(pdfname1, "%sTUnfold_corr_%i.pdf" ,ddtag[idd],iun);sprintf(pdfname2, "%sTUnfold_corr_%i.pdf)" ,ddtag[idd],iun); 
          if(ity==0 && ivar==0 && ipt ==0){cpt5->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt5->Print(pdfname2,"pdf");
          }else{cpt5->Print(pdfname,"pdf");};

       cpt6->cd();
       SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       Set2dHist(Prob[iun][ity][ivar][ipt],lplot_xtitle, lplot_xtitle,"Probability",titoff1, titsize1);
       Prob[iun][ity][ivar][ipt]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       sprintf(pdfname, "%sTUnfold_prob_%i.pdf(" ,ddtag[idd],iun); sprintf(pdfname1, "%sTUnfold_prob_%i.pdf" ,ddtag[idd],iun);sprintf(pdfname2, "%sTUnfold_prob_%i.pdf)" ,ddtag[idd],iun); 
          if(ity==0 && ivar==0 && ipt ==0){cpt6->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt6->Print(pdfname2,"pdf");
          }else{cpt6->Print(pdfname,"pdf");};

       cpt7->cd();
       SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       Set2dHist(Ematrix[iun][ity][ivar][ipt],lplot_xtitle, lplot_xtitle,"",titoff1, titsize1);
       Ematrix[iun][ity][ivar][ipt]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       sprintf(pdfname, "%sTUnfold_COV_%i.pdf(" ,ddtag[idd],iun); sprintf(pdfname1, "%sTUnfold_COV_%i.pdf" ,ddtag[idd],iun);sprintf(pdfname2, "%sTUnfold_COV_%i.pdf)" ,ddtag[idd], iun);
          if(ity==0 && ivar==0 && ipt ==0){cpt7->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt7->Print(pdfname2,"pdf");
          }else{cpt7->Print(pdfname,"pdf");};
	}
      }  //end of phase space cut and variable loop
    }
  }//End of Unfolded plot

//----------------------------------------All Unfold in one plot(Closure Test)
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt = 0 ; ipt <nHLTmx; ipt++){
          
	  TH1D *MyHist  = (TH1D*) Unfold[0][ity][ivar][ipt]->Clone();
    	  cout << MyHist->GetName() << endl;

	  Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:    H_{T,2} > %i %s",itypeN[ity], HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D *unfold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*) Psudo_Data_gen[ity][ivar][ipt]->Clone();
#else
          unfold_input[iout] = (TH1D*) MC_gen[umc][ity][ivar][ipt]->Clone();
#endif
          //unfold_input[iout] = (TH1D*) Unfold[iout][ity][ivar][ipt]->Clone();
          //MC_input[iout]->Rebin(2);
            Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout]; 
          //MCinput_index[iout]= closuretype[iout];
           MCinput_index[iout]= mcnamegen[iout];
      }

      for(int iout = 0 ; iout < unfold_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;
            cout << unfold_input[iout]->GetName()<< endl;
            cout <<  MCinput_index[iout] << endl;
          }
          data_index[0]= closuretype[0];
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.2,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.2,0.8};
          if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%sunfold_plot_%i.pdf(", ddtag[idd],1234); sprintf(pdfname1, "%sunfold_plot_%i.pdf" ,ddtag[idd],1234);sprintf(pdfname2, "%sunfold_plot_%i.pdf)",ddtag[idd],1234); 
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

      }  //end of phase space cut and variable loop
    }
  }


//----------------------------------------All Unfold in one plot(Closure Test) True distibution only for 1D and 2D TUNfolding--------------------------
if(idd>=1){
for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt = 0 ; ipt <nHLTmx; ipt++){

	  TH1D *MyHist  = (TH1D*) Ex_Unfold[0][ity][ivar][ipt]->Clone();    	  cout << MyHist->GetName() << endl;

	  Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:    H_{T,2} > %i %s",itypeN[ity], HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D *unfold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*) Ex_Psudo_Data_gen[ity][ivar][ipt]->Clone();
#else
          unfold_input[iout] = (TH1D*) Ex_MC_gen[umc][ity][ivar][ipt]->Clone();
#endif
          //unfold_input[iout] = (TH1D*) Unfold[iout][ity][ivar][ipt]->Clone();
          //MC_input[iout]->Rebin(2);
            Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout];
          //MCinput_index[iout]= closuretype[iout];

           MCinput_index[iout]= mcnamegen[iout];
      }

      for(int iout = 0 ; iout < unfold_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;
            cout << unfold_input[iout]->GetName()<< endl;
            cout <<  MCinput_index[iout] << endl;
          }

          data_index[0]= closuretype[0];
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.2,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.2,0.8};
          if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname,"%s_Trueunfold_plot_%i.pdf(", ddtag[idd],1234);sprintf(pdfname1,"%s_Trueunfold_plot_%i.pdf",ddtag[idd],1234);sprintf(pdfname2,"%s_Trueunfold_plot_%i.pdf)",ddtag[idd],1234);
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

      }  //end of phase space cut and variable loop
    }
  }
  }//if(idd>=1)

  //-----------------------------------------Refold comparisoin
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	for(int ipt = 0 ; ipt <nHLTmx; ipt++){
	  
	  TH1D *MyHist  = (TH1D*) Refold[iun][ity][ivar][ipt]->Clone();
	  Integralhist(MyHist);
	  //divBybinWidth(MyHist);
	  Myplotset(MyHist,0,0);
	  
	  if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity],  HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
	  else if(ipt==7){ sprintf(Title,"%s:      H_{T,2} > %i %s",itypeN[ity],  HT2range[ipt] ,"GeV/c" );}
	  sprintf(Yaxis," %s" ,Esvlogy[ivar]);
	  MyHist->SetTitle(Title);
	  MyHist->GetXaxis()->SetTitle("");
	  MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  
	  TH1D *MC_input[nmc];
	  const char *MCinput_index[nmc], *data_index[1];
	  for(int iout = 0 ; iout < nmc ; iout++){
	    MC_input[iout] = (TH1D*) MC_reco_fake[iout][ity][ivar][ipt]->Clone();
	   // MC_input[iout]->Rebin(2);
	    Integralhist(MC_input[iout]);
	   // divBybinWidth(MC_input[iout]);
	    MCinput_index[iout]= mcnamerco[iout]; }
	  
	  data_index[0]= RefoldEra[iun];
	  char lplot_xtitle[100];
	  sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
	  //float ratio_range1[2]={1.2,0.9};
	  int num1[2]={nmc,1};
	  float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  cpt0->cd();
	  cpt0->SetBorderSize(0);
	  cpt0->SetRightMargin(0.0);
	  cpt0->SetTopMargin(0.0);
	  cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  sprintf(pdfname, "%sRefold_plot_%i.pdf(" ,ddtag[idd],iun); sprintf(pdfname1, "%sRefold_plot_%i.pdf" ,ddtag[idd],iun);sprintf(pdfname2, "%sRefold_plot_%i.pdf)" ,ddtag[idd],iun); //TRefolded_typ_0_pt0_eta0_3
	  if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
	  }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
	  }else{cpt0->Print(pdfname,"pdf");};
	  
	}
      }  //end of phase space cut and variable loop
    }
  }//End of Refolded plot
  
  
 //Refold in one plot
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt = 0 ; ipt <nHLTmx; ipt++){
          
          TH1D *MyHist  = (TH1D*) MC_reco_fake[umc][ity][ivar][ipt]->Clone();
          Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:      H_{T,2} > %i %s",itypeN[ity] ,  HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);


          TH1D *Refold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
     
	   Refold_input[iout] = (TH1D*) Refold[iout][ity][ivar][ipt]->Clone();
           // MC_input[iout]->Rebin(2);
            Integralhist(Refold_input[iout]);
           // divBybinWidth(Refold_input[iout]);
            MCinput_index[iout]= RefoldEra[iout]; }

          data_index[0]= mcnamerco[umc];
          data_index[1]= "MC";
          data_index[2]= "Refolded";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.2,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%sRefold_plot_%i.pdf(",ddtag[idd], 1234); sprintf(pdfname1, "%sRefold_plot_%i.pdf" ,ddtag[idd],1234);sprintf(pdfname2, "%sRefold_plot_%i.pdf)",ddtag[idd],1234); //TRefolded_typ_0_pt0_eta0_3
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

      }  //end of phase space cut and variable loop
    }
  } 

 //------------------------------Refold True Dist in one plot
 if(idd>=1){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt = 0 ; ipt <nHLTmx; ipt++){

          TH1D *MyHist  = (TH1D*) Ex_Psudo_Data_reco_fake[ity][ivar][ipt]->Clone();
          Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:      H_{T,2} > %i %s",itypeN[ity] ,  HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);


          TH1D *Refold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
            Refold_input[iout] = (TH1D*) Ex_Refold[iout][ity][ivar][ipt]->Clone();
           // MC_input[iout]->Rebin(2);
            Integralhist(Refold_input[iout]);
           // divBybinWidth(Refold_input[iout]);
            MCinput_index[iout]= RefoldEra[iout]; }

          data_index[0]= mcnamerco[umc];
          data_index[1]= "MC";
          data_index[2]= "Refolded";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.2,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%s_TrueRefold_plot_%i.pdf(",ddtag[idd], 1234); sprintf(pdfname1, "%s_TrueRefold_plot_%i.pdf" ,ddtag[idd],1234);sprintf(pdfname2, "%s_TrueRefold_plot_%i.pdf)",ddtag[idd],1234); //TRefolded_typ_0_pt0_eta0_3
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

      }  //end of phase space cut and variable loop
    }
  }
}//if(idd>=1)




} //for(int idd=0; idd <ndd; idd++) 
  
}  // end of main program


//-------------------------------------------Ratio plot function

TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[1]){
  //Nplot[0] = number of MC enetered
  //Nplot[1] = place 1 if upper part is log scale needed
  //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  //plegend[4]= text size
  //plegend[5-6]= ratio plot axis range
  //data = data histogram
  // MC = monte carlo histogram array
  //legendN = name of the legends for mC one by one
  //lowpadx = x axis title of the ratio plot
  
  TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 575,600 );
  canvas->cd();
  canvas->SetRightMargin(0.02);
  canvas->SetTopMargin(0.02);
  
  char ratioXaxis1[100];
  char MCindex[100];
  //float ymax;
  // canvas->SetBottomMargin(0.1); 
  data->GetYaxis()->SetLabelSize(0.03);
  data->GetXaxis()->SetLabelSize(0.03);
  data->GetYaxis()->SetTitleSize(0.045);
  data->GetYaxis()->SetTitleOffset(1.0);
  data->GetXaxis()->SetTitleSize(0.055);
  data->GetXaxis()->SetTitleOffset(0.12);
  data->GetYaxis()->CenterTitle();
  data->GetXaxis()->CenterTitle();
  data->SetLineWidth(2);
  data->SetMarkerStyle(9);
  data->SetMarkerSize(.8);
  int ifont =42;
  data->GetXaxis()->SetTitleFont(ifont);
  data->GetYaxis()->SetTitleFont(ifont);     
  data->SetTitleFont(ifont);
  
  //    ymax = data->GetMaximum();
  //Devide the histogram with bin width
  
  
  TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.30, 1.0, 1.0);
  padfun1->SetBottomMargin(0.01); // Upper and lower plot are joined
  padfun1->SetTopMargin(0.05); // Upper and lowd
  padfun1->SetRightMargin(.04);
  padfun1->Draw();             // Draw the upper pad: pad1
  padfun1->cd();
  data->SetFillColor(kYellow);
  data->SetFillStyle(1111);
  
  //-----------------------------------create the legend of ht2 range 
  char MC_HTindex[100];
  TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);
  HT_range->SetTextSize(0.04);
  sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  HT_range->Draw();
  data->SetTitle("");
  //--------------------------------------end the legend of ht2 range
  data->Draw("e2");
  
  
  //   data->Draw(" ");
  
  
  if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  
  double chi2rat[Nplot[0]];    // for chi2 plot in legend 
  double chi2Ndfrat[Nplot[0]]; // for chi2 plot in legend
  
  int  color[40] = {2,4,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  int   style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
    MC[iup]->SetLineStyle(style[iup]);
    MC[iup]->SetLineColor(color[iup]);
    MC[iup]->SetLineWidth(2);
    //   MC[iup]->Draw("same hist ");// gPad->SetLogy();    
    
    //--------------------------------------------- Addition fo chi2/Ndf  with legend
    int nn = data->GetNbinsX();
    Double_t resl[nn] , chi2l;
    Int_t ndfl ,igoodl;
    TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test
    
    TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");   //for chi square test
    chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);
    
    // chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
    //cout << "Ndf value =" << ndfl << endl;
    //cout << "chi2 value=" << chi2l << endl;
    chi2rat[iup]=chi2l;
    chi2Ndfrat[iup]=chi2l/ndfl;
    //-----------------------------------------------end of chi2/Ndf
    //------------------------------------------------Divide the bin contents/error by bin width
    
    /*int tmpnbn = MC[iup]->GetNbinsX();
      for (int ix=0; ix<tmpnbn; ix++) {
      double awidth = MC[iup]->GetBinWidth(ix+1); // tmpwid;
      MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      double error = MC[iup]->GetBinError(ix+1);
      MC[iup]->SetBinError(ix+1, error/awidth);
      }*/
    //------------------------------------------------end Divide the bin contents/error by bin width
    
    MC[iup]->Draw("same hist e1 ");
    
  }//------------------- end of Montecarlo loop
  
  /*int tmpnbd = data->GetNbinsX();
    for (int ix=0; ix<tmpnbd; ix++) {
    double awidth = data->GetBinWidth(ix+1); // /tmpwid;
    data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    double error = data->GetBinError(ix+1);
    data->SetBinError(ix+1, error/awidth);
    }*/
  
  
  data->Draw(" same e1");
  
  //----------------------------------------------//Maximum Uncertainty with respect to Monash
  // make it off if not needed
  /*
    TH1D *MC_inputerr = (TH1F*)data->Clone("MC_inputerr");
    int MCbinnum = data->GetNbinsX();
    double rel_err[MCbinnum];
    for(int ix =0; ix <MCbinnum ; ix++){rel_err[ix]=0.0; }
    
    for(int iup =0; iup < Nplot[0] ; iup++){
    for (int ix=0; ix<MCbinnum; ix++) {
    double dbin = MC[0]->GetBinContent(ix+1);
    // double dbin = data->GetBinContent(ix+1);
    double ebin = MC[iup]->GetBinContent(ix+1);
    double rel_error=(fabs(dbin-ebin))/dbin;
    MC_inputerr->SetBinContent(ix+1,1);
    
    if(rel_err[ix+1] < rel_error){
    rel_err[ix+1]=rel_error;
    MC_inputerr->SetBinError(ix+1,rel_error);}
    else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);
    }
    } //loop over bin
    }//end of number monte carlo
    
    MC_inputerr->SetFillColorAlpha(30,0.050);
    MC_inputerr->SetMarkerStyle(1);
  */
  //------------------------------------------------------------------end of Uncertainty
  
  
  //  data->Draw("same");
  // TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
  TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
  legendn->SetFillStyle(0);
  legendn->SetBorderSize(0);
  legendn->SetTextSize(plegend[4]);
  legendn->SetTextFont(42);
  
  sprintf(MCindex,"%s" , datanm[0]);      // use if chi2 is not needed in the legend
  legendn->AddEntry(data, MCindex,"lp");
  for(int iup =0 ; iup <  Nplot[0] ; iup++){
    sprintf(MCindex,"%s" , modnam[iup]);      // use if chi2 is not needed in the legend
    //        sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
    legendn->AddEntry(MC[iup], MCindex ,"lp");
  }
  legendn->Draw();
  
  HT_range->Draw();
  
  //----------------------------------------ratio plot pad
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *padfun2 = new TPad("padfun2", "padfun2", 0, 0.02,1.0, 0.30);
  padfun2->SetTopMargin(0);
  padfun2->SetBottomMargin(.4);
  padfun2->SetRightMargin(.04);
  padfun2->SetGridy(); // Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  //   gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    TH1D *rh2;
    //if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    // if(  data->GetBinContent(10) > 0) {
    rh2 = (TH1D*)MC[ilow]->Clone("rh2"); 
    //   rh2->Sumw2();
    rh2->Divide(data);     //MC devide by data
    /* }else{
       rh2 = (TH1D*)data->Clone("rh2");
       rh2->Sumw2();
       rh2->SetLineColor(MC[ilow]->GetLineColor());
       rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       rh2->Divide(MC[ilow]);
       }*/
    rh2->SetTitle("");
    // rh2->SetLineColor(kBlack);
    // rh2->SetMaximum(lowYrange[0]);  // .. range
    // rh2->SetMinimum(lowYrange[1]);  // Define Y ..
    // cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;   
    rh2->SetMinimum(plegend[6]);  // Define Y ..
    rh2->SetMaximum(plegend[5]);  // .. range
    rh2->SetStats(0);      // No statistics on lower plot
    // rh2->Divide(data);     //MC devide by data
    // rh2->SetMarkerStyle(21);
    // rh2->Draw(" same e1");
    rh2->Draw("same e1");
    
    //   MC_inputerr->Draw("same E2");//------------------Draw the Uncertainty in Ratio plot
    
    rh2->GetXaxis()->SetTitleSize(0.13);
    //  rh2->GetXaxis()->SetTitleFont(43);
    rh2->GetXaxis()->SetTitleOffset(1.15);
    rh2->GetXaxis()->SetLabelSize(0.1);
    rh2->GetXaxis()->CenterTitle();
    sprintf(ratioXaxis1," %s" ,lowpadx);
    rh2->GetXaxis()->SetTitle(ratioXaxis1);
    
    rh2->GetYaxis()->SetTitle("MC/Data");
    rh2->GetYaxis()->CenterTitle();
    rh2->GetYaxis()->SetNdivisions(505);
    rh2->GetYaxis()->SetTitleSize(0.12);
    rh2->GetXaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleOffset(0.35);
    //  rh2->GetYaxis()->SetLabelFont(1.0); // Absolute font size in pixel (precision 3)
    rh2->GetYaxis()->SetLabelSize(0.09);
  }
  canvas->Update();
  return canvas;
}//------------------------------------end of  ratio plot function


void Integralhist(TH1D *hist){ hist->Scale(1/(hist->Integral()));}

void divBybinWidth(TH1D *hist){
  int tmpnbn = hist->GetNbinsX();
  for (int ix=0; ix<tmpnbn; ix++) {
    double awidth = hist->GetBinWidth(ix+1); // /tmpwid;
    hist->SetBinContent(ix+1, hist->GetBinContent(ix+1)/awidth);
    double error = hist->GetBinError(ix+1);
    hist->SetBinError(ix+1, error/awidth);
  }
}

void Myplotset(TH1D *MyHist, const char* XTitle, const char* YTitle){
 // int ifornt =102;
  MyHist->SetTitleOffset(0.4);
  //MyHist->SetTitleFont(ifornt);
  MyHist->SetTitleSize(0.02);
  
  MyHist->GetXaxis()->SetLabelSize(0.03);
  MyHist->GetXaxis()->SetTitleSize(0.045);
  MyHist->GetXaxis()->SetTitleOffset(1.0);
  //   MyHist->GetXaxis()->SetTitleFont(ifornt);
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetXaxis()->SetTitle(XTitle);
  
  MyHist->GetYaxis()->SetLabelSize(0.03);
  MyHist->GetYaxis()->SetTitleSize(0.040);
  MyHist->GetYaxis()->SetTitleOffset(1.0);
  MyHist->GetYaxis()->SetTitle(YTitle);
  //   MyHist->GetYaxis()->SetTitleFont(ifornt);     
  MyHist->GetYaxis()->CenterTitle();
    MyHist->SetTitle(""); 
  //  gStyle->SetTitleFontSize(.08);
  MyHist->SetLineWidth(2);
}

void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm){
         cpt->SetBorderSize(bs);
          cpt->SetLeftMargin(lm);
          cpt->SetRightMargin(rm);
          cpt->SetTopMargin(tm); 
          cpt->SetBottomMargin(bm); 
   }

void CTLegend(TLegend *legendn, const char* txt1, const char* txt2){ 
  legendn->SetFillStyle(0);
  legendn->SetBorderSize(0);
  legendn->SetTextSize(0.03);
  legendn->SetTextFont(42);
  legendn->AddEntry((TObject*)0, txt1, "");
  legendn->AddEntry((TObject*)0, txt2, "");
  
}

void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3] ){
  MyHist->SetTitleOffset(0.2);
//  MyHist->SetTitleFont(102);
  MyHist->SetTitleSize(0.02);

  MyHist->GetXaxis()->SetLabelSize(0.03);
  MyHist->GetXaxis()->SetTitleSize(titsize[0]);
  MyHist->GetXaxis()->SetTitleOffset(titoff[0]);
//  MyHist->GetXaxis()->SetTitleFont(102);
//  MyHist->GetXaxis()->CenterTitle();

  MyHist->GetYaxis()->SetLabelSize(0.03);
  MyHist->GetYaxis()->SetTitleSize(titsize[1]);
  MyHist->GetYaxis()->SetTitleOffset(titoff[1]);

  MyHist->GetZaxis()->SetLabelSize(0.03);
  MyHist->GetZaxis()->SetTitleSize(titsize[2]);
  MyHist->GetZaxis()->SetTitleOffset(titoff[2]);

  MyHist->SetTitle("");
  MyHist->GetXaxis()->SetTitle(XTitle);
  MyHist->GetYaxis()->SetTitle(YTitle);
  MyHist->GetZaxis()->SetTitle(ZTitle);
}


TCanvas *ratio_can1(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
  //Nplot[0] = number of MC enetered
  //Nplot[1] = place 1 if upper part is log scale needed
  //Nplot[2] = place 1 if MC/Data or 0 if data/MC
  //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  //plegend[4]= text size
  //plegend[5-6]= ratio plot axis range
  //data = data histogram
  // MC = monte carlo histogram array
  //legendN = name of the legends for mC one by one
  //lowpadx = x axis title of the ratio plot
  //datanm[0] = name of data (or MC in case one MC is used)
  //datanm[1] = name of MC/data Data term
  //datanm[2] = name of Mc/Data MC term
  bool isstat =1;
  TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 575,600 );
  canvas->cd();
  canvas->SetRightMargin(0.02);
  canvas->SetTopMargin(0.02);

  char ratioXaxis1[100]; char ratioYaxis[100]; char MCindex[100];
  //float ymax;
  // canvas->SetBottomMargin(0.1);
  data->GetYaxis()->SetLabelSize(0.03);
  data->GetXaxis()->SetLabelSize(0.03);
  data->GetYaxis()->SetTitleSize(0.045);
  data->GetYaxis()->SetTitleOffset(1.0);
  data->GetXaxis()->SetTitleSize(0.055);
  data->GetXaxis()->SetTitleOffset(0.12);
  data->GetYaxis()->CenterTitle();
  data->GetXaxis()->CenterTitle();
  data->SetLineWidth(2);
  data->SetMarkerStyle(9);
  data->SetMarkerSize(.8);
  int ifont =42;
  data->GetXaxis()->SetTitleFont(ifont);
  data->GetYaxis()->SetTitleFont(ifont);
  data->SetTitleFont(ifont);

  //    ymax = data->GetMaximum();
  //Devide the histogram with bin width


  TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.30, 1.0, 1.0);
  padfun1->SetBottomMargin(0.01); // Upper and lower plot are joined
  padfun1->SetTopMargin(0.05); // Upper and lowd
  padfun1->SetRightMargin(.04);
  padfun1->Draw();             // Draw the upper pad: pad1
  padfun1->cd();
  data->SetFillColor(kYellow);
  data->SetFillStyle(1111);

  //-----------------------------------create the legend of ht2 range
  char MC_HTindex[100];
  TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);
  HT_range->SetTextSize(0.04);
  sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  HT_range->Draw();
  data->SetTitle("");
  //--------------------------------------end the legend of ht2 range
  data->Draw("");


  //   data->Draw(" ");


  if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();

  double chi2rat[Nplot[0]];    // for chi2 plot in legend
  double chi2Ndfrat[Nplot[0]]; // for chi2 plot in legend

  int  color[40] = {2,4,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  int   style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
    MC[iup]->SetLineStyle(style[iup]);
    MC[iup]->SetLineColor(color[iup]);
    MC[iup]->SetLineWidth(2);
    //   MC[iup]->Draw("same hist ");// gPad->SetLogy();

    //--------------------------------------------- Addition fo chi2/Ndf  with legend
    int nn = data->GetNbinsX();
    Double_t resl[nn] , chi2l;
    Int_t ndfl ,igoodl;
    TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test

    TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");   //for chi square test
    chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);

    // chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
    //cout << "Ndf value =" << ndfl << endl;
    //cout << "chi2 value=" << chi2l << endl;
    chi2rat[iup]=chi2l;
    chi2Ndfrat[iup]=chi2l/ndfl;
    //-----------------------------------------------end of chi2/Ndf
    //------------------------------------------------Divide the bin contents/error by bin width

    /*int tmpnbn = MC[iup]->GetNbinsX();
      for (int ix=0; ix<tmpnbn; ix++) {
      double awidth = MC[iup]->GetBinWidth(ix+1); // tmpwid;
      MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      double error = MC[iup]->GetBinError(ix+1);
      MC[iup]->SetBinError(ix+1, error/awidth);
      }*/
    //------------------------------------------------end Divide the bin contents/error by bin width

   // MC[iup]->Draw("same hist e1 ");
    MC[iup]->Draw("same hist");

  }//------------------- end of Montecarlo loop

  /*int tmpnbd = data->GetNbinsX();
    for (int ix=0; ix<tmpnbd; ix++) {
    double awidth = data->GetBinWidth(ix+1); // /tmpwid;
    data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    double error = data->GetBinError(ix+1);
    data->SetBinError(ix+1, error/awidth);
    }*/


  //data->Draw(" same e1");
  data->Draw(" same ");

 //------------------------------------------Stat Error(DATA)
   TH1D *staterr = (TH1D*)data->Clone("staterr");
    staterr->Reset();

    for(int ierr=1; ierr < staterr->GetNbinsX()+1; ierr++){
     double bincont = data->GetBinContent(ierr);
     double serr = data->GetBinError(ierr);

     double error1 =sqrt(pow((serr/bincont),2.));
    staterr->SetBinContent(ierr,1.0);
    staterr->SetBinError(ierr,error1);

    }
    gStyle->SetErrorX(0.5);
    gStyle->IsTransparent();
    staterr->SetFillColor(5);
    staterr->SetMarkerSize(0.1);
    staterr->SetFillStyle(3002);
  //  staterr->SetFillStyle(1111);
 //   staterr->SetFillColor(21);
  //  staterr->SetMarkerSize(0.1);



  //----------------------------------------------//Maximum Uncertainty with respect to Monash
  // make it off if not needed
  /*
    TH1D *MC_inputerr = (TH1F*)data->Clone("MC_inputerr");
    int MCbinnum = data->GetNbinsX();
    double rel_err[MCbinnum];
    for(int ix =0; ix <MCbinnum ; ix++){rel_err[ix]=0.0; }

    for(int iup =0; iup < Nplot[0] ; iup++){
    for (int ix=0; ix<MCbinnum; ix++) {
    double dbin = MC[0]->GetBinContent(ix+1);
    // double dbin = data->GetBinContent(ix+1);
    double ebin = MC[iup]->GetBinContent(ix+1);
    double rel_error=(fabs(dbin-ebin))/dbin;
    MC_inputerr->SetBinContent(ix+1,1);

    if(rel_err[ix+1] < rel_error){
    rel_err[ix+1]=rel_error;
    MC_inputerr->SetBinError(ix+1,rel_error);}
    else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);
    }
    } //loop over bin
    }//end of number monte carlo

    MC_inputerr->SetFillColorAlpha(30,0.050);
    MC_inputerr->SetMarkerStyle(1);
  */
  //------------------------------------------------------------------end of Uncertainty


  //  data->Draw("same");
  // TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
  TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
  legendn->SetFillStyle(0);
  legendn->SetBorderSize(0);
  legendn->SetTextSize(plegend[4]);
  legendn->SetTextFont(42);

  sprintf(MCindex,"%s" , datanm[0]);      // use if chi2 is not needed in the legend
  legendn->AddEntry(data, MCindex,"lp");
  for(int iup =0 ; iup <  Nplot[0] ; iup++){
    sprintf(MCindex,"%s" , modnam[iup]);      // use if chi2 is not needed in the legend
    //        sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
    legendn->AddEntry(MC[iup], MCindex ,"lp");
  }
  if(isstat)legendn->AddEntry(staterr, "Statistical Uncertainty","f");
  
  legendn->Draw();

  HT_range->Draw();



  //----------------------------------------ratio plot pad
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *padfun2 = new TPad("padfun2", "padfun2", 0, 0.02,1.0, 0.30);
  padfun2->SetTopMargin(0);
  padfun2->SetBottomMargin(.4);
  padfun2->SetRightMargin(.04);
  padfun2->SetGridy(); // Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  //   gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    TH1D *rh2;
    //if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    // if(  data->GetBinContent(10) > 0) {
    rh2 = (TH1D*)MC[ilow]->Clone("rh2");
    //   rh2->Sumw2();
    rh2->Divide(data);     //MC devide by data
    /* }else{
       rh2 = (TH1D*)data->Clone("rh2");
       rh2->Sumw2();
       rh2->SetLineColor(MC[ilow]->GetLineColor());
       rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       rh2->Divide(MC[ilow]);
       }*/
    rh2->SetTitle("");
    // rh2->SetLineColor(kBlack);
    // rh2->SetMaximum(lowYrange[0]);  // .. range
    // rh2->SetMinimum(lowYrange[1]);  // Define Y ..
    // cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;
    rh2->SetMinimum(plegend[6]);  // Define Y ..
    rh2->SetMaximum(plegend[5]);  // .. range
    rh2->SetStats(0);      // No statistics on lower plot
//    gStyle->SetStatColor(2);      // No statistics on lower plot
  //  gStyle->SetStatX(0.7);;      // No statistics on lower plot
    // rh2->Divide(data);     //MC devide by data
    // rh2->SetMarkerStyle(21);
    // rh2->Draw(" same e1");
    rh2->Draw("Same P text ");
  
   //-----------------------------------------------------Stat Error
   TH1D *mstaterr = (TH1D*)MC[ilow]->Clone("staterr");
    mstaterr->Reset();

    for(int ierr=1; ierr < mstaterr->GetNbinsX()+1; ierr++){
     double bincont = MC[ilow]->GetBinContent(ierr);
     double serr = MC[ilow]->GetBinError(ierr);

     double error1 =sqrt(pow((serr/bincont),2.));
    mstaterr->SetBinContent(ierr,1.0);
    mstaterr->SetBinError(ierr,error1);

    }
    gStyle->SetErrorX(0.5);
    mstaterr->SetFillColor(5);
    mstaterr->SetMarkerSize(0.1);
    mstaterr->SetFillStyle(1111);
    mstaterr->SetFillStyle(3002);
//    staterr->SetFillColor(21);
//    staterr->SetMarkerSize(0.1);

   //-----------------------------------------------------Stat Error
    
    staterr->Draw("E2 Same");
    //mstaterr->Draw("E2 Same");
    
    //   MC_inputerr->Draw("same E2");//------------------Draw the Uncertainty in Ratio plot

    rh2->GetXaxis()->SetTitleSize(0.13);
    //  rh2->GetXaxis()->SetTitleFont(43);
    rh2->GetXaxis()->SetTitleOffset(1.15);
    rh2->GetXaxis()->SetLabelSize(0.1);
    rh2->GetXaxis()->CenterTitle();
    sprintf(ratioXaxis1," %s" ,lowpadx);
    rh2->GetXaxis()->SetTitle(ratioXaxis1);

    //rh2->GetYaxis()->SetTitle("MC/Data");
    if(Nplot[2]==1){ sprintf(ratioYaxis,"%s/%s" ,datanm[1],datanm[2]);}
    if(Nplot[2]==0){ sprintf(ratioYaxis,"%s/%s" ,datanm[2],datanm[1]);}
    rh2->GetYaxis()->SetTitle(ratioYaxis);
    rh2->GetYaxis()->CenterTitle();
    rh2->GetYaxis()->SetNdivisions(505);
    rh2->GetYaxis()->SetTitleSize(0.1);
    rh2->GetXaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleOffset(0.40);
    //  rh2->GetYaxis()->SetLabelFont(1.0); // Absolute font size in pixel (precision 3)
    rh2->GetYaxis()->SetLabelSize(0.09);
  }
  canvas->Update();
  return canvas;
}//------------------------------------end of  ratio plot function


