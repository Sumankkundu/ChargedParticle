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
//#define CLOUSER
#define TRUEAXIS
void Unfoldplot(){
  
  static const int unfold_ty =1; //Unfold method to be plots 
  static const int clos_ty =4; //Unfold method to be plots 
  static const int nHLTmx=8; //HT2 Range
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =2;
  static const int nmc = 4;
  static const int ndd = 3; //0 : root Hist, 1: 1D , 2: 2D
  static const int umc = 0; //0 for Py8, 1 for MG , 2 for Herwig : which MC have used in Unfolding

  bool isstat =1;  int irbin=1;

  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={2,4,6,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[nusedvar]={3,9,15,18,24};
  Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000};

  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  
  string  Datadir = "Data2D";
  string folddir = "Folded2D";
  string unfdir = "Unfold2D";
  const char*  dirname[nmc]={"Pythia8","Py8Flat","MG8","HW7"};
  string mcdir[nmc]={"Pythia8","Py8Flat","MG8","HW7"};
  
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[8]={"83 < H_{T,2} < 109", "109 < H_{T,2} < 172", "172 < H_{T,2} < 241", "241 < H_{T,2} < 309","309 < H_{T,2} < 377", "377 < H_{T,2} < 462", "462 < H_{T,2} <570","H_{T,2} > 570"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};

  const  char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const  char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};
  const  char* mcname[3]={"Pythia8 CP5 Tune","Madgraph","Herwig++"};
  const  char* mcnamerco[4]={"Pythia8 RECO", "PY8 Flat","Madgraph RECO","Herwig++ RECO"};
  const  char* closuretype[4]={"PY8 by PY8","PY8 Flat By PY8 ", "MG By PY8","HW7 By PY8"};
  const  char* itypeN[ntype]={"Jets","Charged Particles"}; 
  const  char* DataEra[3]={"Data RECO","Data RECO","Data RECO"};
  const  char* UndoldEra[3]={"Unfold 2016","Unfold 2017","Unfold 2018"};
  const  char* RefoldEra[5]={"Refold Pythia8(No Regularisation) ","Refold Pythia8(L-Curve scan)","Refold Pythia8(Scan Tau)","Refold Pythia8(ScanSURE)","Refold Pythia8(Iterative method)"};
  const  char* Unfoldtype[5]={"TUnfold(No Reg) ","TUnfold(L-Curve scan)","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
  const  char* Modelnm[3]={"Pythia8","Madgraph","Herwig"};
  const  char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const  char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  static const int iera = 1;
  int iPeriod = 0;  int iPos=10 ;
 
  TFile *Unf_root[clos_ty];
  //******************************************************************************
//  TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/MC_Unfolded_Result.root");  // Unfolded data 
  TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/Unfolded_Result.root");  // Unfolded data 
  
  Unf_root[0] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/Unfolded_Result_PY8_PY8.root");  // Py- py
  Unf_root[1] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/Unfolded_Result_PY8_PY8Flat.root");  // Py -Py flat 
  Unf_root[2]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/Unfolded_Result_PY8_MG.root");  // Py MG 
  Unf_root[3]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/Unfolded_Result_PY8_Hw7.root");  // Py MG 
  //--------------------------------------Function declearation
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_canV2(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
  void HT2_Normal(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  void HT2_NormalV2(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  TH1D* ReadHist1D(string name,TFile* root, int irbin=1);
  TH2D* ReadHist2D(string name,TFile* root, int irbin=1);
  TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4="");


  TH1D *MC_gen[nmc][ntype][nusedvar];   //MC gen
  TH1D *MC_gen_miss[nmc][ntype][nusedvar]; //MC Gen-miss
  TH1D *Ex_MC_gen[nmc][ntype][nusedvar]; //Truth Dist Gen
  TH1D *Ex_MC_gen_miss[nmc][ntype][nusedvar];  //Truth Gen-miss
  


  TH1D *MC_reco[nmc][ntype][nusedvar];   //MC Reco 
  TH1D *MC_reco_fake[nmc][ntype][nusedvar]; //MC Reco -Fake
  TH1D *Ex_MC_reco[nmc][ntype][nusedvar];   //Truth Dist Reco
  TH1D *Ex_MC_reco_fake[nmc][ntype][nusedvar];  //Truth Reco -Fake
  
  TH2D *MC_Res[nmc][ntype][nusedvar];   //RM

  TH1D* MC_GenHT2[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_RecoHT2[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_FakeHT2[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_missHT2[nmc][ntype][nusedvar][nHLTmx];
#ifdef CLOUSER  
  TH1D *Psudo_Data_gen[ntype][nusedvar];  
  TH1D *Psudo_Data_reco_fake[ntype][nusedvar];  
  TH1D *Ex_Psudo_Data_gen[ntype][nusedvar];
  TH1D *Ex_Psudo_Data_gen_miss[ntype][nusedvar];
  TH1D *Ex_Psudo_Data_reco_fake[ntype][nusedvar];
#endif
  TH1D *Data_reco[ntype][nusedvar];
  TH1D *Ex_Data_reco[ntype][nusedvar];
  
  TH1D *hist_eff[nmc][ntype][nusedvar];
  TH1D *hist_fake[nmc][ntype][nusedvar];
  TH1D *hist_purity[nmc][ntype][nusedvar];
  TH1D *hist_stbl[nmc][ntype][nusedvar];
  TH1D *Ex_hist_eff[nmc][ntype][nusedvar];
  TH1D *Ex_hist_fake[nmc][ntype][nusedvar];
  TH1D *Ex_hist_purity[nmc][ntype][nusedvar];
  TH1D *Ex_hist_stbl[nmc][ntype][nusedvar];

  TH1D *Unfold[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *Refold[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *Ex_Unfold[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *Ex_Refold[clos_ty][unfold_ty][ntype][nusedvar];

  TH2D *Ex_Unfold_HT[clos_ty][unfold_ty][ntype][nusedvar];
  TH2D *Ex_Refold_HT[clos_ty][unfold_ty][ntype][nusedvar];
  
  TH1D* Unf_HT2[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];
  TH1D* Refold_HT2[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];


  TH2D *Ex_MC_reco_HT[nmc][ntype][nusedvar];
  TH2D *Ex_MC_gen_HT[nmc][ntype][nusedvar];

  TH2D *Corr[clos_ty][unfold_ty][ntype][nusedvar];
  TH2D *Prob[clos_ty][unfold_ty][ntype][nusedvar];
  TH2D *Ematrix[clos_ty][unfold_ty][ntype][nusedvar];

  TH1D *BLT_reco[ntype][nusedvar];
  TH1D *BLT_gen[ntype][nusedvar];
  TH2D *Ex_BLT_reco[ntype][nusedvar];
  TH2D *Ex_BLT_gen[ntype][nusedvar];
  TH1D *BLT_reco_HT2[ntype][nusedvar][nHLTmx];
  TH1D *BLT_gen_HT2[ntype][nusedvar][nHLTmx];
  TH1D* BLT_RecoHT2[nHLTmx];

for(int ity=0; ity <ntype; ity++){
      for(int ivar =0;ivar < nusedvar; ivar++){
          Data_reco[ity][ivar] = (TH1D*)ReadHist1D(Datadir+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          BLT_reco[ity][ivar] =(TH1D*)ReadHist1D(unfdir+"/BLTDet_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot); 
#ifdef TRUEAXIS
          Ex_Data_reco[ity][ivar] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unfoldroot);
	  Ex_BLT_reco[ity][ivar] =(TH2D*)ReadHist2D(unfdir+"/EBLTDet_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          HT2_NormalV2(Ex_BLT_reco[ity][ivar], nHLTmx, HT2range, BLT_RecoHT2);
          for(int ipt=0; ipt < nHLTmx; ipt++){
          BLT_reco_HT2[ity][ivar][ipt] = (TH1D*)BLT_RecoHT2[ipt]->Clone();
               }
	   Ex_BLT_gen[ity][ivar] =(TH2D*)ReadHist2D(unfdir+"/EBLTgen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);

      TH1D* BLT_GenHT2[nHLTmx];
      HT2_NormalV2(Ex_BLT_gen[ity][ivar], nHLTmx, HT2range, BLT_GenHT2);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      BLT_gen_HT2[ity][ivar][ipt] = (TH1D*)BLT_GenHT2[ipt]->Clone();
 }

#endif
#ifdef CLOUSER
          Psudo_Data_gen[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot); 
          BLT_gen[ity][ivar] =(TH1D*)ReadHist1D(unfdir+"/BLTgen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot); 
#ifdef TRUEAXIS
          Ex_Psudo_Data_gen[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_Psudo_Data_reco_fake[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/Edd_DataRecominusfake_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot); 
          Ex_Psudo_Data_gen_miss[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/Edd_DataGenminusmiss_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);

          Ex_BLT_gen[ity][ivar] =(TH2D*)ReadHist2D(unfdir+"/EBLTgen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);

      TH1D* BLT_GenHT2[nHLTmx];
      HT2_NormalV2(Ex_BLT_gen[ity][ivar], nHLTmx, HT2range, BLT_GenHT2);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      BLT_gen_HT2[ity][ivar][ipt] = (TH1D*)BLT_GenHT2[ipt]->Clone();
               }
#endif
#endif
      }
    }

  for(int  imc =0; imc < nmc ; imc++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	  MC_gen[imc][ity][ivar] =(TH1D*)ReadHist1D(mcdir[imc]+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);  
	  
	  Ex_MC_gen_HT[imc][ity][ivar] = (TH2D*)ReadHist2D(mcdir[imc]+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
	  MC_reco[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          MC_gen_miss[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Genminusmiss_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          MC_reco_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Recominusfake_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          
	  hist_eff[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_miss_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_fake_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_purity[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Purity_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_stbl[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_stability_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
#ifdef TRUEAXIS
          Ex_MC_gen[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_MC_reco[imc][ity][ivar] =(TH1D*)ReadHist1D(mcdir[imc]+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_MC_gen_miss[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Genminusmiss_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_MC_reco_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Recominusfake_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          
	  Ex_hist_eff[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_miss_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_hist_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_fake_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_hist_purity[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Purity_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
          Ex_hist_stbl[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_stability_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);
#endif
	  //-------------------------------------Read Reso--------------------------------------
	  MC_Res[imc][ity][ivar] = (TH2D*)ReadHist2D(mcdir[imc]+"/dd_corr_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);

      TH1D* gen_var[nHLTmx];    //  TH1D* reco_var[nHLTmx];
      //HT2_NormalV2(Ex_Refold_HTimc][ity][ivar], nHLTmx, HT2range, reco_var);
      HT2_NormalV2(Ex_MC_gen_HT[imc][ity][ivar], nHLTmx, HT2range, gen_var);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      //MC_RecoHT2[imc][ity][ivar][ipt] = (TH1D*) rco_var[ipt]->Clone();
      MC_GenHT2[imc][ity][ivar][ipt] = (TH1D*)gen_var[ipt]->Clone();
       }
    //----------------------------------------------------------------
      }
    }
  }  //end of one MCinput root file  reading
  //---------------------------------------------------------//Read Unfolded
  for(int icl=0; icl < clos_ty; icl++){
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	  Unfold[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
	  Refold[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
	  Ex_Unfold_HT[icl][iun][ity][ivar] = (TH2D*)ReadHist2D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
	  Ex_Refold_HT[icl][iun][ity][ivar] = (TH2D*)ReadHist2D(unfdir+"/Edd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
#ifdef TRUEAXIS
         Ex_Unfold[icl][iun][ity][ivar] =(TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unf_root[icl]); 
	 Ex_Refold[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unf_root[icl]);
#endif
         Corr[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
         Ematrix[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
	 Prob[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_Prob_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);

      TH1D* unf_var[nHLTmx];      TH1D* Refold_var[nHLTmx];
      HT2_NormalV2(Ex_Refold_HT[icl][iun][ity][ivar], nHLTmx, HT2range, Refold_var);
      HT2_NormalV2(Ex_Unfold_HT[icl][iun][ity][ivar], nHLTmx, HT2range, unf_var);
      for(int ipt=0; ipt < nHLTmx; ipt++){
      Unf_HT2[icl][iun][ity][ivar][ipt] = (TH1D*) unf_var[ipt]->Clone();
      Refold_HT2[icl][iun][ity][ivar][ipt] = (TH1D*)Refold_var[ipt]->Clone();
      }
//---------------------------------------------------------
      }
     } 
    }
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
	TH1D *MyHist  = (TH1D*) Data_reco[ity][ivar]->Clone();
	Integralhist(MyHist);
	//divBybinWidth(MyHist);
	Myplotset(MyHist,0,0);
	
	sprintf(Title,"%s:   ",itypeN[ity]);
	sprintf(Yaxis," %s" ,Esvlogy[ivar]);
	MyHist->SetTitle(Title);
	MyHist->GetXaxis()->SetTitle("");
	MyHist->GetYaxis()->SetTitle(Yaxis);
	
	TH1D *MC_input[nmc];
	const char *MCinput_index[nmc+1];
	const char *data_index[1];
	for(int iout = 0 ; iout < nmc ; iout++){
	  MC_input[iout] = (TH1D*) MC_reco[iout][ity][ivar]->Clone(); 
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
	
	sprintf(pdfname, "%sRecoEVS_Plot.pdf(","dd_"); sprintf(pdfname1, "%sRecoEVS_Plot.pdf","dd_"); sprintf(pdfname2, "%sRecoEVS_Plot.pdf)","dd_"); 
	if(ity==0 && ivar==0){cpt0->Print(pdfname,"pdf");
	}else if(ity==1 && ivar==4 ) {cpt0->Print(pdfname2,"pdf");
	}else{  cpt0->Print(pdfname1,"pdf");};
    }  
  }//-------------------end of RECO PLOT    
//Response Matrix------------------------------------------------------
for(int  imc =0; imc < nmc ; imc++){
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
	cpt8->cd();
//	MC_Res[imc][ity][ivar]->RebinY(irbin);
        char lplot_xtitle[100]; char lplot_ytitle[100];
        sprintf(lplot_xtitle, "Reco bin ID ( %s, H_{T,2})",Esvlogx[ivar]); sprintf(lplot_ytitle, "Gen  bin ID ( %s, H_{T,2}) ",Esvlogx[ivar]);	
	double titoff1[3]={1.2,1.3,1.0};
        double titsize1[3] ={0.035,0.035,0.035};
	
	SetMycanvas(cpt8,0,0.1,0.15,0.05,0.1);
        Set2dHist( MC_Res[imc][ity][ivar],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
        gPad->SetLogz();
	MC_Res[imc][ity][ivar]->Draw("colz");
        
         
        sprintf(Title,"%s:   ",itypeN[ity]);
	TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
	CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,itypeN[ity] , "");leg1->SetTextColor(-8);leg1->Draw();
        CMS_lumi( cpt8, iPeriod, iPos ); cpt8->Update();
	sprintf(pdfname, "%sResponse_Mat%i.pdf(","dd_",imc); sprintf(pdfname1, "%sResponse_Mat%i.pdf","dd_",imc); sprintf(pdfname2, "%sResponse_Mat%i.pdf)","dd_",imc);
        if(ity==0 && ivar==0 ){cpt8->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4) {cpt8->Print(pdfname2,"pdf");
        }else{  cpt8->Print(pdfname1,"pdf");};
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

   	TLegend *leg2 = new TLegend(0.4,0.6,0.7,0.9);
	CTLegend(leg2," ","");
	TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        CTLegend(leg1,"", itypeN[ity]); 

        cpt1->cd();
	
        sprintf(Title,"%s:   ",itypeN[ity]);

        for (int i = 1; i <= hist_eff[umc][ity][ivar]->GetNbinsX(); ++i) {
         double content = 1 - hist_eff[umc][ity][ivar]->GetBinContent(i);
         hist_eff[umc][ity][ivar]->SetBinContent(i, content);
       }
         hist_eff[umc][ity][ivar]->SetMinimum(-0.01); hist_eff[umc][ity][ivar]->SetMaximum(1.01);
         hist_fake[umc][ity][ivar]->SetMinimum(-0.01); hist_fake[umc][ity][ivar]->SetMaximum(1.01);

        SetMycanvas(cpt1,0,0.1,0.15,0.05,0.12);
        Myplotset(hist_eff[umc][ity][ivar],Esvlogx[ivar],"Efficiency");
        hist_eff[umc][ity][ivar]->SetLineColor(color[0]);
	hist_eff[umc][ity][ivar]->Draw("same hist "); leg1->Draw();
        leg2->AddEntry(hist_eff[0][ity][ivar], Title ,"lp");

        cpt2->cd();
       //SetMycanvas(cpt2);
        SetMycanvas(cpt2,0,0.1,0.15,0.05,0.12);
 	Myplotset(hist_purity[umc][ity][ivar],Esvlogx[ivar],"Purity");
 	hist_purity[umc][ity][ivar]->SetLineColor(color[0]);
        hist_purity[umc][ity][ivar]->Draw("same hist"); leg1->Draw();
	cpt2->Update();

        cpt3->cd();
        SetMycanvas(cpt3,0,0.1,0.1,0.05,0.12);
        Myplotset(hist_fake[umc][ity][ivar],Esvlogx[ivar],"Fake rate");
        hist_fake[umc][ity][ivar]->SetLineColor(color[0]);
        hist_fake[umc][ity][ivar]->Draw("same hist");leg1->Draw();
         cpt3->Update();
	cpt4->cd();
        SetMycanvas(cpt4,0,0.1,0.1,0.05,0.12);
        Myplotset(hist_stbl[umc][ity][ivar],Esvlogx[ivar],"Stability");
        hist_stbl[umc][ity][ivar]->SetLineColor(color[0]); 
        hist_stbl[umc][ity][ivar]->Draw("same hist"); leg1->Draw();
   	 cpt4->Update();
    sprintf(pdfname, "%seffi_plot.pdf(","dd_"); sprintf(pdfname1, "%seffi_plot.pdf","dd_"); sprintf(pdfname2, "%seffi_plot.pdf)","dd_");
        if(ity==0 && ivar==0 ){cpt1->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt1->Print(pdfname2,"pdf");
        }else{  cpt1->Print(pdfname1,"pdf");};
	cpt1->Clear();
    sprintf(pdfname, "%spuri_plot.pdf(","dd_"); sprintf(pdfname1, "%spuri_plot.pdf","dd_"); sprintf(pdfname2, "%spuri_plot.pdf)","dd_"); 
        if(ity==0 && ivar==0 ){cpt2->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt2->Print(pdfname2,"pdf");
        }else{  cpt2->Print(pdfname1,"pdf");};
	cpt2->Clear();
     sprintf(pdfname, "%sfake_plot.pdf(","dd_"); sprintf(pdfname1, "%sfake_plot.pdf","dd_"); sprintf(pdfname2, "%sfake_plot.pdf)","dd_"); 
        if(ity==0 && ivar==0 ){cpt3->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4) {cpt3->Print(pdfname2,"pdf");
        }else{  cpt3->Print(pdfname1,"pdf");};
	cpt3->Clear();
     sprintf(pdfname, "%sstab_plot.pdf(","dd_"); sprintf(pdfname1, "%sstab_plot.pdf","dd_"); sprintf(pdfname2, "%sstab_plot.pdf)","dd_"); 
        if(ity==0 && ivar==0 ){cpt4->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt4->Print(pdfname2,"pdf");
        }else{  cpt4->Print(pdfname1,"pdf");};
	cpt4->Clear();
    
     }
  }   

cout << "Effi, Stabilty, fake, purity " <<endl;
/*  
  //Unfold comparisoin
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	  
	  TH1D *MyHist  = (TH1D*)Unfold[iun][ity][ivar]->Clone();
	  Integralhist(MyHist);
	  //divBybinWidth(MyHist);
	  Myplotset(MyHist,0,0);
	  sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          sprintf(Title,"%s:   ",itypeN[ity]);
	  MyHist->SetTitle(Title);
	  MyHist->GetXaxis()->SetTitle("");
	  MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  TH1D *MC_input[nmc];
	  const char *MCinput_index[nmc], *data_index[1];
	  for(int iout = 0 ; iout < nmc ; iout++){
	  MC_input[iout] = (TH1D*) MC_gen[iout][ity][ivar]->Clone();
	  Integralhist(MC_input[iout]);
	  //divBybinWidth(MC_input[iout]);
	    
	  MCinput_index[iout]= closuretype[iout]; }
	 
       data_index[0]= Unfoldtype[iun]; 
       //data_index[0]= UndoldEra[1]; 
       char lplot_xtitle[100];
       sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
       //float ratio_range1[2]={1.2,0.9};
       int num1[2]={nmc,1} ;
       float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
	  
       cout <<"OK1" <<endl;
       //cpt0->cd();
       //SetMycanvas(cpt0,0,0,0,0,0);
       cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
           sprintf(pdfname, "%sTUnfold_plot_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_plot_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_plot_%i.pdf)" ,"dd_",iun);
	  if(ity==0 && ivar==0 ){cpt0->Print(pdfname,"pdf");
	 }else if(ity==1 && ivar==4) {cpt0->Print(pdfname2,"pdf");
	  }else{cpt0->Print(pdfname,"pdf");};
	 
       cpt5->cd();
       double titoff1[3]={1.2,1.3,1.0};
       double titsize1[3] ={0.035,0.035,0.035};
       SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       gStyle->SetPaintTextFormat( "4.2f");
       Set2dHist(Corr[iun][ity][ivar],lplot_xtitle, lplot_xtitle,"correlation coefficients", titoff1, titsize1);
       Corr[iun][ity][ivar]->Draw("colz ");
       //Corr[iun][ity][ivar]->Draw("colz text");
       sprintf(Yaxis," %s" ,Esvlogy[ivar]);
       cout <<"OK1" <<endl;
        TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
        CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,itypeN[ity] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       cout <<"OK1" <<endl;
       CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       sprintf(pdfname, "%sTUnfold_corr_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_corr_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_corr_%i.pdf)" ,"dd_",iun); 
          if(ity==0 && ivar==0 ){cpt5->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 ) {cpt5->Print(pdfname2,"pdf");
          }else{cpt5->Print(pdfname,"pdf");};

       cpt6->cd();
       SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       Set2dHist(Prob[iun][ity][ivar],lplot_xtitle, lplot_xtitle,"Probability",titoff1, titsize1);
       gPad->SetLogz();Prob[iun][ity][ivar]->SetMinimum(1e-6); Prob[iun][ity][ivar]->SetMaximum(1);
       Prob[iun][ity][ivar]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       sprintf(pdfname, "%sTUnfold_prob_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_prob_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_prob_%i.pdf)" ,"dd_",iun); 
          if(ity==0 && ivar==0){cpt6->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt6->Print(pdfname2,"pdf");
          }else{cpt6->Print(pdfname,"pdf");};

       cpt7->cd();
       SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       Set2dHist(Ematrix[iun][ity][ivar],lplot_xtitle, lplot_xtitle,"",titoff1, titsize1);
       Ematrix[iun][ity][ivar]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       sprintf(pdfname, "%sTUnfold_COV_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_COV_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_COV_%i.pdf)" ,"dd_", iun);
          if(ity==0 && ivar==0 ){cpt7->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt7->Print(pdfname2,"pdf");
          }else{cpt7->Print(pdfname,"pdf");};
      }  //end of phase space cut and variable loop
    }
  }//End of Unfolded plot
*/
/*
//----------------------------------------All Unfold in one plot(Closure Test)
for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     TLegend *leg2= CTLegendV2(0.4,0.80,0.7,0.95, 0.03, "Closure Test", itypeN[ity]);
     cptr->SetGridy(); // Horizontal grid

     TH1D *histgen[clos_ty];
     TH1D *unf[clos_ty];
     TH1D *ratio[clos_ty];
     cptr->cd();
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     for(int icl=0; icl < clos_ty ; icl++){
     histgen[icl]= (TH1D*) Ex_MC_gen[icl][ity][ivar]->Clone();
     unf[icl]= (TH1D*) Ex_Unfold[icl][0][ity][ivar]->Clone();
     Integralhist(histgen[icl]);     Integralhist(unf[icl]);

     ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     ratio[icl]->SetMinimum(.8); ratio[icl]->SetMaximum(1.2);     
     ratio[icl]->SetLineStyle(1);
     ratio[icl]->SetLineColor(icl+1);
     ratio[icl]->SetLineWidth(2);
     ratio[icl]->Draw("Same hist e1");
     leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     }
     leg2->Draw();
     cptr->Update();
     sprintf(pdfname, "E%sratio_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "E%sratio_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "E%sratio_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 ){cptr->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cptr->Print(pdfname2,"pdf");
          }else{cptr->Print(pdfname,"pdf");};

     cptr->Clear();


      }
   }

//----------------------------------------All Unfold in one plot(Closure Test)
for(int ity=0; ity <ntype; ity++){
  for(int ivar =0 ; ivar < nusedvar ; ivar++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for 
     TLegend *leg2= CTLegendV2(0.4,0.80,0.7,0.95,0.03, "Closure Test",itypeN[ity]);
     cptr->SetGridy(); // Horizontal grid

     TH1D *histgen[clos_ty];
     TH1D *unf[clos_ty];
     TH1D *ratio[clos_ty];
     cptr->cd();
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     for(int icl=0; icl < clos_ty ; icl++){

     histgen[icl]= (TH1D*) MC_gen[icl][ity][ivar]->Clone();
     unf[icl]= (TH1D*) Unfold[icl][0][ity][ivar]->Clone();
     Integralhist(histgen[icl]);     Integralhist(unf[icl]);

     ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     ratio[icl]->SetMinimum(.85); ratio[icl]->SetMaximum(1.15);
     ratio[icl]->SetLineStyle(1);
     ratio[icl]->SetLineColor(icl+1);
     ratio[icl]->SetLineWidth(2);
     ratio[icl]->Draw("Same hist e1");
     leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     }
     leg2->Draw();
     cptr->Update();
     sprintf(pdfname, "%sratio_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "%sratio_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sratio_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 ){cptr->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cptr->Print(pdfname2,"pdf");
          }else{cptr->Print(pdfname,"pdf");};

     cptr->Clear();


      }
   }
*/
//----------------------------------------HT2 BLT 
for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
	 for(int ipt = 0 ; ipt <nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "BLT", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     cptr->SetGridy(); // Horizontal grid

     cptr->cd();
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);

     BLT_reco_HT2[ity][ivar][ipt]->SetMinimum(.85); BLT_reco_HT2[ity][ivar][ipt]->SetMaximum(1.15);
     BLT_reco_HT2[ity][ivar][ipt]->SetLineStyle(1);
     BLT_reco_HT2[ity][ivar][ipt]->SetTitle("");
     BLT_reco_HT2[ity][ivar][ipt]->SetLineColor(1);
     BLT_reco_HT2[ity][ivar][ipt]->SetLineWidth(2);
     BLT_reco_HT2[ity][ivar][ipt]->Draw("Same e1");
     BLT_reco_HT2[ity][ivar][ipt]->GetXaxis()->SetTitle(Esvlogx[ivar]);
     BLT_reco_HT2[ity][ivar][ipt]->GetYaxis()->SetTitle("Data/MC");
     BLT_reco_HT2[ity][ivar][ipt]->GetXaxis()->CenterTitle();
     BLT_reco_HT2[ity][ivar][ipt]->GetYaxis()->CenterTitle();

     BLT_gen_HT2[ity][ivar][ipt]->SetMinimum(.85); BLT_gen_HT2[ity][ivar][ipt]->SetMaximum(1.15);
     BLT_gen_HT2[ity][ivar][ipt]->GetXaxis()->SetTitle(Esvlogx[ivar]);
     BLT_gen_HT2[ity][ivar][ipt]->GetYaxis()->SetTitle("Data/MC");
     BLT_gen_HT2[ity][ivar][ipt]->GetXaxis()->CenterTitle();
     BLT_gen_HT2[ity][ivar][ipt]->GetYaxis()->CenterTitle();

     BLT_gen_HT2[ity][ivar][ipt]->SetLineStyle(1);
     BLT_gen_HT2[ity][ivar][ipt]->SetTitle("");
     BLT_gen_HT2[ity][ivar][ipt]->SetLineColor(2);
     BLT_gen_HT2[ity][ivar][ipt]->SetLineWidth(2);
     BLT_gen_HT2[ity][ivar][ipt]->GetXaxis()->SetTitle(Esvlogx[ivar]);
     BLT_gen_HT2[ity][ivar][ipt]->GetYaxis()->SetTitle("Data/MC");
     BLT_gen_HT2[ity][ivar][ipt]->GetXaxis()->CenterTitle();
     BLT_gen_HT2[ity][ivar][ipt]->GetYaxis()->CenterTitle();

     BLT_gen_HT2[ity][ivar][ipt]->Draw("Same  e1");

     leg2->AddEntry(BLT_reco_HT2[ity][ivar][ipt],"Reco Level","lp");
     leg2->AddEntry(BLT_gen_HT2[ity][ivar][ipt],"Gen Level","lp");
     leg1->Draw(); leg2->Draw();
     cptr->Update();
     sprintf(pdfname, "BLT_HT2%sratio_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "BLT_HT2%sratio_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "BLT_HT2%sratio_plot_%i.pdf)","dd_",1234);
     if(ity==0 && ivar==0 && ipt ==0){cptr->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf");
          }else{  cptr->Print(pdfname1,"pdf");};

     cptr->Clear();


      }
   }
}
/*
//----------------------------------------BLT ratio Plot
 for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "Closure Test", itypeN[ity]);
     TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     cptr->SetGridy(); // Horizontal grid
     
     cptr->cd();
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     BLT_reco[ity][ivar]->SetMinimum(.85); BLT_reco[ity][ivar]->SetMaximum(1.15);
     BLT_reco[ity][ivar]->SetLineStyle(1);
     BLT_reco[ity][ivar]->SetLineColor(1);
     BLT_reco[ity][ivar]->SetLineWidth(2);
     BLT_reco[ity][ivar]->Draw("Same e1");
      BLT_gen[ity][ivar]->SetMinimum(.85); BLT_gen[ity][ivar]->SetMaximum(1.15);
     BLT_gen[ity][ivar]->SetLineStyle(1);
     BLT_gen[ity][ivar]->SetLineColor(2);
     BLT_gen[ity][ivar]->SetLineWidth(2);
     BLT_gen[ity][ivar]->Draw("Same  e1");

        sprintf(pdfname, "%sBLT_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "%sBLT_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sBLT_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 ){cptr->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cptr->Print(pdfname2,"pdf");
          }else{cptr->Print(pdfname,"pdf");};


      
        }
      }
*/
//-----------------------------------------------------------------------------------
for(int ity=0; ity <ntype; ity++){
   for(int ivar =0 ; ivar < nusedvar ; ivar++){
      for(int ipt = 0 ; ipt <nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "Closure Test", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     cptr->SetGridy(); // Horizontal grid

     TH1D *histgen[clos_ty];
     TH1D *unf[clos_ty];
     TH1D *ratio[clos_ty];
     cptr->cd();
     SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     for(int icl=0; icl < clos_ty-1 ; icl++){

     histgen[icl]= (TH1D*) MC_GenHT2[icl][ity][ivar][ipt]->Clone();
     unf[icl]= (TH1D*) Unf_HT2[icl][0][ity][ivar][ipt]->Clone();
     Integralhist(histgen[icl]);     Integralhist(unf[icl]);

     ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     ratio[icl]->SetMinimum(.85); ratio[icl]->SetMaximum(1.15);
     ratio[icl]->SetLineStyle(1);
     ratio[icl]->SetLineColor(color[icl]);
     ratio[icl]->SetLineWidth(2);
     ratio[icl]->Draw("Same hist e1");
     ratio[icl]->GetXaxis()->SetTitle(Esvlogx[ivar]);
     ratio[icl]->SetTitle("");
     ratio[icl]->GetYaxis()->SetTitle("MC/Unfolded");
     ratio[icl]->GetXaxis()->CenterTitle();
     ratio[icl]->GetYaxis()->CenterTitle();


     leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     }
     leg1->Draw(); leg2->Draw();
     cptr->Update();
     sprintf(pdfname, "HT2%sratio_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "HT2%sratio_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "HT2%sratio_plot_%i.pdf)","dd_",1234);
     if(ity==0 && ivar==0 && ipt ==0){cptr->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf");
          }else{  cptr->Print(pdfname1,"pdf");};

     cptr->Clear();


      }
   }
}



//--------------------------------------------------------------------------------------------------
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
          
	  TH1D *MyHist  = (TH1D*) MC_gen[0][ity][ivar]->Clone();    	  cout << MyHist->GetName() << endl;

	  Integralhist(MyHist);
       // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          sprintf(Title,"%s:   ",itypeN[ity]);
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D *unfold_input[clos_ty];
          const char *MCinput_index[clos_ty], *data_index[3];
          for(int iout = 0 ; iout < clos_ty  ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*) Unfold[iout][0][ity][ivar]->Clone();
#endif
          //MC_input[iout]->Rebin(2);
          Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout]; 
          //MCinput_index[iout]= closuretype[iout];
          MCinput_index[iout]= closuretype[iout];
      }
      for(int iout = 0 ; iout < clos_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;     cout << unfold_input[iout]->GetName()<< endl; cout <<  MCinput_index[iout] << endl;
          }
          data_index[0]= closuretype[0];
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, " bin ID ( %s, H_{T,2})",Esvlogx[ivar]);
          //float ratio_range1[2]={1.2,0.9};
          int num1[3]={clos_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
          //if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%sunfold_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "%sunfold_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sunfold_plot_%i.pdf)","dd_",1234); 
          if(ity==0 && ivar==0 ){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};
    }
  }


//----------------------------------------All Unfold in one plot(Closure Test) True distibution only for 1D and 2D TUNfolding--------------------------
   for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){

	  TH1D *MyHist  = (TH1D*) Ex_MC_gen[0][ity][ivar]->Clone();    	  cout <<" unfold extract : " << MyHist->GetName() << endl;

	  Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          sprintf(Title,"%s   ",itypeN[ity]);
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D *unfold_input[clos_ty];
          const char *MCinput_index[clos_ty], *data_index[3];
          for(int iout = 0 ; iout < clos_ty ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*) Ex_Unfold[iout][0][ity][ivar]->Clone(); cout <<unfold_input[iout]->GetName()<< endl;
#endif
          //unfold_input[iout] = (TH1D*) Unfold[iout][ity][ivar]->Clone();
          //MC_input[iout]->Rebin(2);
          Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout];
          //MCinput_index[iout]= closuretype[iout];

          MCinput_index[iout]= closuretype[iout];
      }

      for(int iout = 0 ; iout < clos_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;
            cout << unfold_input[iout]->GetName()<< endl;
            cout <<  MCinput_index[iout] << endl;
          }

          data_index[0]= "PY8 Gen";
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.1,0.9};
          int num1[3]={clos_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
      //    if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname,"%sTrueunfold_plot_%i.pdf(", "dd_",1234);sprintf(pdfname1,"%sTrueunfold_plot_%i.pdf","dd_",1234);sprintf(pdfname2,"%sTrueunfold_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};
     }
   }
/*
  //-----------------------------------------Refold comparisoin
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	  
	  TH1D *MyHist  = (TH1D*) Refold[iun][ity][ivar]->Clone();
	  Integralhist(MyHist);
	  //divBybinWidth(MyHist);
	  Myplotset(MyHist,0,0);
	  
          sprintf(Title,"%s   ",itypeN[ity]);
	  sprintf(Yaxis," %s" ,Esvlogy[ivar]);
	  MyHist->SetTitle(Title);
	  MyHist->GetXaxis()->SetTitle("");
	  MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  
	  TH1D *MC_input[nmc];
	  const char *MCinput_index[nmc], *data_index[1];
	  for(int iout = 0 ; iout < nmc ; iout++){
	    MC_input[iout] = (TH1D*) MC_reco_fake[iout][ity][ivar]->Clone();
	    //MC_input[iout]->Rebin(2);
	    Integralhist(MC_input[iout]);
	    //divBybinWidth(MC_input[iout]);
	    MCinput_index[iout]= mcnamerco[iout]; }
	  
	  data_index[0]= RefoldEra[iun];
	  char lplot_xtitle[100];
	  sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
	  //float ratio_range1[2]={1.1,0.9};
	  int num1[2]={nmc,1};
	  float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
	  
	  cpt0->cd();
	  cpt0->SetBorderSize(0);
	  cpt0->SetRightMargin(0.0);
	  cpt0->SetTopMargin(0.0);
	  cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  sprintf(pdfname, "%sRefold_plot_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sRefold_plot_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sRefold_plot_%i.pdf)" ,"dd_",iun);
	  if(ity==0 && ivar==0 ){cpt0->Print(pdfname,"pdf");
	  }else if(ity==1 && ivar==4) {cpt0->Print(pdfname2,"pdf");
	  }else{cpt0->Print(pdfname,"pdf");};
         } 
       }
     }//End of Refolded plot
  */
/*  
 //Refold in one plot
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
          
          TH1D *MyHist  = (TH1D*) MC_reco_fake[umc][ity][ivar]->Clone();
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
     
	   Refold_input[iout] = (TH1D*) Refold[iout][ity][ivar]->Clone();
           //MC_input[iout]->Rebin(2);
            Integralhist(Refold_input[iout]);
           //divBybinWidth(Refold_input[iout]);
            MCinput_index[iout]= RefoldEra[iout]; }

          data_index[0]= mcnamerco[umc];
          data_index[1]= "MC";
          data_index[2]= "Refolded";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.1,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%sRefold_plot_%i.pdf(","dd_", 1234); sprintf(pdfname1, "%sRefold_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sRefold_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

    }
  } 
*/
  /*
 //------------------------------Refold True Dist in one plot
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){

          TH1D *MyHist  = (TH1D*) Ex_Psudo_Data_reco_fake[ity][ivar]->Clone();
          Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          sprintf(Title,"%s   ",itypeN[ity]);
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);


          TH1D *Refold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
            Refold_input[iout] = (TH1D*) Ex_Refold[iout][ity][ivar]->Clone();
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
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};

          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname, "%sTrueRefold_plot_%i.pdf(","dd_", 1234); sprintf(pdfname1, "%sTrueRefold_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sTrueRefold_plot_%i.pdf)","dd_",1234); 
          if(ity==0 && ivar==0 ){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 ) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

    }
  }
*/
  /*
//////////////////////------------------------------HT2 Extracted---------------------------------
  
TH1D* Unf_HT2[][nHLTmx];        
TH1D* MC_HT2[icl][nHLTmx];        
for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){

for(int icl=0; icl <1; icl++){
TH1D* Unf_HT2[icl][nHLTmx];        
HT2_Normal(Ex_Unfold_HT[icl][0][ity][ivar], nHLTmx, HT2range, Unf_HT2[icl]);

TH1D* MC_HT2[icl][nHLTmx];        
HT2_Normal(Ex_MC_gen_HT[icl][0][ity][ivar], nHLTmx, HT2range, MC_HT2[icl]);
//MC_HT2[0]->Draw();
         for(int ipt = 0 ; ipt <nHLTmx; ipt++){

          TH1D *MyHist  = (TH1D*) Unf_HT2[icl][ipt]->Clone();          cout << MyHist->GetName() << endl;
          Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:    H_{T,2} > %i %s",itypeN[ity], HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D* unfold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*)MC_HT2[ipt]->Clone();
#else
          unfold_input[iout] = (TH1D*)MC_HT2[ipt]->Clone();
#endif
          //unfold_input[iout] = (TH1D*) Unfold[iout][ity][ivar][ipt]->Clone();
          //MC_input[iout]->Rebin(2);
          Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout];
          //MCinput_index[iout]= closuretype[iout];

          MCinput_index[iout]= closuretype[iout];
          //MCinput_index[iout]= "MadGraph Gen";
      }
 
      for(int iout = 0 ; iout < unfold_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;
            cout << unfold_input[iout]->GetName()<< endl;
            cout <<  MCinput_index[iout] << endl;
          }

          data_index[0]= closuretype[0];
          //data_index[0]= "Unfolded MadGraph Using PY8 RM";
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.1,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
          if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}
          cpt0->cd();
          cpt0->SetBorderSize(0);
          cpt0->SetRightMargin(0.0);
          cpt0->SetTopMargin(0.0);
          cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
        //  CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          sprintf(pdfname,"%sHTunfold_plot_%i.pdf(", "dd_",1234);sprintf(pdfname1,"%sHTunfold_plot_%i.pdf","dd_",1234);sprintf(pdfname2,"%sHTunfold_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};

       }  //end of phase space cut and variable loop
     }
   }
}
//////////////////////------------------------------HT2 Extracted All in one---------------------------------
  TCanvas *cptx= new TCanvas("cptx", "canvas0", 700,600 );  //for Reco

      	for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
TH1D* Unf_HT2[nHLTmx];
HT2_Normal(Ex_Unfold_HT[0][ity][ivar], nHLTmx, HT2range, Unf_HT2);

TH1D* MC_HT2[nHLTmx];
HT2_Normal(Ex_MC_gen_HT[0][ity][ivar], nHLTmx, HT2range, MC_HT2);
MC_HT2[0]->Draw();
         for(int ipt = 0 ; ipt <nHLTmx; ipt++){

          TH1D *MyHist  = (TH1D*) Unf_HT2[ipt]->Clone();          cout << MyHist->GetName() << endl;
          Integralhist(MyHist);
         // divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);

          if(ipt<=6){sprintf(Title,"%s:     %i <H_{T,2}< %i %s",itypeN[ity] , HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
          else if(ipt==7){ sprintf(Title,"%s:    H_{T,2} > %i %s",itypeN[ity], HT2range[ipt] ,"GeV/c" );}
          sprintf(Yaxis," %s" ,Esvlogy[ivar]);
          MyHist->SetTitle(Title);
          MyHist->GetXaxis()->SetTitle("");
          MyHist->GetYaxis()->SetTitle(Yaxis);

          TH1D* unfold_input[unfold_ty];
          const char *MCinput_index[unfold_ty], *data_index[3];
          for(int iout = 0 ; iout < unfold_ty ; iout++){
#ifdef CLOUSER
          unfold_input[iout]  = (TH1D*)MC_HT2[ipt]->Clone();
#else
          unfold_input[iout] = (TH1D*)MC_HT2[ipt]->Clone();
#endif
          //unfold_input[iout] = (TH1D*) Unfold[iout][ity][ivar][ipt]->Clone();
          //MC_input[iout]->Rebin(2);
          Integralhist(unfold_input[iout]);
          //divBybinWidth(unfold_input[iout]);
          //MCinput_index[iout]= Unfoldtype[iout];
          //MCinput_index[iout]= closuretype[iout];

          MCinput_index[iout]= closuretype[iout];
          //MCinput_index[iout]= "MadGraph Gen";
      }

      for(int iout = 0 ; iout < unfold_ty ; iout++){
            cout << unfold_input[iout]->GetTitle()<< endl;
            cout << unfold_input[iout]->GetName()<< endl;
            cout <<  MCinput_index[iout] << endl;
          }

          data_index[0]= closuretype[0];
          //data_index[0]= "Unfolded MadGraph Using PY8 RM";
          data_index[1]= "Unfolded";
          data_index[2]= "MC";
          char lplot_xtitle[100];
          sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
          //float ratio_range1[2]={1.1,0.9};
          int num1[3]={unfold_ty,1,0};
          float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.1,0.9};
          if(ity==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}
           cptx->cd();
          //cpt0->SetBorderSize(0);
         // cpt0->SetRightMargin(0.0);
          //cpt0->SetTopMargin(0.0);
          gPad->SetLogy();
          Integralhist(MC_HT2[ipt]);
          Integralhist(unfold_input[0]);
          MC_HT2[ipt]->SetLineColor(color[ipt]);
          // MC_HT2[ipt]->SetMaximum(1);
            MC_HT2[ipt]->SetMinimum(1e-4);

	   unfold_input[0]->Draw("Same hist e2");
	   MC_HT2[ipt]->Draw("Same hist");
	  //cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          //CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       }  //end of phase space cut and variable loop
          sprintf(pdfname,"%sHT_plot_%i.pdf(", "dd_",1234);sprintf(pdfname1,"%sHT_plot_%i.pdf","dd_",1234);sprintf(pdfname2,"%sHT_plot_%i.pdf)","dd_",1234);
          if(ity==0 && ivar==0 ){cptx->Print(pdfname1,"pdf");
          }else if(ity==1 && ivar==4) {cptx->Print(pdfname2,"pdf");
          }else{cptx->Print(pdfname,"pdf");};

	 cptx->Clear();
     }
   }
*/
  
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
  
  //ymax = data->GetMaximum();
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
    cout << "Ndf value =" << ndfl << "chi2 value=" << chi2l << endl;
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
  //MyHist->GetXaxis()->SetTitleFont(ifornt);
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetXaxis()->SetTitle(XTitle);
  
  MyHist->GetYaxis()->SetLabelSize(0.03);
  MyHist->GetYaxis()->SetTitleSize(0.040);
  MyHist->GetYaxis()->SetTitleOffset(1.0);
  MyHist->GetYaxis()->SetTitle(YTitle);
  //MyHist->GetYaxis()->SetTitleFont(ifornt);     
  MyHist->GetYaxis()->CenterTitle();
  MyHist->SetTitle(""); 
  //gStyle->SetTitleFontSize(.08);
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

TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4=""){
  TLegend *leg = new TLegend(x1,y1,x2,y2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(txtsize);
  leg->SetTextFont(42);
  leg->AddEntry((TObject*)0, txt1, "");
  leg->AddEntry((TObject*)0, txt2, "");
  leg->AddEntry((TObject*)0, txt3, "");
  leg->AddEntry((TObject*)0, txt4, "");
  return leg;
}


void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3] ){
  MyHist->SetTitleOffset(0.2);
  //MyHist->SetTitleFont(102);
  MyHist->SetTitleSize(0.02);

  MyHist->GetXaxis()->SetLabelSize(0.03);
  MyHist->GetXaxis()->SetTitleSize(titsize[0]);
  MyHist->GetXaxis()->SetTitleOffset(titoff[0]);
  //MyHist->GetXaxis()->SetTitleFont(102);
  //MyHist->GetXaxis()->CenterTitle();

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


TCanvas *ratio_canV2(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
  //Nplot[0] = number of MC enetered
  //Nplot[1] = place 1 if upper part is log scale needed
  //Nplot[2] = place 1 if MC/Data or 0 if data/MC
  //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  //plegend[4]= text size
  //plegend[5-6]= ratio plot axis range
  //data = data histogram
  //MC = monte carlo histogram array
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
  //canvas->SetBottomMargin(0.1);
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

  //ymax = data->GetMaximum();
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
    //MC[iup]->Draw("same hist ");// gPad->SetLogy();

    //--------------------------------------------- Addition fo chi2/Ndf  with legend
    int nn = data->GetNbinsX();
    Double_t resl[nn] , chi2l;
    Int_t ndfl ,igoodl;
    TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test

    TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");   //for chi square test
    chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);

    //chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
 //   cout << "Ndf value =" << ndfl << "chi2 value=" << chi2l << endl;
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

   //MC[iup]->Draw("same hist e1 ");
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
    //staterr->SetFillStyle(1111);
    //staterr->SetFillColor(21);
    //staterr->SetMarkerSize(0.1);



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
    //sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
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
    //staterr->SetFillColor(21);
    //staterr->SetMarkerSize(0.1);

    //-----------------------------------------------------Stat Error
    
    staterr->Draw("E2 Same");
    //mstaterr->Draw("E2 Same");
    
    //MC_inputerr->Draw("same E2");//------------------Draw the Uncertainty in Ratio plot

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

void HT2_Normal(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]){

TH2D* VarHTin = (TH2D*)VarHT->Clone();
int nbiny= VarHTin->GetNbinsY();
int nbinx= VarHTin->GetNbinsX();

TH1D* h_var = VarHTin->ProjectionX(); h_var->Reset();

for(int iht =1 ; iht <= nbiny; iht++){
   h_var->Reset();
   for(int ivar =1 ; ivar<= nbinx  ; ivar++){
   h_var->SetBinContent(ivar,VarHTin->GetBinContent(ivar,iht));
   h_var->SetBinError(ivar,VarHTin->GetBinError(ivar,iht));   
           }
  Var[iht-1] =(TH1D*)h_var->Clone();
cout << " iht :"<< iht << endl;
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


void HT2_NormalV2(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]){

TH2D* VarHTin = (TH2D*)VarHT->Clone();
int nbiny= VarHTin->GetNbinsY();
int nbinx= VarHTin->GetNbinsX();

//TH1D* h_var = VarHTin->ProjectionX(); h_var->Reset();

for (int iht = 1; iht <= nbiny; ++iht) {
       const char * name_hbin = Form("%s_pt%i",VarHTin->GetName(), iht);
        TH1D* h_var  = VarHTin->ProjectionX(name_hbin, iht, iht);
        Var[iht-1] =(TH1D*)h_var->Clone();        
        //h_var-Write();
      }

  }



