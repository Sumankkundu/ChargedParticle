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

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>

//#define BLT
#define CLOUSER
void Unfoldplot(){
  static const int unfold_ty =1; //Unfold method to be plots 
  static const int clos_ty =4; //Number of Closure Test performed 
  static const int nHLTmx=8; //HT2 Range
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =2;
  static const int nmc = 4;   //Number of MC sample 
  static const int ndd = 3; 
  static const int umc = 0; //0 for Py8, 1 for MG , 2 for Herwig : which MC have used in Unfolding

  bool isstat =1;  int irbin=1;
  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={2,4,6,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[nusedvar]={3,9,15,18,24};
  //Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000}; //2017
  Int_t HT2range[nHLTmx+1]={83, 109, 176, 247, 318, 387, 477, 573, 3000}; //2018

  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  
  string  Datadir = "Data2D";
  string folddir = "Folded2D";
  string unfdir = "Unfold2D";
  const char*  dirname[4]={"Pythia8","MG8","HW7","Py8Flat"};
#ifdef CLOUSER
  string mcdir[4]={"Pythia8","Py8Flat","MG8","HW7"}; //Sequence For Closure Test
#else
  string mcdir[4]={"Pythia8","MG8","HW7","Py8Flat"}; //Sequence For Unfold and Refold plot
#endif
  
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  //const char* htrang[8]={"83 < H_{T,2} < 109", "109 < H_{T,2} < 172", "172 < H_{T,2} < 241", "241 < H_{T,2} < 309","309 < H_{T,2} < 377", "377 < H_{T,2} < 462", "462 < H_{T,2} <570","H_{T,2} > 570"};
  const char* htrang[8]={"83 < H_{T,2} < 109", "109 < H_{T,2} < 176", "176 < H_{T,2} < 247", "247 < H_{T,2} < 318","318 < H_{T,2} < 387", "387 < H_{T,2} < 477", "477 < H_{T,2} <573","H_{T,2} > 573"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};
  const  char* Unit[4]={"GeV","fb^{-1}","Pb","#%"};

  const  char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const  char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};

#ifdef CLOUSER
  const  char* mcname[4]={"Py8-gen","Py8 Flat-gen","MadGraph-gen","Herwig7-gen"};
#else 
  const  char* mcname[4]={"Pythia8","MadGraph","Herwig7","Pythia8 Flat"};
#endif

  const  char* DataEra[3]={"Data","Data","Data"};
  const  char* RunEra[3]={"2016","2017","2018"};
  static const int iera = 2;
  int iPeriod = 0;  int iPos=0 ;
  
  const  char* mcnamerco[4]={"Pythia8 ","MadGraph ","Herwig7","PY8 Flat"};
  const  char* closuretype[4]={"PY8 by PY8","PY8 Flat By PY8 ", "MG By PY8","HW7 By PY8"};
  const  char* itypeN[ntype]={"Jets","Charged Particles"}; 
  const  char* Refoldtype[5]={"Py8 Refold","Py8 Flat Refold ","MadGraph Refold","Herwig7 Refold","Refold Pythia8(Iterative method)"};
#ifdef CLOUSER
  const  char* Unfoldtype[5]={"Py8-unf","Py8 Flat-unf","MadGraph-unf","Herwig7-unf ","Others"}; 
#else 
  const  char* Unfoldtype[5]={"Unfold by Py8 RM","Unfold by MG RM","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
#endif
  const  char* Modelnm[3]={"Pythia8","MadGraph","Herwig7"};
  const  char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const  char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  const  char* data_unf[5]={"Unf by Py8","Unf by MadGraph","Unf by Herwig7","GEN","RECO"};
  const  char* Refoldname[5]={"Py8 Refold","Py8 Flat Refold ","MadGraph Refold","Herwig7 Refold","Refold Pythia8(Iterative method)"};
  const  char* Recomcname[4]={"Py8-reco","Py8 Flat-reco","MadGraph-reco","Herwig7-reco"};
  
  //TFile *inputbinning = TFile::Open("PY8_HT2_83_3000_Binnied.root");
  TFile *inputbinning = TFile::Open("PY8_CP5_Flat_6Jan21.root");
 
  TFile *Unf_root[clos_ty];
  TFile *Unf_dataCT[clos_ty];
  //******************************************************************************
  TFile *Unfoldroot = TFile::Open("../Unfold/Unfolded_Result_Data_PY8.root");  // Unfolded data 
  
  Unf_root[0] = TFile::Open("../Unfold/Unfolded_Result_PY8_PY8.root");  // Py- by Py8
  Unf_root[1] = TFile::Open("../Unfold/Unfolded_Result_Py8_Flat_PY8.root");  // PyFlat - By Py8
  Unf_root[2]= TFile::Open("../Unfold/Unfolded_Result_MG_PY8.root");  // MG by Py8
  Unf_root[3]= TFile::Open("../Unfold/Unfolded_Result_HW_Py8.root");  // HW by Py8 
  
  Unf_dataCT[0]= TFile::Open("../Unfold/Unfolded_Result_Data_PY8.root");  // Py MG 
  Unf_dataCT[1]= TFile::Open("../Unfold/Unfolded_Result_Data_MG.root");  // Py MG 
  Unf_dataCT[2]= TFile::Open("../Unfold/Unfolded_Result_Data_HW.root");  // Py MG 
  //--------------------------------------Function declearation
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_canV2(int Nplot[3],float plegend[8], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
  void HT2_Normal(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  void HT2_NormalV2(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  TH1D* ReadHist1D(string name,TFile* root, int irbin=1);
  TH2D* ReadHist2D(string name,TFile* root, int irbin=1);
  TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4="");
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl);

  TH1D *Data_reco[ntype][nusedvar];
  TH2D *Data_reco2d[ntype][nusedvar];
  TH1D *MC_gen[nmc][ntype][nusedvar];   //MC gen
  TH1D *MC_genM[nmc][ntype][nusedvar];   //MC gen
  TH2D *MC_gen2d[nmc][ntype][nusedvar];   //MC gen
  TH1D *MC_genmiss[nmc][ntype][nusedvar]; //MC Gen-miss
  TH1D *MC_reco[nmc][ntype][nusedvar];   //MC Reco 
  TH2D *MC_reco2d[nmc][ntype][nusedvar];   //MC Reco 
  TH1D *MC_recofake[nmc][ntype][nusedvar]; //MC Reco -Fake
  TH2D *MC_Res[nmc][ntype][nusedvar];   //RM
  
  TH1D *hist_eff[nmc][ntype][nusedvar];
  TH1D *hist_fake[nmc][ntype][nusedvar];
  TH1D *hist_purity[nmc][ntype][nusedvar];
  TH1D *hist_stbl[nmc][ntype][nusedvar];

  TH1D *UnfoldCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *UnfoldCT_M[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *RefoldCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *UnfoldCT_HT[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];
  TH1D *RefoldCT_HT[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];


  TH1D *UnfoldDataCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *UnfoldDataCT_M[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *RefoldDataCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH1D *UnfoldDataCT_HT[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];
  TH1D *RefoldDataCT_HT[clos_ty][unfold_ty][ntype][nusedvar][nHLTmx];


  TH1D *Unfold[ntype][nusedvar];
  TH2D *Unfold2d[ntype][nusedvar];
  TH1D *Refold[ntype][nusedvar];
  TH1D *Unfold_HT[ntype][nusedvar][nHLTmx];
  TH1D *Refold_HT[ntype][nusedvar][nHLTmx];

  //HT2 bin....
  TH1D *Data_recoHT[ntype][nusedvar][nHLTmx];
  TH1D* MC_genHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_recoHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_recofakeHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D* MC_genmissHT[nmc][ntype][nusedvar][nHLTmx];
  
  TH1D *hist_effHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_fakeHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_purityHT[nmc][ntype][nusedvar][nHLTmx];
  TH1D *hist_stblHT[nmc][ntype][nusedvar][nHLTmx];
  
  TH1D *BLTDet[ntype][nusedvar][nHLTmx];
  TH1D *BLTGen[ntype][nusedvar][nHLTmx];
#ifdef CLOUSER  
  TH1D *Psudo_Data_gen[ntype][nusedvar];  
  TH1D *Psudo_Data_genHT[ntype][nusedvar][nHLTmx];  
  TH1D *Psudo_Data_reco_fake[ntype][nusedvar];  
  TH1D *Psudo_Data_gen_miss[ntype][nusedvar];
  TH1D *Psudo_Data_reco_fakeHT[ntype][nusedvar];
#endif

  TH2D *CorrCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH2D *ProbCT[clos_ty][unfold_ty][ntype][nusedvar];
  TH2D *EmatrixCT[clos_ty][unfold_ty][ntype][nusedvar];

  TH2D *Corr[ntype][nusedvar];
  TH2D *Prob[ntype][nusedvar];
  TH2D *Ematrix[ntype][nusedvar];


//---------------------------------------------Read Data
for(int ity=0; ity <ntype; ity++){
      for(int ivar =0;ivar < nusedvar; ivar++){
	  Data_reco[ity][ivar] = (TH1D*)ReadHist1D(Datadir+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          Data_reco2d[ity][ivar] = (TH2D*)ReadHist2D(Datadir+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          for(int ipt=0; ipt < nHLTmx; ipt++){
          Data_recoHT[ity][ivar][ipt] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);             }
#ifdef CLOUSER
   //       Psudo_Data_gen[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot); 
//	  for(int ipt=0; ipt < nHLTmx; ipt++){
  //        Psudo_Data_genHT[ity][ivar][ipt] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);    }
#endif
      }
    }
//--------------------------------------------Read MC
  for(int  imc =0; imc < nmc ; imc++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
	  MC_gen[imc][ity][ivar] =(TH1D*)ReadHist1D(mcdir[imc]+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);  
	  MC_genM[imc][ity][ivar] =(TH1D*)ReadHist1D(mcdir[imc]+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px",Unfoldroot);  
	  MC_gen2d[imc][ity][ivar] =(TH2D*)ReadHist2D(mcdir[imc]+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);  
	  MC_reco[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
	  MC_reco2d[imc][ity][ivar] = (TH2D*)ReadHist2D(mcdir[imc]+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          MC_genmiss[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Genminusmiss_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          MC_recofake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Recominusfake_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          
	  hist_eff[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_miss_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_fake_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_purity[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Purity_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
          hist_stbl[imc][ity][ivar] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_stability_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
	  MC_Res[imc][ity][ivar] = (TH2D*)ReadHist2D(mcdir[imc]+"/dd_corr_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);

         for(int ipt=0; ipt < nHLTmx; ipt++){
          MC_genHT[imc][ity][ivar][ipt] =(TH1D*)ReadHist1D(mcdir[imc]+"/Edd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);  
          MC_recoHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          MC_genmissHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Genminusmiss_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          MC_recofakeHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Recominusfake_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          
          hist_effHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_miss_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          hist_fakeHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_fake_rate_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          hist_purityHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Purity_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
          hist_stblHT[imc][ity][ivar][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_stability_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
	           }
	 //----------------------------------------------------------------
      }
    }
  }//--------------end of one MCinput root file  reading
  //---------------------------------------------------------//Read Unfolded MC for CT
  //
#ifdef CLOUSER  
  cout << "Read The CT  Result : " << endl;
  for(int icl=0; icl < clos_ty; icl++){
  for(int iun=0; iun < unfold_ty; iun++){
     for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	  UnfoldCT[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
	  UnfoldCT_M[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unf_root[icl]);
	  RefoldCT[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
          CorrCT[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
          EmatrixCT[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
          ProbCT[icl][iun][ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_Prob_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_root[icl]);
   
          for(int ipt=0; ipt < nHLTmx; ipt++){
          UnfoldCT_HT[icl][iun][ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unf_root[icl]);
          RefoldCT_HT[icl][iun][ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unf_root[icl]);
         	  }
	       }
            }
         }
      }
//-----------------------------------------------Unfolded data with different MC RM
for(int icl=0; icl < 3; icl++){
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
          UnfoldDataCT[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_dataCT[icl]);
          UnfoldDataCT_M[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_px", Unf_dataCT[icl]);
          RefoldDataCT[icl][iun][ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unf_dataCT[icl]);

          for(int ipt=0; ipt < nHLTmx; ipt++){
          UnfoldDataCT_HT[icl][iun][ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unf_dataCT[icl]);
          RefoldDataCT_HT[icl][iun][ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unf_dataCT[icl]);
             }
          }
         }
      }
     }

#endif
//--------------------------------------------Read Unfolded data  with PY8
  for(int ity=0; ity <ntype; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
          Unfold[ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);
          Unfold2d[ity][ivar] = (TH2D*)ReadHist2D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);
          Refold[ity][ivar] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);
          Corr[ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);
          Ematrix[ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);
          Prob[ity][ivar]= (TH2D*)ReadHist2D(unfdir+"/dd_Prob_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), Unfoldroot);

          for(int ipt=0; ipt < nHLTmx; ipt++){
          Unfold_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unfoldroot);
          Refold_HT[ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unfoldroot);
#ifdef BLT        
	  BLTDet[ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/EBLTDet_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unfoldroot);
          BLTGen[ity][ivar][ipt] = (TH1D*)ReadHist1D(unfdir+"/EBLTgen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), Unfoldroot);
#endif
	  }
          }
        }
  cout <<" Read histogram OK " << endl;

//-----------------------------------------------------Reco comparison 
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
       for(int ipt=0; ipt < nHLTmx; ipt++){
        TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for 
    	TH1D *MyHist  = (TH1D*)Data_recoHT[ity][ivar][ipt]->Clone();  Integralhist(MyHist);	//divBybinWidth(MyHist);
	Myplotset(MyHist,0,0);
	
	//MyHist->SetTitle(Form("%s: %s  ",itypeN[ity],htrang[ipt]));
	MyHist->SetTitle(Form(" %s %s",htrang[ipt],Unit[0]));
	MyHist->GetXaxis()->SetTitle("");
	MyHist->GetYaxis()->SetTitle(Form(" %s" ,Esvlogy[ivar]));
	
	TH1D *MC_input[nmc];
	const char *MCinput_index[nmc+1];
	const char *data_index[1];
	for(int iout = 0 ; iout < nmc ; iout++){
	  MC_input[iout] = (TH1D*) MC_recoHT[iout][ity][ivar][ipt]->Clone();  Integralhist(MC_input[iout]);	 //divBybinWidth(MC_input[iout]);
	  MCinput_index[iout]= mcnamerco[iout]; }
	  data_index[0]= DataEra[iera]; 
	
	char lplot_xtitle[100];
	sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
	//float ratio_range1[2]={1.2,0.9};
	int num1[2]={nmc,1} ;
	float lpos1[7] ={.32,0.2,0.55,0.45, .055, 1.5,0.7};
	
	cpt0->cd();
	SetMycanvas(cpt0,0,0,0,0,0);

	cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
	CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	
	sprintf(pdfname, "%sRecoEVS_%sPlot.pdf(","dd_",RunEra[iera]); sprintf(pdfname1, "%sRecoEVS_%sPlot.pdf","dd_",RunEra[iera]); sprintf(pdfname2, "%sRecoEVS_%sPlot.pdf)","dd_",RunEra[iera]); 
	if(ity==0 && ivar==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 && ipt==7) {cpt0->Print(pdfname2,"pdf");
       }else{cpt0->Print(pdfname1,"pdf");};
      }
    }  
  }//-------------------end of RECO PLOT    
//------------------------------------------------Response Matrix--
for(int  imc =0; imc < nmc ; imc++){
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  cpt0->cd();
       	SetMycanvas(cpt0,0,0.1,0.15,0.05,0.1);
        char lplot_xtitle[100]; char lplot_ytitle[100];
        sprintf(lplot_xtitle, "Reco bin ID ( %s, H_{T,2})",Esvlogx[ivar]); sprintf(lplot_ytitle, "Gen  bin ID ( %s, H_{T,2}) ",Esvlogx[ivar]);	
	double titoff1[3]={1.2,1.3,1.0};
        double titsize1[3] ={0.035,0.035,0.035};
	
        Set2dHist( MC_Res[imc][ity][ivar],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
        gPad->SetLogz();
	MC_Res[imc][ity][ivar]->Draw("colz");
         
        sprintf(Title,"%s:   ",itypeN[ity]);
	TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
	CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,itypeN[ity] , "");leg1->SetTextColor(-8);//leg1->Draw();
        CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	sprintf(pdfname, "%sResponse_Mat%i_%s.pdf(","dd_",imc,RunEra[iera]); sprintf(pdfname1, "%sResponse_Mat%i_%s.pdf","dd_",imc,RunEra[iera]); sprintf(pdfname2, "%sResponse_Mat%i_%s.pdf)","dd_",imc,RunEra[iera]);
        if(ity==0 && ivar==0 ){cpt0->Print(pdfname,"pdf"); }else if(ity==1 && ivar==4) {cpt0->Print(pdfname2,"pdf"); }else{  cpt0->Print(pdfname1,"pdf");};
        }
      }
    }

  cout <<" RECO PLOT OK " << endl;
 
//------------------------------------------------------Efficiency, Purity,Fake rate, stability

  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 600,700 );  //for 
        TCanvas *cpt2 = new TCanvas("cpt2", "canvas2", 600,700 );  //for ESVs
        TCanvas *cpt3 = new TCanvas("cpt3", "canvas3", 600,700 );  //for ESVs
        TCanvas *cpt4 = new TCanvas("cpt4", "canvas4", 600,700 );  //for 

	//void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm){
        SetMycanvas(cpt1,0,0.11,0.02,0.05,0.12);
        SetMycanvas(cpt2,0,0.11,0.03,0.05,0.12);
        SetMycanvas(cpt4,0,0.11,0.03,0.05,0.12);
        SetMycanvas(cpt3,0,0.11,0.03,0.05,0.12);

   	TLegend *leg2 = new TLegend(0.4,0.15,0.7,.55);
	CTLegend(leg2," ","");
	leg2->SetTextFont(132);
	leg2->SetTextSize(0.035);

	TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        CTLegend(leg1,"", itypeN[ity]); 

       for(int ipt=0; ipt < nHLTmx; ipt++){
        cpt1->cd();
	if(ipt<=6){sprintf(Title,"%i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV" );}
        else if(ipt==7){ sprintf(Title," H_{T,2} > %i %s", HT2range[ipt] ,"GeV" );}

        //sprintf(Title,"%s:   ",itypeN[ity]);

        for (int i = 1; i <= hist_effHT[umc][ity][ivar][ipt]->GetNbinsX(); ++i) {
         double content = 1 - hist_effHT[umc][ity][ivar][ipt]->GetBinContent(i);
         hist_effHT[umc][ity][ivar][ipt]->SetBinContent(i, content);
       }
         hist_effHT[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_effHT[umc][ity][ivar][ipt]->SetMaximum(1.1);
         hist_fakeHT[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_fakeHT[umc][ity][ivar][ipt]->SetMaximum(1.01);
         hist_purityHT[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_purityHT[umc][ity][ivar][ipt]->SetMaximum(1.01);
         hist_stblHT[umc][ity][ivar][ipt]->SetMinimum(-0.01); hist_stblHT[umc][ity][ivar][ipt]->SetMaximum(1.01);

        Myplotset(hist_effHT[umc][ity][ivar][ipt],Esvlogx[ivar],"Efficiency");
        hist_effHT[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
	hist_effHT[umc][ity][ivar][ipt]->Draw("same hist "); //leg1->Draw();
        leg2->AddEntry(hist_effHT[umc][ity][ivar][ipt], Title ,"lp");
	if(ipt==7){leg2->Draw();}
        //leg2->AddEntry(hist_effHT[0][ity][ivar][ipt], Title ,"lp");

        cpt2->cd();
       //SetMycanvas(cpt2);
 	Myplotset(hist_purityHT[umc][ity][ivar][ipt],Esvlogx[ivar],"Purity");
 	hist_purityHT[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
        hist_purityHT[umc][ity][ivar][ipt]->Draw("same hist"); //leg1->Draw();
        if(ipt==7){leg2->Draw();}
	cpt2->Update();

        cpt3->cd();
        Myplotset(hist_fakeHT[umc][ity][ivar][ipt],Esvlogx[ivar],"Fake rate");
        hist_fakeHT[umc][ity][ivar][ipt]->SetLineColor(color[ipt]);
        hist_fakeHT[umc][ity][ivar][ipt]->Draw("same hist"); //leg1->Draw();
        if(ipt==7){leg2->Draw();} 
       	cpt3->Update();

	cpt4->cd();
        Myplotset(hist_stblHT[umc][ity][ivar][ipt],Esvlogx[ivar],"Stability");
        hist_stblHT[umc][ity][ivar][ipt]->SetLineColor(color[ipt]); 
        hist_stblHT[umc][ity][ivar][ipt]->Draw("same hist"); //leg1->Draw();
        if(ipt==7){leg2->Draw();}

	cpt4->Update();
        }
         CMS_lumi( cpt1, iPeriod, iPos ); cpt1->Update();
         CMS_lumi( cpt2, iPeriod, iPos ); cpt2->Update();
         CMS_lumi( cpt3, iPeriod, iPos ); cpt3->Update();
         CMS_lumi( cpt4, iPeriod, iPos ); cpt4->Update();

    sprintf(pdfname, "%seffi_%splot.pdf(","dd_",RunEra[iera]); sprintf(pdfname1, "%seffi_%splot.pdf","dd_",RunEra[iera]); sprintf(pdfname2, "%seffi_%splot.pdf)","dd_",RunEra[iera]);
        if(ity==0 && ivar==0){cpt1->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt1->Print(pdfname2,"pdf");
        }else{  cpt1->Print(pdfname1,"pdf");};
    sprintf(pdfname, "%spuri_%splot.pdf(","dd_",RunEra[iera]); sprintf(pdfname1, "%spuri_%splot.pdf","dd_",RunEra[iera]); sprintf(pdfname2, "%spuri_%splot.pdf)","dd_",RunEra[iera]); 
        if(ity==0 && ivar==0 ){cpt2->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt2->Print(pdfname2,"pdf");
        }else{  cpt2->Print(pdfname1,"pdf");};
     sprintf(pdfname, "%sfake_%splot.pdf(","dd_",RunEra[iera]); sprintf(pdfname1, "%sfake_%splot.pdf","dd_",RunEra[iera]); sprintf(pdfname2, "%sfake_%splot.pdf)","dd_",RunEra[iera]); 
        if(ity==0 && ivar==0 ){cpt3->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt3->Print(pdfname2,"pdf");
        }else{  cpt3->Print(pdfname1,"pdf");};
     sprintf(pdfname, "%sstab_%splot.pdf(","dd_",RunEra[iera]); sprintf(pdfname1, "%sstab_%splot.pdf","dd_",RunEra[iera]); sprintf(pdfname2, "%sstab_%splot.pdf)","dd_",RunEra[iera]); 
        if(ity==0 && ivar==0 ){cpt4->Print(pdfname,"pdf");
        }else if(ity==1 && ivar==4 ) {cpt4->Print(pdfname2,"pdf");
        }else{  cpt4->Print(pdfname1,"pdf");};
	cpt1->Clear(); cpt2->Clear(); cpt3->Clear(); cpt4->Clear();
     }
  }   
cout << "Effi, Stabilty, fake, purity " <<endl;


//--------------------------------------------------------------Unfold comparisoin
  for(int iun=0; iun < unfold_ty; iun++){
    for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
       char lplot_xtitle[100];
	 TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,600 );  //for  
	 TCanvas *cpt6 = new TCanvas("cpt6", "canvas6", 600,600 );  //for  
	 TCanvas *cpt7 = new TCanvas("cpt7", "canvas7", 600,600 );  //for  
        for(int ipt=0; ipt < nHLTmx; ipt++){
	 TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,700 );  //for  unfolded plots
	  TH1D *MyHist  = (TH1D*)Unfold_HT[ity][ivar][ipt]->Clone();
	  Integralhist(MyHist);	  divBybinWidth(MyHist); 	  Myplotset(MyHist,0,0);
	  
	  //MyHist->SetTitle(Form("%s: %s %s",itypeN[ity],htrang[ipt],Unit[0]));     
	  if(ipt<=6){MyHist->SetTitle(Form(" %i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
          else if(ipt==7){MyHist->SetTitle(Form(" H_{T,2} > %i %s", HT2range[ipt],Unit[0]));}	  // Ommit Jet or Charged Particle name
	  MyHist->GetXaxis()->SetTitle("");
	  MyHist->GetYaxis()->SetTitle(Form(" %s" ,Esvlogy[ivar]));
	  
	  TH1D *MC_input[nmc];
	  const char *MCinput_index[nmc], *data_index[1];
	  for(int iout = 0 ; iout < nmc ; iout++){
	  MC_input[iout] = (TH1D*) MC_genHT[iout][ity][ivar][ipt]->Clone();
	  Integralhist(MC_input[iout]);
	  divBybinWidth(MC_input[iout]);
	   
	  MCinput_index[iout]= mcname[iout]; }


       data_index[0]= DataEra[iera]; 
       //data_index[0]= UndoldEra[1]; 
       sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
       //float ratio_range1[2]={1.2,0.9};
       int num1[2]={nmc,1} ;
       float lpos1[7] ={.32,0.2,0.55,0.45, .055, 1.45,0.7};
	  
       cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
       CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       sprintf(pdfname, "%sTUnfold_%splot_%i.pdf(" ,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%sTUnfold_%splot_%i.pdf" ,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%sTUnfold_%splot_%i.pdf)" ,"dd_",RunEra[iera],iun);
       if(ity==0 && ivar==0 && ipt==0){cpt0->Print(pdfname,"pdf");
	 }else if(ity==1 && ivar==4 && ipt ==7) {cpt0->Print(pdfname2,"pdf");
	  }else{cpt0->Print(pdfname,"pdf");};
	
	}

       cpt5->cd();
       double titoff1[3]={1.2,1.3,1.0};
       double titsize1[3] ={0.035,0.035,0.035};
       SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       gStyle->SetPaintTextFormat( "4.2f");
       Set2dHist(Corr[ity][ivar],lplot_xtitle, lplot_xtitle,"correlation coefficients", titoff1, titsize1);
       Corr[ity][ivar]->Draw("colz ");
       //Corr[iun][ity][ivar]->Draw("colz text");
       sprintf(Yaxis," %s" ,Esvlogy[ivar]);
       cout <<"OK1" <<endl;
       TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
       CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,itypeN[ity] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       sprintf(pdfname, "%sTUnfold_corr%s_%i.pdf(" ,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%sTUnfold_corr%s_%i.pdf" ,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%sTUnfold_corr%s_%i.pdf)" ,"dd_",RunEra[iera],iun); 
       if(ity==0 && ivar==0 ){cpt5->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4 ) {cpt5->Print(pdfname2,"pdf");
          }else{cpt5->Print(pdfname,"pdf");};

       cpt6->cd();
       SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       Set2dHist(Prob[ity][ivar],lplot_xtitle, lplot_xtitle,"Probability",titoff1, titsize1);
       gPad->SetLogz();Prob[ity][ivar]->SetMinimum(1e-6); Prob[ity][ivar]->SetMaximum(1);
       Prob[ity][ivar]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       sprintf(pdfname, "%sTUnfold_prob%s_%i.pdf(" ,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%sTUnfold_prob%s_%i.pdf" ,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%sTUnfold_prob%s_%i.pdf)" ,"dd_",RunEra[iera],iun); 
       if(ity==0 && ivar==0){cpt6->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt6->Print(pdfname2,"pdf");
          }else{cpt6->Print(pdfname,"pdf");};

       cpt7->cd();
       SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       Set2dHist(Ematrix[ity][ivar],lplot_xtitle, lplot_xtitle,"",titoff1, titsize1);
       Ematrix[ity][ivar]->Draw("colz"); leg1->Draw();
       CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       sprintf(pdfname, "%sTUnfold_COV%s_%i.pdf(" ,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%sTUnfold_COV%s_%i.pdf" ,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%sTUnfold_COV%s_%i.pdf)" ,"dd_",RunEra[iera], iun);
      if(ity==0 && ivar==0 ){cpt7->Print(pdfname,"pdf");
          }else if(ity==1 && ivar==4) {cpt7->Print(pdfname2,"pdf");
          }else{cpt7->Print(pdfname,"pdf");};
      }  //end of phase space cut and variable loop
    }
  }//End of Unfolded plot

#ifdef CLOUSER  

int markersty[4]={3,26,5,6};


//-----------------------------Closure Test Plot with histogram 
for(int icl=0; icl < clos_ty ; icl++){
  for(int ity=0; ity <ntype; ity++){
     for(int ivar =0 ; ivar < nusedvar ; ivar++){
         char lplot_xtitle[100];
	 int iun =0;
         for(int ipt=0; ipt < nHLTmx; ipt++){
         TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for  
          TH1D *MyHist  = (TH1D*)UnfoldCT_HT[icl][0][ity][ivar][ipt]->Clone();; 
          Integralhist(MyHist);          divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);
          MyHist->GetXaxis()->SetTitle("");  
	  if(ipt<=6){MyHist->SetTitle(Form(" %i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
          else if(ipt==7){MyHist->SetTitle(Form(" H_{T,2} > %i %s", HT2range[ipt],Unit[0]));}  
	  MyHist->GetYaxis()->SetTitle(Form(" %s" ,Esvlogy[ivar]));
          
	  TH1D *MC_input[1];
          const char *MCinput_index[1], *data_index[3];
          MC_input[0] = (TH1D*) MC_genHT[icl][ity][ivar][ipt]->Clone();
          Integralhist(MC_input[0]);          divBybinWidth(MC_input[0]);
          MCinput_index[0] = mcname[icl]; 
        
       data_index[0]= Unfoldtype[icl];
       data_index[1]= Unfoldtype[icl]; 
       data_index[2]= mcname[icl]; 
       sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
       //float ratio_range1[2]={1.2,0.9};
       int num1[3]={1,1,0} ;
       float lpos1[8] ={.32,0.2,0.55,0.45, .055, 1.6,0.2,0.1};

       //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
       CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       sprintf(pdfname, "%i_%sClosureHist_%splot_%i.pdf(" ,icl,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%i_%sClosureHist_%splot_%i.pdf" ,icl,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%i_%sClosureHist_%splot_%i.pdf)" ,icl,"dd_",RunEra[iera],iun);
       if(ity==0 && ivar==0 && ipt==0){cpt0->Print(pdfname,"pdf");
         }else if(ity==1 && ivar==4 && ipt ==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};
        }
      }  //end of phase space cut and variable loop
    }
  }//End of Unfolded plot

//-----------------------------Closure same data with different MC
  for(int ity=0; ity <ntype; ity++){
     for(int ivar =0 ; ivar < nusedvar ; ivar++){
         char lplot_xtitle[100];
         int iun =0;
         for(int ipt=0; ipt < nHLTmx; ipt++){
         TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for
          TH1D *MyHist  = (TH1D*)UnfoldDataCT_HT[0][0][ity][ivar][ipt]->Clone();;
          Integralhist(MyHist);          divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);
          MyHist->GetXaxis()->SetTitle("");  
	  if(ipt<=6){MyHist->SetTitle(Form(" %i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
          else if(ipt==7){MyHist->SetTitle(Form(" H_{T,2} > %i %s", HT2range[ipt],Unit[0]));}
	  MyHist->GetYaxis()->SetTitle(Form(" %s" ,Esvlogy[ivar]));

          TH1D *MC_input[2];
          const char *MCinput_index[2], *data_index[3];
  
for(int icl=0; icl < 2 ; icl++){
  	  MC_input[icl] = (TH1D*) UnfoldDataCT_HT[icl+1][0][ity][ivar][ipt]->Clone();
          Integralhist(MC_input[icl]);          divBybinWidth(MC_input[icl]);
          MCinput_index[icl] = data_unf[icl+1];}

       data_index[0]= data_unf[iun];
       data_index[1]= "Unf by other";
       data_index[2]= data_unf[iun];
       sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
       //float ratio_range1[2]={1.2,0.9};
       int num1[3]={2,1,0} ;
       float lpos1[8] ={.32,0.2,0.55,0.45, .05, 1.45,0.7,0.1};

       //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
       CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       sprintf(pdfname, "%sDataCT_%splot_%i.pdf(" ,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%sDataCT_%splot_%i.pdf" ,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%sDataCT_%splot_%i.pdf)" ,"dd_",RunEra[iera],iun);
       if(ity==0 && ivar==0 && ipt==0){cpt0->Print(pdfname,"pdf");
         }else if(ity==1 && ivar==4 && ipt ==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};
        }
      }  //end of phase space cut and variable loop
  }//End of Unfolded plot


//-----------------------------Refold Plot with histogram
for(int icl=0; icl < clos_ty ; icl++){
  for(int ity=0; ity <ntype; ity++){
     for(int ivar =0 ; ivar < nusedvar ; ivar++){
         char lplot_xtitle[100];
         int iun =0;
         for(int ipt=0; ipt < nHLTmx; ipt++){
         TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for
          TH1D *MyHist  = (TH1D*)RefoldCT_HT[icl][0][ity][ivar][ipt]->Clone();;
          Integralhist(MyHist);          divBybinWidth(MyHist);
          Myplotset(MyHist,0,0);
          MyHist->GetXaxis()->SetTitle("");  
	 if(ipt<=6){MyHist->SetTitle(Form(" %i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
          else if(ipt==7){MyHist->SetTitle(Form(" H_{T,2} > %i %s", HT2range[ipt],Unit[0]));} 
	  MyHist->GetYaxis()->SetTitle(Form(" %s" ,Esvlogy[ivar]));

          TH1D *MC_input[1];
          const char *MCinput_index[1], *data_index[3];
          MC_input[0] = (TH1D*) MC_recoHT[icl][ity][ivar][ipt]->Clone();
          Integralhist(MC_input[0]);          divBybinWidth(MC_input[0]);
          MCinput_index[0] = Recomcname[icl];

       data_index[0]= Refoldname[icl];
       data_index[1]= Refoldname[icl];
       data_index[2]= Recomcname[icl];
       sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
       //float ratio_range1[2]={1.2,0.9};
       int num1[3]={1,1,0} ;
       float lpos1[8] ={.32,0.2,0.55,0.45, .055, 1.6,0.2,0.1};

       //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
       CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       sprintf(pdfname, "%i_%sRefoldHist_%splot_%i.pdf(" ,icl,"dd_",RunEra[iera],iun); sprintf(pdfname1, "%i_%sRefoldHist_%splot_%i.pdf" ,icl,"dd_",RunEra[iera],iun);sprintf(pdfname2, "%i_%sRefoldHist_%splot_%i.pdf)" ,icl,"dd_",RunEra[iera],iun);
       if(ity==0 && ivar==0 && ipt==0){cpt0->Print(pdfname,"pdf");
         }else if(ity==1 && ivar==4 && ipt ==7) {cpt0->Print(pdfname2,"pdf");
          }else{cpt0->Print(pdfname,"pdf");};
        }
      }  //end of phase space cut and variable loop
    }
  }//End of Unfolded plot



#endif

  
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
  
  TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 680,750 );
  canvas->cd();
  canvas->SetRightMargin(0.02);
//canvas->SetLeftMargin(0.03);
  canvas->SetTopMargin(0.04);
  
  char ratioXaxis1[100];
  char MCindex[100];
  //float ymax;
  //canvas->SetBottomMargin(0.1); 
  data->GetYaxis()->SetLabelSize(0.055);
  data->GetXaxis()->SetLabelSize(0.0);
  data->GetYaxis()->SetTitleSize(0.065);
  data->GetYaxis()->SetTitleOffset(1.33);
  data->GetXaxis()->SetTitleSize(0.053);
  data->GetXaxis()->SetTitleOffset(0.12);
  data->GetYaxis()->CenterTitle();
  data->GetXaxis()->CenterTitle();
  data->GetXaxis()->SetTickLength(0.05);
  data->GetYaxis()->SetTickLength(0.05);

  data->SetLineWidth(2);
  data->SetMarkerStyle(9);
  data->SetMarkerSize(.9);
  int ifont =42;
  data->GetXaxis()->SetTitleFont(ifont);
  data->GetYaxis()->SetTitleFont(ifont);     
  data->SetTitleFont(ifont);
  
  //ymax = data->GetMaximum();
  //Devide the histogram with bin width
  
  
  TPad *padfun1 = new TPad("padfun1", "padfun1", 0.05, 0.35, 1.0, 1.0);
  padfun1->SetBottomMargin(0.0); // Upper and lower plot are joined
//padfun1->SetTopMargin(0.05); // Upper and lowd
  padfun1->SetLeftMargin(0.17); // Left
  padfun1->SetRightMargin(.02);
  padfun1->Draw();             // Draw the upper pad: pad1
  padfun1->cd();
  data->SetFillColor(kYellow);
  data->SetFillStyle(1111);
  
  //-----------------------------------create the legend of ht2 range 
  char MC_HTindex[100];
  //TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  TLegend *HT_range = new TLegend(plegend[0]-0.08,.06,plegend[2],.1);
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);
  HT_range->SetTextSize(0.06);
  sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  HT_range->Draw();
  data->SetTitle("");
  //--------------------------------------end the legend of ht2 range
  data->Draw("e2");
  
  
  //   data->Draw(" ");
  
  
  if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  
  double chi2rat[Nplot[0]];    // for chi2 plot in legend 
  double chi2Ndfrat[Nplot[0]]; // for chi2 plot in legend
  
  int  color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  int  style[40]={1,1,1,1,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
    MC[iup]->SetLineStyle(style[iup]);
    MC[iup]->SetLineColor(color[iup]);
    MC[iup]->SetLineWidth(3);
    //MC[iup]->Draw("same hist ");// gPad->SetLogy();    
    
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
    
    MC[iup]->Draw("same hist e1 ");
    
  }//------------------- end of Montecarlo loop
  
  
  data->Draw(" same e1");
  
  //------------------------------------------------------------------end of Uncertainty
  
  
  //  data->Draw("same");
  //TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
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
  legendn->Draw();
  
  HT_range->Draw();
  
  //----------------------------------------ratio plot pad
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *padfun2 = new TPad("padfun2", "padfun2", 0.05,0.02,1.0, 0.35);
  padfun2->SetTopMargin(0);
  padfun2->SetBottomMargin(.4);
  padfun2->SetLeftMargin(.17);
  padfun2->SetRightMargin(.02);
  padfun2->SetGridy(); // Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  //gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    TH1D *rh2;
    //if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    // if(data->GetBinContent(10) > 0) {
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
    rh2->SetMarkerStyle(20);
    rh2->SetMarkerStyle(kFullCircle);
    //rh2->Draw(" same hist");
    rh2->Draw(" same ");
    //rh2->Draw("same hist e1");
    
    
    rh2->GetXaxis()->SetTitleSize(0.12);
    //  rh2->GetXaxis()->SetTitleFont(43);
    rh2->GetXaxis()->SetTitleOffset(1.33);
    rh2->GetXaxis()->SetLabelSize(0.13);
    rh2->GetXaxis()->CenterTitle();
    sprintf(ratioXaxis1," %s" ,lowpadx);
    rh2->GetXaxis()->SetTitle(ratioXaxis1);
    
    rh2->GetYaxis()->SetTitle("MC/Data");
    rh2->GetYaxis()->CenterTitle();
    rh2->GetYaxis()->SetNdivisions(505);
    rh2->GetYaxis()->SetTitleSize(0.15);
    rh2->GetXaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleOffset(0.6);
    rh2->GetXaxis()->SetTickLength(0.06);
    rh2->GetYaxis()->SetTickLength(0.06);

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
  MyHist->SetStats(0);
  
  MyHist->GetXaxis()->SetLabelSize(0.03);
  MyHist->GetXaxis()->SetTitleSize(0.050);
  MyHist->GetXaxis()->SetTitleOffset(1.0);
  //MyHist->GetXaxis()->SetTitleFont(ifornt);
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetXaxis()->SetTitle(XTitle);
  
  MyHist->GetYaxis()->SetLabelSize(0.03);
  MyHist->GetYaxis()->SetTitleSize(0.050);
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


TCanvas *ratio_canV2(int Nplot[3],float plegend[8], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
  //Nplot[0] = number of MC enetered
  //Nplot[1] = place 1 if upper part is log scale needed
  //Nplot[2] = place 1 if MC/Data or 0 if data/MC
  //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  //plegend[4]= text size
  //plegend[5-6]= ratio plot axis range
  //plegend[7]= ratio plot Yaxis Title size
  //data = data histogram
  //MC = monte carlo histogram array
  //legendN = name of the legends for mC one by one
  //lowpadx = x axis title of the ratio plot
  //datanm[0] = name of data (or MC in case one MC is used)
  //datanm[1] = name of MC/data Data term
  //datanm[2] = name of Mc/Data MC term
  bool isstat =0;
  TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 680,750 );
  canvas->cd();
  canvas->SetRightMargin(0.02);
  canvas->SetTopMargin(0.04);

  char ratioXaxis1[100]; char ratioYaxis[100]; char MCindex[100];
  //float ymax;
  //canvas->SetBottomMargin(0.1);
  data->GetYaxis()->SetLabelSize(0.055);
  data->GetXaxis()->SetLabelSize(0.0);
  data->GetYaxis()->SetTitleSize(0.065);
  data->GetYaxis()->SetTitleOffset(1.33);
  data->GetXaxis()->SetTitleSize(0.053);
  data->GetXaxis()->SetTitleOffset(0.12);
  data->GetYaxis()->CenterTitle();
  data->GetXaxis()->CenterTitle();
  data->GetXaxis()->SetTickLength(0.05);
  data->GetYaxis()->SetTickLength(0.05);

  data->SetLineWidth(2);
  data->SetMarkerStyle(9);
  data->SetMarkerSize(.9);
  int ifont =42;
  data->GetXaxis()->SetTitleFont(ifont);
  data->GetYaxis()->SetTitleFont(ifont);
  data->SetTitleFont(ifont);

  //ymax = data->GetMaximum();
  //Devide the histogram with bin width

  TPad *padfun1 = new TPad("padfun1", "padfun1", 0.05, 0.35, 1.0, 1.0);
  padfun1->SetBottomMargin(0.0); // Upper and lower plot are joined
 // padfun1->SetTopMargin(0.05); // Upper and lowd
  padfun1->SetLeftMargin(0.17); // Left
  padfun1->SetRightMargin(.02);
  padfun1->Draw();             // Draw the upper pad: pad1
  padfun1->cd();
//  data->SetFillColor(kYellow);
//  data->SetFillStyle(1111);

  //-----------------------------------create the legend of ht2 range
  char MC_HTindex[100];
  //TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  TLegend *HT_range = new TLegend(plegend[0]-0.08,.06,plegend[2],.1);
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);
  HT_range->SetTextSize(0.06);
  sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  HT_range->Draw();
  data->SetTitle("");
  //--------------------------------------end the legend of ht2 range
  data->Draw("e2");

  if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();

  double chi2rat[Nplot[0]];    // for chi2 plot in legend
  double chi2Ndfrat[Nplot[0]]; // for chi2 plot in legend

  int  color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  int   style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
    MC[iup]->SetLineStyle(style[iup]);
    MC[iup]->SetLineColor(color[iup]);
    MC[iup]->SetLineWidth(3);
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

    MC[iup]->Draw("same hist e1 ");
    //MC[iup]->Draw("same hist");

  }//------------------- end of Montecarlo loop



   if(isstat){data->Draw(" same e1");} else
   {data->Draw("same ");}

 //------------------------------------------Stat Error(DATA)
   TH1D *datastaterr = (TH1D*)data->Clone("staterr"); datastaterr->Reset();

    for(int ierr=1; ierr < datastaterr->GetNbinsX()+1; ierr++){
     double bincont = data->GetBinContent(ierr);
     double serr = data->GetBinError(ierr);

     double relerr =sqrt(pow((serr/bincont),2.));
    datastaterr->SetBinContent(ierr,1.0);
    datastaterr->SetBinError(ierr,relerr);

    }
    gStyle->SetErrorX(0.5);
    gStyle->IsTransparent();
    datastaterr->SetFillColor(5);
    datastaterr->SetMarkerSize(0.1);
    datastaterr->SetFillStyle(3002);
    //datastaterr->SetFillStyle(1111);
    //datastaterr->SetFillColor(21);
    //datastaterr->SetMarkerSize(0.1);


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
  if(isstat)legendn->AddEntry(datastaterr, "Statistical Uncertainty","f"); //Stat unc
  legendn->Draw();
  HT_range->Draw();



  //----------------------------------------ratio plot pad
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *padfun2 = new TPad("padfun2", "padfun2", 0.05, 0.02,1.0, 0.35);
  padfun2->SetTopMargin(0);
  padfun2->SetBottomMargin(.4);
  padfun2->SetLeftMargin(.17);
  padfun2->SetRightMargin(.02);
  padfun2->SetGridy(); // Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  //gStyle->SetErrorX(0.5);
 
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop over MC for ratio plot
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
    rh2->SetMinimum(plegend[6]);  // Define Y ..
    rh2->SetMaximum(plegend[5]);  // .. range
    rh2->SetStats(0);      // No statistics on lower plot
    //gStyle->SetStatColor(2);      // No statistics on lower plot
    //gStyle->SetStatX(0.7);;      // No statistics on lower plot
    
    //rh2->SetMarkerStyle(21);
    rh2->SetMarkerStyle(kFullCircle); //For Closure Test
    rh2->Draw("same p e1");
    //rh2->Draw("Same P text ");
    //rh2->Draw("Same P ");
  
   //-----------------------------------------------------Stat Error Monte Carlo
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
    
   if(isstat) datastaterr->Draw("E2 Same");
   if(isstat) mstaterr->Draw("E2 Same");
    
    //MC_inputerr->Draw("same E2");//------------------Draw the Uncertainty in Ratio plot

    rh2->GetXaxis()->SetTitleSize(0.12);
    //  rh2->GetXaxis()->SetTitleFont(43);
    rh2->GetXaxis()->SetTitleOffset(1.33);
    rh2->GetXaxis()->SetLabelSize(0.13);
    rh2->GetXaxis()->CenterTitle();
    sprintf(ratioXaxis1," %s" ,lowpadx);
    rh2->GetXaxis()->SetTitle(ratioXaxis1);

   //rh2->GetYaxis()->SetTitle("#frac{Pythia8-unf}{Pythia8-gen}");
   //if(Nplot[2]==1){ sprintf(ratioYaxis,"%s/%s" ,datanm[1],datanm[2]);}
   //if(Nplot[2]==0){ sprintf(ratioYaxis,"%s/%s" ,datanm[2],datanm[1]);}

    if(Nplot[2]==1){ sprintf(ratioYaxis,"#frac{%s}{%s}" ,datanm[1],datanm[2]);}
    if(Nplot[2]==0){ sprintf(ratioYaxis,"#frac{%s}{%s}" ,datanm[2],datanm[1]);}

    rh2->GetYaxis()->SetTitle(ratioYaxis);
    rh2->GetYaxis()->CenterTitle();
    rh2->GetYaxis()->SetNdivisions(505);
    rh2->GetYaxis()->SetTitleSize(plegend[7]);
    rh2->GetXaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleFont(ifont);
    rh2->GetYaxis()->SetTitleOffset(0.6);
    rh2->GetXaxis()->SetTickLength(0.06);
    rh2->GetYaxis()->SetTickLength(0.06);

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

for (int iht = 1; iht <= nbiny; ++iht) {
       const char * name_hbin = Form("%s_pt%i",VarHTin->GetName(), iht);
        TH1D* h_var  = VarHTin->ProjectionX(name_hbin, iht, iht);
        Var[iht-1] =(TH1D*)h_var->Clone();        
        //h_var-Write();
      }
}

//-----------------------for Ratio plot-------------------------
// xtilte, ytitle, max, min, linewidth, SetLineStyle, SetMarkerStyle, SetMarkerSize, SetLineColor 
void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl){
  
  MyHist->GetXaxis()->SetTitle(XT); MyHist->GetYaxis()->SetTitle(YT);
  MyHist->SetMinimum(min); MyHist->SetMaximum(mx);
  MyHist->SetLineWidth(ilw); MyHist->SetLineStyle(ilsty);
  MyHist->SetMarkerStyle(imsty); MyHist->SetMarkerSize(imstysize); MyHist->SetLineColor(icl);
  MyHist->SetTitle("");
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetYaxis()->CenterTitle();

  MyHist->GetXaxis()->SetLabelSize(0.03);
  MyHist->GetXaxis()->SetTitleSize(0.045);
  MyHist->GetXaxis()->SetTitleOffset(1.0);
  MyHist->GetYaxis()->SetLabelSize(0.03);
  MyHist->GetYaxis()->SetTitleSize(0.040);
  MyHist->GetYaxis()->SetTitleOffset(1.0);
}

void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1){

 TH1D *chidata = (TH1D*)data->Clone("chidata");
 TH1D *chiMC = (TH1D*)MC->Clone("chimc");
chidata->Rebin(rebin);
chiMC->Rebin(rebin);
int n = chidata->GetNbinsX();
Double_t res[n] , chi2;
Int_t ndf ,igood;

 cout << " Numbers of bin = " << n << endl;


 chiMC->Chi2TestX(chidata, chi2, ndf, igood, "WU", res);

cout << "Root chi2/ndf = " << chi2 << " / " << ndf << " = " << (chi2/ndf) << endl;

}

