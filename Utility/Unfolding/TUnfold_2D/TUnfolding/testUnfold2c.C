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

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldIterativeEM.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>



#define NOREG
//#define TIKH_SURE
//#define TIKH_LSCAN
//#define ITER_SURE
//#define ITER
//#define REBIN
//#define UOEx
#define CLOUSER

using namespace std;
static const auto feps = numeric_limits<float>::epsilon();

void testUnfold2c(){
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();
  
  int const nmc=3; //Number of MC -> 0,1,2 PY8, MG, HW7
  int const umc=0; //Number of MC
  int irbin = 1; //Rebin
  int const untype = 1;  // 0 for 2D , 1 for 2D 
  
  TFile *inputbinning = new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  //Input Data and MC histogram
  TFile *inputData=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  //TFile *inputData=new TFile("PY8_UL17_Flat_2D_16Sep20.root");
  
  TFile *inputMC[nmc];
  //TFile *inputMC=new TFile("Test_MC_QCD.root");
  inputMC[0]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  //TFile *inputMC=new TFile("PY8_UL17_Flat_2D_16Sep20.root");
  
  //TFile *inputMC1=new TFile("MG_UL17_binned_6July20.root");
  inputMC[1]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  inputMC[2]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  //TFile *inputMC2=new TFile("HW7_Flat_Binned_8July20.root");
  
  //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  TFile *outputFile=new TFile("Unfolded_Result.root","recreate");
 
  
  ofstream file;
  file.open("Tau_Value.txt");
  file <<"L Curve    "  <<"                   "<< "Scan Sure  "  <<"                "<<" Scan Tau "<<endl;
  
  int const type = 2;         // Jet & Charage particles
  int const itype[type]={0,1};   //{0}--->Jet ; {1}---> Charged Particles
  const  char* itypeN[type]={"Jets","Charged Particles"};
  const char* DirName[3] = {"analyzeBasicPat", "analyzeBasicPat1D", "Binning"};

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
  
  //-------------------------input Histograms, RM 
  TH1D *MC_Reco[nmc][type][nusedvar][njetptmn];  //Reconstructed MC
  TH1D *MC_fake[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_fakerate[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_background[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  
  TH1D *MC_Gen[nmc][type][nusedvar][njetptmn];   //Generator MC
  TH1D *MC_miss[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *MC_missmiss[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *MC_misscorr[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  
  TH1D *PsudoData_Gen[nmc][type][nusedvar][njetptmn];    //Gen Level Psudo Data
  TH2D *h2dGenDetMC[nmc][type][nusedvar][njetptmn];   // MC generator Vs Reco
  TH1D *RMX[nmc][type][nusedvar][njetptmn];  //Reconstructed MC
  TH1D *RMY[nmc][type][nusedvar][njetptmn];  //Reconstructed MC
  
  TH1D *Data_Reco[type][nusedvar][njetptmn];    //Reconstructed Data
  
  TUnfoldBinning *binsRec[type][nusedvar][njetptmn];
  TUnfoldBinning *RecoBinning[type][nusedvar][njetptmn];
  
  TUnfoldBinning *binsGen[type][nusedvar][njetptmn];
  TUnfoldBinning *GenBinning[type][nusedvar][njetptmn];
  
  
  
  TH1D* hist_eff[type][nusedvar][njetptmn];
  TH1D* hist_fake[type][nusedvar][njetptmn];
  TH1D* hist_purity[type][nusedvar][njetptmn];
  TH1D* hist_stbl[type][nusedvar][njetptmn];
  
  
  
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
  
  
  TDirectoryFile *DirData = new TDirectoryFile("Data","Inputs Data");
  
  TDirectoryFile *inputDir[nmc];
  inputDir[0]=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
  inputDir[1]=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
  inputDir[2]=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
  
  
  TDirectoryFile *foldpy8=new TDirectoryFile("fold","folded with Probablility matrix");
  TDirectoryFile *Unfold=new TDirectoryFile("Unfold","Unfolded, Refold, correlation");
  
  int subtract_background1(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
  int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
  void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
  void Condition (TH2 * RM, TH1* miss);
  TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
  TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco,TH1D* MC_gen, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8],int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8] );
  // TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar );
  TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
 
  //----------------------------------------------Read 2D Binning 
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	sprintf(histname, "%s/detector_typ_%i_pt%i_eta0_%i", DirName[2],ity, ipt, var[ivar]); 
	cout << histname <<endl;
	inputbinning->GetObject(histname, binsRec[ity][ivar][ipt]);
	binsRec[ity][ivar][ipt]->PrintStream(cout);
	sprintf(histname, "%s/Generator_typ_%i_pt%i_eta0_%i", DirName[2],  ity, ipt, var[ivar]); 
	inputbinning->GetObject(histname, binsGen[ity][ivar][ipt]);
	binsGen[ity][ivar][ipt]->PrintStream(cout);
      }
    }
  }
  
  
  //Read Input Data MC and Response matrix
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	
	sprintf(histname, "%s/reco_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
	TH1D *RecoData =(TH1D*) inputData->Get(histname);
	//for(int i=1 ; i< RecoData->GetNbinsX()+1; i++){ if(RecoData->GetBinContent(i) == 0) {cout << " Data Bin entry Nil : bin no : "<< i << endl;}}

#ifdef CLOUSER
	sprintf(histname, "%s/gen_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
	TH1D *PsudoDataGen= (TH1D*)inputData->Get(histname);
#else
	TH1D *PsudoDataGen = (TH1D*)inputMC[0]->Get(histname);
#endif
	
for (int imc=0; imc<nmc ; imc++){
	//----------------------------------MC RECO
       sprintf(histname, "%s/reco_typ_%i_pt%i_eta0_%i",DirName[untype],ity, ipt, var[ivar]); 
       TH1D *RecoMC = (TH1D*)inputMC[imc]->Get(histname);      cout << histname ;
      
 //      if(recobins!=Data_Reco[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
       
       //-----------------------------------MC Fake
       sprintf(histname, "%s/fake_reco_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
       TH1D *RecoMC_Fake = (TH1D*)inputMC[imc]->Get(histname);
       
       cout << " Fake= " <<RecoMC_Fake->GetEntries() <<" Reco-fake: " <<(RecoMC->GetEntries() - RecoMC_Fake->GetEntries())<<endl;
       
       //-----------------------------------Gen MC
       sprintf(histname, "%s/gen_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
       TH1D *GenMC = (TH1D*)inputMC[imc]->Get(histname); cout << histname  ;
       
       //for(int i= 1; i < GenMC->GetNbinsX()+1; i++){if(GenMC->GetBinContent(i) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  i  << endl; }}
       
       //----------------------------------MC miss 
       sprintf(histname, "%s/miss_gen_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
       TH1D *GenMC_Miss = (TH1D*)inputMC[imc]->Get(histname);

       cout << " Miss= " << GenMC_Miss->GetEntries() <<" Gen-Miss: " << (GenMC->GetEntries() - GenMC_Miss->GetEntries())<< endl ;
       
       //Response Matrix
       sprintf(histname, "%s/corr_typ_%i_pt%i_eta0_%i", DirName[untype], ity, ipt, var[ivar]); 
       TH2D *RM_RecoGen= (TH2D*) inputMC[imc]->Get(histname);  cout << " Corr = "  << RM_RecoGen->GetEntries() <<endl;

#ifdef UOEx
       TH1D *NewReco; TH1D *NewGen; TH1D *Newmiss; TH1D *Newfake; TH2D *NewRecoGen;
       NewReco = (TH1D*)RecoMC->Clone(); Newfake = (TH1D*)RecoMC_Fake->Clone();
       NewGen  = (TH1D*)GenMC->Clone(); Newmiss = (TH1D*)GenMC_Miss->Clone();
       NewRecoGen = (TH2D*)RM_RecoGen->Clone();

       int NbinxR = NewReco->GetNbinsX(); int NbinxG = NewGen->GetNbinsX(); int NbinMx = NewRecoGen->GetNbinsX(); int NbinMy = NewRecoGen->GetNbinsY();
       NewReco->Reset(); NewGen->Reset(); NewRecoGen->Reset(); Newmiss->Reset(); Newfake->Reset();

       for(int ix=0; ix < NbinxR+2 ; ix++){
         NewReco->SetBinContent(ix,NewReco->GetBinContent(ix));
         NewReco->SetBinError(ix, sqrt(NewReco->GetBinError(ix)* NewReco->GetBinError(ix))); }
       for(int ix=0; ix < NbinxG+2 ; ix++){
         NewGen->SetBinContent(ix, NewGen->GetBinContent(ix));
         NewGen->SetBinError(ix, sqrt(NewGen->GetBinError(ix)* NewGen->GetBinError(ix))); }
       for(int ix=0; ix < NbinMx+2 ; ix++){
         for(int iy=0; iy < NbinMy+2 ; iy++){
           NewRecoGen->SetBinContent(ix, iy, NewRecoGen->GetBinContent(ix,iy));
           NewRecoGen->SetBinError(ix, iy, sqrt(NewRecoGen->GetBinError(ix,iy)*NewRecoGen->GetBinError(ix,iy)));
         }
       }
       for(int ix=0; ix < NbinxR+2 ; ix++){
         Newfake->SetBinContent(ix, Newfake->GetBinContent(ix));
         Newfake->SetBinError(ix, sqrt(Newfake->GetBinError(ix)* Newfake->GetBinError(ix)));}
       for(int ix=0; ix < NbinxG+2 ; ix++){
         Newmiss->SetBinContent(ix, Newmiss->GetBinContent(ix));
         Newmiss->SetBinError(ix, sqrt(Newmiss->GetBinError(ix)* Newmiss->GetBinError(ix))); }
       
       RecoMC = (TH1D*)NewReco->Clone();
       RecoMC_Fake = (TH1D*)Newfake->Clone();
       GenMC = (TH1D*)NewGen->Clone();
       GenMC_Miss = (TH1D*)Newmiss->Clone();
       
       RM_RecoGen = (TH2D*)NewRecoGen->Clone();

#endif 
//--------------------------------------Fake rate and Miss Rate
       TH1D* fakerate = (TH1D*)RecoMC_Fake->Clone();
       sprintf(name,"fake_rate_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); fakerate->SetNameTitle(name,name);
       fakerate->Divide(RecoMC_Fake, RecoMC, 1, 1, "b");
       TH1D* missrate = (TH1D*)GenMC_Miss->Clone();
       sprintf(name,"miss_rate_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); missrate->SetNameTitle(name,name);
       missrate->Divide(missrate, GenMC_Miss, 1, 1, "b");




//--------------------------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20
       TH1* RMx = RM_RecoGen->ProjectionX(); TH1* RMy = RM_RecoGen->ProjectionY();
        
       sprintf(name,"ProjectX_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); RMx->SetNameTitle(name,name);

       sprintf(name,"Recominusfake_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       TH1* RecoFakeCorrect = (TH1D*)RecoMC->Clone(); RecoFakeCorrect->Reset();
       RecoFakeCorrect->SetNameTitle(name,name);

       sprintf(name,"ProjectY_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); RMy->SetNameTitle(name,name);

       sprintf(name,"Genminusmiss_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
       TH1* GenMissCorrect = (TH1D*)GenMC->Clone(); GenMissCorrect->Reset();
       GenMissCorrect->SetNameTitle(name,name);
       
       for (int i = 1; i <= RecoFakeCorrect->GetNbinsX(); ++i) {
         double content = RecoMC->GetBinContent(i); double factor = RecoMC_Fake->GetBinContent(i);
         content -= factor;  RecoFakeCorrect->SetBinContent(i, content);
       }

       for (int i = 1; i <= GenMissCorrect->GetNbinsX(); ++i) {
        double content = GenMC->GetBinContent(i); double factor = GenMC_Miss->GetBinContent(i);
         content -= factor;  GenMissCorrect->SetBinContent(i, content);
       }
       
       inputDir[imc]->cd();
       RMx->Write(); RMy->Write(); RecoFakeCorrect->Write(); GenMissCorrect->Write();

       RecoMC->Write(); RecoMC_Fake->Write(); fakerate->Write();
       GenMC->Write(); GenMC_Miss->Write(); missrate->Write(); 
       RM_RecoGen->Write();
//----------------------------------------------------------------------------
       MC_Reco[imc][ity][ivar][ipt] = (TH1D*)RecoMC->Clone();
       MC_fake[imc][ity][ivar][ipt] = (TH1D*)RecoMC_Fake->Clone();
       MC_fakerate[imc][ity][ivar][ipt] = (TH1D*)fakerate->Clone();
       MC_background[imc][ity][ivar][ipt] = (TH1D*)fakerate->Clone();
       MC_fakerate[imc][ity][ivar][ipt]->SetMinimum(-0.05); MC_fakerate[imc][ity][ivar][ipt]->SetMaximum(1.01);
       

       MC_Gen[imc][ity][ivar][ipt] = (TH1D*)GenMC->Clone();
       MC_miss[imc][ity][ivar][ipt] = (TH1D*)GenMC_Miss->Clone();
       MC_missmiss[imc][ity][ivar][ipt] = (TH1D*)missrate->Clone();
       MC_misscorr[imc][ity][ivar][ipt] = (TH1D*)missrate->Clone();
       MC_missmiss[imc][ity][ivar][ipt]->SetMinimum(-0.05); MC_missmiss[imc][ity][ivar][ipt]->SetMaximum(1.01);

       h2dGenDetMC[imc][ity][ivar][ipt] = (TH2D*)RM_RecoGen->Clone();

} //for (int imc=0;imc<nmc;imc++)
       
//--------------------------------------------A kind of way to remove underflow and overlow bins
#ifdef UOEx
       TH1D *NewData;  NewData = (TH1D*)RecoData->Clone();
       int NbinxD = NewData->GetNbinsX(); NewData->Reset(); 
       for(int ix=1; ix < NbinxD+1 ; ix++){
         NewData->SetBinContent(ix,NewData->GetBinContent(ix));
         NewData->SetBinError(ix, sqrt(NewData->GetBinError(ix)* NewData->GetBinError(ix))); }
       RecoData = (TH1D*)NewData->Clone();
#endif
       DirData->cd();       
       RecoData->Write();
       Data_Reco[ity][ivar][ipt] = (TH1D*)RecoData->Clone();
      }
    }
  }
  cout << "Histogram Read for Pythia8 Done " <<endl;
  
  //------------------Fold check : Patrick 1 Sep 20
  //Get Probability Matrix  
  //(Multiply the gen level by the probability matrix. Of course, don't forget to account for miss and fake entries (if applicable).
  foldpy8->cd();
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	TH2D* RM  = (TH2D*)h2dGenDetMC[umc][ity][ivar][ipt]->Clone();
	TH1D* Reco  = (TH1D*)MC_Reco[umc][ity][ivar][ipt]->Clone();
	TH1D* Gen  = (TH1D*)MC_Gen[umc][ity][ivar][ipt]->Clone();
	TH1D* fake= (TH1D*)MC_fake[umc][ity][ivar][ipt]->Clone();
	TH1D* miss = (TH1D*)MC_miss[umc][ity][ivar][ipt]->Clone();
	
	RM->RebinY(irbin);Gen->Rebin(irbin);miss->Rebin(irbin);
	
	TH1D* Folded = (TH1D*)MC_Reco[umc][ity][ivar][ipt]->Clone(); Folded->Reset();
	
	Fold(RM, Reco, Gen, miss, fake, Folded);
	sprintf(name,"Fold_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
	Folded->SetNameTitle(name,name);
	
	Folded->Write();
      }
    }
  }
  
//----------------------------Condition  number of Probability Matrix
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	TH2D* RM  = (TH2D*)h2dGenDetMC[umc][ity][ivar][ipt]->Clone();
	TH1D* miss = (TH1D*)MC_miss[umc][ity][ivar][ipt]->Clone();
	RM->RebinY(irbin); miss->Rebin(irbin);
	cout <<setw(2) << ity <<setw(5) <<ivar << setw(5) << ipt <<'\n';
	Condition(RM, miss);
      }
    }
  }
  
  
  //--------------------------------------------------------------
  Unfold->cd();
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
        
	sprintf(name,"Purity_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
	hist_purity[ity][ivar][ipt] = (TH1D*)MC_Reco[umc][ity][ivar][ipt]->Clone(); hist_purity[ity][ivar][ipt]->Reset();
	hist_purity[ity][ivar][ipt]->SetNameTitle(name,name);
	sprintf(name,"stability_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);
	hist_stbl[ity][ivar][ipt] =  (TH1D*)MC_Reco[umc][ity][ivar][ipt]->Clone(); hist_stbl[ity][ivar][ipt]->Reset();
	hist_stbl[ity][ivar][ipt]->SetNameTitle(name,name);
	
      }
    }
  }
  
  //cross check efficincy and purity calculation (Tunfold example code)
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){

   /*   	      //----------------------------------------
	for(int binGen=0;binGen<= h2dGenDetMC[umc][ity][ivar][ipt]->GetNbinsY()+1;binGen++) {
	  double sum0=0.;
	  double sum1=0.;
	  for(int binRec=0;binRec<= h2dGenDetMC[umc][ity][ivar][ipt]->GetNbinsX()+1;
	      binRec++) {
	    //double c=  h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	    double c=  h2dGenDetMC[umc][ity][ivar][ipt]->GetBinContent(binRec,binGen);
	    sum0+=c;
	    if((binRec>0)&&(binRec<=h2dGenDetMC[umc][ity][ivar][ipt]->GetNbinsX())) {
	      sum1+=c;
	    }
	  }
	  if(sum0>0.0) {
	    hist_eff[ity][ivar][ipt]->SetBinContent(binGen,sum1/sum0);
	  }
	}
	//---------------------------
	hist_eff1[ity][ivar][ipt]->SetMinimum(0.9); hist_eff1[ity][ivar][ipt]->SetMaximum(1.1);
	hist_eff1[ity][ivar][ipt]->Write();
*/	
	
	for(int binRec=0; binRec<= hist_purity[ity][ivar][ipt]->GetNbinsX()+1; binRec++) {
	  double sum=0.;
	  for(int binGen=0; binGen<=hist_purity[ity][ivar][ipt]->GetNbinsX()+1; binGen++) {
	    //sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	    sum += h2dGenDetMC[umc][ity][ivar][ipt]->GetBinContent(binRec,binGen);
	  }
	  double p=0.;
	  if(sum>0.0) {
	    p = h2dGenDetMC[umc][ity][ivar][ipt]->GetBinContent(binRec,binRec)/sum;
	  }
	  hist_purity[ity][ivar][ipt]->SetBinContent(binRec,p);
	}
	
	hist_purity[ity][ivar][ipt]->SetMinimum(-0.05); hist_purity[ity][ivar][ipt]->SetMaximum(1.01);
	hist_purity[ity][ivar][ipt]->Write();
      }
    }
  }
  
  double bias[5]={0.0,0.0,0.0,0.0,0.0};
  //double bias[5]={1.0,1.0,1.0,1.0,1.0};
  //int bias[5]={2.0,2.0,2.0,2.0,2.0};
  
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	
	//if (ivar==2) continue; 
	cout <<"type "<< ity << " : Variables " << ivar << " HT2 Bin : " << ipt << endl;
	file <<"["<< ity << "," << ivar << "," << ipt <<"] --->" << endl;
	
        
	//Get reco bins
	double rxbins[MC_Reco[umc][ity][ivar][ipt]->GetNbinsX()+1]={0};
	for(int ix=0; ix<MC_Reco[umc][ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	  rxbins[ix] = MC_Reco[umc][ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1); }
	
	//Get Gen bins
	double gxbins[MC_Gen[umc][ity][ivar][ipt]->GetNbinsX()+1]={0};
	for (int ix=0; ix<MC_Gen[umc][ity][ivar][ipt]->GetNbinsX()+1; ix++) {
	  gxbins[ix] = MC_Gen[umc][ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1); }
	
	//Rebin for match the condition of reco vs gen bin 
	h2dGenDetMC[umc][ity][ivar][ipt]->RebinY(irbin);
	MC_Gen[umc][ity][ivar][ipt]->Rebin(irbin);
	//h2dGenDetMC[ity][ivar][ipt]->Rebin(1,2);
	
	
	
	
	TH2* hist_migrationCoarseFine_MC = (TH2D*)h2dGenDetMC[umc][ity][ivar][ipt]->Clone();
	TH1* input = (TH1D*)Data_Reco[ity][ivar][ipt]->Clone();
	TH1* mcgen = (TH1D*)MC_Gen[umc][ity][ivar][ipt]->Clone();
        TH1* mcgenmiss=(TH1D*)MC_misscorr[umc][ity][ivar][ipt]->Clone();	
	
	//correction for Fake as background subtraction : Patrick
	TH1* mcbackground = (TH1D*)MC_background[umc][ity][ivar][ipt]->Clone();
	mcbackground->Multiply(input);
	
	double biasScale =bias[ivar] ;
	const char *REGULARISATION_DISTRIBUTION=0;
	const char *REGULARISATION_AXISSTEERING="*[UOB]";
	char Rhoname[100], Rhotitle[100], Lcursure[100], Lcurtitle[100], TauLsure[100], TuaLtitle[100], probMat[100], probmat_title[100],Rhoname2d[100],
	  Rhotitle2d[100],Ematrix[100],Ematrixtitle[100],foldback[100],foldback_title[100];
	
	// preserve the area
	//TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
	TUnfoldDensity::EConstraint constraintMode= TUnfoldDensity::kEConstraintArea;
	TUnfoldDensity::ERegMode regMode = TUnfoldDensity::kRegModeCurvature;
	TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
	
	//https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
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
	
	//TUnfoldDensity class
	//No  regularisation --------------------------------
	TUnfoldDensity tunfoldNoRegularisation(hist_migrationCoarseFine_MC,
					       TUnfold::kHistMapOutputVert,
					       TUnfoldDensity::kRegModeNone, 
					       TUnfoldDensity::kEConstraintNone,
					       TUnfoldDensity::kDensityModeNone);
					      // binsRec[ity][ivar][ipt],
					      // binsGen[ity][ivar][ipt]);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
	
	// TUnfoldDensity::kDensityModeBinWidthAndUser);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
	
	tunfoldNoRegularisation.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 
	int status = tunfoldNoRegularisation.SetInput(input,biasScale,0,covarianceM);
	
	int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
	cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;
	//tunfoldNoRegularisation.SubtractBackground(mcbackground,"Background", 1.0,0.04); // hist,name,scale, scale error 
	
	//if(>=10000) { std::cout<<"Unfolding result may be wrong\n";  }
	
	
	//the initial bias vector is determined from the response matrix
	//but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
	//tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result
	
	//if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeBins(7,1,4,regMode);}  //Test for Regularization in few  unmatched bins
	//if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeCurvature(10,9,11,1.0,0.6);}
	
	//Choose a value of tau to unfold with tau,0 means no regularization
	tunfoldNoRegularisation.DoUnfold(0.0);//,input,biasScale);
	
         //Binmaps : Thanks to Suman(Tifr)
	  // int binnos = rnbinsx[ity][ivar][ipt];
	  // Int_t *binMap=new Int_t[binnos+2];
	  // for(Int_t i=1;i<=binnos;i++) binMap[i]=i;
	  // binMap[0]=-1;   binMap[binnos+1]=-1;
	
	
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
	TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
	//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,"","*[UO]");//,"signal");
	//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput("hist_PTunfolded_noRegularisation", "P_{T,unfolded} [GeV]","signal");
	
	TH2 *hist_RhoIJ_noRegularisation = tunfoldNoRegularisation.GetRhoIJtotal(NoReg_RhoIJ, NoReg_RhoIJ_tit);//,"signal");
	TH2 *hist_Rho2D_noRegularisation = tunfoldNoRegularisation.GetEmatrixTotal(NoReg_Ematrix, NoReg_Ematrix_tit);//,"signal");
	//tunfoldNoRegularisation.GetEmatrix(COV_Mat_NoReg[ity][ivar][ipt],binMap);//,"signal");
	TH1 *hist_foldedback_NoReg = tunfoldNoRegularisation.GetFoldedOutput(foldback, foldback_title);//,"signal");
	TH2 *hist_prob_Noreg = tunfoldNoRegularisation.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");
	
	// correction for Miss entries  : Partick 
	for (int i = 1; i <= hist_PTunfolded_noRegularisation->GetNbinsX(); ++i) {
	  double content = hist_PTunfolded_noRegularisation->GetBinContent(i);
	  double factor = 1;
	  factor += mcgenmiss->GetBinContent(i);
	  content *= factor;
	  hist_PTunfolded_noRegularisation->SetBinContent(i, content);
	}
	
	  //Quick check for closure Ratio Plots
	   // hist_PTunfolded_noRegularisation->Scale(1/(hist_PTunfolded_noRegularisation->Integral()));
	   // mcgen->Scale(1/mcgen->Integral());
	   // hist_PTunfolded_noRegularisation->Divide(mcgen);
	   // hist_PTunfolded_noRegularisation->SetMinimum(0.85); hist_PTunfolded_noRegularisation->SetMaximum(1.15);
	
	
	hist_PTunfolded_noRegularisation->Write();
	hist_Rho2D_noRegularisation->Write();
	hist_RhoIJ_noRegularisation->Write();
	hist_foldedback_NoReg->Write();
	hist_prob_Noreg->Write();
	
	cout << "Without Regularization finish " << endl;
	//-----------------------------------------------------------------------// L curve Scan
	
	TUnfoldDensity tunfoldTikhonovLCurve
	  (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputVert,regMode, constraintMode,densityFlags,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
	tunfoldTikhonovLCurve.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 
	tunfoldTikhonovLCurve.SetInput(input,biasScale);//,0.,covarianceM);
	// tunfoldTikhonovLCurve.SetBias(mcgen);
	int iBest_TikhonovLCurve=-1;
	double tauBest_TikhonovLCurve=-1.,DF_TikhonovLCurve=-1.;
	int iBest;
	int NPOINT_TikhonovLCurve=50;//200;
	TGraph *graph_LCurve_TikhonovLCurve;
	TSpline *spline_Curvature_TikhonovLCurve;
	taumx = 1e-5; taumi = 1e-9;
	//taumx = 0.0; taumi = 0.0;
	
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
	//TH1 *hist_PTunfolded_TikhonovLCurve = tunfoldTikhonovLCurve.GetOutput(unfoldhist, title);
	TH1 *hist_Rho_TikhonovLCurve = tunfoldTikhonovLCurve.GetRhoItotal(Rhoname, Rhotitle);//,"signal");
	TH2 *hist_RhoIJ_TikhonovLCurve = tunfoldTikhonovLCurve.GetRhoIJtotal(Rhoname2d, Rhotitle2d);//,"signal");
	TH2 *hist_Ematrix_TikhonovLCurve = tunfoldTikhonovLCurve.GetEmatrixTotal(Ematrix, Ematrixtitle);//"*[UO]");//,"signal");
	//TH2 *hist_Ematrix_TikhonovLCurve = tunfoldTikhonovLCurve.GetEmatrix(EmatrixL, EmatrixtitleL);//,"signal");
	TH2 *hist_prob_TikhonovLCurve = tunfoldTikhonovLCurve.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");
	TH1 *hist_foldedback_Lcurve= tunfoldTikhonovLCurve.GetFoldedOutput(foldback, foldback_title);//,"signal");
	
	
	// correction for Miss entries
	for (int i = 1; i <= hist_PTunfolded_TikhonovLCurve->GetNbinsX(); ++i) {
	  double content = hist_PTunfolded_TikhonovLCurve->GetBinContent(i);
	  double factor = 1;
	  factor += mcgenmiss->GetBinContent(i);
	  content *= factor;
	  hist_PTunfolded_TikhonovLCurve->SetBinContent(i, content);
	}
	
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
       
       //-----------------------------------------End of Variables loop
     }
   }
 }
 
 
  delete outputFile;
  file.close();
}




//RooUnfold Subtract Background (Fake, miss ) method : Now replace by patrick code
int subtract_background1(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl) {  
  int nbinx = h2d_correl->GetNbinsX();
  int nbiny = h2d_correl->GetNbinsY();
  const int nbinmx = 100 ;
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};
  
  for (int ix=0; ix<nbinx+1; ix++) {
    for (int iy=0; iy<nbiny+1; iy++) {
      if(ix==0&&iy==0) continue ;
      totalreco[ix] +=h2d_correl->GetBinContent(ix, iy);
      if (iy==0) fakerate[ix-1] = h2d_correl->GetBinContent(ix,iy);      //fake from GEN Underflow bins
      totalgen[iy] +=h2d_correl->GetBinContent(ix, iy);
      if (ix==0) effi[iy-1] =h2d_correl->GetBinContent(ix, iy);         //miss from RECO Underflow bin 
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
    
    // if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
  }//ix
  
  for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
    reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
    reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
  }
  
  for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
    data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
    data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
  }
  
  for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
    // gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]);
    // gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]);
  }
  
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
  return 0;
} //end substract func


//Copy of RooUnfold Subtract Background method with no correction to Reco or Gen --->Can be used  for purity and stability
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
      
    }//iy
  }//ix
  
  for (int iy=0; iy<nbiny; iy++) {
    effi[iy] = (totalgen[iy+1] - effi[iy])/max(1.e-10, totalgen[iy+1]);
    
    // if(iy>10 && effi[iy]>0.1 && effi[iy]<0.9){
    //       effi[iy] = effi[iy-1]-0.001;
    // }
    
    //  if (effi[iy]<1.e-6) effi[iy]=1.e-6;
  } //iy
  
  for (int ix=0; ix<nbinx; ix++){
    if(totalreco[ix+1]>1.e-6) {     fakerate[ix] /= totalreco[ix+1]; }
    else { fakerate[ix] = 0.; }
    
    // if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
  }//ix
  
  /// for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
  //        reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
  //        reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
  // }
  
  // for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
  //       data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
  //        data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
  // }
  
  //  for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
  //  gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]) ;
  //  gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]) ;
  //  }
  
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
  return 0;
} //end substract func



//1D Rebin function:Basically it cut off some bin you  don't needed
TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]){
  TH1D* thout;
  int nxmod2 = -1; double xmod2[300];
  const char *namex; const char *titlex;
  char namey[100], titley[100];
  
  double yvl[200]={0.0};
  double erryvl[200]={0.0};
  int nbinx = thin->GetNbinsX();
  int ifirst = 0;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
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
  for(int ii=0; ii < nxmod2+1 ; ii++){ xmod3[ii]=xmod2[ii]; }
  
  namex = thin->GetName();
  sprintf(namey, "rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "Rebin_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod3);
  
  for (int ix=0; ix<nxmod2+1; ix++) {
    thout->SetBinContent(ix+1, yvl[ix]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix]));     }
  //thout->SetBinContent(0,5);
  return thout;
}   //end of rebin_hist functiomC


//2D Rebin function:Basically it cut off some bin you  don't needed
TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco, TH1D* MC_gen, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8],int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8] ){
  TH2D* thout;
  cout << "Reco rebin:" <<endl;
  double xmod2[200]; double ymod2[200];
  const char* namex; const char* titlex;
  char namey[100],titley[100];
  
  double yvl[200][200]={0.0};
  double erryvl[200][200]={0.0};
  int nbinx = thin->GetNbinsX();
  int nbiny = thin->GetNbinsY();
  cout << "   X bin :=" << nbinx << "   Y bin  : = " <<nbiny <<endl; 
  int ifirst = 0; int jfirst =0;
  int  nreco= -1; int  ngen= -1;
  int xfirst =0;int yfirst =0;

 //Get Reco Level entries
  if (nreco<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
      } else {
        xmod2[ifirst] = MC_reco->GetBinLowEdge(ij);
        xmod2[ifirst+1] = MC_reco->GetBinLowEdge(ij+1);
        ifirst++;
      }
    }
  }
//Get Gen Level entries
  if (ngen<0) {
    for (int ij=0; ij<nbiny+2; ij++) {
      if (ij <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( ij > (nbiny - arrayvar_lastG[itype][ivar][ijetpt])) {
      } else {
        ymod2[jfirst] = MC_gen->GetBinLowEdge(ij);
        ymod2[jfirst+1] = MC_gen->GetBinLowEdge(ij+1);
        jfirst++;
      }
    }
  }
  if (xfirst==0 && yfirst==0) {
    for (int iy=0; iy<nbiny+2; iy++) {
      if (iy <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( iy > (nbiny - arrayvar_lastG[itype][ivar][ijetpt]) ) {
      } else {
	xfirst=0;
      for (int ix=0; ix<nbinx+2; ix++) {
        if (ix <=arrayvar_first[itype][ivar][ijetpt]) {
      } else if (ix > (nbinx - arrayvar_last[itype][ivar][ijetpt])) {
      } else {
        yvl[xfirst][yfirst] = thin->GetBinContent(ix,iy);
        erryvl[xfirst][yfirst] += thin->GetBinError(ix,iy)*thin->GetBinError(ix,iy);
        xfirst++;
        }
      }
      yfirst++;
    }
  }
}

cout << "Reco Bins = " << ifirst <<"  " << xfirst << endl;
cout << "Gen Bins = " << jfirst <<"  " << yfirst << endl;

    nreco = ifirst;
    ngen = jfirst;
double xmod3[nreco+1];
for(int ii=0; ii < nreco+1 ; ii++){ xmod3[ii]=xmod2[ii]; }
double xmod4[ngen+1];
for(int ii=0; ii < ngen+1 ; ii++){  xmod4[ii]=ymod2[ii]; }

  namex = thin->GetName();
  sprintf(namey, "Rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "%s", titlex);
  thout = new TH2D(namey, titley, nreco, xmod3, ngen, xmod4);

  for (int ix=0; ix<nreco+1; ix++) {
  for (int iy=0; iy<ngen+1; iy++) {
    thout->SetBinContent(ix+1,iy+1, yvl[ix][iy]);
    thout->SetBinError(ix+1,iy+1, sqrt(erryvl[ix][iy]));
        }
     }
  return thout;
}   //end of rebin_hist functiomC


//1D Rebin function:Basically it cut off some bin you  don't needed
TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8]){
TH1D* thout;
int nxmod2 = -1; double xmod2[300];

const char *namex; const char *titlex;
char namey[100], titley[100];

  double yvl[200]={0.0};
  double erryvl[200]={0.0};
  int nbinx = thin->GetNbinsX();
  
  int ifirst = 0;
  if (nxmod2<0) {
    for (int ij=0; ij<nbinx+2; ij++) {
      if (ij <=arrayvar_firstG[itype][ivar][ijetpt]) {
      } else if ( ij > (nbinx - arrayvar_lastG[itype][ivar][ijetpt])) {
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
for(int ii=0; ii < nxmod2+1 ; ii++){ xmod3[ii]=xmod2[ii]; }

  namex = thin->GetName();
  sprintf(namey, "rebin_%s", namex);   // Roo means rebinded
  titlex = thin->GetTitle();
  sprintf(titley, "Rebin_%s", titlex);
  thout = new TH1D(namey, titley, nxmod2, xmod3);

  for(int ix=0; ix<nxmod2+1; ix++) {
    thout->SetBinContent(ix+1, yvl[ix]);
    thout->SetBinError(ix+1, sqrt(erryvl[ix]));     }
  //  thout->SetBinContent(0,5);
    return thout;
}   //end of rebin_hist functiomC


void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect){

//calculate Probability Matrix
for(int ij=1; ij<(HistoMatrix->GetNbinsY()+1); ij++){
double row_sum = 0.;
for(int jk=1; jk<(HistoMatrix->GetNbinsX()+1); jk++){
  row_sum+=HistoMatrix->GetBinContent(jk,ij);
}//jk
if(row_sum>1.e-10){
 for(int jk=1; jk<(HistoMatrix->GetNbinsX()+1); jk++){
   HistoMatrix->SetBinContent(jk,ij,HistoMatrix->GetBinContent(jk,ij)*1./row_sum) ; //Probability
  }//jk
}
}//ij

//folding gen level to Reco
for(int i=1;i<HistReco->GetNbinsX()+1;i++){
     double sum=0.; double Err =0.;
       for(int j=1;j<HistGen->GetNbinsX()+1;j++){
       double misscorr = (HistGen->GetBinContent(j))-(miss->GetBinContent(j)); //Miss correction
       sum += HistoMatrix->GetBinContent(i,j)*misscorr;
       Err += (HistoMatrix->GetBinContent(i,j)*HistGen->GetBinError(j))*(HistoMatrix->GetBinContent(i,j)*(HistGen->GetBinError(j)));
           }
      sum = sum +(fake->GetBinContent(i)); //fake correction
      HistoCorrect->SetBinContent(i,sum);
      HistoCorrect->SetBinError(i,sqrt(Err));
    	}
}//end Fold

//Condition number calculation
void Condition (TH2* RM, TH1* miss){
    const int Nx = RM->GetNbinsX(),
              Ny = RM->GetNbinsY();
    cout << Nx <<"  " << Ny << endl;
    if (Ny*2 != Nx) { cout << Nx << ' ' << Ny << endl;  return; }

    TH1D* RMy = RM->ProjectionY("RMx", 0, -1); //Gen Projection

    // normalisation & condition
    TMatrixD m(Ny,Ny);
    for (int i = 1; i <= Ny; ++i) {
        double normalisation = RMy->GetBinContent(i);
               normalisation += miss->GetBinContent(i);
        if (normalisation > 0)
        for (int j = 1; j <= Nx; ++j) {
            double content = RM->GetBinContent(j,i);
            content /= normalisation;
            m((j-1)/2,i-1) += content;
        }
    }
    TDecompSVD svd(m);
    TVectorD v = svd.GetSig();

    double Max = v[0];
    for(int k = 0; k < Ny; ++k) {
        if (abs(v[k]) < feps) break;
        cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/v[k] << '\n';
    }
}



