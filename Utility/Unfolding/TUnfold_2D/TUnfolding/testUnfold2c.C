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
//#define UOEx
#define CLOUSER
#define Unfold1D
#define Unfold2D

using namespace std;
static const auto feps = numeric_limits<float>::epsilon();

void testUnfold2c(){
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();
  
  int const nmc=3; //Number of MC -> 0,1,2 PY8, MG, HW7
  int const umc=0; //Which MC will used for Unfold
  int irbin = 1; //Rebin
  int const untype = 1;  // 0 for 2D , 1 for 2D 
  
  TFile *inputbinning = new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  //Input Data and MC histogram
  TFile *inputData=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  TFile *inputMC[nmc];
  inputMC[0]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  inputMC[1]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  inputMC[2]=new TFile("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/TUnfold_Ready.root");
  
  TFile *outputFile=new TFile("Unfolded_Result.root","recreate");   //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
 
  ofstream file;
  file.open("Tau_Value.txt");
  file <<"L Curve    "  <<"                   "<< "Scan Sure  "  <<"                "<<" Scan Tau "<<endl;
  
  int const type = 2;         // Jet & Charage particles
  int const itype[type]={0,1};   //{0}--->Jet ; {1}---> Charged Particles
  const  char* itypeN[type]={"Jets","Charged Particles"};
  const char* Dirbin[3] = {"","Binning1D","Binning2D"};
  const char* Dirhist[3] = {"analyzeBasicPat", "analyzeBasicPat1D", "analyzeBasicPat2D"};
  
  const char* bintag[3] = {"", "1d", "2d"};
  const char* histtag[3] = {"", "d_", "dd_"};

  char histname[100], name[100], title[100], Axisname[100];
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
  
  
  TDirectoryFile *DirData[3];
  DirData[0] = new TDirectoryFile("Data","Inputs Data");
  DirData[1] = new TDirectoryFile("Data1D","Inputs Data");
  DirData[2] = new TDirectoryFile("Data2D","Inputs Data");
  
  TDirectoryFile *outputDir[nmc];
  outputDir[0]=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
  outputDir[1]=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
  outputDir[2]=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
  
  
  TDirectoryFile *folddir[3];
  folddir[0]=new TDirectoryFile("Folded","folded with Probablility matrix root");
  folddir[1]=new TDirectoryFile("Folded1D","folded with Probablility matrix 1D");
  folddir[2]=new TDirectoryFile("Folded2D","folded with Probablility matrix 2D");
  TDirectoryFile *Unfolddir[3];
  Unfolddir[0]=new TDirectoryFile("Unfold","Unfolded, Refold, correlation root");
  Unfolddir[1]=new TDirectoryFile("Unfold1D","Unfolded, Refold, correlation 1D");
  Unfolddir[2]=new TDirectoryFile("Unfold2D","Unfolded, Refold, correlation 1D");
  
  void setgstyle();
  int subtract_background1(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
  int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
  void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
  void Condition (TH2 * RM, TH1* miss);
  TH1D* rebin1d_hist(TH1D* thin, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
  TH2D* rebin2d_hist(TH2D* thin, TH1D* MC_reco,TH1D* MC_gen, int itype, int ijetpt, int ivar,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8],int arrayvar_firstG[2][5][8],int arrayvar_lastG[2][5][8] );
  //TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar);
  TH1D* rebin1d_hist_gen(TH1D* thin, int itype, int ijetpt, int ivar ,int arrayvar_first[2][5][8],int arrayvar_last[2][5][8]);
  TH1*  GetLocalBinnedHist(TH1* Hist,  TUnfoldBinning* bin, char const* none, char const* axis);

//----------------------------------Different Binning ----------------------------------------------
for(int idd=0; idd <3; idd++){  //0 : Root Hist 1: 1D TunfoldBinning 2: 2D TUnfoldBinning
  TH1D *MC_Reco[nmc][type][nusedvar][njetptmn];  //Reconstructed MC
  TH1D *MC_fake[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_fakerate[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_background[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_reco_fake[nmc][type][nusedvar][njetptmn];  //Fake :  Reco but No Gen

  TH1D *MC_Gen[nmc][type][nusedvar][njetptmn];   //Generator MC
  TH1D *MC_miss[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *MC_Gen_miss[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *MC_missrate[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *MC_misscorr[nmc][type][nusedvar][njetptmn];   //Miss:  No Reco but in Gen

  TH1D *PsudoData_Gen[nmc][type][nusedvar][njetptmn];    //Gen Level Psudo Data Only for Closure
  TH2D *h2dGenDetMC[nmc][type][nusedvar][njetptmn];   // MC Response Matrix
  TH1D *RMX[nmc][type][nusedvar][njetptmn];  //Response Matrix ProjectionX
  TH1D *RMY[nmc][type][nusedvar][njetptmn];  //RM projetcion

  TH1D *Data_Reco[type][nusedvar][njetptmn];    //Reconstructed Data

  TH1D* hist_eff[nmc][type][nusedvar][njetptmn];
  TH1D* hist_fake[nmc][type][nusedvar][njetptmn];
  TH1D* hist_purity[nmc][type][nusedvar][njetptmn];
  TH1D* hist_stbl[nmc][type][nusedvar][njetptmn];

  TUnfoldBinning *binsRec[type][nusedvar][njetptmn]; //Binning Name 
  TUnfoldBinning *RecoBinning[type][nusedvar][njetptmn]; //Node Name

  TUnfoldBinning *binsGen[type][nusedvar][njetptmn];  //Gen Binning
  TUnfoldBinning *GenBinning[type][nusedvar][njetptmn];  //Gen Node 

//---------------------------------------------Binning-------------------------------------
if(idd>0){
for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
        sprintf(histname, "%s/Detector%s_typ_%i_pt%i_eta0_%i", Dirbin[idd], bintag[idd], ity, ipt, var[ivar]);
        cout << histname <<endl;
        inputbinning->GetObject(histname, binsRec[ity][ivar][ipt]);
        binsRec[ity][ivar][ipt]->PrintStream(cout);
        
	sprintf(histname, "%s/Generator%s_typ_%i_pt%i_eta0_%i", Dirbin[idd], bintag[idd], ity, ipt, var[ivar]);
        inputbinning->GetObject(histname, binsGen[ity][ivar][ipt]);
        binsGen[ity][ivar][ipt]->PrintStream(cout);
      }
    }
  }
}
//------------------------------------------
//----------------------------------------------Read Input Data MC and Response matrix
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	
	sprintf(histname, "%s/%sreco_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
	TH1D *RecoData =(TH1D*) inputData->Get(histname);
	//for(int i=1 ; i< RecoData->GetNbinsX()+1; i++){ if(RecoData->GetBinContent(i) == 0) {cout << " Data Bin entry Nil : bin no : "<< i << endl;}}
	sprintf(histname, "%s/%sgen_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
#ifdef CLOUSER
	TH1D *PsudoDataGen= (TH1D*)inputData->Get(histname);
#else
	TH1D *PsudoDataGen = (TH1D*)inputMC[umc]->Get(histname);
#endif

for (int imc=0; imc<nmc ; imc++){
//----------------------------------MC RECO
       sprintf(histname, "%s/%sreco_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
       TH1D *RecoMC = (TH1D*)inputMC[imc]->Get(histname);      cout << histname ;
      
       //if(recobins!=Data_Reco[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
       //-----------------------------------MC Fake
       sprintf(histname, "%s/%sfake_reco_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
       TH1D *RecoMC_Fake = (TH1D*)inputMC[imc]->Get(histname);
       
       cout << " Fake= " <<RecoMC_Fake->GetEntries() <<" Reco-fake: " <<(RecoMC->GetEntries() - RecoMC_Fake->GetEntries())<<endl;
       
       //-----------------------------------Gen MC
       sprintf(histname, "%s/%sgen_typ_%i_pt%i_eta0_%i",  Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
       TH1D *GenMC = (TH1D*)inputMC[imc]->Get(histname); cout << histname  ;
       
       //for(int i= 1; i < GenMC->GetNbinsX()+1; i++){if(GenMC->GetBinContent(i) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  i  << endl; }}
       
       //----------------------------------MC miss 
       sprintf(histname, "%s/%smiss_gen_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
       TH1D *GenMC_Miss = (TH1D*)inputMC[imc]->Get(histname);

       cout << " Miss= " << GenMC_Miss->GetEntries() <<" Gen-Miss: " << (GenMC->GetEntries() - GenMC_Miss->GetEntries()) ;
       
       //Response Matrix
       sprintf(histname, "%s/%scorr_typ_%i_pt%i_eta0_%i", Dirhist[idd],histtag[idd], ity, ipt, var[ivar]); 
       TH2D *RM_RecoGen =(TH2D*)inputMC[imc]->Get(histname);  cout << " Corr = "  << RM_RecoGen->GetEntries() <<endl;

#ifdef UOEx
       TH1D *NewReco; TH1D *NewGen; TH1D *Newmiss; TH1D *Newfake; TH2D *NewRM;
       NewReco = (TH1D*)RecoMC->Clone(); Newfake = (TH1D*)RecoMC_Fake->Clone();
       NewGen  = (TH1D*)GenMC->Clone(); Newmiss = (TH1D*)GenMC_Miss->Clone();
       NewRM = (TH2D*)RM_RecoGen->Clone();

       int NbinxR = NewReco->GetNbinsX(); int NbinxG = NewGen->GetNbinsX(); int NbinMx = NewRM->GetNbinsX(); int NbinMy = NewRecoGen->GetNbinsY();
       NewReco->Reset(); NewGen->Reset(); NewRM->Reset(); Newmiss->Reset(); Newfake->Reset();

       for(int ix=0; ix < NbinxR+2 ; ix++){
         NewReco->SetBinContent(ix,NewReco->GetBinContent(ix));
         NewReco->SetBinError(ix, sqrt(NewReco->GetBinError(ix)* NewReco->GetBinError(ix))); }
       for(int ix=0; ix < NbinxG+2 ; ix++){
         NewGen->SetBinContent(ix, NewGen->GetBinContent(ix));
         NewGen->SetBinError(ix, sqrt(NewGen->GetBinError(ix)* NewGen->GetBinError(ix))); }
       for(int ix=0; ix < NbinMx+2 ; ix++){
         for(int iy=0; iy < NbinMy+2 ; iy++){
           NewRM->SetBinContent(ix, iy, NewRecoGen->GetBinContent(ix,iy));
           NewRM->SetBinError(ix, iy, sqrt(NewRecoGen->GetBinError(ix,iy)*NewRecoGen->GetBinError(ix,iy)));
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
       
       RM_RecoGen = (TH2D*)NewRM->Clone();

#endif 
//--------------------------------------Calculate Fake rate and Miss Rate
       TH1D* fakerate = (TH1D*)RecoMC_Fake->Clone();
       sprintf(name,"%sfake_rate_%i_pt%i_eta0_%i", histtag[idd], ity, ipt, var[ivar]); fakerate->SetNameTitle(name,name);
       fakerate->Divide(RecoMC_Fake, RecoMC, 1, 1, "b");
       TH1D* missrate = (TH1D*)GenMC_Miss->Clone();
       sprintf(name,"%smiss_rate_%i_pt%i_eta0_%i", histtag[idd],ity, ipt, var[ivar]); missrate->SetNameTitle(name,name);
       missrate->Divide(missrate, GenMC_Miss, 1, 1, "b");


//--------------------------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20
       TH1* RMx = RM_RecoGen->ProjectionX(); TH1* RMy = RM_RecoGen->ProjectionY();
        
       sprintf(name,"%sProjectX_%i_pt%i_eta0_%i", histtag[idd], ity, ipt, var[ivar]); RMx->SetNameTitle(name,name);

       sprintf(name,"%sRecominusfake_%i_pt%i_eta0_%i",histtag[idd] ,ity, ipt, var[ivar]);
       TH1* RecoFakeCorrect = (TH1D*)RecoMC->Clone(); RecoFakeCorrect->Reset();
       RecoFakeCorrect->SetNameTitle(name,name);

       sprintf(name,"%sProjectY_%i_pt%i_eta0_%i", histtag[idd],ity, ipt, var[ivar]); RMy->SetNameTitle(name,name);

       sprintf(name,"%sGenminusmiss_%i_pt%i_eta0_%i",histtag[idd] ,ity, ipt, var[ivar]);
       TH1* GenMissCorrect = (TH1D*)GenMC->Clone(); GenMissCorrect->Reset();
       GenMissCorrect->SetNameTitle(name,name);
       
       for (int i = 1; i <= RecoFakeCorrect->GetNbinsX(); ++i) {
         double content = RecoMC->GetBinContent(i); double factor = RecoMC_Fake->GetBinContent(i);
         content -= factor;  RecoFakeCorrect->SetBinContent(i, content);
       }
       MC_reco_fake[imc][ity][ivar][ipt]= (TH1D*)RecoFakeCorrect->Clone();

       for (int i = 1; i <= GenMissCorrect->GetNbinsX(); ++i) {
        double content = GenMC->GetBinContent(i); double factor = GenMC_Miss->GetBinContent(i);
         content -= factor;  GenMissCorrect->SetBinContent(i, content);
       }
      
       MC_Gen_miss[imc][ity][ivar][ipt]= (TH1D*)GenMissCorrect->Clone();
       TH1* RecoFakecorrEx; TH1* GenMisscorrEx; 
     outputDir[imc]->cd(); //..........................MC directory .......................................
       if(idd>=1){
       sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
       sprintf(name,"E%s", MC_reco_fake[imc][ity][ivar][ipt]->GetName());  sprintf(Axisname,"var_%i[UO]",var[ivar]);
       RecoFakecorrEx =binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name,  MC_reco_fake[imc][ity][ivar][ipt], 0, true, Axisname);
     
       sprintf(histname,"Genbin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
       sprintf(name,"E%s", MC_Gen_miss[imc][ity][ivar][ipt]->GetName());   sprintf(Axisname,"var_%i[UO]",var[ivar]);
       GenMisscorrEx =binsGen[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, MC_Gen_miss[imc][ity][ivar][ipt], 0, true, Axisname);
       
       GenMisscorrEx->Write(); RecoFakecorrEx->Write();   
       }
//-------------------------------------------Extact Hist----------------------------- 
     if(idd>=1){
       cout<<" ok  "<<endl;
       sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]); 
       sprintf(name,"E%sreco_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]); 
       sprintf(Axisname,"var_%i[UO]",var[ivar]);
       TH1* RecoExtact =binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, RecoMC, 0, true, Axisname);

       sprintf(histname,"Genbin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
       sprintf(name,"E%sgen_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
       sprintf(Axisname,"var_%i[UO]",var[ivar]);
       TH1* GenExtact =binsGen[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, GenMC , 0, true, Axisname);
      
       RecoExtact->Write(); GenExtact->Write();
     
       //RecoMC =(TH1*) RecoExtact->Clone();
       //GenMC =(TH1*) GenExtact->Clone();
      
       RecoMC->Write(); GenMC->Write();
       }else{
        
       RecoMC->Write(); GenMC->Write();
       
       }

       RMx->Write(); RMy->Write(); RecoFakeCorrect->Write(); GenMissCorrect->Write();
       RecoMC_Fake->Write(); fakerate->Write();
       GenMC_Miss->Write(); missrate->Write(); 
       RM_RecoGen->Write();
//----------------------------------------------------------------------------
       MC_Reco[imc][ity][ivar][ipt] = (TH1D*)RecoMC->Clone();
       MC_fake[imc][ity][ivar][ipt] = (TH1D*)RecoMC_Fake->Clone();
       MC_fakerate[imc][ity][ivar][ipt] = (TH1D*)fakerate->Clone();
       MC_background[imc][ity][ivar][ipt] = (TH1D*)fakerate->Clone();
       MC_fakerate[imc][ity][ivar][ipt]->SetMinimum(-0.05); MC_fakerate[imc][ity][ivar][ipt]->SetMaximum(1.01);
       

       MC_Gen[imc][ity][ivar][ipt] = (TH1D*)GenMC->Clone();
       MC_miss[imc][ity][ivar][ipt] = (TH1D*)GenMC_Miss->Clone();
       MC_missrate[imc][ity][ivar][ipt] = (TH1D*)missrate->Clone();
       MC_misscorr[imc][ity][ivar][ipt] = (TH1D*)missrate->Clone();
       MC_missrate[imc][ity][ivar][ipt]->SetMinimum(-0.05); MC_missrate[imc][ity][ivar][ipt]->SetMaximum(1.01);

       h2dGenDetMC[imc][ity][ivar][ipt] = (TH2D*)RM_RecoGen->Clone();


//------------------------------------------------------Stability and Purity -----------------------------------------
/*
       hist_purity[imc][ity][ivar][ipt] = (TH1D*)MC_fake[imc][ity][ivar][ipt]->Clone(); hist_purity[imc][ity][ivar][ipt]->Reset();
       for(int binRec=0; binRec<= hist_purity[imc][ity][ivar][ipt]->GetNbinsX()+1; binRec++) {
          double sum=0.;
          for(int binGen=0; binGen<=hist_purity[imc][ity][ivar][ipt]->GetNbinsX()+1; binGen++) {
            //sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
            sum += h2dGenDetMC[imc][ity][ivar][ipt]->GetBinContent(binRec,binGen);
          }
          double p=0.;
          if(sum>0.0) {
            p = h2dGenDetMC[imc][ity][ivar][ipt]->GetBinContent(binRec,binRec)/sum;
          }
          hist_purity[imc][ity][ivar][ipt]->SetBinContent(binRec,p);
        }
        hist_purity[imc][ity][ivar][ipt]->SetMinimum(-0.05); hist_purity[imc][ity][ivar][ipt]->SetMaximum(1.01);
        hist_purity[imc][ity][ivar][ipt]->Write();
*/
        TH1D *h_pu = (TH1D*)fakerate->Clone(); h_pu->Reset();
	TH1D *h_st = (TH1D*)fakerate->Clone(); h_st->Reset();
	int ir = h_pu->GetNbinsX(); int ig = missrate->GetNbinsX();
        double fk[ir]; double ef[ig]; double pu[ir]; double st[ir];
        subtract_background(RM_RecoGen, RecoMC, GenMC, RecoData, fk, ef, pu, st);
        
	for(int i =1; i<h_pu->GetNbinsX()+1; i++){
	h_pu->SetBinContent(i,pu[i]);
	h_st->SetBinContent(i,st[i]);
	}

	hist_purity[imc][ity][ivar][ipt]=(TH1D*)h_pu->Clone(); 
	hist_stbl[imc][ity][ivar][ipt]=(TH1D*)h_st->Clone(); 
        hist_purity[imc][ity][ivar][ipt]->SetMinimum(-0.05); hist_purity[imc][ity][ivar][ipt]->SetMaximum(1.01);
        hist_stbl[imc][ity][ivar][ipt]->SetMinimum(-0.05); hist_stbl[imc][ity][ivar][ipt]->SetMaximum(1.01);

	sprintf(name,"%sPurity_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
        hist_purity[imc][ity][ivar][ipt]->SetNameTitle(name,name);
        sprintf(name,"%sstability_%i_pt%i_eta0_%i", histtag[idd], ity, ipt, var[ivar]);
        hist_stbl[imc][ity][ivar][ipt]->SetNameTitle(name,name);

	hist_purity[imc][ity][ivar][ipt]->Write();
	hist_stbl[imc][ity][ivar][ipt]->Write();
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
        DirData[idd]->cd(); 
         if(idd>=1){ 
         sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
	 sprintf(name, "E%sreco_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
         sprintf(Axisname,"var_%i[UO]",var[ivar]);
         TH1* DataRecoExtact = binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, RecoData, 0, true, Axisname);
	 DataRecoExtact->Write();
	}//if(idd>=1)
         RecoData->Write(); 
         Data_Reco[ity][ivar][ipt] = (TH1D*)RecoData->Clone();
#ifdef CLOUSER
	 if(idd>=1){
         sprintf(histname,"Genbin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
         sprintf(name,"E%sgen_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
         sprintf(Axisname,"var_%i[UO]",var[ivar]);
         TH1* GenExtactPsudo =binsGen[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, PsudoDataGen , 0, true, Axisname); 
 	 GenExtactPsudo->Write();
	 }//if(idd>=1)
  	 PsudoDataGen->Write();

       TH1* RecoFakecorrExD; TH1* GenMisscorrExD;
       if(idd>=1){
       sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
       sprintf(name,"E%s", MC_reco_fake[umc][ity][ivar][ipt]->GetName());  sprintf(Axisname,"var_%i[UO]",var[ivar]);
       RecoFakecorrExD =binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name,  MC_reco_fake[umc][ity][ivar][ipt], 0, true, Axisname);

       sprintf(histname,"Genbin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
       sprintf(name,"E%s", MC_Gen_miss[umc][ity][ivar][ipt]->GetName());   sprintf(Axisname,"var_%i[UO]",var[ivar]);
       GenMisscorrExD =binsGen[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, MC_Gen_miss[umc][ity][ivar][ipt], 0, true, Axisname);

       GenMisscorrExD->Write(); RecoFakecorrExD->Write();
       }

        
#endif
      }//for(int ipt = 0 ; ipt < njetptmn ; ipt++)
    }//for(int ivar=0; ivar < nusedvar ; ivar ++)
  }//for(int ity=0; ity <type; ity++)
  cout << "Histogram Read Data and MC Done " <<endl;
  
  //------------------Fold check : Patrick 1 Sep 20
  //Get Probability Matrix  (gen-miss)*probability = (reco-fake). Of course, don't forget to account for miss and fake entries (if applicable).
  folddir[idd]->cd();
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
	
	if(idd>=1){
        sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
        sprintf(name, "E%sFold_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);   sprintf(Axisname,"var_%i[UO]",var[ivar]);
	TH1* TrueFold = binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, Folded, 0, true, Axisname);   TrueFold->Write();
	}
	sprintf(name,"%sFold_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
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
  
  
  
  Unfolddir[idd]->cd();
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
	
	if (idd==0){cout << "Root Hist : "; }
	if (idd==1){cout << "TUnfoldBinning1D : "; }
	if (idd==2){cout << " TUnfoldBinning2D : "; }
	cout <<"typ "<< ity << " : Var " << ivar << " HT2Bin : " << ipt <<"  ";
//	file <<"["<< ity << "," << ivar << "," << ipt <<"] --->" << endl;
	
        
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
//----------------------------------------------------------Define Input---------------------	
	TH2* RMin = (TH2D*)h2dGenDetMC[umc][ity][ivar][ipt]->Clone();
	TH1* input = (TH1D*)Data_Reco[ity][ivar][ipt]->Clone();
	TH1* mcgen = (TH1D*)MC_Gen[umc][ity][ivar][ipt]->Clone();
        TH1* mcgenmiss=(TH1D*)MC_misscorr[umc][ity][ivar][ipt]->Clone();	
	
	//correction for Fake as background subtraction : Patrick
	TH1* mcbackground = (TH1D*)MC_background[umc][ity][ivar][ipt]->Clone();
	mcbackground->Multiply(input);

        double biasScale = 0;
        const char *REGULARISATION_DISTRIBUTION=0;
        const char *REGULARISATION_AXISSTEERING="*[UOB]";
        
	//https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
       sprintf(name,"%sData_covariance_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
       TH2D* covM = new  TH2D(name,name, input->GetNbinsX(), rxbins, input->GetNbinsX(), rxbins);
       covM->Sumw2();
       for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
          double err = input->GetBinError(ix);
          covM->SetBinContent(ix,ix,err*err);
       }

       covM->Write();


       TUnfoldBinning* RecoBin = 0;
       TUnfoldBinning* GenBin = 0;

       if(idd>0){ RecoBin = binsRec[ity][ivar][ipt]; 
                  GenBin =binsGen[ity][ivar][ipt]; };
//----------------------------------------------------------No Reguratization-------------------------------------
	
        TUnfoldDensity density(RMin,TUnfold::kHistMapOutputVert,
                                               TUnfoldDensity::kRegModeNone,
                                               TUnfoldDensity::kEConstraintNone,
                                               TUnfoldDensity::kDensityModeNone,
                                               GenBin,RecoBin);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);

        density.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 
        int status = density.SetInput(input,biasScale,0,covM);

        int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
        cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;

        //but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
        //tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result

        density.DoUnfold(0.0);//,input,biasScale);  // tau 0.0 means no regularization

        sprintf(histname, "%sTUnfold_NoReg_typ_%i_pt%i_eta0_%i", histtag[idd], ity, ipt, var[ivar]); //unfolded_typ_0_pt2_eta0_3
        sprintf(title, "Unfolded No Reg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
        bool AxisBin =true;
	if(idd>=1){AxisBin = false;}//Keep or Not ?
	TH1 *Unfolded = density.GetOutput(histname,title,0,"*[UO]", AxisBin);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");

        // correction for Miss entries  : Partick
        for (int i = 1; i <= Unfolded->GetNbinsX(); ++i) {
          double content = Unfolded->GetBinContent(i);
          double factor = 1;
          factor += mcgenmiss->GetBinContent(i);
          content *= factor;
          Unfolded->SetBinContent(i, content);
        }

        sprintf(histname, "%scorr_NoReg_typ_%i_pt%i_eta0_%i", histtag[idd], ity, ipt, var[ivar]);
        sprintf(title, "RhoIJtotal No Reg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
        TH2 *Hist_RhoIJ = density.GetRhoIJtotal(histname, title);//,"signal");

        sprintf(histname, "%sEmat_NoReg_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
        sprintf(title, "Ematrix No Reg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
        TH2 *hist_Emat = density.GetEmatrixTotal(histname,title);//,"signal");

        sprintf(histname, "%sRefold_NoReg_typ_%i_pt%i_eta0_%i", histtag[idd],ity, ipt, var[ivar]);
        sprintf(title, "Back folded No Reg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
        TH1 *hist_folded = density.GetFoldedOutput(histname,title,0,"*[UO]",AxisBin);//,"signal");

        sprintf(histname, "%sProb_NoReg_typ_%i_pt%i_eta0_%i",histtag[idd], ity, ipt, var[ivar]);
        sprintf(title, "ProbabilityMatrix No Reg %s %i 2.4 %s ", itypeN[ity], int(leadingPtThreshold[ipt]), vartitle[var[ivar]]);
        TH2 *hist_prob = density.GetProbabilityMatrix(histname,title,TUnfold::kHistMapOutputVert);//,"signal");

        Unfolded->Write(); Hist_RhoIJ->Write(); hist_Emat->Write(); hist_folded->Write(); hist_prob->Write();

if(idd>=1){
         sprintf(histname,"Recobin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
         sprintf(name, "E%s",hist_folded->GetName()); sprintf(Axisname,"var_%i[UO]",var[ivar]);
         TH1* foldedExtact =binsRec[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, hist_folded, 0, true, Axisname);
         
	 sprintf(histname,"Genbin%s_typ_%i_pt%i_eta0_%i", bintag[idd], ity, ipt, var[ivar]);
         sprintf(name,"E%s",Unfolded->GetName()); sprintf(Axisname,"var_%i[UO]",var[ivar]);
         TH1* UnfoldedExtact =binsGen[ity][ivar][ipt]->FindNode(histname)->ExtractHistogram(name, Unfolded , 0, true, Axisname);
          
	 foldedExtact->Write(); UnfoldedExtact->Write();


          }
       //-----------------------------------------End of Variables loop
     }
   }
 }

if(idd==0){ cout << "Root OK " <<endl;};
if(idd==1){ cout << "TUnfoldBinning 1D OK " <<endl;};
if(idd==2){ cout << "TUnfoldBinning 2D OK " <<endl;};


} //idd
 
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
  const int nbinmx = 6000 ;
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


void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect){

TH2D* Histprob = (TH2D*) HistoMatrix->Clone(); Histprob->Reset();
//calculate Probability Matrix
for(int ij=0; ij<(HistoMatrix->GetNbinsY()+2); ij++){
double row_sum = 0.;
for(int jk=0; jk<(HistoMatrix->GetNbinsX()+2); jk++){
  row_sum+=HistoMatrix->GetBinContent(jk,ij);
}//jk
if(row_sum>1.e-10){
 for(int jk=0; jk<(HistoMatrix->GetNbinsX()+2); jk++){
   Histprob->SetBinContent(jk,ij,(HistoMatrix->GetBinContent(jk,ij)*1./row_sum)) ; //Probability
  }//jk
}
}//ij

//folding gen level to Reco
for(int i=0;i<Histprob->GetNbinsX()+2;i++){
     double sum=0.; double Err =0.;
       for(int j=0;j<HistGen->GetNbinsX()+2;j++){
       double misscorr = (HistGen->GetBinContent(j))-(miss->GetBinContent(j)); //Miss correction
       sum += Histprob->GetBinContent(i,j)*misscorr;
       Err += (Histprob->GetBinContent(i,j)*HistGen->GetBinError(j))*(Histprob->GetBinContent(i,j)*(HistGen->GetBinError(j)));
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

//-----------------------------Set gstyle (Copy from TUnfold code)
void setgstyle(){     
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
}

TH1* GetLocalBinnedHist(TH1* Hist,  TUnfoldBinning* bin, const char* node, const char* axis){
    const char* name = Hist->GetName();
    TH1* localhist = bin->FindNode(node)->ExtractHistogram(name, Hist, 0, true, axis);
    localhist->SetNameTitle(name,name);
return localhist;
}
