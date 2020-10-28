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

#define CLOUSER
#define BLTest

using namespace std;
static const auto feps = numeric_limits<float>::epsilon();

void testUnfold2c(){
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();
  
  int const nmc=4; //Number of MC -> 0,1,2 PY8, MG, HW7
  int const umc=0; //Which MC will used for Unfold
  int irbin = 1; //Rebin
  int const untype = 1;  // 0 for 2D , 1 for 2D 
  
  //const TString Pyinput = "PY8_UL17_2D_ALLHT2.root";
  const TString Pyinput = "PY8_UL17_Bin17_JecV5_JERV2_2D_17Oct20.root";
  //const TString Pyinput = "MG_UL17_JECV5_JEERV2_2D_17oct20.root";
  //const TString Pyinput = "PY8_UL17_Flat_JECV5_JERV2_17oct20.root";
  //const TString Pyinput = "Herwig_UL17_Flat_17Oct20.root";
  
  //const TString datainput = "PY8_UL17_Flat_JECV5_JERV2_17oct20.root"; 
  //const TString datainput = "PY8_UL17_Bin17_JecV5_JERV2_2D_17Oct20.root"; 
  //const TString datainput = "Data_UL17v1_JRV2_JECV4_2D_14Oct20.root"; 
  const TString datainput = "MG_UL17_JECV5_JEERV2_2D_17oct20.root"; 
  //const TString datainput = "Herwig_UL17_Flat_17Oct20.root"; 
  
  TFile *inputbinning = new TFile("PY8_HT2_83_3000_Binnied.root");
  
  //Input Data and MC histogram
  TFile *inputData=new TFile(datainput);
  TFile *RMinput=new TFile(Pyinput);
  
  TFile *inputMC[nmc];
  inputMC[0]=new TFile("PY8_UL17_Bin17_JecV5_JERV2_2D_17Oct20.root");
  inputMC[1]=new TFile("PY8_UL17_Flat_JECV5_JERV2_17oct20.root");
  inputMC[2]=new TFile("MG_UL17_JECV5_JEERV2_2D_17oct20.root");
  inputMC[3]=new TFile("Herwig_UL17_Flat_17Oct20.root");
  
  TFile *outputFile=new TFile("Unfolded_Result.root","recreate");   //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  
  int const type = 2;         // Jet & Charage particles
  int const itype[type]={0,1};   //{0}--->Jet ; {1}---> Charged Particles
  const  char* itypeN[type]={"Jets","Charged Particles"};
  const char* BinDir = "Binning2D";
  string HistDir = "analyzeBasicPat2D";
  
  const char* Dimtag = "2d";
  const char* Histtag = "dd_";
  
  char histname[100], name[100], title[100], Axisname[100];
  const int nHLTmx=8; //HT2 Range
  const int njetetamn=1;  //eta value used 2.4
  //double etarange[njetetamn] ={2.4}; //
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass
  static const int nhist=10; //We need 8 But define 10
  const int njetptmn = nHLTmx;
  Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000};
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
  TDirectoryFile *DirData = new TDirectoryFile("Data2D","Inputs Data");
  
  TDirectoryFile *outputDir[nmc];
  outputDir[0]=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
  outputDir[1]=new TDirectoryFile("Py8Flat"," Pythia8 Flat sample , MC and Probability Matrix");
  outputDir[2]=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
  outputDir[3]=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
  
  TDirectoryFile *DirRMinput = new TDirectoryFile("DirRMinput","MC used in Unfolding");
  TDirectoryFile *folddir = new TDirectoryFile("Folded2D","Gen fold with Probablility matrix 2D");
  TDirectoryFile *Unfolddir = new TDirectoryFile("Unfold2D","Unfolded, Refold, correlation 1D");
  
  void setgstyle();
  void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
  void Condition (TH2 * RM, TH1* miss);
  void BLT (TH1 * dataDist, TH2 * dataCov, TH1 * MC, int rebin = 1);
  void ConditionV2(TH2* prob_mat);
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void Integralhist(TH1 *hist);
  double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/);
  void Normalise(TH1* h, TH2* covmat);
  void Extract(TH1* global2d, TUnfoldBinning* Bin , char* Axisname, bool iw=0);
  TH1D* ReadHist1D(string name,TFile* root);
  TH2D* ReadHist2D(string name,TFile* root);
  void HT2_NormalV3(TH1* global2d, TUnfoldBinning* Bin, char* Axisname, int nht,  bool iw=0);
  
  //---------------------------------Different Binning ----------------------------------------------
  TH1D *MC_Reco[nmc][type][nusedvar];  //Reconstructed MC
  TH1D *MC_fake[nmc][type][nusedvar];  //Fake :  Reco but No Gen
  TH1D *MC_fakerate[nmc][type][nusedvar];  //Fake rate
  TH1D *MC_reco_fake[nmc][type][nusedvar];  //reco -fake
  
  TH1D *MC_Gen[nmc][type][nusedvar];   //Generator MC
  TH1D *MC_miss[nmc][type][nusedvar];   //Miss:  No Reco but in Gen
  TH1D *MC_Gen_miss[nmc][type][nusedvar];   //Gen -miss
  TH1D *MC_missrate[nmc][type][nusedvar];   //Missrate
  
  TH2D *h2dGenDetMC[nmc][type][nusedvar];   // MC Response Matrix
  TH1D *RMX[nmc][type][nusedvar];  //Response Matrix ProjectionX
  TH1D *RMY[nmc][type][nusedvar];  //RM projetcion
  
  TH1D* hist_eff[nmc][type][nusedvar];
  TH1D* hist_fake[nmc][type][nusedvar];
  TH1D* hist_purity[nmc][type][nusedvar];
  TH1D* hist_stbl[nmc][type][nusedvar];
  
  TH1D *Data_Reco[type][nusedvar];    //Reconstructed Data
#ifdef CLOUSER
  TH1D *Data_Gen[type][nusedvar];    //     For Closure 
  TH1D *Data_fakerate[type][nusedvar];    //Reconstructed Data
  TH1D *Data_fakerateInv[type][nusedvar];    //Reconstructed Data
  TH1D *Data_missrateInv[type][nusedvar];    //Reconstructed Data
  TH1D *Data_missrate[type][nusedvar];    //Reconstructed Data
  TH1D *Data_reco_fake[type][nusedvar];    //Reconstructed Data
  TH1D *Data_gen_miss[type][nusedvar];    //Reconstructed Data
#endif 
  TH1D *RMinput_Reco[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_Gen[type][nusedvar];    //     For Closure
  TH1D *RMinput_fakerate[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_fakerateInv[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_missrateInv[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_missrate[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_reco_fake[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_fake[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_gen_miss[type][nusedvar];    //Reconstructed RMinput
  TH1D *RMinput_miss[type][nusedvar];    //Reconstructed RMinput
  TH2D *RMinput_RM[type][nusedvar];    //Reconstructed RMinput
  
  TUnfoldBinning *binsRec[type][nusedvar]; //Binning Name 
  TUnfoldBinning *RecoBinning[type][nusedvar]; //Node Name
  TUnfoldBinning *binsGen[type][nusedvar];  //Gen Binning
  TUnfoldBinning *GenBinning[type][nusedvar];  //Gen Node 
  //---------------------------------------------Binning-------------------------------------
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      sprintf(histname, "%s/Detector%s_typ_%i_eta0_%i", BinDir, Dimtag, ity,  var[ivar]);	 
      inputbinning->GetObject(histname, binsRec[ity][ivar]);  cout << histname <<endl;
      //binsRec[ity][ivar][ipt]->PrintStream(cout);        
      sprintf(histname, "%s/Generator%s_typ_%i_eta0_%i", BinDir, Dimtag, ity, var[ivar]);
      inputbinning->GetObject(histname, binsGen[ity][ivar]);  cout << histname <<endl;
      //binsGen[ity][ivar][ipt]->PrintStream(cout);
    }
  }
  //----------------------------------------------Read Input Data MC and Response matrix
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
        sprintf(Axisname,"var_%i[UO]",var[ivar]);
	DirData->cd();
	Data_Reco[ity][ivar]= (TH1D*)ReadHist1D(HistDir+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),inputData)->Clone();
        //for(int i=1 ; i< Data_Reco[ity][ivar]->GetNbinsX()+1; i++){ if(Data_Reco[ity][ivar]->GetBinContent(i) == 0) {cout << " Data Bin entry Nil : bin no : "<< i << endl;}}
	HT2_NormalV3(Data_Reco[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
#ifdef CLOUSER
        Data_Gen[ity][ivar]=(TH1D*)ReadHist1D( HistDir+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),inputData)->Clone();
	HT2_NormalV3(Data_Gen[ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
	
        TH1D *Datafake= (TH1D*)ReadHist1D( HistDir+"/dd_fake_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),inputData)->Clone();
        TH1D *Datamiss= (TH1D*)ReadHist1D( HistDir+"/dd_miss_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),inputData)->Clone();
	
        Data_fakerate[ity][ivar] = (TH1D*)Datafake->Clone();   Data_fakerate[ity][ivar]->Divide(Data_fakerate[ity][ivar], Data_Reco[ity][ivar], 1, 1, "b");
        Data_fakerate[ity][ivar]->SetNameTitle(Form("%sData_fakerate_%i_eta0_%i", Histtag, ity, var[ivar]),Form("%sData_fakerate_%i_eta0_%i", Histtag, ity, var[ivar]));
	
	Data_fakerateInv[ity][ivar]=(TH1D*)Data_fakerate[ity][ivar]->Clone();  Data_fakerateInv[ity][ivar]->Reset();
        for (int i = 1; i <= Data_fakerate[ity][ivar]->GetNbinsX(); ++i) {
	  double factor = Data_fakerate[ity][ivar]->GetBinContent(i); 
	  double content =1-factor;  Data_fakerateInv[ity][ivar]->SetBinContent(i, content);
	}	
	Data_missrate[ity][ivar]= (TH1D*)Datamiss->Clone();  Data_missrate[ity][ivar]->Divide(Data_missrate[ity][ivar], Data_Gen[ity][ivar], 1, 1, "b");
        sprintf(name,"%sDatamiss_rate_%i_eta0_%i", Histtag,ity, var[ivar]); Data_missrate[ity][ivar]->SetNameTitle(name,name);
	
	Data_missrateInv[ity][ivar]=(TH1D*)Data_missrate[ity][ivar]->Clone();  Data_missrateInv[ity][ivar]->Reset();
        for (int i = 1; i <= Data_missrate[ity][ivar]->GetNbinsX(); ++i) {
	  double factor = Data_missrate[ity][ivar]->GetBinContent(i);
	  double content =1-factor; Data_missrateInv[ity][ivar]->SetBinContent(i, content);
	}
        sprintf(name,"%sDataRecominusfake_%i_eta0_%i",Histtag ,ity, var[ivar]);
        Data_reco_fake[ity][ivar] = (TH1D*)Data_fakerateInv[ity][ivar]->Clone();
        Data_reco_fake[ity][ivar]->SetNameTitle(name,name); 
	Data_reco_fake[ity][ivar]->Multiply(Data_Reco[ity][ivar]);
	Extract(Data_reco_fake[ity][ivar], binsRec[ity][ivar], Axisname);

        sprintf(name,"%sDataGenminusmiss_%i_eta0_%i",Histtag ,ity, var[ivar]);
        Data_gen_miss[ity][ivar] = (TH1D*)Data_missrateInv[ity][ivar]->Clone(); 
        Data_gen_miss[ity][ivar]->SetNameTitle(name,name);
	Data_gen_miss[ity][ivar]->Multiply(Data_Gen[ity][ivar]);
        Extract(Data_gen_miss[ity][ivar], binsGen[ity][ivar], Axisname);
#endif
	//--------------------------------------------MC for Read RM------------------------------------
	DirRMinput->cd();
        RMinput_Reco[ity][ivar]= (TH1D*)ReadHist1D(HistDir+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),RMinput)->Clone();
        for(int i=1 ; i< RMinput_Reco[ity][ivar]->GetNbinsX()+1; i++){ if(RMinput_Reco[ity][ivar]->GetBinContent(i) == 0) {cout << " RMinput Bin entry Nil : bin no : "<< i << endl;}}
	HT2_NormalV3(RMinput_Reco[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
        
	RMinput_Gen[ity][ivar]=(TH1D*)ReadHist1D( HistDir+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),RMinput)->Clone();
	HT2_NormalV3(RMinput_Gen[ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
	
        RMinput_RM[ity][ivar] = (TH2D*)ReadHist2D(HistDir+"/dd_corr_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),RMinput)->Clone();
	
        RMinput_fake[ity][ivar]= (TH1D*)ReadHist1D( HistDir+"/dd_fake_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),RMinput)->Clone();
        RMinput_miss[ity][ivar]= (TH1D*)ReadHist1D( HistDir+"/dd_miss_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),RMinput)->Clone();
	
        RMinput_fakerate[ity][ivar] = (TH1D*) RMinput_fake[ity][ivar]->Clone();   RMinput_fakerate[ity][ivar]->Divide(RMinput_fakerate[ity][ivar], RMinput_Reco[ity][ivar], 1, 1, "b");
        RMinput_fakerate[ity][ivar]->SetNameTitle(Form("%sRMinputfake_rate_%i_eta0_%i", Histtag, ity, var[ivar]),Form("%s fakerate %i eta0 %i", Histtag, ity, var[ivar]));
	HT2_NormalV3(RMinput_fakerate[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
	
        RMinput_fakerateInv[ity][ivar]=(TH1D*)RMinput_fakerate[ity][ivar]->Clone();  RMinput_fakerateInv[ity][ivar]->Reset();
        for (int i = 1; i <= RMinput_fakerate[ity][ivar]->GetNbinsX(); ++i) {
          double factor = RMinput_fakerate[ity][ivar]->GetBinContent(i);
          double content =1-factor;  RMinput_fakerateInv[ity][ivar]->SetBinContent(i, content);
        }
        RMinput_missrate[ity][ivar]= (TH1D*)RMinput_miss[ity][ivar]->Clone();  RMinput_missrate[ity][ivar]->Divide(RMinput_missrate[ity][ivar], RMinput_Gen[ity][ivar], 1, 1, "b");
        RMinput_missrate[ity][ivar]->SetNameTitle(Form("%sRMinputmiss_rate_%i_eta0_%i", Histtag,ity, var[ivar]),Form("%s missrate%i eta0 %i", Histtag,ity, var[ivar]));
	HT2_NormalV3(RMinput_missrate[ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
	
        RMinput_missrateInv[ity][ivar]=(TH1D*)RMinput_missrate[ity][ivar]->Clone();  RMinput_missrateInv[ity][ivar]->Reset();
        for (int i = 1; i <= RMinput_missrate[ity][ivar]->GetNbinsX(); ++i) {
          double factor = RMinput_missrate[ity][ivar]->GetBinContent(i);
          double content =1-factor; RMinput_missrateInv[ity][ivar]->SetBinContent(i, content);
        }
        RMinput_reco_fake[ity][ivar] = (TH1D*)RMinput_fakerateInv[ity][ivar]->Clone();
        RMinput_reco_fake[ity][ivar]->SetNameTitle(Form("%sRMinputRecominusfake_%i_eta0_%i",Histtag ,ity, var[ivar]),Form("Reco - fake %i eta0 %i",Histtag ,ity, var[ivar]));
        RMinput_reco_fake[ity][ivar]->Multiply(RMinput_Reco[ity][ivar]);
	HT2_NormalV3(RMinput_reco_fake[ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
	
        RMinput_gen_miss[ity][ivar] = (TH1D*)RMinput_missrateInv[ity][ivar]->Clone();
        RMinput_gen_miss[ity][ivar]->SetNameTitle(Form("%sRMinputGenminusmiss_%i_eta0_%i",Histtag ,ity, var[ivar]),Form("Gen - miss %i eta0 %i",Histtag ,ity, var[ivar]));
        RMinput_gen_miss[ity][ivar]->Multiply(RMinput_Gen[ity][ivar]);
	HT2_NormalV3(RMinput_gen_miss[ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
//--------------------------------------------Read MC---------------------------------------------	
	for (int imc=0; imc<nmc ; imc++){
	  cout <<  " MC number : " << imc << endl;
	  outputDir[imc]->cd(); //..........................MC directory .......................................
	  MC_Reco[imc][ity][ivar] = (TH1D*)ReadHist1D(HistDir+"/dd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), inputMC[imc])->Clone();
          MC_Gen[imc][ity][ivar] = (TH1D*)ReadHist1D(HistDir+"/dd_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), inputMC[imc])->Clone();
          MC_fake[imc][ity][ivar] = (TH1D*)ReadHist1D(HistDir+"/dd_fake_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), inputMC[imc])->Clone();
          MC_miss[imc][ity][ivar] = (TH1D*)ReadHist1D(HistDir+"/dd_miss_gen_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), inputMC[imc])->Clone();
	  
	  HT2_NormalV3(MC_Reco[imc][ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
	  HT2_NormalV3(MC_Gen[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
          HT2_NormalV3(MC_fake[imc][ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx);
	  HT2_NormalV3(MC_miss[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx);
	  
	  cout << " Fake= " <<MC_fake[imc][ity][ivar]->GetEntries() <<" Reco-fake: " <<(MC_Reco[imc][ity][ivar]->GetEntries() -  MC_fake[imc][ity][ivar]->GetEntries())<<endl;
	  cout << " Miss= " <<MC_miss[imc][ity][ivar]->GetEntries() <<" Gen-Miss: " << (MC_Gen[imc][ity][ivar]->GetEntries() - MC_miss[imc][ity][ivar]->GetEntries()) ;
	  //Response Matrix
	  h2dGenDetMC[imc][ity][ivar] = (TH2D*)ReadHist2D(HistDir+"/dd_corr_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar]), inputMC[imc])->Clone();
	  cout << " Corr = "  << h2dGenDetMC[imc][ity][ivar]->GetEntries() <<endl;
	  //--------------------------------------Calculate Fake rate and Miss Rate
	  TH1D* fakerate = (TH1D*)MC_fake[imc][ity][ivar]->Clone();
	  fakerate->Divide(fakerate,MC_Reco[imc][ity][ivar], 1, 1, "b");
	  sprintf(name,"%sfake_rate_%i_eta0_%i", Histtag, ity, var[ivar]); fakerate->SetNameTitle(name,Form("%s fake rate %i eta0 %i", Histtag, ity, var[ivar]));
	  TH1D* missrate = (TH1D*)MC_miss[imc][ity][ivar]->Clone();
	  missrate->Divide(missrate, MC_Gen[imc][ity][ivar], 1, 1, "b");
	  sprintf(name,"%smiss_rate_%i_eta0_%i", Histtag,ity, var[ivar]); missrate->SetNameTitle(name,Form("%s miss rate %i eta0 %i", Histtag,ity, var[ivar]));
	  fakerate->SetMinimum(-0.05); fakerate->SetMaximum(1.01);  missrate->SetMinimum(-0.05); missrate->SetMaximum(1.01);
	  MC_fakerate[imc][ity][ivar] = (TH1D*)fakerate->Clone();  MC_missrate[imc][ity][ivar] = (TH1D*)missrate->Clone();
	  HT2_NormalV3(MC_fakerate[imc][ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx,1);
	  HT2_NormalV3(MC_missrate[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx,1);
	  //--------------------------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20
	  TH1* RMx = h2dGenDetMC[imc][ity][ivar]->ProjectionX(); TH1* RMy = h2dGenDetMC[imc][ity][ivar]->ProjectionY();
	  sprintf(name,"%sProjectX_%i_eta0_%i", Histtag, ity, var[ivar]); RMx->SetNameTitle(name,name);
	  sprintf(name,"%sProjectY_%i_eta0_%i", Histtag,ity, var[ivar]); RMy->SetNameTitle(name,name);
	  
	  sprintf(name,"%sRecominusfake_%i_eta0_%i",Histtag ,ity, var[ivar]);
	  MC_reco_fake[imc][ity][ivar] = (TH1D*)MC_Reco[imc][ity][ivar]->Clone(); MC_reco_fake[imc][ity][ivar]->Reset();
	  MC_reco_fake[imc][ity][ivar]->SetNameTitle(name,name);
	  
	  sprintf(name,"%sGenminusmiss_%i_eta0_%i",Histtag ,ity, var[ivar]);
	  MC_Gen_miss[imc][ity][ivar] = (TH1D*)MC_Gen[imc][ity][ivar]->Clone();  MC_Gen_miss[imc][ity][ivar]->Reset();
	  MC_Gen_miss[imc][ity][ivar]->SetNameTitle(name,name);
	  
	  for (int i = 1; i <= MC_reco_fake[imc][ity][ivar]->GetNbinsX(); ++i) {
	    double content =  MC_Reco[imc][ity][ivar]->GetBinContent(i); double factor =  MC_fake[imc][ity][ivar]->GetBinContent(i);
	    content -= factor;  MC_reco_fake[imc][ity][ivar]->SetBinContent(i, content);
	    MC_reco_fake[imc][ity][ivar]->SetBinError(i, MC_fake[imc][ity][ivar]->GetBinError(i)+ MC_Reco[imc][ity][ivar]->GetBinError(i));
	  }
	  
	  for (int i = 1; i <=  MC_Gen_miss[imc][ity][ivar]->GetNbinsX(); ++i) {
	    double content = MC_Gen[imc][ity][ivar]->GetBinContent(i); double factor = MC_miss[imc][ity][ivar]->GetBinContent(i);
	    content -= factor;   MC_Gen_miss[imc][ity][ivar]->SetBinContent(i, content);
	    MC_Gen_miss[imc][ity][ivar]->SetBinError(i, MC_Gen_miss[imc][ity][ivar]->GetBinError(i)+ MC_Gen[imc][ity][ivar]->GetBinError(i));
	  }
	  MC_reco_fake[imc][ity][ivar]->Write();  MC_Gen_miss[imc][ity][ivar]->Write(); RMx->Write(); RMy->Write(); 
          HT2_NormalV3(MC_reco_fake[imc][ity][ivar], binsRec[ity][ivar], Axisname,nHLTmx,1);
          HT2_NormalV3(MC_Gen_miss[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx,1);
	  
	  //------------------------------------------------------Stability and Purity -----------------------------------------
	  TH1* hist_pu= (TH1D*)MC_miss[imc][ity][ivar]->Clone(); hist_pu->Reset();
	  TH1* hist_st= (TH1D*)MC_miss[imc][ity][ivar]->Clone(); hist_st->Reset();
	  TH2* RMcopy = (TH2D*)h2dGenDetMC[imc][ity][ivar]->Clone(); RMcopy->RebinX(2);
	  for(int binRec=1; binRec<= hist_pu->GetNbinsX(); binRec++) {
	    double sum=0.;
	    for(int binGen=1; binGen<=hist_pu->GetNbinsX(); binGen++) {
	      //sum += RMcopy->GetBinContent(binGen,binRec);
	      sum += RMcopy->GetBinContent(binRec,binGen);
	   }
	  double p=0.;
	   if(sum>0.0) {
	      p = RMcopy->GetBinContent(binRec,binRec)/sum;
	    }
	    hist_pu->SetBinContent(binRec,p);
	  }
	  hist_pu->SetMinimum(-0.07); hist_pu->SetMaximum(1.01);
	  //hist_pu->Write();
	  //TH2* RMcopy = (TH2D*)RM_RecoGen->Clone(); RMcopy->RebinX(2);
	  TH1* RMxcopy= (TH1D*)RMx->Clone(); RMxcopy->Rebin(2);
	  //for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_pu->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
	  for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_st->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
	  //hist_pu->Divide(hist_pu,RMxcopy, 1, 1, "b");
	  hist_st->Divide(hist_st,RMy, 1, 1, "b");
	  
	  //TH1* RMMx= (TH1D*)RMx->Clone(); RMMx->Rebin(2);
	  //hist_pu->Divide(RMy,RMMx, 1, 1, "b");
	  
	  TH1D *h_pu = (TH1D*)fakerate->Clone(); h_pu->Reset();
	  TH1D *h_st = (TH1D*)fakerate->Clone(); h_st->Reset();
	  int ir = h_pu->GetNbinsX(); int ig = missrate->GetNbinsX();
	  double fk[ir]; double ef[ig]; double pu[ir]; double st[ir];
	  
	  for(int i =1; i<h_pu->GetNbinsX()+1; i++){
	    h_pu->SetBinContent(i,pu[i]);
	    h_st->SetBinContent(i,st[i]); 
	  }
	  hist_purity[imc][ity][ivar]=(TH1D*)hist_pu->Clone(); 
	  hist_stbl[imc][ity][ivar]=(TH1D*)hist_st->Clone(); 

	  sprintf(name,"%sPurity_%i_eta0_%i",Histtag, ity, var[ivar]);
	  hist_purity[imc][ity][ivar]->SetNameTitle(name,name);
	  sprintf(name,"%sstability_%i_eta0_%i", Histtag, ity, var[ivar]);
	  hist_stbl[imc][ity][ivar]->SetNameTitle(name,name);
	  
          HT2_NormalV3(hist_purity[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx,1);
          HT2_NormalV3(hist_stbl[imc][ity][ivar], binsGen[ity][ivar], Axisname,nHLTmx,1);
	} //for (int imc=0;imc<nmc;imc++)
      }//for(int ivar=0; ivar < nusedvar ; ivar ++)
    }//for(int ity=0; ity <type; ity++)
    cout << "Histogram Read Data and MC Done " <<endl;
    
    //------------------Fold check : Patrick 1 Sep 20
    //Get Probability Matrix  (gen-miss)*probability = (reco-fake). Of course, don't forget to account for miss and fake entries (if applicable).
    folddir->cd();
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	TH2D* RM  = (TH2D*)h2dGenDetMC[umc][ity][ivar]->Clone();
	TH1D* Reco  = (TH1D*)MC_Reco[umc][ity][ivar]->Clone();
	TH1D* Gen  = (TH1D*)MC_Gen[umc][ity][ivar]->Clone();
	TH1D* fake= (TH1D*)MC_fake[umc][ity][ivar]->Clone();
	TH1D* miss = (TH1D*)MC_miss[umc][ity][ivar]->Clone();
	
	RM->RebinY(irbin);Gen->Rebin(irbin);miss->Rebin(irbin);
	
	TH1D* Folded = (TH1D*)MC_Reco[umc][ity][ivar]->Clone(); Folded->Reset();
	Fold(RM, Reco, Gen, miss, fake, Folded);
	
	sprintf(name,"%sFold_%i_eta0_%i",Histtag, ity, var[ivar]);
	Folded->SetNameTitle(name,name);
        Extract(Folded, binsRec[ity][ivar], Axisname);	
	Folded->Write();
      }
    }
    //----------------------------Condition  number of Probability Matrix
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	TH2D* RM  = (TH2D*)RMinput_RM[ity][ivar]->Clone();
	TH1D* miss = (TH1D*)RMinput_miss[ity][ivar]->Clone();
	//TH1D* miss = (TH1D*)MC_missrate[umc][ity][ivar][ipt]->Clone();
	RM->RebinY(irbin); miss->Rebin(irbin);
	cout <<setw(2) << ity <<setw(5) <<ivar << setw(5) <<'\n';
	Condition(RM, miss);
      }
    }
    
    Unfolddir->cd();
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	cout << " TUnfoldBinning2D : "<< " typ : "<< ity << " : Var " << ivar << " ";
        
	double rxbins[MC_Reco[umc][ity][ivar]->GetNbinsX()+1]={0}; //get Reco bin array
	for(int ix=0; ix<MC_Reco[umc][ity][ivar]->GetNbinsX()+1; ix++){ rxbins[ix] = MC_Reco[umc][ity][ivar]->GetXaxis()->GetBinLowEdge(ix+1); }
	double gxbins[MC_Gen[umc][ity][ivar]->GetNbinsX()+1]={0}; //Get Gen Bin array
	for (int ix=0; ix<MC_Gen[umc][ity][ivar]->GetNbinsX()+1; ix++) {gxbins[ix] = MC_Gen[umc][ity][ivar]->GetXaxis()->GetBinLowEdge(ix+1); }
	
	//Rebin for match the condition of reco vs gen bin 
	RMinput_RM[ity][ivar]->RebinY(irbin); 	RMinput_Gen[ity][ivar]->Rebin(irbin);
	//----------------------------------------------------------Define Input---------------------	
        TUnfoldBinning* RecoBin = binsRec[ity][ivar];
        TUnfoldBinning* GenBin = binsGen[ity][ivar];
	
	TH2* RMin = (TH2D*)RMinput_RM[ity][ivar]->Clone();
	TH1* input = (TH1D*)Data_Reco[ity][ivar]->Clone();
	TH1* mc_missrate = (TH1D*)RMinput_missrate[ity][ivar]->Clone();
	TH1* mc_miss = (TH1D*)RMinput_miss[ity][ivar]->Clone();
        TH1* mcbackground = (TH1D*)RMinput_fakerate[ity][ivar]->Clone();
	TH1* mcgen = (TH1D*)RMinput_Gen[ity][ivar]->Clone();
	mcbackground->Multiply(input); //correction for Fake as background subtraction : Patrick
	
        double biasScale = 0;
        const char *REGULARISATION_DISTRIBUTION=0;
        const char *REGULARISATION_AXISSTEERING="*[UO]";
        
        //https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
        TH2D* covM = RecoBin->CreateErrorMatrixHistogram(Form("%sData_covariance_%i_eta0_%i",Histtag, ity, var[ivar]),false); covM->Sumw2();
	for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
	  double err = input->GetBinError(ix);
	  covM->SetBinContent(ix,ix,err*err);
	}
	
	covM->Write();
	//----------------------------------------------------------No Reguratization-------------------------------------	
        TUnfoldDensity density(RMin,TUnfold::kHistMapOutputVert,
			       TUnfoldDensity::kRegModeNone,
			       TUnfoldDensity::kEConstraintNone,
			       TUnfoldDensity::kDensityModeNone,
			       GenBin,RecoBin, REGULARISATION_DISTRIBUTION, REGULARISATION_AXISSTEERING);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);
	
        density.SubtractBackground(mcbackground, "fake", 1.0, 0.00); // hist,name,scale, scale error 
        int status = density.SetInput(input,biasScale,0,covM);
	
        int nBadErrors = status%10000, nUnconstrOutBins = status/10000; cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;
	
        //tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result
        density.DoUnfold(0.0);//,input,biasScale);  // tau 0.0 means no regularization

        sprintf(histname, "%sTUnfold_NoReg_typ_%i_eta0_%i", Histtag, ity, var[ivar]); //unfolded_typ_0_pt2_eta0_3
        sprintf(title, "Unfolded No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
	bool AxisBin = false;
	
	TH1 *Unfolded = density.GetOutput(histname,title,0,"*[UO]", AxisBin);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
        
        sprintf(title, "Input Ematrix No Reg %s 2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *hist_In_Emat = density.GetEmatrixInput(Form("%sInEmat_NoReg_typ_%i_eta0_%i",Histtag, ity, var[ivar]),title);//,"signal");
	
        sprintf(title, "RhoIJtotal No Reg %s  2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *Hist_RhoIJ = density.GetRhoIJtotal(Form("%scorr_NoReg_typ_%i_eta0_%i", Histtag, ity, var[ivar]), title);//,"signal");
	
        sprintf(title, "Ematrix No Reg %s 2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *hist_Emat = density.GetEmatrixTotal(Form("%sEmat_NoReg_typ_%i_eta0_%i",Histtag, ity, var[ivar]),title);//,"signal");
	
        sprintf(title, "Back folded No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
        TH1 *hist_folded = density.GetFoldedOutput(Form("%sRefold_NoReg_typ_%i_eta0_%i", Histtag,ity, var[ivar]),title,0,"*[UO]",AxisBin);//,"signal");
	
        sprintf(title, "ProbabilityMatrix No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
        TH2 *hist_prob = density.GetProbabilityMatrix(Form("%sProb_NoReg_typ_%i_eta0_%i",Histtag, ity, var[ivar]),title,TUnfold::kHistMapOutputVert);//,"signal");

        TH1 * BLTUnf = (TH1*)Unfolded->Clone();  //for BLT test
	//correction for Miss entries  : Partick
        for (int i = 1; i <= Unfolded->GetNbinsX(); ++i) {
          double content = Unfolded->GetBinContent(i);
          double factor = 1;
        //factor += (mc_missrate->GetBinContent(i)/(1-mc_missrate->GetBinContent(i)));
          factor += (mc_miss->GetBinContent(i)/(mcgen->GetBinContent(i) - mc_miss->GetBinContent(i)));
          content *= factor;
          Unfolded->SetBinContent(i, content);  }
	for (int i = 1; i <= Unfolded->GetNbinsX(); ++i) {cout<<" "<<Unfolded->GetBinContent(i)/mcgen->GetBinContent(i); } ; cout << endl;
	Unfolded->Write(); Hist_RhoIJ->Write(); hist_Emat->Write(); hist_folded->Write(); hist_prob->Write(); hist_In_Emat->Write();
	ConditionV2(hist_prob);
        //----------------------Extract The True Distrubuntions
	HT2_NormalV3(hist_folded,binsRec[ity][ivar], Axisname,nHLTmx,1);
	HT2_NormalV3(Unfolded, binsGen[ity][ivar], Axisname,nHLTmx,1);
	
#ifdef BLTest  
        int ib =0;	
        TH1 * BLTdata = (TH1*)Data_Reco[ity][ivar]->Clone(); //Data Reco
        TH1 * BLTreco = (TH1*)MC_Reco[ib][ity][ivar]->Clone(); //Data Reco
        TH1 * BLTdata1 = (TH1*)RMinput_fakerateInv[ity][ivar]->Clone(); BLTdata1->Multiply(Data_Reco[ity][ivar]);  //Data -background
	
        TH1 * BLTgen = (TH1*)MC_Gen[ib][ity][ivar]->Clone(); //Data Reco
        TH1 * BLTunf = (TH1*)Unfolded->Clone(); //Data Reco

        sprintf(name,"%sBLT_Data_covariance_%i_eta0_%i",Histtag, ity, var[ivar]);
        TH2D* covMBLT= RecoBin->CreateErrorMatrixHistogram(name,false); covMBLT->Sumw2();
        for (int ix=1; ix<BLTdata1->GetNbinsX()+1; ix++) {
          double err = BLTdata1->GetBinError(ix);
          covMBLT->SetBinContent(ix,ix,err*err);
       }
	//https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
        sprintf(name,"%sMC_covariance_%i_eta0_%i",Histtag, ity, var[ivar]);
        TH2D* covUf = GenBin->CreateErrorMatrixHistogram(name,false); covUf->Sumw2();
        for (int ix=1; ix<covUf->GetNbinsX()+1; ix++) {
          double err = Unfolded->GetBinError(ix);
          covUf->SetBinContent(ix,ix,err*err);  
	}
	cout << " From Density Chi2A() : " << density.GetChi2A() << " " << density.GetChi2L()  << " " << density.GetNdf() <<" Chi2/NDf : "  <<density.GetChi2A()/(density.GetNdf()+1) <<endl;
        cout << " chi2 det level 1: "; BLT(BLTdata, covM, MC_Reco[ib][ity][ivar], 1);
        cout << " chi2 det level 2: "; BLT(BLTdata1, covM, MC_reco_fake[ib][ity][ivar], 1); 
        cout << " chi2 det level 3: "; BLT(hist_folded, covMBLT, MC_reco_fake[ib][ity][ivar], 1) ; cout << endl;
	
        cout << " Unfold before Miss input uncert: "; BLT(BLTUnf, hist_In_Emat, MC_Gen_miss[ib][ity][ivar]);
        cout << " Unfold before Miss Tolat: "; BLT(BLTUnf, hist_Emat, MC_Gen_miss[ib][ity][ivar]); cout << endl;
	
	cout << " Unfold after Miss 2: "; BLT(Unfolded, hist_In_Emat,MC_Gen[ib][ity][ivar]);	
	cout << " Unfold after Miss 1: "; BLT(Unfolded, hist_Emat, MC_Gen[ib][ity][ivar]);

	BLTdata->Rebin(2);
	BLTreco->Rebin(2);
	Integralhist(BLTdata);
        Integralhist(BLTreco);
	TH1D* BLTDetRatio= (TH1D*)BLTreco->Clone(); BLTDetRatio->SetNameTitle(Form("BLTDet_typ_%i_eta0_%i",ity,var[ivar]),Form("BLTDet typ%i eta0 %i",ity,var[ivar]));
        BLTDetRatio->Divide(BLTdata);
        HT2_NormalV3(BLTDetRatio, binsGen[ity][ivar], Axisname,nHLTmx,1);

        Integralhist(Unfolded);
        Integralhist(MC_Gen[ib][ity][ivar]);
        TH1D* BLTGenRatio= (TH1D*)MC_Gen[ib][ity][ivar]->Clone(); BLTGenRatio->SetNameTitle(Form("BLTgen_typ_%i_eta0_%i",ity,var[ivar]),Form("BLTgen typ%i eta0 %i",ity,var[ivar]));
        BLTGenRatio->Divide(Unfolded);
	HT2_NormalV3(BLTGenRatio, binsGen[ity][ivar], Axisname,nHLTmx,1);
#endif
	//-----------------------------------------End of Variables loop	
	cout << endl ; 
      }
    }
   cout << "TUnfoldBinning 2D OK " <<endl;
  
    delete outputFile;
}

//-----------------------------------------------------------------------

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


//Condition number calculation------------------Using Patrick code 
void Condition (TH2* RM, TH1* miss){
    const int Nx = RM->GetNbinsX(),
              Ny = RM->GetNbinsY();
    cout << Nx <<"  " << Ny << endl;
    if (Ny*2 != Nx) { cout << Nx << ' ' << Ny << endl;  return; }

    TH1D* RMy = RM->ProjectionY("RMy", 0, -1); //Gen Projection

    // normalisation & condition
    TMatrixD m(Ny,Ny);
    for (int i = 1; i <= Ny; ++i) {
        double normalisation = RMy->GetBinContent(i);
      //         normalisation += miss->GetBinContent(i);
        if (normalisation > 0)
        for (int j = 1; j <= Nx; ++j) {
            double content = RM->GetBinContent(j,i);
           content /= normalisation;
            m((j-1)/2,i-1) += content;
        }
    }
    TDecompSVD svd(m);
    TVectorD v = svd.GetSig();
    cout <<" Condition :" << svd.Condition() << endl;
    double Max = v[0];
    for(int k = 0; k < Ny; ++k) {
        if (abs(v[k]) < feps) break;
    //    cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/v[k] << '\n';
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


//Bottom Line test --------------
void BLT (TH1 * dataDistX, TH2 * dataCovX, TH1 * MCX, int rebin = 1)
{

    TH1 * dataDist = (TH1*)dataDistX->Clone();
    TH2 * dataCov = (TH2*)dataCovX->Clone();
    TH1 * MC = (TH1*)MCX->Clone();

    dataDist->Rebin(rebin);
    MC->Rebin(rebin);
    dataCov->Rebin2D(rebin, rebin);
    vector<int> indices;
    for (int i = 1; i <= dataCov->GetNbinsX(); ++i) 
        if (dataCov->GetBinContent(i,i) > 0) indices.push_back(i);
    int ndf = indices.size();

    TVectorD v(ndf);
    TMatrixD M(ndf,ndf);
    for (int i = 0; i < ndf; ++i) {
        int index = indices.at(i);
        double iData = dataDist->GetBinContent(index),
               iMC   = MC      ->GetBinContent(index);
        v(i) = iData-iMC;
        for (int j = 0; j < ndf; ++j) {
            int jndex = indices.at(j);
            M(i,j) = dataCov->GetBinContent(index,jndex);
        }
    }
//    TMatrixD vT = v;
  //  vT.T();

    M.Invert();
    
    double chi2 = v*(M*v);

    cout << "chi2/ndf = " << chi2 << " / " << ndf << " = " << (chi2/ndf) << endl;
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
//-----------------------------------Using Ashlye code
void ConditionV2(TH2* prob_mat){
       TH2* prob = (TH2*)prob_mat->Clone();
       prob->RebinY(2);

int nbiny =prob->GetNbinsX(); int nbinx =prob->GetNbinsY();
cout << "Bin " << nbiny << " " << nbinx ;
TMatrixD  matd(nbinx, nbiny);

for( int i =0; i <nbiny ; i++){
    for( int j =0; j <nbinx ; j++){
matd[i][j] = prob->GetBinContent( j+1, i+1);
    }
}

TDecompSVD svd(nbinx,nbiny);
svd.SetMatrix(matd);

cout <<" Condition :" << svd.Condition() << endl;
svd.Decompose();

TVectorD sig = svd.GetSig();
cout << " SigmaMax : " << sig[0] <<endl;
for(int i =0; i<nbiny; i++){
//cout << i << " condition : " << sig[0]/sig[i] << endl;
}

}

void Integralhist(TH1 *hist){ hist->Scale(1/(hist->Integral()));}

double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/)
{
  //hGen->Print("all");
  //covmat->Print("all");
  //hData->Print("all");
  //throw;
  int n = hData->GetNbinsX();
  if(skip)
    n -= 1;
  TMatrixD res(1, n);
  for(int i = 0; i < n; i++)
    res(0, i) = hData->GetBinContent(i + 1) - ((hGen) ? hGen->GetBinContent(i + 1) : 0.0);
  //res.Print("all");
  //hData->Print("all");
  //hGen->Print("all");
  TMatrixD resT = res;
  resT.T();
  TMatrixD cov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      cov(i, j) = covmat->GetBinContent(i + 1, j + 1);
  cov.Invert();
  double chi2 = (res * cov * resT)(0, 0);
  //hData->Print("all");
  //hGen->Print("all");
  //printf("CHI2 %.6f\n", chi2);
  return chi2;
}


void Normalise(TH1* h, TH2* covmat)
{
  int n = h->GetNbinsX();
  if(n != covmat->GetNbinsX() || n != covmat->GetNbinsY())
  {
    printf("Error in Normalise() incosistent input h->GetNbinsX() = %d covmat->GetNbinsX() = %d covmat->GetNbinsY() = %d\n", h->GetNbinsX(), covmat->GetNbinsX(), covmat->GetNbinsY());
    throw;
  }
  double integral = h->Integral();
  double integral2 = integral * integral;

  TMatrixD matrixG(n, n);
  TMatrixD matrixCov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
    {
      matrixCov(i, j) = covmat->GetBinContent(i + 1, j + 1);
      if(i == j)
        matrixG(i, i) = (integral - h->GetBinContent(i + 1)) / integral2;
      else
        matrixG(i, j) = -1 * h->GetBinContent(i + 1) / integral2;
    }
  TMatrixD matrixGT = matrixG;
  matrixGT.T();
  TMatrixD res = matrixG * matrixCov * matrixGT;

  for(int i = 0; i < n; i++)
  {
    h->SetBinContent(i + 1, h->GetBinContent(i + 1) / integral);
    h->SetBinError(i + 1, TMath::Sqrt(res(i, i)));
    for(int j = 0; j < n; j++)
      covmat->SetBinContent(i + 1, j + 1, res(i, j));
  }
}
//---------------------------------------------------------------------------------
void Extract(TH1* global2d, TUnfoldBinning* Bin, char* Axisname, bool iw=0){
	if(iw){global2d->Write();}
	char  name[100], title[100];
	sprintf(name,"E%s", global2d->GetName()); sprintf(title,"True bin %s", global2d->GetName());
	auto extract=  dynamic_cast<TH2*>(Bin->ExtractHistogram(name, global2d, 0, true, Axisname));
	extract->SetTitle(title);

        auto extractpx = extract->ProjectionX(); 
	sprintf(title,"Projection %s", global2d->GetName());
	extractpx->SetTitle(title);

	extract->Write();
	extractpx->Write();
}
//----------------------------------------------------------------------
TH1D* ReadHist1D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH1D* hist=(TH1D*)root->Get(histname);
hist->Write();
return hist;
}
//-----------------------------------------------------------
TH2D* ReadHist2D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH2D* hist=(TH2D*)root->Get(histname);
hist->Write();
return hist;
}
//-------------------------------------------------
void HT2_NormalV3(TH1* global2d, TUnfoldBinning* Bin, char* Axisname, int nht, bool iw=0){
    if(iw){global2d->Write();}
    char  name[100], title[100];
    sprintf(name,"E%s", global2d->GetName()); sprintf(title,"True bin %s", global2d->GetName());
    auto extract=  dynamic_cast<TH2*>(Bin->ExtractHistogram(name, global2d, 0, true, Axisname));
    extract->SetTitle(title);
    extract->Write();
int nbiny= extract->GetNbinsY();  int nbinx= extract->GetNbinsX();

for (int iht = 1; iht <= nbiny; ++iht) {
       const char * name_hbin = Form("%s_pt%i",extract->GetName(), iht-1);
        TH1D* h_var  = extract->ProjectionX(name_hbin, iht, iht);
	h_var->SetTitle(Form("HT2 %s pt%i",global2d->GetTitle(),iht-1));
        h_var->Write();
      }
  }


