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
  
  const TString Pyinput = "PY8_UL17_2D_ALLHT2.root";
  //const TString datainput = "PY8_UL17_2D_ALLHT2.root"; 
  //const TString datainput = "PY8_UL17_Flat_Hist_HTmarge.root"; 
  const TString datainput = "MG_UL17_HT2_83_3000.root"; 
  //const TString Pyinput = "PY8_UL17_Binned_TUnfold.root";
  TFile *inputbinning = new TFile("PY8_HT2_83_3000_Binnied.root");
  
  //Input Data and MC histogram
  TFile *inputData=new TFile(datainput);
  
  TFile *inputMC[nmc];
  inputMC[0]=new TFile(Pyinput);
  
  inputMC[1]=new TFile("MG_UL17_HT2_83_3000.root");
  
  inputMC[2]=new TFile("PY8_UL17_2D_ALLHT2.root");
  
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
  TH1*  GetLocalBinnedHist(TH1* Hist,  TUnfoldBinning* bin, char const* none, char const* axis);
  void BLT (TH1 * dataDist, TH2 * dataCov, TH1 * MC, int rebin = 1);
  void ConditionV2(TH2* prob_mat);
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void Integralhist(TH1 *hist);
  double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/);
  void Normalise(TH1* h, TH2* covmat);
  //----------------------------------Different Binning ----------------------------------------------
  for(int idd=2; idd <3; idd++){  //0 : Root Hist 1: 1D TunfoldBinning 2: 2D TUnfoldBinning
    
    TH1D *MC_Reco[nmc][type][nusedvar];  //Reconstructed MC
    TH1D *MC_fake[nmc][type][nusedvar];  //Fake :  Reco but No Gen
    TH1D *MC_fakerate[nmc][type][nusedvar];  //Fake :  Reco but No Gen
    TH1D *MC_reco_fake[nmc][type][nusedvar];  //Fake :  Reco but No Gen
    
    TH1D *MC_Gen[nmc][type][nusedvar];   //Generator MC
    TH1D *MC_miss[nmc][type][nusedvar];   //Miss:  No Reco but in Gen
    TH1D *MC_Gen_miss[nmc][type][nusedvar];   //Miss:  No Reco but in Gen
    TH1D *MC_missrate[nmc][type][nusedvar];   //Miss:  No Reco but in Gen
    
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
    
    TUnfoldBinning *binsRec[type][nusedvar]; //Binning Name 
    TUnfoldBinning *RecoBinning[type][nusedvar]; //Node Name
    
    TUnfoldBinning *binsGen[type][nusedvar];  //Gen Binning
    TUnfoldBinning *GenBinning[type][nusedvar];  //Gen Node 
    
    //---------------------------------------------Binning-------------------------------------
    if(idd>0){
      for(int ity=0; ity <type; ity++){
	for(int ivar=0; ivar < nusedvar ; ivar ++){
	  sprintf(histname, "%s/Detector%s_typ_%i_eta0_%i", Dirbin[idd], bintag[idd], ity,  var[ivar]);
	  cout << histname <<endl;
	  inputbinning->GetObject(histname, binsRec[ity][ivar]);
	  //binsRec[ity][ivar][ipt]->PrintStream(cout);        
	  sprintf(histname, "%s/Generator%s_typ_%i_eta0_%i", Dirbin[idd], bintag[idd], ity, var[ivar]);
	  inputbinning->GetObject(histname, binsGen[ity][ivar]);
	  //binsGen[ity][ivar][ipt]->PrintStream(cout);
	}
      }
    }
    //------------------------------------------
    //----------------------------------------------Read Input Data MC and Response matrix
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
     	
	sprintf(histname, "%s/%sreco_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]); 
	TH1D *DataReco =(TH1D*) inputData->Get(histname);
	for(int i=1 ; i< DataReco->GetNbinsX()+1; i++){ if(DataReco->GetBinContent(i) == 0) {cout << " Data Bin entry Nil : bin no : "<< i << endl;}}
#ifdef CLOUSER
	sprintf(histname, "%s/%sgen_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]); 
	TH1D *DataGen= (TH1D*)inputData->Get(histname);
        Data_Gen[ity][ivar]=(TH1D*)DataGen->Clone(); 
	
	sprintf(histname, "%s/%sfake_reco_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]);
        TH1D *Datafake= (TH1D*)inputData->Get(histname);
	
	sprintf(histname, "%s/%smiss_gen_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]);
        TH1D *Datamiss= (TH1D*)inputData->Get(histname);
	
        TH1D* Datafakerate = (TH1D*)Datafake->Clone();
        Datafakerate->Divide(Datafakerate, DataReco, 1, 1, "b");
        sprintf(name,"%sDatafake_rate_%i_eta0_%i", histtag[idd], ity, var[ivar]); Datafakerate->SetNameTitle(name,name);
        Data_fakerate[ity][ivar]= (TH1D*)Datafakerate->Clone();
	
	Data_fakerateInv[ity][ivar]=(TH1D*)Datafakerate->Clone();  Data_fakerateInv[ity][ivar]->Reset();
        for (int i = 1; i <= Datafakerate->GetNbinsX(); ++i) {
	  double factor = Datafakerate->GetBinContent(i); 
	  double content =1-factor;  Data_fakerateInv[ity][ivar]->SetBinContent(i, content);
	}
	
	
	TH1D* Datamissrate = (TH1D*)Datamiss->Clone();
        Datamissrate->Divide(Datamissrate, DataGen, 1, 1, "b");
        sprintf(name,"%sDatamiss_rate_%i_eta0_%i", histtag[idd],ity, var[ivar]); Datamissrate->SetNameTitle(name,name);
        Data_missrate[ity][ivar]= (TH1D*)Datamissrate->Clone();
        
	Data_missrateInv[ity][ivar]=(TH1D*)Datamissrate->Clone();  Data_missrateInv[ity][ivar]->Reset();
        for (int i = 1; i <= Datamissrate->GetNbinsX(); ++i) {
	  double factor = Datamissrate->GetBinContent(i);
	  double content =1-factor; Data_missrateInv[ity][ivar]->SetBinContent(i, content);
	}
        sprintf(name,"%sDataRecominusfake_%i_eta0_%i",histtag[idd] ,ity, var[ivar]);
        TH1* DataFakeCorrect = (TH1D*)DataReco->Clone(); DataFakeCorrect->Reset();
        DataFakeCorrect->SetNameTitle(name,name);
	
        sprintf(name,"%sDataGenminusmiss_%i_eta0_%i",histtag[idd] ,ity, var[ivar]);
        TH1* DataMissCorrect = (TH1D*)DataGen->Clone(); DataMissCorrect->Reset();
        DataMissCorrect->SetNameTitle(name,name);
	
	for (int i = 1; i <= DataFakeCorrect->GetNbinsX(); ++i) {
	  double content = DataReco->GetBinContent(i); double factor = Datafake->GetBinContent(i);
	  content -= factor;  DataFakeCorrect->SetBinContent(i, content);
	  DataFakeCorrect->SetBinError(i, Datafake->GetBinError(i)+ DataReco->GetBinError(i));
	}
	Data_reco_fake[ity][ivar]= (TH1D*)DataFakeCorrect->Clone();
	
	for (int i = 1; i <= DataMissCorrect->GetNbinsX(); ++i) {
	  double content = DataGen->GetBinContent(i); double factor = Datamiss->GetBinContent(i);
	  content -= factor;  DataMissCorrect->SetBinContent(i, content);
	  DataMissCorrect->SetBinError(i, Datamiss->GetBinError(i)+ DataGen->GetBinError(i));
	}
	Data_gen_miss[ity][ivar]= (TH1D*)DataMissCorrect->Clone();
	
#endif
	
	for (int imc=0; imc<nmc ; imc++){
	  //----------------------------------MC RECO
	  sprintf(histname, "%s/%sreco_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]); 
	  TH1D *RecoMC = (TH1D*)inputMC[imc]->Get(histname);      cout << histname ;
	  
	  //if(recobins!=Data_Reco[ity][ivar][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
	  //-----------------------------------MC Fake
	  sprintf(histname, "%s/%sfake_reco_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity,  var[ivar]); 
	  TH1D *RecoMC_Fake = (TH1D*)inputMC[imc]->Get(histname);
	  
	  cout << " Fake= " <<RecoMC_Fake->GetEntries() <<" Reco-fake: " <<(RecoMC->GetEntries() - RecoMC_Fake->GetEntries())<<endl;
	  
	  //-----------------------------------Gen MC
	  sprintf(histname, "%s/%sgen_typ_%i_eta0_%i",  Dirhist[idd],histtag[idd], ity, var[ivar]); 
	  TH1D *GenMC = (TH1D*)inputMC[imc]->Get(histname); cout << histname  ;
	  
	  //for(int i= 1; i < GenMC->GetNbinsX()+1; i++){if(GenMC->GetBinContent(i) == 0) { cout << " MC gen Bin is Zero for bin number :********** "<<  i  << endl; }}
	  
	  //----------------------------------MC miss 
	  sprintf(histname, "%s/%smiss_gen_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]); 
	  TH1D *GenMC_Miss = (TH1D*)inputMC[imc]->Get(histname);
	  
	  cout << " Miss= " << GenMC_Miss->GetEntries() <<" Gen-Miss: " << (GenMC->GetEntries() - GenMC_Miss->GetEntries()) ;
	  
	  //Response Matrix
	  sprintf(histname, "%s/%scorr_typ_%i_eta0_%i", Dirhist[idd],histtag[idd], ity, var[ivar]); 
	  TH2D *RM_RecoGen =(TH2D*)inputMC[imc]->Get(histname);  cout << " Corr = "  << RM_RecoGen->GetEntries() <<endl;
	  
	  //--------------------------------------Calculate Fake rate and Miss Rate
	  TH1D* fakerate = (TH1D*)RecoMC_Fake->Clone();
	  fakerate->Divide(fakerate, RecoMC, 1, 1, "b");
	  sprintf(name,"%sfake_rate_%i_eta0_%i", histtag[idd], ity, var[ivar]); fakerate->SetNameTitle(name,name);
	  TH1D* missrate = (TH1D*)GenMC_Miss->Clone();
	  missrate->Divide(missrate, GenMC, 1, 1, "b");
	  sprintf(name,"%smiss_rate_%i_eta0_%i", histtag[idd],ity, var[ivar]); missrate->SetNameTitle(name,name);
	  fakerate->SetMinimum(-0.05); fakerate->SetMaximum(1.01);
	  missrate->SetMinimum(-0.05); missrate->SetMaximum(1.01);
	  
	  
	  //--------------------------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20
	  TH1* RMx = RM_RecoGen->ProjectionX(); TH1* RMy = RM_RecoGen->ProjectionY();
	  sprintf(name,"%sProjectX_%i_eta0_%i", histtag[idd], ity, var[ivar]); RMx->SetNameTitle(name,name);
	  sprintf(name,"%sProjectY_%i_eta0_%i", histtag[idd],ity, var[ivar]); RMy->SetNameTitle(name,name);
	  
	  sprintf(name,"%sRecominusfake_%i_eta0_%i",histtag[idd] ,ity, var[ivar]);
	  TH1* RecoFakeCorrect = (TH1D*)RecoMC->Clone(); RecoFakeCorrect->Reset();
	  RecoFakeCorrect->SetNameTitle(name,name);
	  
	  sprintf(name,"%sGenminusmiss_%i_eta0_%i",histtag[idd] ,ity, var[ivar]);
	  TH1* GenMissCorrect = (TH1D*)GenMC->Clone(); GenMissCorrect->Reset();
	  GenMissCorrect->SetNameTitle(name,name);
	  
	  for (int i = 1; i <= RecoFakeCorrect->GetNbinsX(); ++i) {
	    double content = RecoMC->GetBinContent(i); double factor = RecoMC_Fake->GetBinContent(i);
	    content -= factor;  RecoFakeCorrect->SetBinContent(i, content);
	    RecoFakeCorrect->SetBinError(i, RecoMC_Fake->GetBinError(i)+ RecoMC->GetBinError(i));
	  }
	  MC_reco_fake[imc][ity][ivar]= (TH1D*)RecoFakeCorrect->Clone();
	  
	  for (int i = 1; i <= GenMissCorrect->GetNbinsX(); ++i) {
	    double content = GenMC->GetBinContent(i); double factor = GenMC_Miss->GetBinContent(i);
	    content -= factor;  GenMissCorrect->SetBinContent(i, content);
	    GenMissCorrect->SetBinError(i, GenMC_Miss->GetBinError(i)+ GenMC->GetBinError(i));
	  }
	  
	  MC_Gen_miss[imc][ity][ivar]= (TH1D*)GenMissCorrect->Clone();
	  TH1* RecoFakecorrEx; TH1* GenMisscorrEx; 
	  outputDir[imc]->cd(); //..........................MC directory .......................................
	    sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	    sprintf(name,"E%s", MC_reco_fake[imc][ity][ivar]->GetName());  sprintf(Axisname,"var_%i[UO]",var[ivar]);
	    RecoFakecorrEx =binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name,  MC_reco_fake[imc][ity][ivar], 0, true, Axisname);
	    
	    sprintf(histname,"Genbin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	    sprintf(name,"E%s", MC_Gen_miss[imc][ity][ivar]->GetName());   sprintf(Axisname,"var_%i[UO]",var[ivar]);
	    GenMisscorrEx =binsGen[ity][ivar]->FindNode(histname)->ExtractHistogram(name, MC_Gen_miss[imc][ity][ivar], 0, true, Axisname);
	    
	    GenMisscorrEx->Write(); RecoFakecorrEx->Write();   
	  //-------------------------------------------Extact Hist----------------------------- 
	    
	    sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd],ity ,var[ivar]); 
	    sprintf(name,"E%sreco_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]); 
	    sprintf(Axisname,"var_%i[UO]",var[ivar]);
	    TH1* RecoExtact =binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name, RecoMC, 0, true, Axisname);
	    
	    sprintf(histname,"Genbin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	    sprintf(name,"E%sgen_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	    sprintf(Axisname,"var_%i[UO]",var[ivar]);
	    TH1* GenExtact =binsGen[ity][ivar]->FindNode(histname)->ExtractHistogram(name, GenMC , 0, true, Axisname);
	    
	    RecoExtact->Write(); GenExtact->Write();
	    
	    //RecoMC =(TH1*) RecoExtact->Clone();
	    //GenMC =(TH1*) GenExtact->Clone();
	    
	    RecoMC->Write(); GenMC->Write();
	  
	  RMx->Write(); RMy->Write(); RecoFakeCorrect->Write(); GenMissCorrect->Write();
	  RecoMC_Fake->Write(); fakerate->Write();
	  GenMC_Miss->Write(); missrate->Write(); 
	  RM_RecoGen->Write();
	  //----------------------------------------------------------------------------
	  MC_Reco[imc][ity][ivar] = (TH1D*)RecoMC->Clone();
	  MC_fake[imc][ity][ivar] = (TH1D*)RecoMC_Fake->Clone();
	  MC_fakerate[imc][ity][ivar] = (TH1D*)fakerate->Clone();
	  MC_reco_fake[imc][ity][ivar] = (TH1D*)RecoFakeCorrect->Clone();
	  MC_fakerate[imc][ity][ivar]->SetMinimum(-0.05); MC_fakerate[imc][ity][ivar]->SetMaximum(1.01);
	  
	  
	  MC_Gen[imc][ity][ivar] = (TH1D*)GenMC->Clone();
	  MC_miss[imc][ity][ivar] = (TH1D*)GenMC_Miss->Clone();
	  MC_missrate[imc][ity][ivar] = (TH1D*)missrate->Clone();
	  MC_Gen_miss[imc][ity][ivar] = (TH1D*)GenMissCorrect->Clone();
	  MC_missrate[imc][ity][ivar]->SetMinimum(-0.05); MC_missrate[imc][ity][ivar]->SetMaximum(1.01);
	  
	  h2dGenDetMC[imc][ity][ivar] = (TH2D*)RM_RecoGen->Clone();
	  
	  
	  //------------------------------------------------------Stability and Purity -----------------------------------------
	  
	  //TH1* hist_pu= (TH1D*)MC_fake[imc][ity][ivar][ipt]->Clone(); hist_pu->Reset();
	  TH1* hist_pu= (TH1D*)MC_miss[imc][ity][ivar]->Clone(); hist_pu->Reset();
	  TH1* hist_st= (TH1D*)MC_miss[imc][ity][ivar]->Clone(); hist_st->Reset();
	  TH2* RMcopy = (TH2D*)RM_RecoGen->Clone(); RMcopy->RebinX(2);
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
	  // hist_pu->Write();
	  //    TH2* RMcopy = (TH2D*)RM_RecoGen->Clone(); RMcopy->RebinX(2);
	  TH1* RMxcopy= (TH1D*)RMx->Clone(); RMxcopy->Rebin(2);
	  //	for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_pu->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
	  for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_st->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
	  // hist_pu->Divide(hist_pu,RMxcopy, 1, 1, "b");
	  hist_st->Divide(hist_st,RMy, 1, 1, "b");
	  
	  
	  //       TH1* RMMx= (TH1D*)RMx->Clone(); RMMx->Rebin(2);
	  //      hist_pu->Divide(RMy,RMMx, 1, 1, "b");
	  
	  TH1D *h_pu = (TH1D*)fakerate->Clone(); h_pu->Reset();
	  TH1D *h_st = (TH1D*)fakerate->Clone(); h_st->Reset();
	  int ir = h_pu->GetNbinsX(); int ig = missrate->GetNbinsX();
	  double fk[ir]; double ef[ig]; double pu[ir]; double st[ir];
	  subtract_background(RM_RecoGen, RecoMC, GenMC, DataReco, fk, ef, pu, st);
	  
	  for(int i =1; i<h_pu->GetNbinsX()+1; i++){
	    h_pu->SetBinContent(i,pu[i]);
	    h_st->SetBinContent(i,st[i]);
	  }
	  
	  hist_purity[imc][ity][ivar]=(TH1D*)hist_pu->Clone(); 
	  //hist_purity[imc][ity][ivar][ipt]=(TH1D*)h_pu->Clone(); 
	  //hist_stbl[imc][ity][ivar][ipt]=(TH1D*)h_st->Clone(); 
	  hist_stbl[imc][ity][ivar]=(TH1D*)hist_st->Clone(); 
	  hist_purity[imc][ity][ivar]->SetMinimum(-0.05); hist_purity[imc][ity][ivar]->SetMaximum(1.01);
	  hist_stbl[imc][ity][ivar]->SetMinimum(-0.05); hist_stbl[imc][ity][ivar]->SetMaximum(1.01);

	  sprintf(name,"%sPurity_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	  hist_purity[imc][ity][ivar]->SetNameTitle(name,name);
	  sprintf(name,"%sstability_%i_eta0_%i", histtag[idd], ity, var[ivar]);
	  hist_stbl[imc][ity][ivar]->SetNameTitle(name,name);
	  
	  hist_purity[imc][ity][ivar]->Write();
	  hist_stbl[imc][ity][ivar]->Write();
	} //for (int imc=0;imc<nmc;imc++)
	
        DirData[idd]->cd(); 
        
	sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	sprintf(name, "E%sreco_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	sprintf(Axisname,"var_%i[UO]",var[ivar]);
	TH1* DataRecoExtact = binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name, DataReco, 0, true, Axisname);
	DataRecoExtact->Write();
	DataReco->Write(); 
	Data_Reco[ity][ivar] = (TH1D*)DataReco->Clone();
#ifdef CLOUSER
	sprintf(histname,"Genbin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	sprintf(name,"E%sgen_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	sprintf(Axisname,"var_%i[UO]",var[ivar]);
	TH1* GenExtactPsudo =binsGen[ity][ivar]->FindNode(histname)->ExtractHistogram(name, DataGen , 0, true, Axisname); 
	GenExtactPsudo->Write();
	DataGen->Write();
	
	TH1* RecoFakecorrExD; TH1* GenMisscorrExD;
	sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	sprintf(name,"E%s", MC_reco_fake[umc][ity][ivar]->GetName());  sprintf(Axisname,"var_%i[UO]",var[ivar]);
	RecoFakecorrExD =binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name, Data_reco_fake[ity][ivar], 0, true, Axisname);
	
	sprintf(histname,"Genbin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	sprintf(name,"E%s", MC_Gen_miss[umc][ity][ivar]->GetName());   sprintf(Axisname,"var_%i[UO]",var[ivar]);
	GenMisscorrExD =binsGen[ity][ivar]->FindNode(histname)->ExtractHistogram(name, Data_gen_miss[ity][ivar], 0, true, Axisname);
	
	GenMisscorrExD->Write(); RecoFakecorrExD->Write();
	
	Data_reco_fake[ity][ivar]->Write();     Data_gen_miss[ity][ivar]->Write();
#endif
      }//for(int ivar=0; ivar < nusedvar ; ivar ++)
    }//for(int ity=0; ity <type; ity++)
    cout << "Histogram Read Data and MC Done " <<endl;
    
    //------------------Fold check : Patrick 1 Sep 20
    //Get Probability Matrix  (gen-miss)*probability = (reco-fake). Of course, don't forget to account for miss and fake entries (if applicable).
    folddir[idd]->cd();
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
	
	if(idd>=1){
	  sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	  sprintf(name, "E%sFold_%i_eta0_%i",histtag[idd], ity, var[ivar]);   sprintf(Axisname,"var_%i[UO]",var[ivar]);
	  TH1* TrueFold = binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name, Folded, 0, true, Axisname);   TrueFold->Write();
	}
	sprintf(name,"%sFold_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	Folded->SetNameTitle(name,name);
	
	Folded->Write();
      }
    }
    //----------------------------Condition  number of Probability Matrix
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	TH2D* RM  = (TH2D*)h2dGenDetMC[umc][ity][ivar]->Clone();
	TH1D* miss = (TH1D*)MC_miss[umc][ity][ivar]->Clone();
	//TH1D* miss = (TH1D*)MC_missrate[umc][ity][ivar][ipt]->Clone();
	RM->RebinY(irbin); miss->Rebin(irbin);
	cout <<setw(2) << ity <<setw(5) <<ivar << setw(5) <<'\n';
	Condition(RM, miss);
      }
    }
    
    Unfolddir[idd]->cd();
    for(int ity=0; ity <type; ity++){
      for(int ivar=0; ivar < nusedvar ; ivar ++){
	if (idd==0){cout << "Root Hist : "; }
	if (idd==1){cout << "TUnfoldBinning1D : "; }
	if (idd==2){cout << " TUnfoldBinning2D : "; }
	cout <<"typ "<< ity << " : Var " << ivar << " ";
        //file <<"["<< ity << "," << ivar << "," << ipt <<"] --->" << endl;
	
        
	//Get reco bins
	double rxbins[MC_Reco[umc][ity][ivar]->GetNbinsX()+1]={0};
	for(int ix=0; ix<MC_Reco[umc][ity][ivar]->GetNbinsX()+1; ix++) {
	  rxbins[ix] = MC_Reco[umc][ity][ivar]->GetXaxis()->GetBinLowEdge(ix+1); }
	
	//Get Gen bins
	double gxbins[MC_Gen[umc][ity][ivar]->GetNbinsX()+1]={0};
	for (int ix=0; ix<MC_Gen[umc][ity][ivar]->GetNbinsX()+1; ix++) {
	  gxbins[ix] = MC_Gen[umc][ity][ivar]->GetXaxis()->GetBinLowEdge(ix+1); }
	
	//Rebin for match the condition of reco vs gen bin 
	h2dGenDetMC[umc][ity][ivar]->RebinY(irbin);
	MC_Gen[umc][ity][ivar]->Rebin(irbin);
	//----------------------------------------------------------Define Input---------------------	
        TUnfoldBinning* RecoBin = binsRec[ity][ivar];
        TUnfoldBinning* GenBin = binsGen[ity][ivar];
	
	
	TH2* RMin = (TH2D*)h2dGenDetMC[umc][ity][ivar]->Clone();
	TH1* input = (TH1D*)Data_Reco[ity][ivar]->Clone();
	TH1* mcgen = (TH1D*)MC_Gen[umc][ity][ivar]->Clone();
	TH1* mcreco = (TH1D*)MC_Reco[umc][ity][ivar]->Clone();
	TH1* mcgen_miss = (TH1D*)MC_Gen_miss[umc][ity][ivar]->Clone();
	TH1* mc_miss = (TH1D*)MC_miss[umc][ity][ivar]->Clone();
	TH1* mc_missrate = (TH1D*)MC_missrate[umc][ity][ivar]->Clone();
	
	//correction for Fake as background subtraction : Patrick
        TH1* mcbackground = (TH1D*)MC_fakerate[umc][ity][ivar]->Clone();
	mcbackground->Multiply(input);
	
        double biasScale = 0;
        const char *REGULARISATION_DISTRIBUTION=0;
        const char *REGULARISATION_AXISSTEERING="*[UO]";
        
        //https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
        sprintf(name,"%sData_covariance_%i_eta0_%i",histtag[idd], ity, var[ivar]);
        TH2D* covM = RecoBin->CreateErrorMatrixHistogram(name,false); covM->Sumw2();
	for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
	  double err = input->GetBinError(ix);
	  covM->SetBinContent(ix,ix,err*err);
	}
	/*       for (int ix=1; ix<covM->GetNbinsX()+1; ix++) {
		 for (int iy=1; iy<covM->GetNbinsY()+1; iy++) {
		 double cTmp = input->GetBinContent(ix)*input->GetBinContent(iy);
		 if(ix==iy)covM->SetBinContent(ix,iy,cTmp);
		 }
		 }*/
      
	covM->Write();
	
	//TH1 * BLTdata = (TH1*)Data_reco_fake[ity][ivar]->Clone();//data - back ground
	TH1 * BLTdata = (TH1*)Data_fakerateInv[ity][ivar]->Clone(); BLTdata->Multiply(input);
	
	//TH1 * BLTgen= (TH1*)Data_gen_miss[ity][ivar]->Clone();
	TH1 * BLTgen= (TH1*)Data_missrateInv[ity][ivar]->Clone();  BLTgen->Multiply(Data_Gen[ity][ivar]);
	sprintf(name,"%sBLT_Data_covariance_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	TH2D* covMBLT= RecoBin->CreateErrorMatrixHistogram(name,false); covMBLT->Sumw2();
	for (int ix=1; ix<BLTdata->GetNbinsX()+1; ix++) {
	  double err = BLTdata->GetBinError(ix);
	  covMBLT->SetBinContent(ix,ix,err*err);
       }
	//BLT------------------
	//----------------------------------------------------------No Reguratization-------------------------------------	
        TUnfoldDensity density(RMin,TUnfold::kHistMapOutputVert,
			       TUnfoldDensity::kRegModeNone,
			       TUnfoldDensity::kEConstraintNone,
			       TUnfoldDensity::kDensityModeNone,
			       GenBin,RecoBin, REGULARISATION_DISTRIBUTION, REGULARISATION_AXISSTEERING);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);
	
        density.SubtractBackground(mcbackground, "fake", 1.0, 0.00); // hist,name,scale, scale error 
        int status = density.SetInput(input,biasScale,0,covM);
	
        int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
        cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;
	
        //but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
        //tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result
        density.DoUnfold(0.0);//,input,biasScale);  // tau 0.0 means no regularization

        sprintf(histname, "%sTUnfold_NoReg_typ_%i_eta0_%i", histtag[idd], ity, var[ivar]); //unfolded_typ_0_pt2_eta0_3
        sprintf(title, "Unfolded No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
        bool AxisBin =true;
	if(idd>=1){AxisBin = false;}//Keep or Not ?
	TH1 *Unfolded = density.GetOutput(histname,title,0,"*[UO]", AxisBin);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
        
	sprintf(histname, "%sInEmat_NoReg_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
        sprintf(title, "Input Ematrix No Reg %s 2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *hist_In_Emat = density.GetEmatrixInput(histname,title);//,"signal");
	
	sprintf(histname, "%scorr_NoReg_typ_%i_eta0_%i", histtag[idd], ity, var[ivar]);
        sprintf(title, "RhoIJtotal No Reg %s  2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *Hist_RhoIJ = density.GetRhoIJtotal(histname, title);//,"signal");
	
        sprintf(histname, "%sEmat_NoReg_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
        sprintf(title, "Ematrix No Reg %s 2.4 %s ", itypeN[ity],  vartitle[var[ivar]]);
        TH2 *hist_Emat = density.GetEmatrixTotal(histname,title);//,"signal");
	
        sprintf(histname, "%sRefold_NoReg_typ_%i_eta0_%i", histtag[idd],ity, var[ivar]);
        sprintf(title, "Back folded No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
        TH1 *hist_folded = density.GetFoldedOutput(histname,title,0,"*[UO]",AxisBin);//,"signal");
	
        sprintf(histname, "%sProb_NoReg_typ_%i_eta0_%i",histtag[idd], ity, var[ivar]);
        sprintf(title, "ProbabilityMatrix No Reg %s 2.4 %s ", itypeN[ity], vartitle[var[ivar]]);
        TH2 *hist_prob = density.GetProbabilityMatrix(histname,title,TUnfold::kHistMapOutputVert);//,"signal");
	
	
	//Calculate Unfolded covarince 
	//https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
	sprintf(name,"%sMC_covariance_%i_eta0_%i",histtag[idd], ity, var[ivar]);
	TH2D* covUf = GenBin->CreateErrorMatrixHistogram(name,false); covUf->Sumw2();
	for (int ix=1; ix<covUf->GetNbinsX()+1; ix++) {
	  //double errmc = mcgen->GetBinError(ix);
	  double err = Unfolded->GetBinError(ix);
	  covUf->SetBinContent(ix,ix,err*err);
	  //covUf->SetBinContent(ix,ix,errmc*errmc);
	  //covUf->SetBinContent(ix,ix,err*err+errmc*errmc);
	}
	//covMC->Add(hist_Emat);
        TH2 *hist_Emat_Norm = (TH2*)hist_Emat->Clone();
        cout << " From Density Chi2A() : " << density.GetChi2A() << " " << density.GetChi2L()  << " " << density.GetNdf() <<" Chi2/NDf : "  <<density.GetChi2A()/(density.GetNdf()+1) <<endl;
        cout << " Smeared detector level 1: ";BLT(BLTdata, covM, hist_folded, 1);
        cout << " Smeared detector level 2: ";BLT(BLTdata, covMBLT, hist_folded, 1) ; cout<< endl;
	
	cout << " Unfold before Miss correction 1 : "; BLT(Unfolded, covUf, BLTgen);
	cout << " Unfold before Miss 2: "; BLT(Unfolded, hist_Emat, BLTgen);
	cout << " Unfold before Miss 3: "; BLT(Unfolded, hist_In_Emat, BLTgen) ;cout << endl;
	
	// correction for Miss entries  : Partick
        for (int i = 1; i <= Unfolded->GetNbinsX(); ++i) {
          double content = Unfolded->GetBinContent(i);
          double factor = 1;
          //factor += (mc_missrate->GetBinContent(i)/(1-mc_missrate->GetBinContent(i)));
          factor += (mc_miss->GetBinContent(i)/(mcgen->GetBinContent(i) - mc_miss->GetBinContent(i)));
          content *= factor;
          Unfolded->SetBinContent(i, content);  }
          //for (int i = 1; i <= Unfolded->GetNbinsX(); ++i) {cout<<" "<<Unfolded->GetBinContent(i)/mcgen->GetBinContent(i); }
	
	cout << " Unfold after Miss 3: "; BLT(Unfolded, covUf,Data_Gen[ity][ivar]);
	cout << " Unfold after Miss 3: "; BLT(Unfolded, hist_Emat,Data_Gen[ity][ivar]);
	cout << " Unfold after Miss 3: "; BLT(Unfolded, hist_In_Emat,Data_Gen[ity][ivar]);
	
	Unfolded->Write(); Hist_RhoIJ->Write(); hist_Emat->Write(); hist_folded->Write(); hist_prob->Write(); hist_In_Emat->Write();
        ConditionV2(hist_prob); 
	
	
	//----------------------Extract The True Distrubuntions 
	if(idd>=1){
	  sprintf(histname,"Recobin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	  sprintf(name, "E%s",hist_folded->GetName()); sprintf(Axisname,"var_%i[UO]",var[ivar]);
	  TH1* foldedExtact =binsRec[ity][ivar]->FindNode(histname)->ExtractHistogram(name, hist_folded, 0, true, Axisname);
	  
	  sprintf(histname,"Genbin%s_typ_%i_eta0_%i", bintag[idd], ity, var[ivar]);
	  sprintf(name,"E%s",Unfolded->GetName()); sprintf(Axisname,"var_%i[UO]",var[ivar]);
	  TH1* UnfoldedExtact =binsGen[ity][ivar]->FindNode(histname)->ExtractHistogram(name, Unfolded , 0, true, Axisname);
          
	  foldedExtact->Write(); UnfoldedExtact->Write();
	}//if(idd>=1)
	//-----------------------------------------End of Variables loop
	
	cout << endl ; 
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

