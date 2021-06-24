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

#include "CMS_lumi.C"

using namespace std;

void unc_plot(){
  bool percen=false;
  int const njer = 2;
  int const nmc=4; //Number of MC -> 0,1,2 PY8, MG, HW7
  int const umc=0; //Which MC will used for Unfold

  int irbin = 1; //Rebin
  int const untype = 1;  // 0 for 2D , 1 for 2D

  const int nHLTmx=8; //HT2 Range
  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =2;
  Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass
  static const int nhist=10; //We need 8 But define 10
  const int njetptmn = nHLTmx;
  const  char* Unit[4]={"GeV","fb^{-1}","Pb","#%"};
  const  char* DataEra[3]={"Data","Data","Data"};
  const  char* RunEra[3]={"2016","2017","2018"};
  static const int iera = 2;
  int iPeriod = 0;  int iPos=0 ;


  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={2,4,6,28,66,46,1,28,38,42};  // define the color for different histograms

  //Int_t HT2range[nHLTmx+1]={84, 111, 172, 241, 309, 376, 462, 568, 3000}; // 2017 1June21
  //Int_t HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000}; // 2017
  Int_t HT2range[nHLTmx+1]={83, 109, 176, 247, 318, 387, 477, 573, 3000}; //2018

  const  char* itypeN[ntype]={"Jets","Charged Particles"};
  //const  char* unc_name[8]={"JES Unc.","JER Unc.","PDF Unc.", "PU Unc.", "Model Unc.","Track Unc.", "Stat Unc. ","Unf Unc."};
  const  char* unc_name[8]={"JES ","JER ","PDF ", "PU ", "Model ","Track ", "Stat ","Unf "};
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[8]={"83 < H_{T,2} < 109", "109 < H_{T,2} < 172", "172 < H_{T,2} < 241", "241 < H_{T,2} < 309","309 < H_{T,2} < 377", "377 < H_{T,2} < 462", "462 < H_{T,2} <570","H_{T,2} > 570"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};


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
  TLegend* CTLegendV3(float x1, float y1, float x2, float y2, float txtsize );
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl);
  void MyplotsetV2(TH1D *MyHist, const char* XTitle, const char* YTitle);
TCanvas *GridPlot(int nplot, TH1D *Hist[nplot] ,const char* Legend[nplot]);
TCanvas *GridPlotV2(int nplot[2],     vector<vector<TH1D*>> Hist ,const char* Legend[nplot[2]], const char* unc_leg[nplot[0]] );

  

  TH1D *JEC_Unc[ntype][nusedvar][nHLTmx];
  TH1D *JER_Unc[ntype][nusedvar][nHLTmx];
  TH1D *PDF_Unc[ntype][nusedvar][nHLTmx];
  TH1D *stat_Unc[ntype][nusedvar][nHLTmx];
  TH1D *PU_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Unfold_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Unfstat_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Track_Unc[ntype][nusedvar][nHLTmx];
  TH1D *Total_Unc[ntype][nusedvar][nHLTmx];

  TFile *JEC_root = TFile::Open("JEC_Result.root");
  TFile *JER_root = TFile::Open("JER_Result.root");
  TFile *PDF_root = TFile::Open("PDF_Result.root");
  TFile *stat_root = TFile::Open("stat_Result.root");
  TFile *PU_root = TFile::Open("PU_Result.root");
  TFile *unf_root = TFile::Open("Unf_Result.root");
  TFile *unfstat_root = TFile::Open("Unf_stat_Result.root");
  TFile *track_root = TFile::Open("Track_Result.root");

 TFile *Total_unc =new TFile("Total_unc.root","recreate");

  for(int ity=0; ity <ntype; ity++){
      for(int ivar =0;ivar < nusedvar; ivar++){
          for(int ipt=0; ipt < nHLTmx; ipt++){
          JEC_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("jec_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),JEC_root);
          JER_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("jer_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),JER_root);
          PDF_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("PDF_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),PDF_root);
          //stat_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("stat_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),stat_root);  //Direct from unfold hist
          stat_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("Inputstat_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),unfstat_root);  //only from data input
          PU_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("pu_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), PU_root);
          Unfold_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("unfold_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unf_root);
          Unfstat_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("Unfstat_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), unfstat_root);  //only from MC RM
          Track_Unc[ity][ivar][ipt] = (TH1D*)ReadHist1D("track_erro_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt), track_root);
 	  
	  //Total
	  Total_Unc[ity][ivar][ipt]= (TH1D*)JEC_Unc[ity][ivar][ipt]->Clone();
          Total_Unc[ity][ivar][ipt]->Reset();
          Total_Unc[ity][ivar][ipt]->SetNameTitle(Form("total_erro_%i_eta0_%i_pt%i", ity, var[ivar], ipt),Form("Total Relative %i_eta0_%i_pt%i", ity, var[ivar], ipt));
	  
	  }
       }
    }
//------------------------------------------------------------Total
for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt=0; ipt < nHLTmx; ipt++){
        int nbins = JEC_Unc[ity][ivar][ipt]->GetNbinsX();

	cout <<" " <<ity << " " <<ivar << " " <<ipt <<endl;
        for (int ix=1; ix<nbins+1; ix++) {


        double  jec =  JEC_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  jer =  JER_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  pdf =  PDF_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  stat =  stat_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  pileup =  PU_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  unf =  Unfold_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  Unfstat =  Unfstat_Unc[ity][ivar][ipt]->GetBinContent(ix);
        double  track =  Track_Unc[ity][ivar][ipt]->GetBinContent(ix);

        double toterr = sqrt( jec*jec + jer*jer + pdf*pdf + stat*stat + pileup*pileup + unf*unf+ Unfstat*Unfstat );//+ track*track); 
	cout << jec <<" "<< jer << " " << pdf << " " << stat <<" " <<pileup << " " << unf << " " << track << " " << toterr <<" " << endl;  
	Total_Unc[ity][ivar][ipt]->SetBinContent(ix, toterr);
          }
	
	Total_Unc[ity][ivar][ipt]->Write();;
       }
       
    }
}



for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt=0; ipt < nHLTmx; ipt++){

        int nbins = JEC_Unc[ity][ivar][ipt]->GetNbinsX();
        for (int ix=0; ix<nbins+2; ix++) {
        if(JEC_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "JEC" << JEC_Unc[ity][ivar][ipt]->GetBinContent(ix) <<endl ; JEC_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(JER_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "JER";JER_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(PDF_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "PDF";PDF_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(stat_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "STAT";stat_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(PU_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "PU";PU_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(Unfold_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "Unf";Unfold_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(Unfstat_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "corr";Unfstat_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(Track_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "Track";Track_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
        if(Total_Unc[ity][ivar][ipt]->GetBinContent(ix)<=0){cout << "Total "; Total_Unc[ity][ivar][ipt]->SetBinContent(ix,1.e-6);}
	
	}
      
       }

    }
}
if(percen){
for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt=0; ipt < nHLTmx; ipt++){

        int nbins = JEC_Unc[ity][ivar][ipt]->GetNbinsX();
        for (int ix=0; ix<nbins+2; ix++) {
        JEC_Unc[ity][ivar][ipt]->SetBinContent(ix,JEC_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        JER_Unc[ity][ivar][ipt]->SetBinContent(ix, JER_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        PDF_Unc[ity][ivar][ipt]->SetBinContent(ix,PDF_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        stat_Unc[ity][ivar][ipt]->SetBinContent(ix,stat_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        PU_Unc[ity][ivar][ipt]->SetBinContent(ix,PU_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        Unfold_Unc[ity][ivar][ipt]->SetBinContent(ix,Unfold_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        Unfstat_Unc[ity][ivar][ipt]->SetBinContent(ix,Unfstat_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        Track_Unc[ity][ivar][ipt]->SetBinContent(ix,Track_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);
        Total_Unc[ity][ivar][ipt]->SetBinContent(ix,Total_Unc[ity][ivar][ipt]->GetBinContent(ix)*100);

        }

       }

    }
}
}





//cout << "ok" <<endl; 
  //---------------------------------------------------------- ALL
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
    TCanvas *cpgrid = new TCanvas("cpgrid", "cpgrid", 1000,600 );    SetMycanvas(cpgrid,0,0.0,0.00,0.00,0.0);
     cpgrid->SetGridy();
     cpgrid->Divide(4,2);
    gStyle->SetPadGridY(3);
    for(int ipt=0; ipt < nHLTmx; ipt++){
     TLegend *leg1;
    TCanvas *cptr = new TCanvas("cptr", "cptr", 500,400 );    SetMycanvas(cptr,0,0.15,0.02,0.06,0.17);
    // TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "Sytematic Uncertainty", itypeN[ity], htrang[ipt]);
     if(ipt==7){leg1= CTLegendV2(0.1,0.70,0.5,0.90,0.053, Form(" H_{T,2} > %i %s",  HT2range[ipt],Unit[0]));}else
     {leg1= CTLegendV2(0.1,0.70,0.5,0.90,0.053, Form(" %i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
     //TLegend *leg1= CTLegendV2(0.1,0.70,0.5,0.90,0.06, htrang[ipt]);
     TLegend *leg2= CTLegendV3(0.6,0.55,0.8,.93,0.044);
     
     char Drawm[100];
     if(ipt==0 && ipt==4){sprintf(Drawm,"SAME");} else{sprintf(Drawm,"A SAME");};
             

//     gPad->SetLogy();
     cptr->SetGridy(); cptr->cd();
//cout << "Ok 0" << endl;
     TH1D *jec_unc = (TH1D*)JEC_Unc[ity][ivar][ipt]->Clone();
     jec_unc->SetMinimum(1.e-6); jec_unc->SetMaximum(0.30);
     Myplotset(jec_unc, Esvlogx[ivar] ,"Relative Unc.");
     leg2->AddEntry(jec_unc,unc_name[0],"lp");
     jec_unc->SetLineColor(color[0]);
     cptr->cd();
     jec_unc->Draw();

     cpgrid->cd(ipt+1);
     jec_unc->Draw(Drawm);
//cout << "Ok 1" << endl;
     cptr->SetGridy(); cptr->cd();
     TH1D *jer_unc = (TH1D*)JER_Unc[ity][ivar][ipt]->Clone();
     jer_unc->SetMinimum(1.e-6); jer_unc->SetMaximum(0.30);
     Myplotset(jer_unc, Esvlogx[ivar] ,"Relative Unc.");
     jer_unc->SetLineColor(color[1]);
     leg2->AddEntry(jer_unc,unc_name[1],"lp");
     jer_unc->Draw("  SAME");
//cout << "Ok 2" << endl;
     
     cpgrid->cd(ipt+1);
     jer_unc->Draw(Drawm);

     cptr->SetGridy(); cptr->cd();
     TH1D *pdf_unc = (TH1D*)PDF_Unc[ity][ivar][ipt]->Clone();
     pdf_unc->SetMinimum(1.e-6); pdf_unc->SetMaximum(0.30);
     Myplotset(pdf_unc, Esvlogx[ivar] ,"Relative Unc.");
     pdf_unc->SetLineColor(color[2]);
     leg2->AddEntry(pdf_unc,unc_name[2],"lp");
     pdf_unc->Draw("SAME");
//cout << "Ok 3" << endl;

     cpgrid->cd(ipt+1);
     pdf_unc->Draw(Drawm);

     cptr->SetGridy(); cptr->cd();
     TH1D *pu_unc = (TH1D*)PU_Unc[ity][ivar][ipt]->Clone();
     pu_unc->SetMinimum(1.e-6); pu_unc->SetMaximum(0.30);
     Myplotset(pu_unc, Esvlogx[ivar] ,"Relative Unc.");
     pu_unc->SetLineColor(color[3]);
     leg2->AddEntry(pu_unc,unc_name[3],"lp");
     pu_unc->Draw("SAME");
//cout << "Ok 4" << endl;

     cpgrid->cd(ipt+1);
     pu_unc->Draw(Drawm);

     cptr->SetGridy(); cptr->cd();
     TH1D *unf_unc = (TH1D*)Unfold_Unc[ity][ivar][ipt]->Clone();
     unf_unc->SetMinimum(1.e-6); unf_unc->SetMaximum(0.30);
     Myplotset(unf_unc, Esvlogx[ivar] ,"Relative Unc.");
     unf_unc->SetLineColor(color[4]);
     leg2->AddEntry(unf_unc,unc_name[4],"lp");
     unf_unc->Draw("SAME");
//cout << "Ok 5" << endl;
     cpgrid->cd(ipt+1);
     unf_unc->Draw(Drawm);
/*

     cptr->SetGridy(); cptr->cd();
     TH1D *track_unc = (TH1D*)Track_Unc[ity][ivar][ipt]->Clone();
     track_unc->SetMinimum(-0.02); track_unc->SetMaximum(0.30);
     Myplotset(track_unc, Esvlogx[ivar] ,"Relative Unc.");
     track_unc->SetLineColor(color[5]);
     leg2->AddEntry(track_unc,unc_name[5],"lp");
     track_unc->Draw("SAME");
*/

     cptr->SetGridy(); cptr->cd();
     TH1D *stat_unc = (TH1D*)stat_Unc[ity][ivar][ipt]->Clone();
     stat_unc->SetMinimum(1.e-6); stat_unc->SetMaximum(0.30);
     Myplotset(stat_unc, Esvlogx[ivar] ,"Relative Unc.");
     stat_unc->SetLineColor(color[6]);
     leg2->AddEntry(stat_unc,unc_name[6],"lp");
     stat_unc->Draw(" SAME");
//cout << "Ok 6" << endl;

     cpgrid->cd(ipt+1);
     stat_unc->Draw(Drawm);
     
     cptr->SetGridy(); cptr->cd();
     TH1D *unfstat_unc = (TH1D*)Unfstat_Unc[ity][ivar][ipt]->Clone();
     unfstat_unc->SetMinimum(1.e-1); unfstat_unc->SetMaximum(0.30);
     Myplotset(unfstat_unc, Esvlogx[ivar] ,"Relative Unc.");
     unfstat_unc->SetLineColor(color[7]);
     leg2->AddEntry(unfstat_unc,unc_name[7],"lp");
     unfstat_unc->Draw("SAME");
//cout << "Ok 7" << endl;


     cpgrid->cd(ipt+1);
     unfstat_unc->Draw(Drawm);
//cout << "Ok 7" << endl;
     
     cptr->SetGridy(); cptr->cd();
     TH1D *tot_unc = (TH1D*)Total_Unc[ity][ivar][ipt]->Clone();
     tot_unc->SetMinimum(1.e-6); //tot_unc->SetMaximum(0.30);
     Myplotset(tot_unc, Esvlogx[ivar] ,"Relative Unc.");
     tot_unc->SetLineColor(1);
     tot_unc->SetLineStyle(1);
     leg2->AddEntry(tot_unc,"Total Unc.","lp");
     tot_unc->Draw("SAME");
//cout << "Ok 8" << endl;

     cpgrid->cd(ipt+1);
  //   tot_unc->GetXaxis()->SetLabelOffset(999);
//     tot_unc->GetXaxis()->SetLabelSize(0);
     tot_unc->Draw(Drawm);

     
     //     cptr->SetLogy();
 //   gPad->SetLogy();
     cptr->SetGridy(); cptr->cd();
     leg2->Draw();leg1->Draw();
     CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "Unc_plot_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "Unc_plot_%s.pdf",RunEra[iera]);sprintf(pdfname2, "Unc_plot_%s.pdf)",RunEra[iera]);
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
//     cptr->Clear();
        }

      }

   }

 //---------------------------------------------------------- Total
   const  char* unc_source[8]={"JES ","JER ", "PU","Unf","Model", "PDF","Stat", "Total "};

  for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
    TCanvas *cpgrid = new TCanvas("cpgrid", "cpgrid", 1000,600 ); 
    vector<vector<TH1D*>> Hist( 8 , vector<TH1D*> (8));
    int nplot[2]={8,8};// 8 number of unc , 8 number of Ht range
    const char* Htrange[nHLTmx];
    for(int ipt=0; ipt < nHLTmx; ipt++){
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "JEC Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");

     if(ipt==7){Htrange[ipt]= Form(" H_{T,2} > %i",  HT2range[ipt]);}else{ Htrange[ipt]= Form(" %i < H_{T,2} < %i ", HT2range[ipt] , HT2range[ipt+1]);}
     Hist[0][ipt] = (TH1D*)JEC_Unc[ity][ivar][ipt]->Clone();
     Hist[1][ipt] = (TH1D*)JER_Unc[ity][ivar][ipt]->Clone();
     Hist[2][ipt] = (TH1D*)PU_Unc[ity][ivar][ipt]->Clone();
     Hist[3][ipt] = (TH1D*)Unfstat_Unc[ity][ivar][ipt]->Clone();
     Hist[4][ipt] = (TH1D*)Unfold_Unc[ity][ivar][ipt]->Clone();
     Hist[5][ipt] = (TH1D*)PDF_Unc[ity][ivar][ipt]->Clone();
     Hist[6][ipt] = (TH1D*)stat_Unc[ity][ivar][ipt]->Clone();
     Hist[7][ipt] = (TH1D*)Total_Unc[ity][ivar][ipt]->Clone();
        }

        for(int jj=0; jj < nHLTmx; jj++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
         Myplotset( Hist[jj][ipt], Esvlogx[ivar] ,"Relative Unc");
          Hist[jj][ipt]->SetMinimum(0.0001); Hist[jj][ipt]->SetMaximum(.19);
//          if(percen)Hist[jj][ipt]->SetMinimum(0.0001); Hist[jj][ipt]->SetMaximum(100);
          }
	}
     cpgrid =(TCanvas*)(GridPlotV2(nplot, Hist,Htrange,unc_source));

     SetMycanvas(cpgrid,0,0.1,0.02,0.05,0.12);
     CMS_lumi( cpgrid, iPeriod, iPos ); cpgrid->Update();
     sprintf(pdfname, "Total_unc_plot_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "Total_unc_plot_%s.pdf",RunEra[iera]);sprintf(pdfname2, "Total_unc_plot_%s.pdf)",RunEra[iera]);
     if(ity==0  && ivar==0){cpgrid->Print(pdfname,"pdf");}else if(ity==1 && ivar==4) {cpgrid->Print(pdfname2,"pdf"); }else{cpgrid->Print(pdfname,"pdf");};
      }
   }

 //---------------------------------------------------------- JEC
  for(int ity=0; ity <ntype; ity++){
      for(int ivar =0 ; ivar < nusedvar ; ivar++){
    TCanvas *cpgrid = new TCanvas("cpgrid", "cpgrid", 1000,600 );    SetMycanvas(cpgrid,0,0.0,0.00,0.00,0.0);
    TH1D *Hist[nHLTmx];
    const char* Htrange[nHLTmx];
    for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "JEC Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *jec_unc = (TH1D*)JEC_Unc[ity][ivar][ipt]->Clone();
     jec_unc->SetMinimum(0.0); jec_unc->SetMaximum(0.20);
     Myplotset(jec_unc, Esvlogx[ivar] ,"Relative Unc");
     jec_unc->Draw();

     leg2->Draw();leg1->Draw();

     if(ipt==7){Htrange[ipt]= Form(" H_{T,2} > %i",  HT2range[ipt]);}else{ Htrange[ipt]= Form(" %i < H_{T,2} < %i ", HT2range[ipt] , HT2range[ipt+1]);}


     cpgrid->cd(ipt+1);
     jec_unc->Draw();
     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "JEC_plot.pdf("); sprintf(pdfname1, "JEC_plot.pdf");sprintf(pdfname2, "JEC_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
    // cptr->Clear();

    Hist[ipt]= (TH1D*)jec_unc->Clone();
        }
     cpgrid =(TCanvas*)(GridPlot(nHLTmx, Hist,Htrange));

     SetMycanvas(cpgrid,0,0.1,0.02,0.04,0.12);
     CMS_lumi( cpgrid, iPeriod, iPos ); cpgrid->Update();
     sprintf(pdfname, "JEC_plot_2018.pdf("); sprintf(pdfname1, "JEC_plot_2018.pdf");sprintf(pdfname2, "JEC_plot_2018.pdf)");
     if(ity==0  && ivar==0){cpgrid->Print(pdfname,"pdf");}else if(ity==1 && ivar==4) {cpgrid->Print(pdfname2,"pdf"); }else{cpgrid->Print(pdfname,"pdf");};
      }
   }


  //---------------------------------------------------------- JER
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "JER Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *jer_unc = (TH1D*)JER_Unc[ity][ivar][ipt]->Clone();
     jer_unc->SetMinimum(0.0); jer_unc->SetMaximum(0.30);
     Myplotset(jer_unc, Esvlogx[ivar] ,"Relative Unc");
     jer_unc->Draw();

     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "JER_plot.pdf("); sprintf(pdfname1, "JER_plot.pdf");sprintf(pdfname2, "JER_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
  //   cptr->Clear();
        }
      }
   }


  //---------------------------------------------------------- PDF
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "PDF Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *pdf_unc = (TH1D*)PDF_Unc[ity][ivar][ipt]->Clone();
     pdf_unc->SetMinimum(0.0); pdf_unc->SetMaximum(0.20);
     Myplotset(pdf_unc, Esvlogx[ivar] ,"Relative Unc");
     pdf_unc->Draw();

     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "PDF_plot.pdf("); sprintf(pdfname1, "PDF_plot.pdf");sprintf(pdfname2, "PDF_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
    // cptr->Clear();
        }
      }
   }

   //---------------------------------------------------------- Unfold
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "unf Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *unf_unc = (TH1D*)Unfold_Unc[ity][ivar][ipt]->Clone();
     unf_unc->SetMinimum(0.0); unf_unc->SetMaximum(0.20);
     Myplotset(unf_unc, Esvlogx[ivar] ,"Relative Unc");
     unf_unc->Draw();

     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "unf_plot.pdf("); sprintf(pdfname1, "unf_plot.pdf");sprintf(pdfname2, "unf_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
  //   cptr->Clear();
        }
      }
   }


//---------------------------------------------------------- Track unc
  for(int ity=0; ity <ntype; ity++){
    for(int ivar =0 ; ivar < nusedvar ; ivar++){
        for(int ipt=0; ipt < nHLTmx; ipt++){
     TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "Track Uncertainty", itypeN[ity], htrang[ipt]);
     TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     cptr->SetGridy(); cptr->cd();
     TH1D *track_unc = (TH1D*)Track_Unc[ity][ivar][ipt]->Clone();
     track_unc->SetMinimum(0.0); track_unc->SetMaximum(0.20);
     Myplotset(track_unc, Esvlogx[ivar] ,"Relative Unc");
     track_unc->Draw();

     leg2->Draw();leg1->Draw();
   //  CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     sprintf(pdfname, "track_plot.pdf("); sprintf(pdfname1, "track_plot.pdf");sprintf(pdfname2, "track_plot.pdf)");
     if(ity==0 && ivar==0 && ipt==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{cptr->Print(pdfname,"pdf");};
   //  cptr->Clear();
        }
      }
   }

delete Total_unc;

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

void Myplotset(TH1D *MyHist, const char* XTitle, const char* YTitle){
  int ifornt =42;
  MyHist->SetTitleOffset(0.4);
  MyHist->SetTitleFont(ifornt);
  MyHist->SetTitleSize(0.02);
  MyHist->SetStats(0);

  MyHist->GetXaxis()->SetLabelSize(0.05);
  MyHist->GetXaxis()->SetTitleSize(0.06);
  MyHist->GetXaxis()->SetTitleOffset(1.0);
  MyHist->GetXaxis()->SetTitleFont(ifornt);
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetXaxis()->SetTitle(XTitle);

  MyHist->GetYaxis()->SetLabelSize(0.05);
  MyHist->GetYaxis()->SetTitleSize(0.06);
  MyHist->GetYaxis()->SetTitleOffset(1.2);
  MyHist->GetYaxis()->SetTitle(YTitle);
  MyHist->GetYaxis()->SetTitleFont(ifornt);
  MyHist->GetYaxis()->CenterTitle();
  MyHist->SetTitle("");
  //gStyle->SetTitleFontSize(.08);
  MyHist->SetLineStyle(2);
  MyHist->SetLineWidth(4);
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

TLegend* CTLegendV3(float x1, float y1, float x2, float y2, float txtsize ){
  TLegend *leg = new TLegend(x1,y1,x2,y2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(txtsize);
  leg->SetTextFont(42);
  return leg;
}

void MyplotsetV2(TH1D *MyHist, const char* XTitle, const char* YTitle){
  int ifornt =42;
  MyHist->SetTitleOffset(0.0);
  MyHist->SetTitleFont(ifornt);
  MyHist->SetTitleSize(0.0);
  MyHist->SetStats(0);

  MyHist->GetXaxis()->SetLabelSize(0.05);
  MyHist->GetXaxis()->SetTitleSize(0.06);
  MyHist->GetXaxis()->SetTitleOffset(1.0);
  MyHist->GetXaxis()->SetTitleFont(ifornt);
  MyHist->GetXaxis()->CenterTitle();
  MyHist->GetXaxis()->SetTitle(XTitle);
 
  MyHist->GetYaxis()->SetLabelOffset(999);
  MyHist->GetYaxis()->SetLabelSize(0.0);
  MyHist->GetYaxis()->SetTitleSize(0.0);
  MyHist->GetYaxis()->SetTitleOffset(0);
  MyHist->GetYaxis()->SetTitle(YTitle);
  MyHist->GetYaxis()->SetTitleFont(ifornt);
  MyHist->GetYaxis()->CenterTitle();
  MyHist->SetTitle("");
  //gStyle->SetTitleFontSize(.08);
  MyHist->SetLineStyle(2);
  MyHist->SetLineWidth(4);
}
//------------------------------------------------------------------------
TCanvas *GridPlot(int nplot, TH1D *Hist[nplot] ,const char* Legend[nplot]){
//nplot Total Number of Float (Should Be in even number) 
//if odd Put a value of next even number
//Hist[nplot] Number of plot

TPad *padfun[nplot];
TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 1000,750 );
  canvas->cd();
  canvas->SetRightMargin(0.0);
  canvas->SetLeftMargin(0.0);
//  canvas->SetTopMargin(0.045);

float padx1[8]={0.0,.29,.53,0.77,0.0, .29,.53,.77};
float padx2[8]={.29,.53,0.77, 1.0,.29,.53,0.77,1.0};
float pady1[8]={0.55,.55,.55,.55, 0.0,0.0,0.0,0.0};
float pady2[8]={.96,.96,.96,.96, 0.54,.54,0.54,0.54};

float LeftM[8]={.22,.0,.0,.0, 0.22,.0,0.0,0.0};
//float RightM[8]={.95,.95,.95,.95, 0.5,.5,0.5,0.5};
float BottomM[8]={0,0.0,0,0, 0.25,.25,0.25,0.25};

  Hist[0]->GetYaxis()->SetTitleSize(.095);
  Hist[4]->GetYaxis()->SetTitleSize(.091);
  Hist[0]->GetYaxis()->SetLabelSize(0.07);
  Hist[4]->GetYaxis()->SetLabelSize(0.06);


for(int i =0 ; i < nplot ; i++){  //loop for ratio plot
padfun[i] = new TPad("padfun1", "padfun1", padx1[i], pady1[i], padx2[i], pady2[i]);
  TLegend *HT_range;
  if(i==0 || i==4){HT_range = new TLegend(0.25,.8,.5,.9);}else{HT_range = new TLegend(0.01,.8,.5,.9);}
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);
  HT_range->SetTextSize(0.09);


  padfun[i]->SetTopMargin(0);
  padfun[i]->SetBottomMargin(BottomM[i]);
  padfun[i]->SetLeftMargin(LeftM[i]);
  padfun[i]->SetRightMargin(0);
  padfun[i]->SetGridy(); // Horizontal grid
  padfun[i]->Draw();
  padfun[i]->cd();
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  //HT_range->AddEntry(Hist[i], Legend[i],"l");
  HT_range->AddEntry((TObject*)0, Legend[i], "");
//  gPad->SetLogy();
  Hist[i]->SetLineWidth(2);
  Hist[i]->SetLineStyle(2);
  Hist[i]->GetXaxis()->SetLabelSize(0.07);
  Hist[i]->GetXaxis()->SetTitleSize(.095);
  if(i!=7)Hist[i]->GetXaxis()->SetTitleSize(0);
  Hist[i]->Draw(" SAME ");
  HT_range->Draw(" SAME");
  canvas->cd();

}

  return canvas;
}







//------------------------------------------------------------------------
TCanvas *GridPlotV2(int nplot[2],     vector<vector<TH1D*>> Hist ,const char* Legend[nplot[2]], const char* unc_leg[nplot[0]] ){
//nplot Total Number of Float (Should Be in even number) 
//if odd Put a value of next even number
//Hist[nplot] Number of plot

int ifront =42;
Int_t color[10] ={2,4,6,28,66,46,30,1,38,42};  // define the color for different histograms
TPad *padfun[nplot[1]];
TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 1000,650 );
  canvas->cd();
  canvas->SetRightMargin(0.0);
  canvas->SetLeftMargin(0.0);
  canvas->SetTopMargin(0.045);

float padx1[8]={0.0,.29,.53,0.77,0.0, .29,.53,.77};
float padx2[8]={.29,.53,0.77, 1.0,.29,.53,0.77,1.0};
float pady1[8]={0.55,.55,.55,.55, 0.0,0.0,0.0,0.0};
float pady2[8]={.95,.95,.95,.95, 0.54,.54,0.54,0.54};

float LeftM[8]={.22,.0,.0,.0, 0.22,.0,0.0,0.0};
//float RightM[8]={.95,.95,.95,.95, 0.5,.5,0.5,0.5};
float BottomM[8]={0,0.0,0,0, 0.25,.25,0.25,0.25};



  TLegend *unc_label = new TLegend(0.52,0.6,0.85,0.98);
  unc_label->SetFillStyle(0);
  unc_label->SetBorderSize(0);
  unc_label->SetTextSize(0.09);
  unc_label->SetTextFont(42);

for(int i =0 ; i < nplot[1] ; i++){  //loop for ratio plot
  padfun[i] = new TPad("padfun1", "padfun1", padx1[i], pady1[i], padx2[i], pady2[i]);
  for(int j =0 ; j < nplot[0] ; j++){  //loop for ratio plot
  //TLegend *unc_label = CTLegendV3(CTLegendV30.5,.3,1.0,1.0,0.1);
  TLegend *HT_range;
  if(i==0 || i==4){HT_range = new TLegend(0.25,.8,.5,.9); HT_range->SetTextSize(0.082); }else{HT_range = new TLegend(0.001,.8,.5,.9);HT_range->SetTextSize(0.09);}
  HT_range->SetFillStyle(0);
  HT_range->SetBorderSize(0);



  Hist[j][0]->GetYaxis()->SetTitleSize(.095);
  Hist[j][4]->GetYaxis()->SetTitleSize(.085);
  Hist[j][0]->GetYaxis()->SetLabelSize(0.067);
  Hist[j][4]->GetYaxis()->SetLabelSize(0.063);
  Hist[j][4]->GetXaxis()->SetLabelSize(0.06);
  if(i==4) HT_range->SetTextSize(0.075);



  padfun[i]->SetTopMargin(0);
  padfun[i]->SetBottomMargin(BottomM[i]);
  padfun[i]->SetLeftMargin(LeftM[i]);
  padfun[i]->SetRightMargin(0);
  padfun[i]->SetGridy(1); // Horizontal grid
  padfun[i]->Draw();
  padfun[i]->cd();
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  //HT_range->AddEntry(Hist[i], Legend[i],"l");
  HT_range->AddEntry((TObject*)0, Legend[i], "");
  if(i==7)unc_label->AddEntry(Hist[j][i], unc_leg[j],"l");
//  gPad->SetLogy();
  Hist[j][i]->GetXaxis()->SetLabelSize(0.07);
  Hist[j][i]->SetMarkerSize(1);
  Hist[j][i]->SetMarkerStyle(1);
  Hist[j][i]->SetLineWidth(2);
  Hist[j][i]->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
  Hist[j][i]->GetYaxis()->SetNdivisions(6,5,0, kTRUE);
  Hist[j][i]->GetYaxis()->SetTickLength(0.08);
//  Hist[j][i]->GetYaxis()->SetTickLength(0.06);


  Hist[j][i]->SetLineColor(color[j]);
  if(j==7)Hist[j][i]->SetLineStyle(1);
  Hist[j][i]->GetXaxis()->SetTitleSize(.095);
  if(i!=7)Hist[j][i]->GetXaxis()->SetTitleSize(0);
  if(i!=7)Hist[j][i]->GetXaxis()->SetTitleFont(ifront);
 
 
  Hist[j][i]->Draw(" SAME ");
  if(j==7) HT_range->Draw(" SAME");
  if(i==7&& j==7)unc_label->Draw("SAME");
  canvas->cd();

   }
}

  return canvas;
}
