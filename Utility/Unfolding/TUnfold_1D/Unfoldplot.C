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
#include <TH2F.h>
#include <TTree.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"

void Unfoldplot(){

TH1F *MC_hist[4][2][5][8];
TH1F *datahist[5][2][5][8];

  char histname1[100];
  char histname2[100];
  char Title[100];
  char Xaxis[100];
  char Yaxis[100];
  char ratioXaxis[100];
  char ratioYaxis[100];
  
  bool Reco, Gen;
  Int_t color[10] ={1,2,4,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[5]={3,9,15,18,24};
  Int_t HT2range[9]={83, 109, 172, 241, 309, 377, 462, 570, 3000};

  const int nHLTmx=8; //HT2 Range
  const int njetetamn=1;  //eta value used 2.4
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  static const int ntype =3;  
  static const int unfold_ty =3;  
  static const int nmc =1;
  
  const  char* regN[4]={"Tunfolded_Noreg_typ_","Tunfolded_TikhLscan_typ_","Tunfolded_TikhSURE_typ_","Tunfolded_ScanTau_typ_"}
  const  char* mcname[3]={"Pythia8","Madgraph","Herwig++"}
  const  char* itypeN[ntype]={"Jets","Charged Particles"};
 
  //******************************************************************************
   TFile *Unfoldroot = TFile::Open("testunfold2c_unfolded.root");  // Unfolded data 
  
 for(int  imc =0; imc < nmc <<imc++){ 
  for(int ity=0; ity <ntype; ity++){
   for(int ivar =0 ; ivar < nusedvar ; ivar++){
     for(Int_t ipt =0; ipt < nHLTmx ; ipt++){      
	     sprintf(histname2, "RebinMC/rebin_gen_typ_0_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt6_eta0_15
	     MC_hist[imc][ity][ivar][ipt] = (TH1F*) Unfoldroot->Get(histname2);
	     cout << histname2 << endl;
	     MC_hist[imc][ity][ivar][ipt]->Scale(1/(MC_hist[imc][ity][ivar][ipt]->Integral()));

	     //devision of bin width
	     int tmpnbn = MC_hist[imc][ity][ivar][ipt]->GetNbinsX();
	     for (int ix=0; ix<tmpnbn; ix++) {
		     double awidth = MC_hist[imc][ity][ivar][ipt]->GetBinWidth(ix+1); // tmpwid;
		     MC_hist[imc][ity][ivar][ipt]->SetBinContent(ix+1, MC_hist[imc][ivar][ipt]->GetBinContent(ix+1)/awidth);
		     double error = MC_hist[imc][ity][ivar][ipt]->GetBinError(ix+1);
		     MC_hist[imc][ity][ivar][ipt]->SetBinError(ix+1, error/awidth);
	        } 
            }
          }
       }
    }  //end of one MCinput root file  reading
 

  for(int iun=0; iun < unfold_ty; iun++){
	  for(int ity=0; ity <type; ity++){
		  for(int ivar=0; ivar < nusedvar ; ivar ++){
			  for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
				  sprintf(histname1, "Unfold/%s%i_pt%i_eta0_%i",regN[iun],ity,ipt, var[ivar]); //Tunfolded_typ_0_pt0_eta0_3
				  datahist[iun][ity][ivar][ipt] = (TH1F*) Unfoldroot->Get(histname1);
				  cout << histname1 << endl;
				  
				  datahist[iun][ity][ivar][ipt]->Scale(1/(datahist[iun][ity][ivar][ipt]->Integral()));

				  //Devide the histogram with bin width
				  int tmpnbn = datahist[iun][ity][ivar][ipt]->GetNbinsX();
				  for (int ix=0; ix<tmpnbn; ix++) {
					  double awidth = datahist[iun][ity][ivar][ipt]->GetBinWidth(ix+1); // /tmpwid;
					  datahist[iun][ity][ivar][ipt]->SetBinContent(ix+1, datahist[iun][ity][ivar][ipt]->GetBinContent(ix+1)/awidth);
					  double error = datahist[iun][ity][ivar][ipt]->GetBinError(ix+1);
					  datahist[iun][ity][ivar][ipt]->SetBinError(ix+1, error/awidth);
				  } 
////////////////////////////////////
			  }
		  }
	  } //unfolded hist
	  
   //****************************************************************************************Declearation of the Functions
   TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1F* data, TH1F* MC[Nplot[0]], char* lowpadx);

    TCanvas *cpt0 = new TCanvas("cpt0", "canvas", 575,600 );  //for ESVs

    const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
    const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
    const char* htrang[8]={"73 < H_{T,2} < 93", "93 < H_{T,2} < 165", "165 < H_{T,2} < 225", "225 < H_{T,2} < 298", "298 < H_{T,2} < 365", "365 < H_{T,2} <452", "452 < H_{T2} <557","H_{T,2} >557"};
    const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
    const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};

  
 for(int ity=0; ity <type; ity++){
  for(int ivar =0 ; ivar < 5 ; ivar++){
    for(int iplot = 0 ; iplot <8 ; iplot++){
     int ifornt =102; 
     datahist[ity][ivar][iplot]->SetTitleOffset(0.4);
     datahist[ity][ivar][iplot]->SetTitleFont(ifornt);
      datahist[ity][ivar][iplot]->SetTitleSize(0.02);
 //   datahist[ity][ivar][iplot]->SetTitleFont(ifornt);
     
      datahist[ity][ivar][iplot]->GetXaxis()->SetLabelSize(0.03);
      datahist[ity][ivar][iplot]->GetXaxis()->SetTitleSize(0.045);
      datahist[ity][ivar][iplot]->GetXaxis()->SetTitleOffset(0.7);
 //   datahist[ity][ivar][iplot]->GetXaxis()->SetTitleFont(ifornt);
      datahist[ity][ivar][iplot]->GetXaxis()->CenterTitle();    
      
      datahist[ity][ivar][iplot]->GetYaxis()->SetLabelSize(0.03);
      datahist[ity][ivar][iplot]->GetYaxis()->SetTitleSize(0.040);
      datahist[ity][ivar][iplot]->GetYaxis()->SetTitleOffset(1.0);
 //   datahist[ity][ivar][iplot]->GetYaxis()->SetTitleFont(ifornt);     
      datahist[ity][ivar][iplot]->GetYaxis()->CenterTitle();
  
      //  gStyle->SetTitleFontSize(.08);
      datahist[ity][ivar][iplot]->SetLineWidth(2);
  
 if(iplot<=6){sprintf(Title,"%s: %s     %i <H_{T,2}< %i %s",itypeN[ity], Esvsym[ivar] , HT2range[iplot] , HT2range[iplot+1] ,"GeV/c" );}
 else if(iplot==7){ sprintf(Title,"%s:  %s      H_{T,2} > %i %s",itypeN[ity], Esvsym[ivar] ,  HT2range[iplot] ,"GeV/c" );}
   sprintf(Yaxis," %s" ,Esvlogy[ivar]);
   datahist[ivar][ity][iplot]->SetTitle(Title);
   datahist[ivar][ity][iplot]->GetXaxis()->SetTitle("");
   datahist[ivar][ity][iplot]->GetYaxis()->SetTitle(Yaxis);
 
   
   //**********************************************************************//calcuation of histogram means
 
   for(int iout = 0; iout < nmc ;iout++){
//   MC_hist[iout][ivar][iplot]->SetLineColor(color[iout]);
      MC_hist[iout][ity][ivar][iplot]->SetLineWidth(1);
   
      meanMC[iout][ity][ivar]->SetBinContent(iplot+1 , MC_hist[iout][ivar][iplot]->GetMean());
      meanData[iout][ity][ivar]->SetBinContent(iplot+1, datahist[ivar][iplot]->GetMean());
     //meanMC[iout][ity][ivar]->SetLineColor(color[iout]);

    sprintf(Title,"%s",Esvname[ivar] );
   //sprintf(Xaxis," %s" ,Esvlogx[ivar]);
   sprintf(Yaxis," %s" ,Esvlogx[ivar]);
   meanMC[iout][ity][ivar]->SetTitle(Title);
   meanData[iout][ity][ivar]->SetTitle(Title);
   meanData[iout][ity][ivar]->SetTitleFont(122);
   meanData[iout][ity][ivar]->GetXaxis()->SetTitle("");
   meanData[iout][ity][ivar]->GetYaxis()->SetTitle(Yaxis);
  }
   
 TH1F *MC_input[nmc];
 const char *MCinput_index[nmc];
 for(int iout = 0 ; iout < outnum ; iout++){
   MC_input[iout]=MC_hist[iout][ity][ivar][iplot];
  }
char lplot_xtitle[100];
  sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
//  float ratio_range1[2]={1.2,0.9};
  int num1[2]={nout,1} ;
  float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

  cpt0->cd();
cpt0->SetBorderSize(0);
cpt0->SetRightMargin(0.0);
cpt0->SetTopMargin(0.0);

  cpt0 =(TCanvas*)(ratio_can(num1, lpos1, datahist[ivar][iplot], MC_input, lplot_xtitle));

//Nplot[0] = number of MC enetered
//Nplot[1] = place 1 if upper part is log scale needed
//plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
//plegend[4]= text size
//plegend[5-6]= ratio plot axis range
//data = data histogram
// MC = monte carlo histogram array
//legendN = name of the legends for mC one by one
//lowpadx = x axis title of the ratio plot



   if(ity==0 && ivar==0 && iplot ==0){cpt0->Print("Unfolded_plot.pdf(","pdf");
   }else if(ity=1 && ivar==4 && iplot==7) {cpt0->Print("Unfolded_plot.pdf)","pdf");
   }else{
     cpt0->Print("Unfolded_plot.pdf.pdf","pdf");};

     }
   }  //end of phase space cut and variable loop

 } //End of 




}  // end of main program


//-------------------------------------------Ratio plot function

TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1F* data, TH1F* MC[Nplot[0]], char* lowpadx){
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
   TH1F *chidatal = (TH1F*)data->Clone("chidatal");   //for chi square test

   TH1F *chiMCl = (TH1F*)MC[iup]->Clone("chiMCl");   //for chi square test
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
   
MC[iup]->Draw("same hist  ");

    }//------------------- end of Montecarlo loop

/*int tmpnbd = data->GetNbinsX();
      for (int ix=0; ix<tmpnbd; ix++) {
        double awidth = data->GetBinWidth(ix+1); // /tmpwid;
        data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
        double error = data->GetBinError(ix+1);
        data->SetBinError(ix+1, error/awidth);
      }*/


 data->Draw(" same e2");

//----------------------------------------------//Maximum Uncertainty with respect to Monash
// make it off if not needed
/*
TH1F *MC_inputerr = (TH1F*)data->Clone("MC_inputerr");
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
   legendn->AddEntry(data, "Unfold Data 2017","lp");
  const char* legendN[40] ={"MC I","MC II","#alpha_{s} -4%","#alpha_{s} -6%","#alpha_{s} -8%","#alpha_{s} -10%","#alpha_{s} -12%","#alpha_{s} -14%","#alpha_{s} -16%","#alpha_{s} -18%","#alpha_{s} -20%","pythia12","pythia13","pythia14","pythia15","pythia16","pythia17","pythia18","pythia19","pythia20","Monash","#alpha_{s} -2%","#alpha_{s} -4%","#alpha_{s} -6%","#alpha_{s} -8%","#alpha_{s} -10%","#alpha_{s} -12%","#alpha_{s} -14%","#alpha_{s} -16%","#alpha_{s} -18%","#alpha_{s} -20%","pythia12","pythia13","pythia14","pythia15","pythia16","pythia17","pythia18","pythia19","pythia20"}; 
   for(int iup =0 ; iup <  Nplot[0] ; iup++){
    sprintf(MCindex,"%s" ,legendN[iup]);      // use if chi2 is not needed in the legend
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
   TH1F *rh2;
  //if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
 // if(  data->GetBinContent(10) > 0) {
   rh2 = (TH1F*)MC[ilow]->Clone("rh2"); 
   rh2->Sumw2();
   rh2->Divide(data);     //MC devide by data
  /* }else{
   rh2 = (TH1F*)data->Clone("rh2");
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
   rh2->Draw("hist same e1 ");

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



