// DataMC_Hist.C

//Read the Data root file and the MC root file and plot them 
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

void DataMC_Hist_Jet(){

TH1F *MC_hist[30][10][10];
//TFile *MC_root[30];

  char histname1[100];
  char histname2[100];
  char Title[100];
  char Xaxis[100];
  char Yaxis[100];
  char ratioXaxis[100];
  
Int_t color[10] ={1,2,4,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t var[5]={3,9,15,18,24};
//  Int_t HT2range[9]={62, 83, 145, 208, 268, 329, 413, 515, 3000};
  Int_t HT2range[9]={83, 109, 172, 241, 309, 377, 462, 570, 3000};

//  Int_t iter[4][8]={{4,4,4,4,5,4,4,3},{4,4,4,4,4,4,4,4},{4,7,7,4,5,5,5,4},{4,5,4,5,4,5,5,5}};

//******************************************************************************
Int_t outnum =0;    //initialize the number of root file readed 

char line[256];
  ifstream myfile ("list1.txt");
  if (myfile.is_open())
  {
    while (myfile)
    {                                                     //Start the loop using output number
       myfile.getline(line,256); 
   if  (line[0] != '\0') {
    cout << "Root file name = "<< line << endl;
    TFile *MC_root = TFile::Open(line);
    

    for(int ivar =0 ; ivar < 5 ; ivar++){
    for(Int_t ipt =0; ipt < 8 ; ipt++){

	sprintf(histname2, "analyzeBasicPat/reco_typ_0_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt6_eta0_15
//       sprintf(histname2, "analyzeBasicPat/reco_typ_1_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt6_eta0_15
//      sprintf(histname2, "Roo_Esv_var%i_pt%i", var[ivar],ipt);
      MC_hist[outnum][ivar][ipt] = (TH1F*) MC_root->Get(histname2);
      cout << histname2 << endl;

 MC_hist[outnum][ivar][ipt]->Scale(1/(MC_hist[outnum][ivar][ipt]->Integral()));

//devision of bin width
//////////////////////////// 
int tmpnbn = MC_hist[outnum][ivar][ipt]->GetNbinsX();
     for (int ix=0; ix<tmpnbn; ix++) {
        double awidth = MC_hist[outnum][ivar][ipt]->GetBinWidth(ix+1); // tmpwid;
        MC_hist[outnum][ivar][ipt]->SetBinContent(ix+1, MC_hist[outnum][ivar][ipt]->GetBinContent(ix+1)/awidth);
        double error = MC_hist[outnum][ivar][ipt]->GetBinError(ix+1);
        MC_hist[outnum][ivar][ipt]->SetBinError(ix+1, error/awidth);
      } 

    }
  }  //end of one MCinput root file  reading
outnum++;
}          

  }
myfile.close();

}//end of loop using output root file number

cout <<"number of root file present in that directory = " <<outnum << endl;
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//    TFile *file1 = TFile::Open("/home/suman/Paradox/ESV_charged/QCD_ESV_CMSSW/root_output/PLOTS_V1/ESV_MINIAOD17_BCDEF_Data_un15.root");  // data root file
   //TFile *file1 = TFile::Open("Data2017_C5TuneBinned_4april.root");  // data root file
 TFile *file1 = TFile::Open("JetHT_UL_JECv4_JRv2.root");  // data root file
//TFile *file2 = TFile::Open("/home/suman/FPythia2/CLHEPESV/DataMCcomp/Esv_rebined_testing.root");
  
//Data ******************************************************************
  TH1F *datahist[10][10];
  
  for(int ivar=0; ivar < 5 ; ivar ++){
    for(int ipt = 0 ; ipt < 8 ; ipt++){
      
    //  sprintf(histname1, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_typ_0_pt%i_eta0_%i", iter[ivar][ipt],ipt, var[ivar]); //bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt5_eta1_9
//        sprintf(histname1, "analyzeBasicPat/reco_typ_1_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
      sprintf(histname1, "analyzeBasicPat/reco_typ_0_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
	    datahist[ivar][ipt] = (TH1F*) file1->Get(histname1);
      cout << histname1 << endl;
/////////////////////////////////////////////
//normalized before deviding with bin width
datahist[ivar][ipt]->Scale(1/(datahist[ivar][ipt]->Integral()));

//Devide the histogram with bin width
int tmpnbn = datahist[ivar][ipt]->GetNbinsX();
      for (int ix=0; ix<tmpnbn; ix++) {
        double awidth = datahist[ivar][ipt]->GetBinWidth(ix+1); // /tmpwid;
        datahist[ivar][ipt]->SetBinContent(ix+1, datahist[ivar][ipt]->GetBinContent(ix+1)/awidth);
        double error = datahist[ivar][ipt]->GetBinError(ix+1);
        datahist[ivar][ipt]->SetBinError(ix+1, error/awidth);
      } 
////////////////////////////////////
    }
  }//end of data histogram reading
 
//////////////////////////////////////////////////// mean plot

double meanplotx[9]={73, 93, 165, 225, 298, 365, 452, 557, 800};
TH1F *meanMC[outnum][5];
TH1F *meanData[outnum][5];
char meanMC_name[200];
char meandata_name[200];

for(int iout =0 ; iout < outnum ; iout++){

    for(int imean =0; imean < 5 ; imean++){
    sprintf(meanMC_name, "Mean_of%i_%i", imean,iout);
    meanMC[iout][imean]=new TH1F(meanMC_name,"mean_hist_MC", 8 , meanplotx);
    meanMC[iout][imean]->Sumw2();
    sprintf(meandata_name, "Mean_ofData%i_%i", imean,iout);
    meanData[iout][imean]=new TH1F(meandata_name,"mean_hist_Data", 8 , meanplotx);
    meanData[iout][imean]->Sumw2();
          
 }
}//end of mean plot histogram define loop

  
//****************************************************************************************Declearation of the Functions
TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1F* data, TH1F* MC[Nplot[0]], char* lowpadx);
TCanvas* table_can(int nvar[2], Float_t range[3],Double_t var[4][8] , const char* varnam[nvar[1]]);
  

  TCanvas *cpt0 = new TCanvas("cpt0", "canvas", 575,600 );  //for ESVs
  TCanvas *cpt1 = new TCanvas("cpt1", "canvas_mean", 900,1000 ); //for mean of eventshape
  TCanvas *cpt2 = new TCanvas("cpt2", "canvas_chi", 900,1000 ); //for chi2 table
  TCanvas *cpt3 = new TCanvas("cpt3", "canvas_chi_para", 500,500 ); //for chi2 graph

  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[8]={"73 < H_{T,2} < 93", "93 < H_{T,2} < 165", "165 < H_{T,2} < 225", "225 < H_{T,2} < 298", "298 < H_{T,2} < 365", "365 < H_{T,2} <452", "452 < H_{T2} <557","H_{T,2} >557"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};

  
  for(int ivar =0 ; ivar < 5 ; ivar++){
    for(int iplot = 0 ; iplot <8 ; iplot++){
     int ifornt =102; 
      datahist[ivar][iplot]->SetTitleOffset(0.4);
  //    datahist[ivar][iplot]->SetTitleFont(ifornt);
      
      datahist[ivar][iplot]->SetTitleSize(0.02);
      datahist[ivar][iplot]->GetYaxis()->SetLabelSize(0.03);
      datahist[ivar][iplot]->GetXaxis()->SetLabelSize(0.03);
      
      datahist[ivar][iplot]->GetYaxis()->SetTitleSize(0.040);
      datahist[ivar][iplot]->GetYaxis()->SetTitleOffset(1.0);
      datahist[ivar][iplot]->GetXaxis()->SetTitleSize(0.045);
      datahist[ivar][iplot]->GetXaxis()->SetTitleOffset(0.7);
//      datahist[ivar][iplot]->GetXaxis()->SetTitleFont(ifornt);
//      datahist[ivar][iplot]->GetYaxis()->SetTitleFont(ifornt);     
  //    datahist[ivar][iplot]->SetTitleFont(ifornt);
      datahist[ivar][iplot]->GetYaxis()->CenterTitle();
      datahist[ivar][iplot]->GetXaxis()->CenterTitle();    
  //   gStyle->SetTitleFontSize(.08);
   datahist[ivar][iplot]->SetLineWidth(2);
   
//   sprintf(Title,"%s    %s %s",Esvname[ivar] ,htrang[iplot] ,"GeV" );
//   sprintf(Xaxis," %s" ,Esvlogx[ivar]);
//   sprintf(Yaxis," %s" ,Esvlogy[ivar]);
//   datahist[ivar][iplot]->SetTitle(Title);
//   datahist[ivar][iplot]->GetXaxis()->SetTitle(Xaxis);
//   datahist[ivar][iplot]->GetYaxis()->SetTitle(Yaxis);
   
   //**********************************************************************//calcuation of histogram means
   for(int iout = 0; iout < outnum ;iout++){
//   MC_hist[iout][ivar][iplot]->SetLineColor(color[iout]);
   MC_hist[iout][ivar][iplot]->SetLineWidth(1);
   
meanMC[iout][ivar]->SetBinContent(iplot+1 , MC_hist[iout][ivar][iplot]->GetMean());
meanData[iout][ivar]->SetBinContent(iplot+1, datahist[ivar][iplot]->GetMean());
//meanMC[iout][ivar]->SetLineColor(color[iout]);

sprintf(Title,"%s",Esvname[ivar] );
  // sprintf(Xaxis," %s" ,Esvlogx[ivar]);
   sprintf(Yaxis," %s" ,Esvlogx[ivar]);
   meanMC[iout][ivar]->SetTitle(Title);
   meanData[iout][ivar]->SetTitle(Title);
   meanData[iout][ivar]->SetTitleFont(122);
   meanData[iout][ivar]->GetXaxis()->SetTitle("");
   meanData[iout][ivar]->GetYaxis()->SetTitle(Yaxis);


  }
   
   //**********************************************************//
 //sprintf(Title,"%s    %s %s",Esvname[ivar] ,htrang[iplot] ,"GeV/c" );

//   sprintf(Title,"  %s %s",htrang[iplot] ," GeV" );
 // sprintf(Xaxis," %s" ,Esvlogx[ivar]);
 if(iplot<=6){ sprintf(Title,"Reco Jets:  %s      %i <H_{T,2}< %i %s",Esvsym[ivar] , HT2range[iplot] , HT2range[iplot+1] ,"GeV/c" );}
 else if(iplot==7){ sprintf(Title,"Reco Jets :  %s      H_{T,2} > %i %s",Esvsym[ivar] ,  HT2range[iplot] ,"GeV/c" );}

   sprintf(Yaxis," %s" ,Esvlogy[ivar]);
   datahist[ivar][iplot]->SetTitle(Title);
   datahist[ivar][iplot]->GetXaxis()->SetTitle("");
   datahist[ivar][iplot]->GetYaxis()->SetTitle(Yaxis);


//*******************************************************************Ratio plot  for Eventshape  Variables

  TH1F *MC_input[outnum];
 const char *MCinput_index[outnum];
 for(int iout = 0 ; iout < outnum ; iout++){
   MC_input[iout]=MC_hist[iout][ivar][iplot];
  }
char lplot_xtitle[100];
  sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
//  float ratio_range1[2]={1.2,0.9};
  int num1[2]={outnum,1} ;
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



   if(ivar==0 && iplot ==0){cpt0->Print("Esv_dataMC_Jets.pdf(","pdf");
   }else if(ivar==4 && iplot==7) {cpt0->Print("Esv_dataMC_Jets.pdf)","pdf");
   }else{
     cpt0->Print("Esv_dataMC_Jets.pdf","pdf");};

    }
  }  //end of phase space cut and variable loop


//-------------------------------------------------------//chi2 value calculation
double chi2value[outnum][5][8];
double chi2Y[5][8][outnum];
double chi2X[outnum];

double chi2Ndfvalue[outnum][5][8];
double chiNdf2Y[5][8][outnum];



for(int iout =0; iout < outnum ;iout++){
for(int ivar =0 ; ivar < 5 ; ivar++){
    for(int iplot = 0 ; iplot <8 ; iplot++){
  int n = datahist[ivar][iplot]->GetNbinsX();
  cout << " Numbers of bin = " << n << endl;
  Double_t res[n] , chi2;
  Int_t ndf ,igood;
   TH1F *chidata = (TH1F*)datahist[ivar][iplot]->Clone("chidata");   //for chi square test

   TH1F *chiMC = (TH1F*)MC_hist[iout][ivar][iplot]->Clone("chiMC");   //for chi square test
   chiMC->Chi2TestX(chidata, chi2, ndf, igood, "WU", res);

cout << "chi2 value=" << chi2 << endl;
chi2value[iout][ivar][iplot]=chi2;

chi2Y[ivar][iplot][iout]=chi2;
chi2X[iout]=iout;

chi2Ndfvalue[iout][ivar][iplot]=chi2/ndf;
chiNdf2Y[ivar][iplot][iout]=chi2/ndf;


cout <<"Chi2 Values" << chi2Y[ivar][iplot][iout] << endl;
   }
 }
}//end of chi2 calculation

//---------------------------------------------------------//chi2 value Plot
//double chi2X_edit[11]={0,-2,-4,-6,-8,-10,-12,-14,-16,-18,-20};
//for(int iout =0; iout < outnum ;iout++){
//double chi2value[4][8];

/*
for(int ivar =0 ; ivar < 5 ; ivar++){
    for(int iplot = 0 ; iplot <8 ; iplot++){
cpt3->cd();
TGraph* chi_graph = new TGraph(outnum,chi2X, chi2Y[ivar][iplot]);
//TGraph* chi_graph = new TGraph(outnum,chi2X_edit, chi2Y[ivar][iplot]);
chi_graph->Draw("AC*");
chi_graph->SetLineColor(2);
sprintf(Title,"#chi_{2} value for %s    %s %s",Esvname[ivar] ,htrang[iplot] ,"GeV/c" );
chi_graph->SetTitle(Title);
chi_graph->GetYaxis()->SetTitle("#chi^{2}");
chi_graph->GetXaxis()->SetTitle("% of ");
chi_graph->GetYaxis()->SetTitleOffset(1.4);
cpt3->Update();
if(ivar==0 && iplot==0){cpt3->Print("chi2_vs_Par.pdf(","pdf");
   }else if(ivar==4 && iplot==7) {cpt3->Print("chi2_vs_Par.pdf)","pdf");
    }else{
     cpt3->Print("chi2_vs_Par.pdf","pdf");
      };
cpt3->Clear();  
 
  }
 }
*/

//-----------------------------------------------------------------------Graph for chi2/ndf
/*
//double chi2X_edit[11]={0,-2,-4,-6,-8,-10,-12,-14,-16,-18,-20};
//double chi2value[4][8];
for(int ivar =0 ; ivar < 4 ; ivar++){
    for(int iplot = 0 ; iplot <8 ; iplot++){
cpt3->cd();
TGraph* chiNdf_graph = new TGraph(outnum,chi2X, chiNdf2Y[ivar][iplot]);
//TGraph* chi_graph = new TGraph(outnum,chi2X_edit, chi2Y[ivar][iplot]);
chiNdf_graph->Draw("AC*");
chiNdf_graph->SetLineColor(2);
sprintf(Title,"#chi^{2}/Ndf value for %s    %s %s",Esvname[ivar] ,htrang[iplot] ,"GeV/c" );
chiNdf_graph->SetTitle(Title);
chiNdf_graph->GetYaxis()->SetTitle("#chi^{2}/Ndf");
chiNdf_graph->GetXaxis()->SetTitle("Weight variation");
chiNdf_graph->GetYaxis()->SetTitleOffset(1.4);
cpt3->Update();
if(ivar==0 && iplot==0){cpt3->Print("chi2Ndf_vs_Par.pdf(","pdf");
   }else if(ivar==3 && iplot==7) {cpt3->Print("chi2Ndf_vs_Par.pdf)","pdf");
    }else{
     cpt3->Print("chi2Ndf_vs_Par.pdf","pdf");
      };
cpt3->Clear();

  }
 }
*/

//----------------------------------------------chi2 of different energy range in the same plot
/*
int chi_cl[8]={1,2,3,4,5,6,7,8};

for(int ivar =0 ; ivar < 4 ; ivar++){
cpt3->cd();
auto g = new TMultiGraph();

for(int iplot = 0 ; iplot <8 ; iplot++){
//TGraph* chi_graph1 = new TGraph(outnum,chi2X_edit, chi2Y[ivar][iplot]);
TGraph* chi_graph1 = new TGraph(outnum,chi2X, chi2Y[ivar][iplot]);
chi_graph1->SetTitle(htrang[iplot]);
chi_graph1->SetLineColor(chi_cl[iplot]);

g->Add(chi_graph1);

//if(iplot==0){chi_graph1->Draw("apl");}
//if(iplot>0){chi_graph1->Draw("pl SAME");}

cpt3->Modified();
cpt3->Update();
  }
g->Draw("AC");
sprintf(Title,"#chi^{2} variation of %s  ",Esvname[ivar] );
g->SetTitle(Title);
g->GetYaxis()->SetTitle("#chi^{2}");
g->GetXaxis()->SetTitle("% of Change for ISR #alpha_{s}");


cpt3->BuildLegend(0.7,.7,.95,.95,"","lp");
if(ivar==0){cpt3->Print("chi_all_ph.pdf(","pdf");
      }else if(ivar==3){cpt3->Print("chi_all_ph.pdf)","pdf");
   }else{ cpt3->Print("chi_all_ph.pdf","pdf");
}

cpt3->Clear();

}
*/
//-----------------------------------------------------chi2/Ndf for all HT2 range
/*
int chi_cl[8]={1,2,3,4,5,6,7,8};

for(int ivar =0 ; ivar < 5 ; ivar++){
cpt3->cd();
auto g = new TMultiGraph();

for(int iplot = 0 ; iplot <8 ; iplot++){
//TGraph* chi_graph1 = new TGraph(outnum,chi2X_edit, chi2Y[ivar][iplot]);
TGraph* chi_graph1 = new TGraph(outnum,chi2X, chiNdf2Y[ivar][iplot]);
chi_graph1->SetLineWidth(2);
chi_graph1->SetTitle(htrang[iplot]);
chi_graph1->SetLineColor(chi_cl[iplot]);

g->Add(chi_graph1);

//if(iplot==0){chi_graph1->Draw("apl");}
//if(iplot>0){chi_graph1->Draw("pl SAME");}

cpt3->Modified();
cpt3->Update();
  }
g->Draw("AC");
sprintf(Title,"#chi^{2}/Ndf Variation of %s  ",Esvname[ivar] );
g->SetTitle(Title);
g->GetYaxis()->SetTitle("#chi^{2}/Ndf");
g->GetXaxis()->SetTitle("% of change Space shower #alpha_{s}");


cpt3->BuildLegend(0.8,0.65,.98,.90,"","lp");
if(ivar==0){cpt3->Print("chi2Ndf.pdf(","pdf");
      }else if(ivar==3){cpt3->Print("chi2Ndf.pdf)","pdf");
   }else{ cpt3->Print("chi2Ndf.pdf","pdf");
}

cpt3->Clear();

}


*/


//--------------------------------------------------Chi2 value in table
/*for(int iout =0 ; iout<outnum ; iout++){
  Float_t xrange = 30;
  Float_t yrange = 40;
   Int_t w = 650;
   Int_t h = w*yrange/xrange;

  // TCanvas *c1 = new TCanvas("c1","c1",200,10,w,h);
   int nvar1[2]={4,8};
   float range1[3] = {0.0,30.0,40.0};
   cpt2=(TCanvas*)(table_can(nvar1, range1, chi2value[iout] ,Esvsym));
//   table(0.5*xrange+0.5,xrange-0.5,yrange,t,symbol2,0);
   TText *tlabel = new TText(0,0,"a");
   tlabel->SetTextFont(72);
   tlabel->SetTextSize(0.018);
   tlabel->SetTextAlign(22);
//   tlabel->DrawText(0.5*xrange,1.3,  "Input characters are standard keyboard characters");
   cpt2->Modified();
   cpt2->Update();
//   c1->Print("pstable1.pdf");

if(iout==0){cpt2->Print("chi2table.pdf(","pdf");
   }else if(iout==outnum-1) {cpt2->Print("chi2table.pdf)","pdf");
   }else{
     cpt2->Print("chi2table.pdf","pdf");};

} //end of chi2test loop

*/

//---------------------------------------------------Ratio plot for mean value
/*
for(int ivar =0 ; ivar <5 ;ivar++){
 
 TH1F *mean_MCinput[outnum];
for(int iout =0 ; iout < outnum ; iout++){
 mean_MCinput[iout]= meanMC[iout][ivar];
 }
 char lplot_xtitle[100];
 sprintf(lplot_xtitle, "%s","H_{T,2} (GeV)");
 int num2[2]={outnum,0} ;
 float lpos2[7] ={.65,0.40,0.87,0.68, .02, 1.1,.9};

  cpt1 =(TCanvas*)(ratio_can(num2,lpos2,meanData[0][ivar],mean_MCinput, lplot_xtitle));

if(ivar==0){cpt1->Print("Esv_mean.pdf(","pdf");
   }else if(ivar==4) {cpt1->Print("Esv_mean.pdf)","pdf");
   }else{
     cpt1->Print("Esv_mean.pdf","pdf");};
   } //end of mean plot

*/
char meanYaxis[100];
char meantitle[100];


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
      data->SetMarkerStyle(8);
      data->SetMarkerSize(1);
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
//   data->Draw("e2");


   data->Draw(" ");


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
cout << "Ndf value =" << ndfl << endl;
cout << "chi2 value=" << chi2l << endl;
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


// data->Draw(" same e2");

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
   legendn->AddEntry(data, "Data 2017","lp");
  const char* legendN[40] ={"PY8 CP5 Tune","Pythia8(Prof2-tune) ","#alpha_{s} -4%","#alpha_{s} -6%","#alpha_{s} -8%","#alpha_{s} -10%","#alpha_{s} -12%","#alpha_{s} -14%","#alpha_{s} -16%","#alpha_{s} -18%","#alpha_{s} -20%","pythia12","pythia13","pythia14","pythia15","pythia16","pythia17","pythia18","pythia19","pythia20","Monash","#alpha_{s} -2%","#alpha_{s} -4%","#alpha_{s} -6%","#alpha_{s} -8%","#alpha_{s} -10%","#alpha_{s} -12%","#alpha_{s} -14%","#alpha_{s} -16%","#alpha_{s} -18%","#alpha_{s} -20%","pythia12","pythia13","pythia14","pythia15","pythia16","pythia17","pythia18","pythia19","pythia20"}; 
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


//------------------------------------------//Function to store chi2 value in table

TCanvas* table_can(int nvar[2], Float_t range[3],Double_t var[4][8], const char* varnam[nvar[1]]){
   Float_t xrange = 30;
   Float_t yrange = 40;
   Int_t w = 650;
   Int_t h = w*yrange/xrange;

   TLatex *t = new TLatex(0,0,"a");
   t->SetTextSize(0.02);
   t->SetTextFont(62);
   t->SetTextAlign(22);

   TCanvas *c1 = new TCanvas("c1","c1",200,10, w, h);
   c1->Range(0,0,xrange,yrange);


   float x1 =range[0];
   float x2 = range[1];
   Float_t y1  = 2.5;
   Float_t y2  = range[2] - 0.5;
   Float_t dx  = (x2-x1)/(nvar[1]+1);
   Float_t dy  = 2.5;
   Float_t yy   = y2 -dy;
   Float_t xx = x1  + 0.5*dx;

   TLine *line = new TLine();

   line->DrawLine(x1,y1,x1,y2);
   line->DrawLine(x1,y1,x2,y1);
   line->DrawLine(x1,y2,x2,y2);
   line->DrawLine(x2,y1,x2,y2);
   line->DrawLine(x1,y2-2.0,x2,y2-2.0);

     int xt=xx;
     int yt = yy-dy/2;

   for(int ichi=0; ichi < nvar[0] ; ichi++){
   char vartext[100];
   TLatex *vartit = new TLatex(0,0,"a");
   vartit->SetTextSize(0.015);
   vartit->SetTextFont(72);
   vartit->SetTextAlign(22);
   sprintf(vartext,"%s",varnam[ichi]);
   vartit->DrawLatex(xt,yt+dy/2,vartext);
   line->DrawLine(x1, yt, x2, yt);
    int xv=xx;
       for(int iipt =0 ; iipt < nvar[1];iipt++){
   int yv = yt+dy/2;

   line->DrawLine(xv+dx/2, y1, xv+dx/2, y2);
   TLatex *tit = new  TLatex(0,0,"a");
   tit->SetTextSize(0.015);
   tit->SetTextFont(72);
   tit->SetTextAlign(22);
   char range[100];
  const char* phrang[8]={"73<H_{T,2}<93", "93<H_{T,2}<165", "165<H_{T,2}<225", "225<H_{T,2}<298", "298<H_{T,2}<365", "365<H_{T,2}<452", "452<H_{T,2}<557","H_{T,2}>557"};

   sprintf(range,"%s",phrang[iipt]);

   tit->DrawLatex(xv+dx,y2-0.6,range);

   char text[12];
      sprintf(text,"%.2f",var[ichi][iipt]);
      t->DrawLatex(xv+dx,yt+dy/2,text);
      xv += dx;
      }
yt -=dy;

  }
return c1;
} //end of chi2 function


