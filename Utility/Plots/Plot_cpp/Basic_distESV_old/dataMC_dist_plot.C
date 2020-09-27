// dataMC_dist_plot.C

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

void dataMC_dist_plot(){

int numplot = 11;
char Title[100];
char Xaxis[100];
char Yaxis[100];

char histname2[100];  // For monte carlo histogram reading

//Check the folder for all Monte carlo file for plots
 TH1F *MC_hist[30][11];  //maximum number of Monte carlo can be analysis in 30
 Int_t color[10] ={1,2,4,5,6,46,3,28,38,42};  // define the color for different histograms
 
 //--------------------------------------------------------------------------------------
 TH1F *Numjet0= new TH1F("numjet0","Number_of_jets",30,0.5,30.5);
 Numjet0->Sumw2();
 
 TH1F *Numjet1= new TH1F("numjet1","Number_of_jets",30,0.5,30.5);
 Numjet1->Sumw2();
 
 TH1F *Numjet2= new TH1F("numjet2","Number_of_jets",30,0.5,30.5);
 Numjet2->Sumw2();
 TH1F *Numjet3= new TH1F("numjet3","Number_of_jets",30,0.5,30.5);
 Numjet3->Sumw2();
 
 TH1F *Numjetd= new TH1F("numjetd","Number_of_jets",30,0.5,30.5);
 Numjetd->Sumw2();
 
 
 Int_t outnum =0;
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
	   
	   
	   
	   for(int ivar=0; ivar < numplot ; ivar ++){
	     if(ivar==0){sprintf(histname2, "analyzeBasicPat/recojet1_pt_0");};
	     if(ivar==1){sprintf(histname2, "analyzeBasicPat/recojet2_pt_0");};
	     if(ivar==2){sprintf(histname2, "analyzeBasicPat/recojetallave_pt_0");};
	     if(ivar==3){sprintf(histname2, "analyzeBasicPat/hjetdpt_0");};
	     if(ivar==4){sprintf(histname2, "analyzeBasicPat/hjetptbypl_0");};
	     if(ivar==5){sprintf(histname2, "analyzeBasicPat/recojet1_eta");};
	     if(ivar==6){sprintf(histname2, "analyzeBasicPat/recojet2_eta");};
	     if(ivar==7){sprintf(histname2, "analyzeBasicPat/recojt_phi");};
	     if(ivar==8){sprintf(histname2, "analyzeBasicPat/hjetdphi_0");};
	     if(ivar==9){sprintf(histname2, "analyzeBasicPat/njets_0");};
	     if(ivar==10){sprintf(histname2, "analyzeBasicPat/ncharges_0");};
	     /*
	     if(ivar==0){sprintf(histname2, "recojet1_pt_0");};
             if(ivar==1){sprintf(histname2, "recojet2_pt_0");};
             if(ivar==2){sprintf(histname2, "recojetallave_pt_0");};
             if(ivar==3){sprintf(histname2, "hjetdpt_0");};
             if(ivar==4){sprintf(histname2, "hjetptbypl_0");};
             if(ivar==5){sprintf(histname2, "recojet1_eta");};
             if(ivar==6){sprintf(histname2, "recojet2_eta");};
             if(ivar==7){sprintf(histname2, "recojt_phi");};
             if(ivar==8){sprintf(histname2, "hjetdphi_0");};
             if(ivar==9){sprintf(histname2, "njets_0");};
             if(ivar==10){sprintf(histname2, "ncharges_0");};
	     */
	     MC_hist[outnum][ivar]= (TH1F*) MC_root->Get(histname2);
	     cout << histname2 << endl;
	     //int MCbin = MChist[ivar]->GetNbinsX();
	     //int MCbinedge = MChist[ivar]->GetBinLowEdge(0);
	     //cout << "binedge = " <<MCbinedge << endl;
	     //for(int ibin =0 ; ibin <= MCbin ; ibin++ ){
	     //cout <<"bincontent["<<ibin <<"] =" << MChist[ivar]->GetBinContent(ibin) << endl;
	     //          }
	     
	     //if(ivar==9){ MC_hist[outnum][ivar]->Rebin(2);}
	     if(ivar==10){ MC_hist[outnum][ivar]->Rebin(4);}
	     if(ivar==8){MC_hist[outnum][ivar]->GetXaxis()->SetRangeUser(0,3);}
	     if(ivar==10){MC_hist[outnum][ivar]->GetXaxis()->SetRangeUser(0, 140);}
	     if(ivar==5 ||ivar==6){MC_hist[outnum][ivar]->GetXaxis()->SetRangeUser(-2.4,2.4);}
	     
	     
	     //Numjet->Reset();
	     if(ivar==9 && outnum==0){
	       for(int i=1;i<=30;i++){
		 double entries = MC_hist[outnum][ivar]->GetBinContent(i*2+1);
		 double error = MC_hist[outnum][ivar]->GetBinError(i*2+1);
		 Numjet0->SetBinContent(i,entries);
		 Numjet0->SetBinError(i,error);
	       }
	       
	       //MC_hist[outnum][ivar]= (TH1F)Numjet->Clone();
	       MC_hist[outnum][ivar]= Numjet0;
	     }
	     
	     if(ivar==9 && outnum==1){
	       for(int i=1;i<=30;i++){
		 double entries = MC_hist[outnum][ivar]->GetBinContent(i*2+1);
		 double error = MC_hist[outnum][ivar]->GetBinError(i*2+1);
		 Numjet1->SetBinContent(i,entries);
		 Numjet1->SetBinError(i,error);
	       }
	       
	       //MC_hist[outnum][ivar]= (TH1F)Numjet->Clone();
	       MC_hist[outnum][ivar]= Numjet1;
	     }
	     
	     if(ivar==9 && outnum==2){
	       for(int i=1;i<=30;i++){
		 double entries = MC_hist[outnum][ivar]->GetBinContent(i*2+1);
		 double error = MC_hist[outnum][ivar]->GetBinError(i*2+1);
		 Numjet2->SetBinContent(i,entries);
		 Numjet2->SetBinError(i,error);
	       }
	       
	       //MC_hist[outnum][ivar]= (TH1F)Numjet->Clone();
	       MC_hist[outnum][ivar]= Numjet2;
	     }
	     
	     if(ivar==9 && outnum==3){
	       for(int i=1;i<=30;i++){
		 double entries = MC_hist[outnum][ivar]->GetBinContent(i*2+1);
		 double error = MC_hist[outnum][ivar]->GetBinError(i*2+1);
		 Numjet3->SetBinContent(i,entries);
		 Numjet3->SetBinError(i,error);
	       }
	       
	       //MC_hist[outnum][ivar]= (TH1F)Numjet->Clone();
	       MC_hist[outnum][ivar]= Numjet3;
	     }
	     
	     if(ivar==9){MC_hist[outnum][ivar]->GetXaxis()->SetRangeUser(1,10);}
	     
	     
	     MC_hist[outnum][ivar]->Scale(1/(MC_hist[outnum][ivar]->Integral()));
	     
	     //devision of bin width
	     //////////////////////////// 
	     /*
	       int tmpnbn = MC_hist[outnum][ivar]->GetNbinsX();
	       for (int ix=0; ix<tmpnbn; ix++) {
	       double awidth = MC_hist[outnum][ivar]->GetBinWidth(ix+1); // tmpwid;
	       MC_hist[outnum][ivar]->SetBinContent(ix+1, MC_hist[outnum][ivar]->GetBinContent(ix+1)/awidth);
	       double error = MC_hist[outnum][ivar]->GetBinError(ix+1);
	       MC_hist[outnum][ivar]->SetBinError(ix+1, error/awidth);
	       }*/
	   } //end of one MCinput root file  reading
	   outnum++;
	 } //end of the file list2.txt
       }
     myfile.close();
     
   }//end of loop using output root file number
 
 
 cout <<"number of root file present in that directory = " <<outnum << endl;
 
 //Input root files for Data
 //TFile *file1 = TFile::Open("");  // data root file
 //TFile *file1 = TFile::Open("JETHT2017_31Mar18DEF_TrigV2_6Mar20.root");  // data root file
 //TFile *file1 = TFile::Open("Test_Data2017_31Mar2018.root");  // data root file
 //TFile *file1 = TFile::Open("Data2017_GT102_3sigmaBin_16April.root");  // data root file
 TFile *file1 = TFile::Open("JetHT_UL_JECv4_JRv2.root");  // data root file
//TFile *file1 = TFile::Open("Data2017_GT102_3sigmaBin_16April.root");  // data root file
 
 char histname1[100];
 
 TH1F *datahist[11];
 
 
 
 for(int ivar=0; ivar < numplot ; ivar ++){
   if(ivar==0){sprintf(histname1, "analyzeBasicPat/recojet1_pt_0");};
   if(ivar==1){sprintf(histname1, "analyzeBasicPat/recojet2_pt_0");};
   if(ivar==2){sprintf(histname1, "analyzeBasicPat/recojetallave_pt_0");};
   if(ivar==3){sprintf(histname1, "analyzeBasicPat/hjetdpt_0");};
   if(ivar==4){sprintf(histname1, "analyzeBasicPat/hjetptbypl_0");};
   if(ivar==5){sprintf(histname1, "analyzeBasicPat/recojet1_eta");};
   if(ivar==6){sprintf(histname1, "analyzeBasicPat/recojet2_eta");};
   if(ivar==7){sprintf(histname1, "analyzeBasicPat/recojt_phi");};
   if(ivar==8){sprintf(histname1, "analyzeBasicPat/hjetdphi_0");};
   if(ivar==9){sprintf(histname1, "analyzeBasicPat/njets_0");};
   if(ivar==10){sprintf(histname1, "analyzeBasicPat/ncharges_0");};
   
   datahist[ivar]= (TH1F*) file1->Get(histname1);
   cout << histname1 << endl;
   
   
   //if(ivar==9){datahist[ivar]->Rebin(2);}
   if(ivar==10){datahist[ivar]->Rebin(4);}
   if(ivar==8){datahist[ivar]->GetXaxis()->SetRangeUser(0,3);}
   if(ivar==10){datahist[ivar]->GetXaxis()->SetRangeUser(0,140);}
   if(ivar==5 ||ivar==6){datahist[ivar]->GetXaxis()->SetRangeUser(-2.4,2.4);}
   datahist[ivar]->Scale(1/(datahist[ivar]->Integral()));
   
   if(ivar==9){
     for(int i=1;i<=30;i++){
       double entries = datahist[ivar]->GetBinContent(i*2+1);
       double error = datahist[ivar]->GetBinError(i*2+1);
       Numjetd->SetBinContent(i,entries);
       Numjetd->SetBinError(i,error);
     }
     
     //datahist[ivar]= (TH1F)Numjet->Clone();
       datahist[ivar]= Numjetd;
   }
   
   
   if(ivar==9){datahist[ivar]->GetXaxis()->SetRangeUser(1,10);}
   
   /*
     int tmpnbn = datahist[ivar]->GetNbinsX();
     for (int ix=0; ix<tmpnbn; ix++) {
     double awidth = datahist[ivar]->GetBinWidth(ix+1); // /tmpwid;
     datahist[ivar]->SetBinContent(ix+1, datahist[ivar]->GetBinContent(ix+1)/awidth);
     double error = datahist[ivar]->GetBinError(ix+1);
     datahist[ivar]->SetBinError(ix+1, error/awidth);
     }*/
 }
 
 
 //-------------------------------------------------------------------------------------------------------
 const char* var_name[11] ={ "Pt of leading jet (GeV/c)", "Pt of second leading jet (GeV/c)","H_{T2} (GeV/c)","#Delta Pt of two leading jets (GeV/c)","Pt2 x sin( #Delta #phi )/Pt1","#eta of leading jet", "#eta of second leading jet", "#phi of leading jet","#Delta#phi of Jets","No. of jet"," No. of Charged particles"};
 const char* varlogx[11] ={"Pt of leading jet (GeV/c)", "Pt of second leading jet (GeV/c)","H_{T2} (GeV/c)","#Delta Pt of two leading jets (GeV/c)","Pt2 x sin( #Delta #phi )/Pt1","#eta of leading jet", "#eta of second leading jet", "#phi of leading jet","#Delta#phi of Jets","No. of jet","No. of Charged particles"};
 const char* varlogy[11] ={"1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d"};
 
 
 //for(int ivar=0; ivar < 10 ; ivar ++){ // loop for variables
 TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1F* data, TH1F* MC[Nplot[0]], char* lowpadx);  // define the ratio plot canvas
 
 TCanvas *cpt0 = new TCanvas("cpt0", "canvas", 900,1000 );  
 //TCanvas *cpt3 = new TCanvas("cpt3", "canvas_chi_para", 500,500 ); //for chi2 graph
 
 
 for(int ivar=0; ivar < 11 ; ivar ++){ // loop for variables
   
   datahist[ivar]->SetTitleOffset(0.4);
   datahist[ivar]->SetTitleSize(0.02);
   datahist[ivar]->GetYaxis()->SetLabelSize(0.03);
   datahist[ivar]->GetXaxis()->SetLabelSize(0.03);
   datahist[ivar]->GetYaxis()->SetTitleSize(0.040);
   datahist[ivar]->GetYaxis()->SetTitleOffset(1.0);
   datahist[ivar]->GetXaxis()->SetTitleSize(0.045);
   datahist[ivar]->GetXaxis()->SetTitleOffset(0.7);
   datahist[ivar]->GetYaxis()->CenterTitle();
   datahist[ivar]->GetXaxis()->CenterTitle();
   datahist[ivar]->SetLineWidth(1);
   
   sprintf(Title,"%s  ",var_name[ivar] );
   sprintf(Xaxis," %s" ,varlogx[ivar]);
   sprintf(Yaxis," %s%s" ,varlogy[ivar],var_name[ivar]);
   datahist[ivar]->SetTitle(Title);
   //   datahist[ivar]->GetXaxis()->SetTitle(Xaxis);
   datahist[ivar]->GetYaxis()->SetTitle(Yaxis);
   //for(int iout =0 ; iout < outnum ; iout++){    // loop on the file number
   
   //define arguments for ration plot function
   TH1F *MC_input[outnum];
   // const char *MCinput_index[outnum];
   for(int iout = 0 ; iout < outnum ; iout++){
     MC_input[iout]=MC_hist[iout][ivar];
   }
   char lplot_xtitle[100];
   sprintf(lplot_xtitle, "%s",varlogx[ivar]);  //have to change
   //  float ratio_range1[2]={1.2,0.9};
  int num1[2]={outnum,1} ;
  float lpos1[7] ={.6,0.7,0.9,0.88, .033, 2.2,.20};
  cout << "variable =" << ivar << endl;
  cpt0 =(TCanvas*)(ratio_can(num1, lpos1, datahist[ivar], MC_input, lplot_xtitle));
  if(ivar==0 ){cpt0->Print("Basic_dist_comp.pdf(","pdf");
  }else if(ivar==10) {cpt0->Print("Basic_dist_comp.pdf)","pdf");
  }else{
    cpt0->Print("Basic_dist_comp.pdf","pdf");};
  
  // end of file loop
  //cpt0->Clear();
  //------------------------------------------------------------
 } // end of variable loop
 
 //------------------------------------------------------------------------
 //chi2 value
 double chi2value[outnum][10];
 double chi2Y[10][outnum];
 double chi2X[outnum];
 
 double chi2Ndfvalue[outnum][10];
 double chiNdf2Y[10][outnum];
 /*
   
   for(int iout =0; iout < outnum ;iout++){
   //double chi2value[4][8];
   for(int ivar =0 ; ivar < 10 ; ivar++){
   
   //apply the chi2 test and retrieve the residuals
   int n = datahist[ivar]->GetNbinsX();
   //  int hist_low =datahist[ivar][iplot]->GetBinLowEdge(1);
   //  int hist_high = (datahist[ivar][iplot]->GetBinLowEdge(n)+datahist[ivar][iplot]->GetBinWidth(n));
   Double_t res[n] , chi2;
   Int_t ndf ,igood;
   //    for (int ibin =0 ; ibin< n ; ibin++){
   //    x[ibin] = (datahist[ivar][iplot]->GetBinLowEdge(ibin+1) + datahist[ivar][iplot]->GetBinLowEdge(ibin+2))/2;
   //           }
   TH1F *chidata = (TH1F*)datahist[ivar]->Clone("chidata");   //for chi square test
   
   TH1F *chiMC = (TH1F*)MC_hist[iout][ivar]->Clone("chiMC");   //for chi square test
   chiMC->Chi2TestX(chidata, chi2, ndf, igood, "WU", res);
   
   cout << "chi2 value=" << chi2 << endl;
   chi2value[iout][ivar]=chi2;
   
   chi2Y[ivar][iout]=chi2;
   chi2X[iout]=iout;
   chi2Ndfvalue[iout][ivar]=chi2/ndf;
   chiNdf2Y[ivar][iout]=chi2/ndf;
   
   cout <<"TGraph=" << chi2Y[ivar][iout] << endl;
   
   }
   
   }//end of chi2 calculation
 */
 //------------------------------------------------------------chi2/Ndf plot
 
 /*double chi2X_edit[21]={-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20};
   
 //Chi2 plot for with parameter change
 for(int ivar =0 ; ivar < 10 ; ivar++){
 cpt3->cd();
 //TGraph* chi_graph = new TGraph(outnum,chi2X, chi2Y[ivar]);
 TGraph* chi_graph = new TGraph(outnum,chi2X_edit, chiNdf2Y[ivar]);
 chi_graph->Draw("AC*");
 chi_graph->SetLineColor(2);
 sprintf(Title,"#chi^{2}/Ndf value for %s ",var_name[ivar] );
 chi_graph->SetTitle(Title);
 chi_graph->GetYaxis()->SetTitle("#chi^{2}/Ndf");
 chi_graph->GetXaxis()->SetTitle("% of Change for both ISR  #alpha_{s}");
 chi_graph->GetYaxis()->SetTitleOffset(1.4);
 cpt3->Update();
 if(ivar==0){cpt3->Print("chi2Ndf.pdf(","pdf");
 }else if(ivar==9) {cpt3->Print("chi2Ndf.pdf)","pdf");
 }else{
 cpt3->Print("chi2Ndf.pdf","pdf");
 };
 cpt3->Clear();
 }//--------------------------------------------------------chi2/Ndf plot
 */
 
}   // end of main program


//Ratio plot function
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
  
  TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 900,1000 );
  canvas->cd();
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
  data->SetTitle("");
  //    ymax = data->GetMaximum();
  
  TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.35, 1.0, 1.0);
  padfun1->SetBottomMargin(0.01); // Upper and lower plot are joined
  padfun1->SetTopMargin(0.1); // Upper and lowd
  padfun1->Draw();             // Draw the upper pad: pad1
  padfun1->cd();
  data->Draw(" ");
  if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.1);
  
  
  
  
  int  color[21] = {2,6,4,8,46,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,9};
  int   style[21]={1,1,1,1,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
    MC[iup]->SetLineStyle(style[iup]);
    MC[iup]->SetLineColor(color[iup]);
    MC[iup]->SetLineWidth(2);
    MC[iup]->Draw("same hist");// gPad->SetLogy();    
    //  if(fabs(MC[iup]->GetMaximum())>fabs(ymax)){ymax =MC[iup]->GetMaximum();}
    //  MC[iup]->SetMaximum(ymax);
  }
  //  data->Draw("same");
  // TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
  TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
  legendn->SetFillStyle(0);
  legendn->SetBorderSize(0);
  legendn->SetTextSize(plegend[4]);
  legendn->AddEntry(data, "Data 2017","lp");
  const char* legendN[21] ={"Pythia8 ","Madgraph ","Herwig","PY8 CUET PU HLT true","Herwig CUET","#alpha_{s} -10%","#alpha_{s} -8%","#alpha_{s} -6%","#alpha_{s} -4%","#alpha_{s} -2%","Monash tune","#alpha_{s} +2%","#alpha_{s} +4%","#alpha_{s} +6%","#alpha_{s} +8%","#alpha_{s} +10%","#alpha_{s} +12%","#alpha_{s} +14%","#alpha_{s} +16%","#alpha_{s} +18%","#alpha_{s} +20%"}; 
  for(int iup =0 ; iup <  Nplot[0] ; iup++){
    sprintf(MCindex,"%s" ,legendN[iup]);
    legendn->AddEntry(MC[iup], MCindex ,"lp");
  }
  legendn->Draw();
  
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *padfun2 = new TPad("padfun2", "padfun2", 0, 0.1,1.0, 0.35);
  padfun2->SetTopMargin(0);
  padfun2->SetBottomMargin(.35);
  padfun2->SetGridy(); // Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetRightMargin(0.03);
  //   gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    TH1F *rh2;
    if(data->GetBinContent(1) >= 0){
      rh2 = (TH1F*)MC[ilow]->Clone("rh2"); 
      // rh2->Sumw2();
      rh2->Divide(data);     //MC devide by data
    }else{
      rh2 = (TH1F*)data->Clone("rh2");
      //   rh2->Sumw2();
      rh2->SetLineColor(MC[ilow]->GetLineColor());
      rh2->SetLineWidth(MC[ilow]->GetLineWidth());
      rh2->SetLineStyle(MC[ilow]->GetLineStyle());
      rh2->Divide(MC[ilow]);
    }
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
    rh2->Draw("hist same");
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
    //  rh2->GetYaxis()->SetTitleFont(3);
    rh2->GetYaxis()->SetTitleOffset(0.35);
    //  rh2->GetYaxis()->SetLabelFont(1.0); // Absolute font size in pixel (precision 3)
    rh2->GetYaxis()->SetLabelSize(0.09);
  }
  canvas->Update();
  return canvas;
}//end of  ratio plot function


