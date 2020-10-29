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
#include "CMS_lumi.h"
//#include "/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/anal_jetqcd_run2.C"
void mean()
{
//TFile *file0 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/unfold/outfile13Tev_10_3_Madgh_1_c0_e13_130529_kcov.root");
//TFile *file0 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/unfold/outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.root");
//TFile *file0 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/RECO/qcdevt_byPythia8_outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov_v3.root");
TFile *file0 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/RECO/qcd_evt_by_pythia8_eta0_eta1_outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.root");
TFile *file1 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/RECO/qcd_evt_by_pythia8_eta0_eta1_outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.root");
   if (file0->IsOpen()){
     cout<<"File1 opened successfully\n"<<endl;
   }

 //sprintf(namey, "outfile%s.root", titley);
  TFile* fileOut = new TFile("mean.root", "recreate");

  int strtiter[5][3][8]={{{4,4,4,4,5,4,4,3},
                        {3,4,5,5,4,5,5,7},
                        {4,4,4,4,4,5,5,4}},// thrustc
                       {{4,4,4,4,4,4,4,4},
                        {4,4,4,4,4,6,4,5},
                        {5,6,6,5,6,6,6,5}},//t3mass
                       {{4,7,7,4,5,5,5,4},
                        {3,5,5,6,7,5,4,7},
                        {6,6,7,7,6,5,8,5}},//broadt
                       {{4,5,4,5,4,5,5,5},
                        {4,5,5,7,8,8,8,7},
                        {5,5,4,5,4,5,5,4}}};//ttmass

 const int nvarx=4;
  const int varindex[nvarx+1]={3,9,18,24};
 const int ntyp=4;
  const int njetptmn=8;
  const int njetetamn=2;
  double mean[100][100][100];
  double mean1[100][100][100];
  double mean2[100][100][100];
  double mean3[100][100][100];
  double poserr[100];
  double negerr[100];
  char histname[100];
  char histname1[100];
  char histname2[100];
 char namey[100];
 const char* titlex= "madgh_weight_logfile";
 char titley[100];
const int n = 9;
 double x[n]={73,93,165,225,298,365,452,557,800};
  TH1F* fdfx[ntyp][njetetamn][nvarx];
  TH1F* fdfx1[ntyp][njetetamn][nvarx];
  TH1F* fdfx2[ntyp][njetetamn][nvarx];
  TH1F* fdfx3[ntyp][njetetamn][nvarx];
  char name[100];
  char title[100];
 for (int ityp=0; ityp<ntyp; ityp++) {
    //for (int ipt=0; ipt<njetptmn; ipt++) {
      for (int iet=1; iet<njetetamn; iet++) {
                                for (int ij=0; ij<nvarx; ij++) {
                                                        sprintf(name, "gen_mean_unolded_bayes_typ_%i_pt_eta%i_%i", ityp, iet, varindex[ij]);
                                                        sprintf(title, "genmean unolded %i %i %i", ityp, iet, varindex[ij]);
                                                        //fdfx[ityp][iet][ij] = new TH1F(name, title, 8, 100,1500);
                                                        fdfx[ityp][iet][ij] = new TH1F(name, title, n-1,x);
                                                        fdfx[ityp][iet][ij]->Sumw2();
                                                        
                                                        sprintf(name, "gen_mean_Pythia8_typ_%i_pt_eta%i_%i", ityp, iet, varindex[ij]);
                                                        sprintf(title, "genmean Pythia8 %i %i %i", ityp, iet, varindex[ij]);
                                                        //fdfx[ityp][iet][ij] = new TH1F(name, title, 8, 100,1500);
                                                        fdfx1[ityp][iet][ij] = new TH1F(name, title, n-1,x);
                                                        fdfx1[ityp][iet][ij]->Sumw2(); 

                                                        sprintf(name, "gen_mean_Madgh_typ_%i_pt_eta%i_%i", ityp, iet, varindex[ij]);
                                                        sprintf(title, "genmean Madgh %i %i %i", ityp, iet, varindex[ij]);
                                                        //fdfx[ityp][iet][ij] = new TH1F(name, title, 8, 100,1500);
                                                        fdfx2[ityp][iet][ij] = new TH1F(name, title, n-1,x);
                                                        fdfx2[ityp][iet][ij]->Sumw2();

                                                        sprintf(name, "gen_mean_Herwigpp_typ_%i_pt_eta%i_%i", ityp, iet, varindex[ij]);
                                                        sprintf(title, "genmean Herwigpp %i %i %i", ityp, iet, varindex[ij]);
                                                        //fdfx[ityp][iet][ij] = new TH1F(name, title, 8, 100,1500);
                                                        fdfx3[ityp][iet][ij] = new TH1F(name, title, n-1,x);
                                                        fdfx3[ityp][iet][ij]->Sumw2(); 


                                                } 
         }
      // }
     }
  char a;
  //Double_t x[n]={73,93,165,225,298,365,452,557};

/*Roo_Pythia8_gen_typ_2_pt7_eta0_24;1	Roo_Pythia8_Gen Charged Particles 557 3 #rho^{T}_{Tot,C} 
  KEY: TH1D	Roo_Madgh_gen_typ_2_pt7_eta0_24;1	Roo_Madgh_Gen Charged Particles 557 3 #rho^{T}_{Tot,C} 
  KEY: TH1D	Roo_Herwigpp_gen_typ_2_pt7_eta0_24
*/

for (int klx=0; klx<nvarx; klx++) {
//for (int klx=2; klx<3; klx++) {
    int kl=varindex[klx];
    for (int ij=0; ij<ntyp-1; ij++){
    //for (int ij=0; ij<1; ij++){
      for(int jk=0; jk<njetptmn; jk++){
      //for(int jk=6; jk<7; jk++){
          for (int mn=1;mn<njetetamn;mn++) {
            fileOut->cd();
            sprintf(histname,"bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_typ_%i_pt%i_eta%i_%i", strtiter[klx][ij][jk], ij, jk, mn, kl);
            cout << "Name = " << histname << endl;
            if(kl==18 && ij==0 && jk==6) TH1F* default_hist =(TH1F*) file1->Get(histname);
            else TH1F* default_hist =(TH1F*) file0->Get(histname);
           // cout << " ht = " << jk << " Mean = " <<default_hist->GetNbinsX() <<endl;
            mean[klx][ij][jk]=(default_hist->GetMean());


           sprintf(histname,"Roo_Pythia8_gen_typ_%i_pt%i_eta%i_%i", ij, jk, mn, kl);
            //if(kl==18 && ij==0 && jk==6) TH1F* default_hist =(TH1F*) file1->Get(histname);
            TH1F* default_hist1 =(TH1F*) file0->Get(histname);
           // cout << default_hist1->GetNbinsX() <<endl;
            mean1[klx][ij][jk]=(default_hist1->GetMean());
  
           sprintf(histname,"Roo_Madgh_gen_typ_%i_pt%i_eta%i_%i", ij, jk, mn, kl);
            //if(kl==18 && ij==0 && jk==6) TH1F* default_hist =(TH1F*) file1->Get(histname);
            TH1F* default_hist2 =(TH1F*) file0->Get(histname);
            //cout << default_hist2->GetNbinsX() <<endl;
            mean2[klx][ij][jk]=(default_hist2->GetMean());

           sprintf(histname,"Roo_Herwigpp_gen_typ_%i_pt%i_eta%i_%i", ij, jk, mn, kl);
            //if(kl==18 && ij==0 && jk==6) TH1F* default_hist =(TH1F*) file1->Get(histname);
            TH1F* default_hist3 =(TH1F*) file0->Get(histname);
            //cout << default_hist3->GetNbinsX() <<endl;
            mean3[klx][ij][jk]=(default_hist3->GetMean());
            //default_hist->Write();
          }
        }
      }
    }

/*for(ij=0; ij<8; ij++){
 cout << mean[ij] <<endl;
}*/

for (int klx=0; klx<nvarx; klx++) {
//for (int klx=2; klx<3; klx++) {
    int kl=varindex[klx];
    for (int ij=0; ij<ntyp-1; ij++){
    //for (int ij=0; ij<1; ij++){
      for(int jk=0; jk<njetptmn; jk++){
      //for(int jk=6; jk<7; jk++){
          for (int mn=1;mn<njetetamn;mn++) {
          cout << "jk = " <<jk << " ;Mean="<<mean[klx][ij][jk] << endl;
          //fdfx[ij][mn][klx]->SetBinContent(mean[klx][ij][jk], x[jk]);
          //fdfx[ij][mn][klx]->SetBinContent(jk, sqrt(mean[klx][ij][jk]*mean[klx][ij][jk]));
          fdfx[ij][mn][klx]->SetBinContent(jk+1, mean[klx][ij][jk]);
          if(jk==7) fdfx[ij][mn][klx]->Write();

          fdfx1[ij][mn][klx]->SetBinContent(jk+1, mean1[klx][ij][jk]);
          if(jk==7) fdfx1[ij][mn][klx]->Write();
 
         fdfx2[ij][mn][klx]->SetBinContent(jk+1, mean2[klx][ij][jk]);
          if(jk==7) fdfx2[ij][mn][klx]->Write();
 
         fdfx3[ij][mn][klx]->SetBinContent(jk+1, mean3[klx][ij][jk]);
          if(jk==7) fdfx3[ij][mn][klx]->Write();

         }
        }
      }
    }
delete file0;

fileOut->Close();
}


void mean_plot(const char* title="thrustc", const char* etapt="c3_e0", const int typindex=1 )
{
TFile *file1 = TFile::Open("mean.root");
char histname[100];
  char histname1[100];
  char histname2[100];
char typname[100];
char fname[100];
  int iPeriod=4;
  int iPos=0;
  TH1F* fdfx;
  TH1F* fdfx1;
  TH1F* fdfx2;
  TH1F* fdfx3;
const char* typ[3]={"Jet", "All Particles", "Charged Particles"};
  sprintf(histname, "gen_mean_unolded_bayes_%s", etapt);
  default_hist =(TH1F*) file1->Get(histname);
          sprintf(histname, "gen_mean_Pythia8_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
  default_hist1 =(TH1F*) file1->Get(histname);
          sprintf(histname, "gen_mean_Madgh_%s", etapt);
  default_hist2 =(TH1F*) file1->Get(histname);
         sprintf(histname, "gen_mean_Herwigpp_%s", etapt);
  default_hist3 =(TH1F*) file1->Get(histname);

  // Define the Canvas
   TCanvas *c1 = new TCanvas("c1", "canvas", 600, 800);

// Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 0.92);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetTopMargin(0); // Upper and lower plot are joined
   pad1->SetLeftMargin(0.15);
   //pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   default_hist->SetLineColor(1);
   default_hist->SetTitle(" ");
   default_hist->GetYaxis()->CenterTitle();
   if (strstr(title,"thrustc")) default_hist->GetYaxis()->SetTitle(" <log_{e} (#tau_{_{#perp} _{   ,C}})>");
   if (strstr(title,"rhotot")) default_hist->GetYaxis()->SetTitle("<log_{e} (#rho_{Tot,C})>");
   if (strstr(title,"broadt")) default_hist->GetYaxis()->SetTitle("<log_{e} (B_{ T,C})>");
   if (strstr(title,"rhottot")) default_hist->GetYaxis()->SetTitle("<log_{e} (#rho^{T}_{Tot,C})>");

   default_hist->GetYaxis()->SetLabelSize(0.05);
   default_hist->GetYaxis()->SetTitleSize(0.060);
   default_hist->GetYaxis()->SetTitleOffset(1.1);
   default_hist1->SetLineColor(2);
   default_hist2->SetLineColor(3);
   default_hist3->SetLineColor(4);
  
   default_hist1->SetLineStyle(2);
   default_hist2->SetLineStyle(2);
   default_hist3->SetLineStyle(2);

   default_hist->SetLineWidth(3);
   default_hist1->SetLineWidth(3);
   default_hist2->SetLineWidth(3);
   default_hist3->SetLineWidth(3);
 
   //default_hist->Scale(1./(default_hist->Integral()));
/*   default_hist1->Scale(1./(default_hist->Integral()));
   default_hist2->Scale(1./(default_hist->Integral()));
   default_hist3->Scale(1./(default_hist->Integral()));
    
   default_hist->Scale(1./(default_hist->Integral()));*/
   default_hist->SetStats(0);          // No statistics on upper plot
   default_hist->Draw();               // Draw h1
   default_hist1->Draw("same");         // Draw h2 on top of h1
   default_hist2->Draw("same");         // Draw h2 on top of h1
   default_hist3->Draw("same");         // Draw h2 on top of h1
TLegend *legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->AddEntry(default_hist, "Data","lp");
  legend->AddEntry(default_hist1, "Pythia8","lp");
  legend->AddEntry(default_hist2, "Madgraph","lp");
  legend->AddEntry(default_hist3, "Herwigpp","lp");
  legend->Draw();
  if(typindex==0) sprintf(typname,"%s",typ[0]);
    if(typindex==1) sprintf(typname,"%s",typ[1]);
    if(typindex==2) sprintf(typname,"%s",typ[2]);

    TPaveText *ptstt = new TPaveText(.7,0.4,0.8,0.5,"brNDC");
      ptstt->SetFillColor(10);
    ptstt->SetBorderSize(0);
    ptstt->SetTextFont(42);

    TText* text = ptstt->AddText(typname);
      text->SetTextSize(0.06);
    ptstt->Draw();
CMS_lumi(c1,iPeriod,iPos);
   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   //default_hist->GetYaxis()->SetLabelSize(0.);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(15);
   //axis->Draw();

   // lower plot will be in pad
   c1->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0.02);
   pad2->SetBottomMargin(0.24);
   pad2->SetLeftMargin(0.15);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

   // Define the ratio plot
   TH1F *h = (TH1F*)default_hist->Clone("h");
   TH1F *h1 = (TH1F*)default_hist1->Clone("h1");
   TH1F *h2 = (TH1F*)default_hist2->Clone("h2");
   TH1F *h3 = (TH1F*)default_hist3->Clone("h3");
   //h->SetLineColor(kBlack);
   h1->SetMinimum(0.5);  // Define Y ..
   h1->SetMaximum(1.5); // .. range
   //h->Sumw2();
   //h->SetMarkerStyle(20);
   /*h1->SetMarkerColor(2);
   h1->SetMarkerStyle(20);
   h2->SetMarkerColor(3);
   h2->SetMarkerStyle(20);
   h3->SetMarkerColor(4);
   h3->SetMarkerStyle(20);*/
    

  int yyv=1;
  int alow=73;
  int ahg=800; 
  TLine lnn(alow, yyv, ahg+.5, yyv);
//        ratHist[ij]->GetXaxis()->GetXmin();

          TLine lnn();
          lnn.SetLineColor(1);
          lnn.SetLineWidth(1.2);
          lnn.SetLineStyle(0);

          //lnn.DrawLine(alow, 0, ahg+.5, 0);
                //lnn.DrawLine(h1->GetXaxis()->GetXmin(), yyv, h1->GetXaxis()->GetXmax(), yyv);
             cout << "RATTT HIST X MIN ======" << h1->GetXaxis()->GetXmin() << endl;


   TH1F *div = (TH1F*)default_hist->Clone("h");
   h1->SetStats(0);      // No statistics on lower plot
   div->SetLineStyle(2);
   div->SetLineWidth(3);
  // h1->Divide(div);
   div->Divide(h1);
   div->SetMinimum(0.9);  // Define Y ..
   div->SetMaximum(1.12);
   div->SetLineColor(2);
   //h1->Draw("");       // Draw the ratio plot
   div->SetStats(0);      // No statistics on lower plot
   div->Draw(" ");       // Draw the ratio plot
   //div->Reset();
   TH1F *div1 = (TH1F*)default_hist->Clone("h");
   div1->SetLineStyle(2);
   div1->SetLineWidth(3);
   div1->Divide(h2);
   div1->SetLineColor(3);
   div1->SetStats(0);      // No statistics on lower plot
  // h2->Divide(h);
   //h2->Draw("same");       // Draw the ratio plot
   div1->Draw("same");       // Draw the ratio plot
   //div->Reset();
   TH1F *div2 = (TH1F*)default_hist->Clone("h");
   div2->SetLineStyle(2);
   div2->SetLineWidth(3);
   div2->SetStats(0);      // No statistics on lower plot
   div2->Divide(h3);
   div2->SetLineColor(4);
  // h3->Divide(h);
   div2->Draw("same");
   //h3->Draw("same");       // Draw the ratio plot
   //h->SetMarkerStyle(21);
   //h3->Draw("ep");       // Draw the ratio plot
   // h1 settings
  // h1->SetLineColor(kBlue+1);
   //h1->SetLineWidth(2);

                lnn.DrawLine(h1->GetXaxis()->GetXmin(), yyv, h1->GetXaxis()->GetXmax(), yyv);
   // Y axis h1 plot settings
   h->GetYaxis()->SetTitleSize(20);
   h->GetYaxis()->SetTitleFont(43);
   h->GetYaxis()->SetTitleOffset(1.75);

   // h2 settings
   //h2->SetLineColor(kRed);
   //h2->SetLineWidth(2);


   // Ratio plot (h3) settings
   div->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   div->GetYaxis()->SetTitle("MC/Data");
   div->GetYaxis()->CenterTitle();
   div->GetYaxis()->SetNdivisions(505);
   div->GetYaxis()->SetTitleSize(20);
   div->GetYaxis()->SetTitleFont(43);
   div->GetYaxis()->SetTitleOffset(2.3);
   div->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   div->GetYaxis()->SetLabelSize(25);

   // X axis ratio plot settings
   div->GetXaxis()->SetTitle("H_{T,2}");
   div->GetXaxis()->CenterTitle();
   div->GetXaxis()->SetTitleSize(22);
   div->GetXaxis()->SetTitleFont(43);
   div->GetXaxis()->SetTitleOffset(3.5);
   div->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   div->GetXaxis()->SetLabelSize(20);
   sprintf(fname, "%s.pdf",etapt);
      c1->SaveAs(fname); 

}

void lumitest(){
TCanvas *canv = new TCanvas("c1", "runfile", 700., 400.);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
  int iPeriod=4;
  int iPos=33;
CMS_lumi( canv, iPeriod, iPos );

}
void CMS_lumi(TPad* pad, int iPeriod, int iPosX){
  bool outOfFrame    = false;
  if( iPosX/10==0 )
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();
TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 )
    {
      lumiText = lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
else if ( iPeriod==7 )
    {
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }
  else if ( iPeriod==0 )
    {
      lumiText += lumi_sqrtS;
    }

  std::cout << lumiText << endl;
TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(lumiTextSize*t);
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      cout << "OUT OF FRAME ====== "<< endl;
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11);
      latex.SetTextFont(42);
  latex.SetTextAlign(31);
      latex.SetTextSize(cmsTextSize*t);
      //latex.DrawLatex(0.47-l,1-t+lumiTextOffset*t,cmsText); //all other plot
      latex.DrawLatex(0.36-l,1-t+lumiTextOffset*t,cmsText); //only for delta
     // latex.DrawLatex(1-r,1-t+lumiTextOffset*t,cmsText);
    }

  pad->cd();
float posX_=0;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
        {
          posX_ =   l + 0.045*(1-l-r)*W/H;
          posY_ = 1-t - 0.045*(1-t-b);
          float xl_0 = posX_;
          float yl_0 = posY_ - 0.15;
          float xl_1 = posX_ + 0.15*H/W;
          float yl_1 = posY_;
          TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
          TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
          pad_logo->Draw();
          pad_logo->cd();
          CMS_logo->Draw("X");
          pad_logo->Modified();
          pad->cd();
        }
else
        {
          latex.SetTextFont(cmsTextFont);
          latex.SetTextSize(cmsTextSize*t);
          latex.SetTextAlign(align_);
          latex.DrawLatex(posX_, posY_, cmsText);
          if( writeExtraText )
            {
              latex.SetTextFont(extraTextFont);
              latex.SetTextAlign(align_);
              latex.SetTextSize(extraTextSize*t);
              latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
            }
        }
    }
  else if( writeExtraText )
    {
      if( iPosX==0)
        {
          posX_ =   l +  relPosX*(1-l-r);
          posY_ =   1-t+lumiTextOffset*t;
        }
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
     // latex.DrawLatex(0.51-posX_, posY_, extraText);  //on for all plots  
      latex.DrawLatex(0.40-posX_, posY_, extraText);   // only for delta    
    }
  return;
}
