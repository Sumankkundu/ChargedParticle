#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPoisson.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooBifurGauss.h"

using namespace RooFit;
void effi_fit(){
  
  const int nDiJetHLTmx=8;
  const int njetetamn=1;                     // one eta space is choosen
  
  
  double l1Pt[nDiJetHLTmx] = {0,35,60,90,120,170,170,170};
  double jethlt_thr[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
  double leadingPtThreshold[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
  const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
  double etarange[njetetamn] ={2.4};
  
  const char* dijethlt_name[nDiJetHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};
  
  
  TH1F* Triggered_fired[nDiJetHLTmx][njetetamn];
  TH1F* Total_event[nDiJetHLTmx][njetetamn];
  TH1F* TrigEff[nDiJetHLTmx][njetetamn];
  
  RooDataHist *Roofired[nDiJetHLTmx][njetetamn];
  RooDataHist *RooTotal[nDiJetHLTmx][njetetamn];
  //  RooDataHist *RooEff[nDiJetHLTmx][njetetamn];
  
  char name1[200];
  char name2[200];
  char name3[200];
  
  
  
  
  //TFile *file1 = TFile::Open("Test_QCD_data.root");  // data root file
//  TFile *file1 = TFile::Open("Triggereffi_31Mar2018.root");  // data root file
  TFile *file1 = TFile::Open("Data_2018_Trig_effi_12April.root");  // data root file
  
  for(int jk = 0; jk < njetetamn ; jk++){
    for (int ij=0; ij<nDiJetHLTmx; ij++) {
      
      sprintf(name1, "analyzeBasicPat/hlt_dijet_effi_%i_%i", ij+1, jk); //, jk);
      Triggered_fired[ij][jk] = (TH1F*)file1->Get(name1);
      TrigEff[ij][jk] = (TH1F*)Triggered_fired[ij][jk]->Clone();
      TrigEff[ij][jk]->Rebin(2);
      //  TrigEff[ij][jk]->Sumw2();
      
      sprintf(name2, "analyzeBasicPat/hlt_dijet_all_evt_%i_%i", ij+1, jk);//, jk);
      Total_event[ij][jk] = (TH1F*)file1->Get(name2);
      Total_event[ij][jk]->Rebin(2);
      //    Total_event[ij][jk]->Sumw2();
    }
  }
  
  for(int jk = 0; jk < njetetamn ; jk++){
    for (int ij=0; ij<nDiJetHLTmx; ij++) {
      TrigEff[ij][jk]->Divide(Total_event[ij][jk]);
      
      
      TrigEff[ij][jk]->SetMarkerStyle(20);
      TrigEff[ij][jk]->SetMarkerSize(.55);
      TrigEff[ij][jk]->SetMarkerColor(2);
    }
  }
  
  
  
  // C r e a t e   p d f   f o r   s a m p l i n g
   // ---------------------------------------------
  
    RooPlot* frame[nDiJetHLTmx];
    
    
    double turnon[nDiJetHLTmx];
  //  double minx [nDiJetHLTmx]={24,32,56,80,104,128,160,200};
    double minx [nDiJetHLTmx]={40,40,40,40,40,40,40,40};
//    double maxx [nDiJetHLTmx]={150,200,350,500,650,800,1000,1250};
    double maxx [nDiJetHLTmx]={1240,1240,1240,1240,1240,1240,1240,1240};

    for(int iet = 0; iet < njetetamn ; iet++){
      for (int iht=0; iht<nDiJetHLTmx; iht++) {
	sprintf(name3, "#font[82]{%s  |#eta| #leq 2.4}",dijethlt_name[iht]); //, jk);
	
	
//	RooRealVar integral( "integral", "", 0, 20000 );
	RooRealVar integral( "integral", "", minx[iht], maxx[iht] );
	TH1F* hist = dynamic_cast<TH1F*> (TrigEff[iht][iet]->Clone());
	RooDataHist RooEff("data", "Data from histogram", integral, hist );
	
	RooHistPdf histpdf("histpdf", "histpdf", integral, RooEff, 3);

	// Plot unbinned data and histogram pdf overlaid
	frame[iht] = integral.frame(Title(name3), Bins(100));
        	
	//  RooEff.DrawOption(MarkerStyle(2))
	
	RooEff.plotOn(frame[iht], MarkerStyle(20), MarkerSize(.90), MarkerColor(2)	);
	histpdf.plotOn(frame[iht], LineWidth(2), LineColor(4),MarkerSize(.60));
	
	RooCurve *func = (RooCurve*)frame[iht]->getObject(0);

	double nbinx = hist->GetNbinsX();
	double xmin = hist->GetBinLowEdge(0); 
	double xmax = hist->GetBinLowEdge(nbinx)+hist->GetBinWidth(nbinx);
	
	cout <<" Plot range "  <<xmin << " " << xmax << endl;
	double tlnc = 1;
	double  tnon = 0.0;
	double  fmax = 0.0;
	for(int ix= (int) xmin ; ix <=xmax ; ix++){
	  double value1 = func->Eval(ix);
		 if(value1 > fmax){fmax=value1;} 
	} 
	cout << " Func max : " << fmax << endl;
	
	for(int ix= (int) xmin ; ; ix ++){
	double ixx=.01*ix ; 
	double value1 = func->Eval(ixx);
	//	 cout <<"HT : "<< ixx << " Fun : " << value1 <<endl;
	  if(value1 >= .99*fmax){tnon= ixx;
	    break;}
	if(ixx >= xmax) break;
	}
	
	turnon[iht] = tnon;
//	cout<< "Tunr on for " << jethlt_thr[iht] <<" is  " <<  turnon[iht] << endl;
      }
      
    }

     for (int ij=0; ij<nDiJetHLTmx; ij++) {
            cout << "turn on for " << jethlt_thr[ij] << " is  "  << turnon[ij] << endl;
    }

 //    int turnonI[nDiJetHLTmx]
//for (int ij=0; ij<nDiJetHLTmx; ij++) {
  //          turnonI[ij]= turnon[ij];
  //  }


    
    TCanvas *c1 ;
    TCanvas *c3 ;
    char canvname[100] ;
    char plotname[100] ;
    
    gStyle->SetPadGridX(3);
    gStyle->SetPadGridY(3);
    gStyle->SetGridStyle(3);
    gStyle->SetStatX(.9);
    gStyle->SetStatY(.9);
    //    gStyle->SetTitleW(0.5);
//    gStyle->SetTitleH(0.2);
//    gStyle->SetTitleX(0.55);
//    gStyle->SetTitleY(0.4);
    
    double pminx[nDiJetHLTmx]={40,50,110,170,230,280,350,450};
    double pmaxx[nDiJetHLTmx]={300,300,500,500,500,800,800,800};

    
    sprintf(canvname,"Trig_Eff_with_rootfit") ;
    c1= new TCanvas(canvname,canvname,1600,1600) ;
 //   c2= new TCanvas("canvname1","canvname1",1200,900) ;
    c3= new TCanvas("canvname3","canvname3",1200,900) ;
    c1->Divide(2,4); 
 //   c2->Divide(2,4); 
    
    for (int iht=0; iht<nDiJetHLTmx; iht++) {
      c1->cd(iht+1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.3);
      frame[iht]->GetXaxis()->SetTickLength(0.03);
      frame[iht]->GetYaxis()->SetTitleOffset(1);
      frame[iht]->GetXaxis()->SetTitleOffset(2.5);
      frame[iht]->GetYaxis()->SetLabelSize(0.04);
      frame[iht]->GetXaxis()->SetLabelSize(0.04);
      frame[iht]->GetYaxis()->SetTitleOffset(1);
      frame[iht]->GetYaxis()->SetTitle("Trigger Efficiency");
      frame[iht]->GetXaxis()->SetTitle("H_{T,2} in GeV");   
      frame[iht]->GetXaxis()->CenterTitle();
      frame[iht]->GetYaxis()->CenterTitle();
      frame[iht]->GetXaxis()->SetTitleSize(.03);
      frame[iht]->GetYaxis()->SetTitleSize(.045);
   //hpxpy->GetXaxis()->SetRange(18,19);

      frame[iht]->Draw();
      
    }
   
    TCanvas *c2 ;
    char canvname1[100] ;
    sprintf(canvname1,"Trig_Eff_page") ;
    c2= new TCanvas("canvname1","canvname1",1200,900) ;
    double ndiv[nDiJetHLTmx]={512,512,1010,1012,1212,1212,1212,1212};
    

 for (int iht=0; iht<nDiJetHLTmx; iht++) {
      
	 c2->Clear();
	 c2->cd();
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      frame[iht]->GetXaxis()->SetTickLength(0.03);
      frame[iht]->GetYaxis()->SetTitleOffset(1);
      frame[iht]->GetXaxis()->SetTitleOffset(1.2);
      frame[iht]->GetYaxis()->SetLabelSize(0.04);
      frame[iht]->GetXaxis()->SetLabelSize(0.04);
      frame[iht]->GetYaxis()->SetTitleOffset(1);
      frame[iht]->GetYaxis()->SetTitle("Trigger Efficiency");
      frame[iht]->GetXaxis()->SetTitle("H_{T,2} in GeV");
      frame[iht]->GetXaxis()->CenterTitle();
      frame[iht]->GetYaxis()->CenterTitle();
      frame[iht]->GetXaxis()->SetTitleSize(.045);
      frame[iht]->GetYaxis()->SetTitleSize(.045);
     
      frame[iht]->GetXaxis()->SetNdivisions(ndiv[iht]);
//  if(iht==0){frame[iht]->GetXaxis()->SetRange(10,10000);}
 if(iht==0){frame[iht]->GetXaxis()->SetLimits(40,120);}
//    frame[iht]->GetXaxis()->SetLimits(40,1200);
     	// frame[iht]->Draw("SAME");
     	 frame[iht]->Draw("");

    char ledtxt1[100];
    char ledtxt2[100];
    sprintf(ledtxt1,"HLT_DiPFJetAve_%0.0f_v", jethlt_thr[iht]);
    sprintf(ledtxt2,"Trun on point : %0.2f GeV ", roundf(turnon[iht])) ;

    TLegend *legendn = new TLegend(.4,0.18,0.62,0.30);
    legendn->SetFillStyle(0);
    legendn->SetBorderSize(0);
    legendn->SetTextSize(.04);
//    legendn->AddEntry("Turn on");
    legendn->AddEntry((TObject*)0, ledtxt1, "");
    legendn->AddEntry((TObject*)0, ledtxt2, "");
    legendn->Draw();

    
    sprintf(plotname,"HLT_DiPFJetAve_%0.0f.pdf", jethlt_thr[iht]);
    c2->Print(plotname,"pdf");
    if(iht==0){c2->Print("DiJetHLT_effi.pdf(","pdf");
    }else if(iht==7) {c2->Print("DiJetHLT_effi.pdf)","pdf");
    }else{
      c2->Print("DiJetHLT_effi.pdf","pdf");};

 }
 
/*    for (int iht=0; iht<nDiJetHLTmx; iht++) {

      gPad->SetLeftMargin(0.15);
      frame[iht]->GetYaxis()->SetTitleOffset(1.2);
      frame[iht]->GetXaxis()->SetTitleOffset(1.2);
      frame[iht]->GetYaxis()->SetLabelSize(0.03);
      frame[iht]->GetXaxis()->SetLabelSize(0.03);
      frame[iht]->SetTitleOffset(1.8);
      //c2->cd(iht+1);
      
      frame[iht]->Draw("SAME");
      
    }*/
 
 TLegend *legendn3 = new TLegend(.7,0.18,0.8,0.5); 
 TLegend *year = new TLegend(.7,.7,0.8,0.5); 
 year->AddEntry((TObject*)0, "2018 Data", "");
 year->SetFillStyle(0);
 year->SetBorderSize(0);
 year->SetTextSize(.03);

 int mkcl[8]={2,3,4,46,6,30,91,1};
for(int jk = 0; jk < njetetamn ; jk++){
    for (int ij=0; ij<nDiJetHLTmx; ij++) {
      gStyle->SetOptStat(0);
      gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      gPad->SetLeftMargin(0.1);
       gPad->SetTickx();
      gPad->SetTicky();
      TrigEff[ij][jk]->GetYaxis()->SetTitleOffset(1.2);
      TrigEff[ij][jk]->GetXaxis()->SetTitleOffset(1.2);
      TrigEff[ij][jk]->GetYaxis()->SetLabelSize(0.03);
      TrigEff[ij][jk]->GetXaxis()->SetLabelSize(0.03);
     // TrigEff[ij][jk]->GetXaxis()->SetLimits(40,700);
      TrigEff[ij][jk]->GetYaxis()->SetTitle("Trigger Efficiency");
      TrigEff[ij][jk]->GetXaxis()->SetTitle("H_{T,2} in GeV");
            TrigEff[ij][jk]->SetMarkerSize(.70);
            TrigEff[ij][jk]->SetMarkerColor(mkcl[ij]);
      TrigEff[ij][jk]->SetTitleOffset(0);
   //   TrigEff[ij][jk]->Rebin(2);
       TrigEff[ij][jk]->SetTitle("");
  //     gStyle->SetErrorY(0.00001);
      c3->cd();

      char ledtxt3[100];
    sprintf(ledtxt3,"HLT_DiPFJetAve_%0.0f", jethlt_thr[ij]);
       legendn3->SetFillStyle(0);
    legendn3->SetBorderSize(0);
    legendn3->SetTextSize(.03);
//    legendn3->SetMarkerSize(20);
//    legendn->AddEntry("Turn on");
    legendn3->AddEntry(TrigEff[ij][jk], ledtxt3, "lp");


       TrigEff[ij][jk]->Draw("same ");

    }
    gStyle->SetMarkerSize(1.3);
    year->Draw("same");
    legendn3->Draw( "same");
}
    
}
