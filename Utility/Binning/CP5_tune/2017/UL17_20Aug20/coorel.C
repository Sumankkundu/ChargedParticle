//To plots The correlation plots for the MC and data

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

void coorel(){
  
  
  TH2F *coorel_hist_jet[5][8];
  TH2F *coorel_hist_char[5][8];
  Int_t var[5]={3,9,15,18,24};
  
  //  TFile *file1 = TFile::Open("/home/suman/Paradox/ESV_PYTHIA/DATA_4v/qcd_evt_v4.root");  // data root file
  
  char histname1[100];
  char histname2[100];
  TCanvas *cpt0 = new TCanvas("cpt0", "canvas", 1500,900 );  //for ESVs
  TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 1500,900 );  //for ESVs
  cpt0->SetBorderSize(0);
  cpt0->SetRightMargin(0.1);
  cpt1->SetRightMargin(0.1);
  cpt0->SetLeftMargin(0.08);
  cpt1->SetLeftMargin(0.08);
//  cpt0->SetTopMargin(0.0);
 
  gStyle->SetPaintTextFormat( "4.2f");
  //for Jets
  for(int ivar =0 ; ivar < 5 ; ivar++){
    for(Int_t ipt =0; ipt < 8 ; ipt++){   
      sprintf(histname1, "analyzeBasicPat/corr_typ_0_pt%i_eta0_%i",ipt,var[ivar]);  //corr_typ_1_pt7_eta0_24
      cout << histname1 << endl;
      coorel_hist_jet[ivar][ipt] = (TH2F*)gDirectory->Get(histname1);
      coorel_hist_jet[ivar][ipt]->SetStats(0);
      // coorel_hist_jet[ivar][ipt] = (TH1F*) coorel_hist_jet->Get(histname1);
      coorel_hist_jet[ivar][ipt]->GetZaxis()->SetTitleOffset(1.0);
      coorel_hist_jet[ivar][ipt]->GetXaxis()->SetTitle("Reco");
      coorel_hist_jet[ivar][ipt]->GetYaxis()->SetTitle("Gen");
        cpt0->cd();
      coorel_hist_jet[ivar][ipt]->Draw("colz");      
      //coorel_hist_jet[ivar][ipt]->Draw("colz text");      
     // coorel_hist_jet[ivar][ipt]->Draw("");      
      cpt0->Update();
      
      if(ivar==0 && ipt ==0){cpt0->Print("coorelation_Jet.pdf(","pdf");
      }else if(ivar==4 && ipt==7) {cpt0->Print("coorelation_Jet.pdf)","pdf");
      }else{
	cpt0->Print("coorelation_Jet.pdf","pdf");};
//	cpt0->Print("coorelation_Jet1.pdf","pdf");};
      
      cpt0->Clear();
      
    }
  }
  
  //For Charge Particles
  for(int ivar =0 ; ivar < 5 ; ivar++){
    for(Int_t ipt =0; ipt < 8 ; ipt++){                
      sprintf(histname2, "analyzeBasicPat/corr_typ_1_pt%i_eta0_%i",ipt,var[ivar]);  //corr_typ_1_pt7_eta0_24
      coorel_hist_char[ivar][ipt] = (TH2F*)gDirectory->Get(histname2);
      coorel_hist_char[ivar][ipt]->SetStats(0);
      coorel_hist_char[ivar][ipt]->GetXaxis()->SetTitle("Reco");
      coorel_hist_char[ivar][ipt]->GetYaxis()->SetTitle("Gen");
      
      cpt1->cd();
    //  coorel_hist_char[ivar][ipt]->Draw("colz text");
     coorel_hist_char[ivar][ipt]->Draw("colz");
     //  coorel_hist_char[ivar][ipt]->Draw("");
      
      cpt1->Update();
      
      if(ivar==0 && ipt ==0){cpt1->Print("coorelation_char.pdf(","pdf");
      }else if(ivar==4 && ipt==7) {cpt1->Print("coorelation_char.pdf)","pdf");
      }else{
//	cpt1->Print("coorelation_char1.pdf","pdf");};
//
	cpt1->Print("coorelation_char.pdf","pdf");};
      
      cpt1->Clear();
      
    }
  }
  
  
  
  
}
