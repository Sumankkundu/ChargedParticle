#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
//#include "LHAPDF/LHAPDF.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>

#include <TRandom.h>

#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include <string.h>
#include <stdlib.h>
#include <Riostream.h>
#include <stddef.h>
#include "TROOT.h"
#include "TApplication.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"



using namespace std;
//using namespace CLHEP;



int main() { //int main(argc, char**argv)
  gSystem->Load("libTree");
  const int njetetamn=3;
  const int ntype=4;
  const int njetptmn=5;
  const int nvar=32;
  
  char outfile[100];
  char outfilx[100];
  char infile[200];
  char rootfiles[100];
  char datafile[100];
  char name[200];
  char title[100];
  char title2[100];

  ofstream file_out("output_test.txt");

  cout <<"Give the input file name"<<endl;;
  cin>> rootfiles;
  
  double xbins[400]={0};
  double ybins[400]={0};
  int len = strlen(rootfiles);
  strncpy(outfilx, rootfiles, len-4);
  outfilx[len-4]='\0';
  sprintf (outfile,"%s.root",outfilx);

  TFile* fileOut = new TFile(outfile, "recreate");


  TH1F* h_recoevtvar[17];

  TH1F* h_genevtvar[17];

  ifstream file_db;
  file_db.open(rootfiles);   
  float weight2=1.0, weight=1.0;
  int nfile = 0;
  const TH1F* hist1d;
  const TH2F* hist2d;
  
  while(!(file_db.eof())) {
    file_db >> datafile >> weight2; // >> qlow>> qhigh;
    if (strstr(datafile,"#")) continue;
    if(file_db.eof()) break;
    sprintf(infile, "%s", datafile);
    nfile++;
    TFile* fileIn = new TFile(infile, "read");
  for(int iplot =0; iplot <27; iplot++){     
    //for (int iplot=0; iet<njetetamn-2; iet++) {
      //cout <<"datafile "<< datafile<< ""<<weight2<<" "<<iplot<<endl;
      fileIn->cd();
      if(iplot == 0) {sprintf(name, "analyzeBasicPat/recojetallave_pt_0");}
      if(iplot == 1) {sprintf(name, "analyzeBasicPat/recojet1_pt_0");}
      if(iplot == 2) {sprintf(name, "analyzeBasicPat/recojet2_pt_0");}
      if(iplot == 3) {sprintf(name, "analyzeBasicPat/recojet1_eta");}
      if(iplot == 4) {sprintf(name, "analyzeBasicPat/recojet2_eta");}
      if(iplot == 5) {sprintf(name, "analyzeBasicPat/recojet1_phi");}
      if(iplot == 6) {sprintf(name, "analyzeBasicPat/recojet2_phi");}
      if(iplot == 7) {sprintf(name, "analyzeBasicPat/recojt_eta");}
      if(iplot == 8) {sprintf(name, "analyzeBasicPat/recojt_phi");}
      if(iplot == 9) {sprintf(name, "analyzeBasicPat/hjetdphi_0");}
      if(iplot == 10) {sprintf(name, "analyzeBasicPat/hjetptbypl_0");}
      if(iplot == 11) {sprintf(name, "analyzeBasicPat/hjetdpt_0");}
      if(iplot == 12) {sprintf(name, "analyzeBasicPat/njets_0");}
      if(iplot == 13) {sprintf(name, "analyzeBasicPat/recojetHT2_0");}


      if(iplot == 14) {sprintf(name, "analyzeBasicPat/genjetallave_pt_1");}
      if(iplot == 15) {sprintf(name, "analyzeBasicPat/genjet1_pt_1");}
      if(iplot == 16) {sprintf(name, "analyzeBasicPat/genjet2_pt_1");}
      if(iplot == 17) {sprintf(name, "analyzeBasicPat/genjet1_eta");}
      if(iplot == 18) {sprintf(name, "analyzeBasicPat/genjet2_eta");}
      if(iplot == 19) {sprintf(name, "analyzeBasicPat/genjet1_phi");}
      if(iplot == 20) {sprintf(name, "analyzeBasicPat/genjet2_phi");}
      if(iplot == 21) {sprintf(name, "analyzeBasicPat/genjt_eta");}
      if(iplot == 22) {sprintf(name, "analyzeBasicPat/genjt_phi");}
      if(iplot == 23) {sprintf(name, "analyzeBasicPat/genjetdphi_1");}
      if(iplot == 24) {sprintf(name, "analyzeBasicPat/genjetptbypl_1");}
      if(iplot == 25) {sprintf(name, "analyzeBasicPat/genjetdpt_1");}
      if(iplot == 26) {sprintf(name, "analyzeBasicPat/gennjets_1");}
    //  if(iplot == 13) {sprintf(name, "analyzeBasicPat/recojetHT2_0");}
  hist1d = (TH1F*)fileIn->Get(name);
      fileOut->cd();
      cout<< "NAME= "<< name <<endl;
      if (nfile==1) {
	//cout<< "NAME1= "<< name <<endl;
	for (int ix=0; ix<hist1d->GetNbinsX()+1; ix++) {
	  xbins[ix] = hist1d->GetXaxis()->GetBinLowEdge(ix+1);
	}
	
	h_recoevtvar[iplot] = new TH1F(
				     hist1d->GetName(),
				     hist1d->GetTitle(),
				     hist1d->GetNbinsX(),  xbins);
	h_recoevtvar[iplot]->Sumw2();
      }
      for (int ix=0; ix<hist1d->GetNbinsX()+2; ix++) {
	double aent = hist1d->GetBinContent(ix);
//	if (iplot==0) {cout <<"bince "<<name<<" "<< ix<<" "<< hist1d->GetNbinsX() <<" "<<hist1d->GetXaxis()->GetBinCenter(ix)<<" "<<aent<<endl;}
	//cout<< "NAME2= "<< name <<endl;
	for (int iy=0; iy<aent; iy++) {
	  h_recoevtvar[iplot]->Fill(hist1d->GetXaxis()->GetBinCenter(ix), weight2);
	}
      }
      //cout<<"ok1====="<<endl; 
    //}
  }
    //cout<<"ok3====="<<endl;
    fileIn->cd();
    delete fileIn;
    //cout<<"ok2====="<<endl;
    
 } 

  fileOut->cd();
  fileOut->Write();
}
