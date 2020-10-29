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
  void jec()
  {

 char histname[100];
 char histname1[100];
 char histname2[100];
TLegend *tleg = new TLegend(.5,0.2,0.75,0.5,"","brNDC");
 tleg->SetFillColor(10);
 tleg->SetTextSize(0.06);
 tleg->SetBorderSize(0);
 tleg->SetTextFont(42);

 ofstream file_out("jecunc_thrustc_v3_test.inc");


 //TFile *file1 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/unfold_jec/check/outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.root");
 TFile *file1 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/re-combine/unfold_jec/jec_unc_Madgh_1_c0_e13_130529_kcov.root");
   if (file1->IsOpen()){
     cout<<"File1 opened successfully\n"<<endl;
   }
 //TFile *file2 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/combined/unfold/qcd_evt_byMadgh_outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.root");
 TFile *file2 = TFile::Open("/home/tanmay/QCD/13TeV/miniaod/re-combine/unfold/data_unfolded_by_pythia8.root");
  

   //sprintf(histname,"analyzeBasicPat/reco_typ_0_pt0_eta0_53");

   //TH1F* default_hist = file1->Get(histname);
   // sprintf(histname,"analyzeBasicPat/reco_typ_53_pt0_eta0_3");

   // TH1F* new  = file1->Get(histname);

   //mean_pt->Add(up_pt,-1);
   //mean_pt->Draw()
/*;
int strtiter[5][3][8]={{{4,4,4,4,5,4,4,3},
                        {3,4,5,5,4,5,5,7},
                        {4,4,4,4,4,5,5,4}},// thrustc
                       {{4,4,4,4,4,4,4,4},
                        {4,4,4,4,4,6,4,5},
                        {5,6,6,5,6,6,6,5}},//t3mass
                       {{6,5,8,4,4,4,4,4},
                        {6,5,8,4,4,4,4,4},
                        {6,5,8,4,4,4,4,4}},//y3c
                       {{4,7,7,4,5,5,5,4},
                        {3,5,5,6,7,5,4,7},
                        {6,6,7,7,6,5,8,5}},//broadt
                       {{4,5,4,5,4,5,5,5},
                        {4,5,5,7,8,8,8,7},
                        {5,5,4,5,4,5,5,4}}};//ttmass
*/
int strtiter[4][3][8]={{{4,4,4,4,5,4,4,3},
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
//   int nvar=34;
const int ntyp=4;
const int njetptmn=8;
const int njetetamn=2;
 double mean[100];
 double poserr[100];
 double negerr[100];
 TH1F* fdfz;
 TH1F* fdfx;
 TH1F* fdfy;
 //TH1F* default_hist;
 //TH1F* up_hist;
 //TH1F* down_hist;
 TPostScript ps("outfile",111);
 ps.Range(20,28);
 TCanvas *c1 = new TCanvas("c1", " Fitres", 800, 800);
 ps.NewPage();

  TH1F* fdfx[ntyp][njetptmn][njetetamn][nvarx];
 //Merge all up dn hist
/*  for (int klx=0; klx<nvarx; klx++) {
   int kl=varindex[klx];
   for (int ij=0; ij<ntyp-1; ij++){
   for(int jk=0; jk<njetptmn; jk++){
   //for (int ij=0; ij<ntyp-1; ij++){   
   for (int mn=1;mn<njetetamn;mn++) {
    //  sprintf(histname,"analyzeBasicPat/reco_typ_53_pt%i_eta0_%i", jk, kl );
      sprintf(histname,"bayes_Data_Jets_%i_unfold_Roo_Madgh_reco_typ_%i_pt%i_eta%i_%i", strtiter[klx][ij][jk], ij, jk, mn, kl);
   // sprintf(histname,"analyzeBasicPat/reco_typ_53_pt0_eta0_3" ); 
      fdfx[ij][jk][mn][kl] =(TH1F*) file2->Get(histname);
      }
     }
   }
 }*/
 //

 for (int klx=0; klx<nvarx; klx++) {
 //for (int klx=0; klx<1; klx++) {
   int kl=varindex[klx];
   double meanerr[100]={0};
   for (int ij=0; ij<1; ij++){   
   for(int jk=0; jk<njetptmn; jk++){
   //for(int jk=0; jk<1; jk++){
   //for (int ij=0; ij<ntyp-1; ij++){   
   for (int mn=0;mn<1;mn++) { 
    //if (mn!=0) continue;
     double xpos[100]={0};
     double xerr[100]={0};
     double mean[100]={0};

     double nnmean[100]={0};
     double nnrms[100]={0};

     double poserr[100]={0};
     double negerr[100]={0};

    //  sprintf(histname,"analyzeBasicPat/reco_typ_53_pt%i_eta0_%i", jk, kl );
      sprintf(histname,"bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_typ_%i_pt%i_eta0_%i", strtiter[klx][ij][jk], ij, jk, kl);     
   // sprintf(histname,"analyzeBasicPat/reco_typ_53_pt0_eta0_3" ); 
      TH1F* default_hist =(TH1F*) file2->Get(histname);
     cout<<"histname =  " << histname << endl;
      //default_hist = (TH1F*)fdfz->Clone();
     double histmean=abs(default_hist->GetMean()); 
     //cout<<"histname =  " << histname << endl;
      // TH1F* default_hist = file1->Get(histname);  
      int nbins = default_hist->GetNbinsX();
     //cout<<"histname =  " << histname << endl;
       double binwid = default_hist->GetBinWidth(nbins/2);
     if (binwid<0) continue;
     if (nbins>99) nbins=99;
     for (int ix=0; ix<nbins+1; ix++) {
       double thewid = binwid/default_hist->GetBinWidth(ix);

       mean[ix] = default_hist->GetBinContent(ix);



     for (int iy=1; iy<53; iy +=2){
      // sprintf(histname1,"analyzeBasicPat/reco_typ_%i_pt%i_eta%i_%i", ij, jk,mn,kl);
      // 	sprintf(histname1,"bayes_Data_Jets_4_unfold_Roo_Pythia8_recoreso_typ_0_pt3_eta0_3_1");  
      sprintf(histname1,"bayes_Data_Jets_%i_%i_unfold_Roo_Madgh_reco_typ_%i_pt%i_eta0_%i", strtiter[klx][ij][jk], iy, ij, jk, kl);     
      //cout<<"histname1 =  " << histname1 << endl;
       //TH1F* up_hist = file1->Get(histname1);

       TH1F* up_hist = (TH1F*) file1->Get(histname1);
      // up_hist = (TH1F*)fdfz->Clone();
       double histmean1=abs(up_hist->GetMean());

      // sprintf(histname2,"analyzeBasicPat/reco_typ_%i_pt%i_eta%i_%i", ij+1, jk,mn,kl);
//       sprintf(histname2,"bayes_Data_Jets_4_unfold_Roo_Pythia8_recoreso_typ_0_pt3_eta0_3_2");       
      sprintf(histname2,"bayes_Data_Jets_%i_%i_unfold_Roo_Madgh_reco_typ_%i_pt%i_eta0_%i", strtiter[klx][ij][jk], iy+1, ij, jk, kl);     
 //cout<<"histname2 =  " << histname2 << endl;
       // TH1F* down_hist = file1->Get(histname1);

       down_hist = (TH1F*) file1->Get(histname2);
       //down_hist = (TH1F*)fdfz->Clone();
       double histmean2=abs(down_hist->GetMean());

       double val1 =  up_hist->GetBinContent(ix);
       double val2 =  down_hist->GetBinContent(ix);

       cout<<"val1= "<< val1 << "val2= " << val2 <<"; "<<mean[ix] << endl;

       double dif1 = val1 - mean[ix];
       double dif2 = val2 - mean[ix];

       double meandiff3=0;
       cout<< histmean1 << " Mean up "<< histmean1 <<" "<< histmean2<<endl; 
       double meandiff1=abs(histmean-histmean1);
       double meandiff2=abs(histmean-histmean2);
       if (meandiff1>meandiff2) meandiff3=meandiff1;
       else meandiff3=meandiff2;
       meanerr[jk]+= meandiff3*meandiff3; 
      cout << " Mear err sum = "<<  meanerr[jk]<<endl;
       //meanerr[jk]=sqrt(meanerr[jk]);
              cout<< "dif1 = " <<dif1 << "; " <<dif2 << ";" << mean[ix] << endl;

       if (dif1*dif2<0) {
	 if (dif1>0) {
	   poserr[ix] += dif1*dif1; negerr[ix] += dif2*dif2;
	 } else {
	   poserr[ix] += dif2*dif2; negerr[ix] += dif1*dif1;
	 }
       } else {
	 if (dif1>0) {
	   poserr[ix] += pow(TMath::Max(dif1, dif2),2);
	 } else {
	   negerr[ix] += pow(TMath::Min(dif1, dif2),2);
	 }
       }
       //            }
     } //for (int iy=1; iy < ntyp; iy +=2) 
     } //for (int ix=0; ix<h_pdfweightvar[ij][mn][kl][0]->GetNbinsx(); ix++)

     for (int ix=0; ix<nbins+1; ix++) {
       xpos[ix] = default_hist->GetBinCenter(ix);
       poserr[ix] = sqrt(poserr[ix]);
       negerr[ix] = sqrt(negerr[ix]);
     
     }

    file_out<<"double jecxpos_"<<ij<<"_"<<jk<<"_"<<mn<<"_"<<kl<<"["<<nbins+1<<"] ={";
        for (int ix=1; ix<nbins+1; ix++) {file_out <<xpos[ix]<<", ";}
        file_out<<xpos[nbins]<<"};"<<endl;

        file_out<<"double jecvar_"<<ij<<"_"<<jk<<"_"<<mn<<"_"<<kl<<"["<<nbins+1<<"] ={";
        for (int ix=1; ix<nbins+1; ix++) {file_out <<mean[ix]<<", ";}
        file_out<<mean[nbins]<<"};"<<endl;

        file_out<<"double jecvarpe_"<<ij<<"_"<<jk<<"_"<<mn<<"_"<<kl<<"["<<nbins+1<<"] ={";
        for (int ix=1; ix<nbins+1; ix++) {file_out <<poserr[ix]/max(1.e-6,mean[ix])<<", ";}
        file_out<<poserr[nbins]/max(1.e-6,mean[nbins])<<"};"<<endl;
        cout<<"jecposerr=" <<poserr[nbins]/max(1.e-6,mean[nbins]) << endl;

        file_out<<"double jecvarne_"<<ij<<"_"<<jk<<"_"<<mn<<"_"<<kl<<"["<<nbins+1<<"] ={";
        for (int ix=1; ix<nbins+1; ix++) {file_out <<negerr[ix]/max(1.e-6,mean[ix])<<", ";}
        file_out<<negerr[nbins]/max(1.e-6,mean[nbins])<<"};"<<endl;
        cout<< "jecnegerr= " <<negerr[nbins]/max(1.e-6,mean[nbins]) << endl;

       int naddbin=nbins+1;
       
       TGraphAsymmErrors* gry = new TGraphAsymmErrors(naddbin, xpos, mean,  xerr, xerr, negerr, poserr);
       //  sprintf(title2, "graph_%s", h_pdfweightvar[ij][mn][kl][0]->GetName());
       gry->SetMarkerColor(2);
       gry->SetMarkerStyle(24);
       gry->SetMarkerSize(0.7);
       gry->SetLineColor(4);
       // gry->SetTitle(h_pdfweightvar[ij][mn][kl][0]->GetTitle());
       // gry->GetXaxis()->SetTitle(h_pdfweightvar[ij][mn][kl][0]->GetTitle());
       gry->GetXaxis()->CenterTitle();
       gry->Draw("AP");
       //  gry->Write(title2);
       c1->Update();
       ps.NewPage();
       if (gry) { delete gry; gry=0;}
       //extract through 
       //TGraph *gr1 = (TGraph*)f->Get("name of gry");
       //gr1->Draw("alp"); //or any other option#ifndef PDFSTORE
       //deleting all pdf fil
       
       } // Fpr mn
   }   
    cout << "Mean Error = " << " ; " << meanerr[0]  << endl;
   file_out<<"double jecmeanhist_"<<ij<<"_"<<kl<<"["<<8<<"] ={";
    for(int ixx=0; ixx<8; ixx++) {
    cout<< "Mean error= "<< meanerr[ixx] << endl;
     if(ixx<7) file_out <<meanerr[ixx]<<", ";
     else file_out <<meanerr[ixx]<<"};"<<endl;
     }
   
  } //for (int ij=0; ij<int(njetptall); ij++)
    
    
   } //for (int kl=0;kl<nvar;kl++)

     //  ps.Close();
       
 }
