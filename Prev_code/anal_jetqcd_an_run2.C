//#include "tdrstyle.C"
//#include "CMS_lumi.C"
//#include "SetStyle.C"
#include "CMS_lumi.h"
#include <iostream>
#include "TStyle.h"
#include "TH1.h"
#include "TH1F.h"
#include <iostream>       // std::cout
#include <string>
Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kTRUE));
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return 1000;
}

TH1F* rebin_hist(TH1F* thin, int nxbin, double* xbins) {
  //Unfolded distribution has n+1 bins, including underflow bin
  const int nbn=thin->GetNbinsX();
  double yval[100]={0};
  double yerr[100]={0};
  for (int ij=0; ij<=nbn; ij++) {
    int ix = getbinid(thin->GetBinCenter(ij), nxbin, xbins);
    if (ix>=0 && ix<nxbin) {
      yval[ix] +=thin->GetBinContent(ij);
      double err = thin->GetBinError(ij);
      yerr[ix] += err*err;
    }
  }
  
  TH1F* thout = new TH1F(thin->GetName(), thin->GetTitle(), nxbin, xbins); 
  for (int ix=0; ix<nxbin; ix++) {
    thout->SetBinContent(ix+1, yval[ix]);
    thout->SetBinError(ix+1, sqrt(yerr[ix]));
  }
  return thout;
}

TH1F* rebin_histpdf(TH1F* thin, int nxbin, double* xbins) {
  //without underflow bin
  char name[100];
  const int nbn=thin->GetNbinsX();
  double yval[100]={0};
  double xwid[100]={0};
  double yerr[100]={0};
  double total=0;
  for (int ij=0; ij<=nbn; ij++) {
    int ix = getbinid(thin->GetBinCenter(ij), nxbin, xbins);
    //    cout <<" IJ "<< ij<<" "<< ix<<" "<< thin->GetBinCenter(ij)<<endl;
    if (ix>=0 && ix<nxbin) {
      if (ij>0) {total +=thin->GetBinContent(ij);}
      yval[ix] +=thin->GetBinContent(ij);
      xwid[ix] +=thin->GetBinWidth(ij);
      double err = thin->GetBinError(ij);
      yerr[ix] += 0.000000001; //err*err;
    }
  }
  
  for (int ix=0; ix<=nxbin; ix++) {
    yval[ix] /=total;
  }

  sprintf(name, "pdf_%s", thin->GetName());
  TH1F* thout = new TH1F(name, thin->GetTitle(), nxbin, xbins); 
  for (int ix=0; ix<nxbin; ix++) {
    thout->SetBinContent(ix+1, yval[ix]/xwid[ix]);
    thout->SetBinError(ix+1, sqrt(yerr[ix]));
  }
  return thout;
}



TH1F* rebin_histwithoutlast(TH1F* thin, int nxbin, double* xbins) {
  //Unfolded distribution has n+1 bins, including underflow bin
  const int nbn=thin->GetNbinsX();
  double yval[100]={0};
  double yerr[100]={0};
  int isft=4;
  int tmpbn=0;
  for (int ij=isft; ij<=nbn; ij++) {
    
    yval[tmpbn] =thin->GetBinContent(ij+1);
    yerr[tmpbn] = thin->GetBinError(ij+1);
    xbins[tmpbn] = thin->GetBinLowEdge(ij+1);
    cout <<"ij "<< ij<<" "<<yval[ij]<<" "<<yerr[ij] <<" "<<xbins[ij]<<endl; 
    tmpbn++;
  }
  
  //  nxbin = (xbins[nbn]>=0) ? nbn-1-isft : nbn-isft;
  //  if (xbins[nbn-1-isft]<=0) nxbin--;
  
  nxbin = nbn-isft;

  TH1F* thout = new TH1F(thin->GetName(), thin->GetTitle(), nxbin, xbins); 
  for (int ix=0; ix<nxbin; ix++) {
    thout->SetBinContent(ix+1, yval[ix]);
    thout->SetBinError(ix+1, sqrt(yerr[ix]));
  }
  return thout;
}



Double_t turnon(Double_t* x, Double_t* par){
  double xx = x[0] - par[0];
  //  if (xx<0) {
  //    return 0.5*exp(-0.5*(xx*xx/(par[1]*par[1])));
  //  } else {
  //    return (1.0 - 0.5*exp(-0.5*(xx*xx/(par[1]*par[1]))));
  //  }

  if (xx<0) {
    return 0.5*exp(xx/par[1]);
  } else {
    return (1.0 - 0.5*exp(-xx/par[2]));
  }

}

//const char* subttlx[15]={"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)","(m)","(n)", "(o)"};

  const int nvar=34;
  const char* varname[nvar]={"thrust", "minor", "y3", "thrustc", "thrustce", "thrustcr",
			     "minorc", "minorce", "minorcr", "t3mass", "t3masse", "t3massr",
			     "h3mass", "h3masse", "h3massr", "y3c", "y3ce", "y3cr",
			     "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
			     "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
			     "sphericity", "cparameter", "p3byp12", "p3byp12c"};
  const char* vartitle[nvar]={"#tau_{#perp} ", "T_{ m}", "Y_{23} ", 
			      "#tau_{_{#perp} _{   ,C}}", "#tau_{_{#perp} _{   ,E}}", "#tau_{_{#perp} _{   ,R}}",
			      "T_{ m,C}", "T_{ m,E}", "T_{ m,R}",
			      "#rho_{Tot,C}", "#rho_{Tot,E}", "#rho_{Tot,R}",
			      "#rho_{H,C}", "#rho_{H,E}", "#rho_{H,R}",
			      "Y_{23,C}", "Y_{23,E}", "Y_{23,R}",
			      "B_{ T,C}", "B_{ T,E}", "B_{ T,R}", 
			      "B_{ W,C}", "B_{ W,E}", "B_{ W,R}",
			      "#rho^{T}_{Tot,C}", "#rho^{T}_{Tot,E}", "#rho^{T}_{Tot,R}",
				"#rho^{T}_{H,C}", "#rho^{T}_{H,E}", "#rho^{T}_{H,R}",
			      "S_{_{#perp} _{   ,C}}", "C_{_{#perp} _{   ,C}}", 
			      "2#times P_{T3}/(P_{T2} + P_{T3})", "2#times P_{T3}/(P_{T2} + P_{T3})c"};

void d2plots(const char* name1="thrustc", const char* name2="minorc") {
  //  d2plots( "tmass", "hmass");
  // d2plots( "ttmass", "htmass");
  // d2plots( "broadt", "broadw");
  // d2plots( "sphericity", "cparameter");
  // d2plots( "y3c", "p3byp12c");

  gStyle->SetTitle(0);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(.22);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.09);
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // align right

  TH2F* histx[10];

  char title[100];
  sprintf(title, "%s_PFjet_Gen_c5_e1", name1);
  cout <<title<<endl;
  histx[0] = (TH2F*)gDirectory->Get(title);
  sprintf(title, "%s_PFjet_Gen_c5_e1", name2);
  cout <<title<<endl;
  histx[1] = (TH2F*)gDirectory->Get(title);

  //  histx[0] = (TH2F*)gDirectory->Get("thrustc_PFjet_Gen_c5_e1");
  //  histx[1] = (TH2F*)gDirectory->Get("minorc_PFjet_Gen_c5_e1");

  TCanvas *c1 = new TCanvas("c1", "runfile", 700., 400.);
  c1->Divide(2,2);
  for (int ij=0; ij<2; ij++) {
    c1->cd(2*ij+1);
    histx[ij]->GetYaxis()->SetTitle("Gen Level Object");
    histx[ij]->GetYaxis()->CenterTitle();
    histx[ij]->GetYaxis()->SetTitleSize(0.095);
    histx[ij]->GetYaxis()->SetTitleOffset(0.65);
    histx[ij]->GetYaxis()->SetLabelSize(0.09);
    histx[ij]->GetXaxis()->SetTickLength(0.05);
    histx[ij]->GetXaxis()->SetLabelSize(0.09);
    histx[ij]->GetXaxis()->SetTitle((ij==0)?"Reco log #tau_{_{#perp} _{   ,C}} ":"Reco log T_{ m,C} ");
    histx[ij]->GetXaxis()->SetTitleColor(1);
    histx[ij]->GetXaxis()->CenterTitle();
    histx[ij]->GetXaxis()->SetTitleSize(0.095);
    histx[ij]->GetXaxis()->SetTitleOffset(1.0);
    histx[ij]->GetXaxis()->SetLabelFont(42);
    histx[ij]->GetXaxis()->SetTitleFont(42);
    histx[ij]->GetYaxis()->SetLabelFont(42);
    histx[ij]->GetYaxis()->SetTitleFont(42); 
    
    histx[ij]->Draw("colz");
    //    latex.DrawLatex(.24, .91, subttlx[2*ij+1]);
    histx[ij]->SetMarkerSize(.1);
    histx[ij]->SetMarkerStyle(24);
    histx[ij]->SetMarkerColor(2);

    c1->cd(2*ij+2); histx[ij]->Draw(""); //latex.DrawLatex(.24, .91, subttlx[2*ij+2]);
  }
}
 
void dphibias() {
  gStyle->SetTitle(0);
 gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(.22);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.11);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.09);
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // align right

  TH2F* histx[10];



  histx[0] = (TH2F*)gDirectory->Get("thrustc_GenPhi_c5_e1");
  histx[1] = (TH2F*)gDirectory->Get("thrustc_GenForPhi_c5_e1");
  histx[2] = (TH2F*)gDirectory->Get("browdt_GenPhi_c5_e1");
  histx[3] = (TH2F*)gDirectory->Get("broadt_GenForPhi_c5_e1");


  TCanvas *c1 = new TCanvas("c1", "runfile", 800., 500.);
  c1->Divide(2,2);
  for (int ij=0; ij<4; ij++) {
    c1->cd(ij+1);
    histx[ij]->GetYaxis()->CenterTitle();
    histx[ij]->GetYaxis()->SetTitleSize(0.075);
    histx[ij]->GetYaxis()->SetTitleOffset(0.75);
    histx[ij]->GetYaxis()->SetLabelSize(0.07);
    histx[ij]->GetXaxis()->SetTickLength(0.05);
    histx[ij]->GetXaxis()->SetLabelSize(0.07);

    switch(ij) { 
    case 0 : 
      histx[ij]->GetXaxis()->SetTitle("|#phi_{Jet1} - #phi_{Jet2}|");
      histx[ij]->GetYaxis()->SetTitle("#Delta log #tau_{_{#perp} _{   ,C}}");
      break;
    case 1 :
      histx[ij]->GetYaxis()->SetTitle("log #tau_{_{#perp} _{   ,C}} (Gen)");
       histx[ij]->GetXaxis()->SetTitle("log #tau_{_{#perp} _{   ,C}} (Reco)");
      break;
    case 2 : 
      histx[ij]->GetXaxis()->SetTitle("|#phi_{Jet1} - #phi_{Jet2}|");
      histx[ij]->GetYaxis()->SetTitle("#Delta log T_{ m,C} ");
      break;
    case 3 :
      histx[ij]->GetYaxis()->SetTitle("log T_{ m,C} (Gen)");
       histx[ij]->GetXaxis()->SetTitle("log T_{ m,C} (Reco)");
      break;
    default : break;
      break;

    }
    histx[ij]->GetXaxis()->SetTitleColor(1);
    histx[ij]->GetXaxis()->CenterTitle();
    histx[ij]->GetXaxis()->SetTitleSize(0.075);
    histx[ij]->GetXaxis()->SetTitleOffset(1.20);
    histx[ij]->GetXaxis()->SetLabelFont(42);
    histx[ij]->GetXaxis()->SetTitleFont(42);
    histx[ij]->GetYaxis()->SetLabelFont(42);
    histx[ij]->GetYaxis()->SetTitleFont(42);    
    histx[ij]->Draw("colz");
    latex.DrawLatex(.24, .91, subttlx[ij]);
    //    histx[ij]->SetMarkerSize(.1);
    //    histx[ij]->SetMarkerStyle(24);
    //    histx[ij]->SetMarkerColor(2);

    //    c1->cd(2*ij+2); histx[ij]->Draw("");latex.DrawLatex(.24, .91, subttlx[2*ij+2]);
  }
}
 

const int nvar=34;
const char* varname[nvar]={"thrust", "minor", "y3", "thrustc", "thrustce", "thrustcr",
			   "minorc", "minorce", "minorcr", "t3mass", "t3masse", "t3massr",
			   "h3mass", "h3masse", "h3massr", "y3c", "y3ce", "y3cr",
			   "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
			   "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
			   "sphericity", "cparameter", "p3byp12", "p3byp12c"};
const char* vartitle[nvar]={"#tau_{#perp} ", "T_{ m}", "Y_{23} ", 
			    "log(#tau_{_{#perp} _{  ,C}})", "#tau_{_{#perp} _{   ,E}}", "#tau_{_{#perp} _{   ,R}}",
			    "T_{ m,C}", "T_{ m,E}", "T_{ m,R}",
			    "log(#rho_{Tot,C})", "#rho_{Tot,E}", "#rho_{Tot,R}",
			    "#rho_{H,C}", "#rho_{H,E}", "#rho_{H,R}",
			    "Y_{23,C}", "Y_{23,E}", "Y_{23,R}",
			    "log(B_{ T,C})", "B_{ T,E}", "B_{ T,R}", 
			    "B_{ W,C}", "B_{ W,E}", "B_{ W,R}",
			    "log(#rho^{T}_{Tot,C})", "#rho^{T}_{Tot,E}", "#rho^{T}_{Tot,R}",
			    "#rho^{T}_{H,C}", "#rho^{T}_{H,E}", "#rho^{T}_{H,R}",
			    "S_{_{#perp} _{   ,C}}", "C_{_{#perp} _{   ,C}}", 
			    "2#times P_{T3}/(P_{T2} + P_{T3})", "2#times P_{T3}/(P_{T2} + P_{T3})c"};


const char* basic_vartitle[16]= {"H_{T2} (GeV/c)", "Pt of leading jet (GeV/c)", "Pt of second leading jet (GeV/c)", 
                             "#eta of leading jet", "#eta of second leading jet", "#phi of leading jet", 
                             "#phi of second leading jet","#eta jet", "#phi jet", "#Delta#phi of Jets",
                              "Pt2 x sin( #Delta #phi )/Pt1", "#Delta Pt of two leading jets (GeV/c)",
                              "No. of jet"};


const char* typ[3]={"Jet", "All Particles", "Charged Particles"};
const char* htrang[8]={"73<H_{T,2}<93", "93<H_{T,2}<165", "165<H_{T,2}<225", "225<H_{T,2}<298", "298<H_{T,2}<365", "365<H_{T,2}<452", "452<H_{T2}<557","H_{T,2}>557"};
//#include "jec_systemetic2015.inc"
//#include "jer_systemetic2015.inc"
//#include "unfolding_systemetic2015.inc"
//#include "track_systemetic2015.inc"
//#include "pdf_systemetic2015.inc"
//#include "/home/tanmay/QCD/13TeV/miniaod/combined/systematics/total_systemetic2015.inc"
//#include "/home/tanmay/QCD/13TeV/miniaod/combined/systematics/v3_eta1/total_jets_systemetic2015_v2.inc" //GMA need to change only for test Broad
#include "/home/tanmay/QCD/13TeV/miniaod/combined/systematics/v3_eta1/total_all_systemetic2015_v3.inc" //GMA need to change only for test Broad
#include "/home/tanmay/QCD/13TeV/miniaod/combined/systematics/CUETP8M1/CUETP8M1_eta1.inc" //GMA need to change only for test Broad
//#include "/home/tanmay/QCD/13TeV/miniaod/combined/unfold/only_typ0_pt6_18outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.inc"
//#include "/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/RECO/outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov_v3.inc"
#include "/home/tanmay/QCD/13TeV/miniaod/combined/rootfileplot/RECO/by_pythia8_eta1_outfile13Tev_10_3_Pythia8_1_c0_e13_130529_kcov.inc"
/*
************* including all covariance matrix for eta < 2.4 ***************
#include "outfile_10_3_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_6_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_9_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_12_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_15_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_18_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_21_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_24_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_27_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_30_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_31_Pythia6_1_c0_e24_130529.inc"
#include "outfile_10_33_Pythia6_1_c0_e24_130529.inc"

*/



//void compare1(const char* title="hThrust_a0", int ipad=0, ofstream& fileout="file", int idx=-1, int isame2=-1, int ipteta=0, const char* unf="bayes", const char* etapt="c3_e0", int istrt=4, int incre=-1) {
//void compare1(const char* title="hThrust_a0", int ipad=0, ofstream& fileout="file", int idx=-1, int isame2=-1, int ipteta=0, const char* unf="bayes", const char* etapt="c3_e0", int istrt=4, int incre=-1, int titleindex=1) { //to append all in outfile using next line

void compare1(const char* title="hThrust_a0", int ipad=0, const char* namefile="file", int idx=-1, int isame2=-1, int ipteta=0, const char* unf="bayes", const char* etapt="c3_e0", int istrt=4, int incre=-1, int titleindex=1, int typindex=-1) {

  ofstream fileout;
  fileout.open (namefile, ios::app); //open olly when want add everything 
  //fileout.open (namefile);
  
  // input title : Name of histogramme for isame=-1, else part of it
  //       ipad  : ==0 Only one pad
  //               ==1 Divide canvas in four pad and put histogram in first pad
  //               ==2-4 put histogram in that pad
  //       fileout : ascii output file name
  //       idx   : >=0 for ascii output, -ve for no output 
  //       isame2 : <0 different root files, one by one  (started with this)
  //             : =0, Same root file (_file(ij-1)->cd()
  //             : >0 : individual choices (by hand, no meaning of istrt and incre)
  //             : 1-100 : First file _file0
  //             : >100  : Second file _file1
  //       ipteta : valid of isame=0 only
  //              : 0 for variation of pt threshold
  //              : 1 for variation of eta threshold 
  //              : 2 for variation of vertex bins
  //              : 3 for variation of individual jets

  //              : 8 Comparison of Bayesian and SVD
  //              : 9 unfold Set0 Comparison of Pythia two sets/iteration (same as 10)
  //              : 10 unfold Set0 Comparison of Pythia two sets/iteration
  //              : 11 unfold Set1 Same resposne matrix on different MC
  //              : 12 unfold Set2 Diferent response matrix on data etapt=Data(1), Pythia6(2), Pythia8(3), Madgraph(3) 
  //              : 13 Unfolded data and Genrator MC (combination of reco and gen
  //              : 14 unfold Set4 Unfolded data and Genrator MC (all from reco sample)
  //              : 15 : mixeture of unfolded and standard gen level MC
  //              : 16 Comparison of unfolded and reco data/MC
  //              : 17 same as 16, but Gen level is the statring point
  //              : 18 Uncertain in jet energ ysmearing
  //              : 19 Uncertainty in jet energy scale
  //       istrt : startgin value of _c%i (integer)
  //               for unfolding, # of iteration/regularisation paramters (12-19) & 8
  //               sets 0:"Pythia6", 1:"Pythia8", 2:"Herwig", 2:"Madgraph", 4:"Data" (9-10)
  //       incre : increment of _c%i (integer) variable from isame
  //               for unfolding, this is unfolded(1) /reco distribution (2)
  
  int isame = isame2;
  int ifile = 0;
  if (isame2>0) { 
    isame = isame2%100;
    ifile = int(isame2/100.);
  }

  int basic_title = 0;
  int typnamerange=1;
  int isUnfolding = 0; //1 for unfolding plots, error is imported from last iteration
  // for unfolded data
  if (incre>0) isUnfolding=incre;

  //only variables name in title;
  int isitForPaper=2; // 0, nothing, 1 final plots for note, 2 final plots for paper/pas

  int unfoldederror=0; //Put error by hand if error is less than 1.5\times sqrt{N}

  int mRatio=1; // 1 : plot both superimpose distributions and ratios
                // 0 : only superimposed plots
                // -1 : plot only ratio
  int isLargerDeno = 1; //0; // For calcualtion of chi^2, invert the ratio if it is more than 1 (one), make
                        // symmetric such that order of root file done not matter
  int isAssym = 0; // Asymmetry or ratio plot

  int istatoth = -1; // all(1), only first(0) none(-1); //Statistics for remaining histogramme
  int isrootchi = 1; // 1; //rootfit comparison of different histogrammes or not

  int nzone4= 1 ; //0; //zone 2,2 (otherwise 1,2)
  int isttl2=1; //Use of ttl2 for Tlegent

  int irebin = 2; // rebinning of histogrammes with this
  if (isUnfolding==1) irebin = 1;
  double anorm =1.4; //2.1; //2.1; //1.1; // SetMaximum scale for 1st histogramme;
  //gStyle->SetOptLogy(1);

  double statxpos=0.42; //0.99;
  
  int isArbitrary = 0; //Arbitrary number etc
  int isDoublerat = 0; //formting output for double ratio plots
  int isNote=0; // only for note
  int isSlide=1; // only for slides
  int isSepRatio = 1;//0; //1; // All ratio in same plot (1) / individual(0)  
  int isExternalError=0; // 0 : none
                         // 1 : Calcuate and store it first
                         // 2 : External error for final plot (add all other errors
    //if (isExternalError==2) isSepRatio=0;
   int mcsys=0;
   int isstat=0;
   int sysonmc=0;
  int isCovariance_matrix =0;//0;//1; //0; //1;//1; // if we want to calculate chiq2 using covarince matrix
  int isUnitScale = 1;//1;//0; //Normalised to unit (1), or normalised wrt to 1st histogramme(0) 
  // should be zero for isSepratio==0

  int idError = -1;
  int isIndex = 1; // Put a, b .. in plots (for note)
  const char* subttl[12]={"a", "b", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"};
  const char* ptrange[8]={"110 < p_{T,1} < 170 GeV/c", "170 < p_{T,1} < 250 GeV/c", "250 < p_{T,1} < 320 GeV/c", "320 < p_{T,1} < 390 GeV/c", "p_{T,1} > 390 GeV/c", "p_{T,1} > 110 GeV/c", "p_{T,1} > 170 GeV/c", "p_{T,1} > 250 GeV/c"};
  bool   noYtitle=0; //No title in Y-axis, reduce leftmargin size
  bool noStatOnMC = 0; // for mixed sample do not use sqrt(N) for statisticl error, by hand put it zero
  bool isAida = false; //true;


  //  const double relwt[10]={1.0, 0.00197631, 0.00109191, 0.00188911, 1., 0.00470994, 1., 1., 1., 1.};
  const double relwt[10]={1.0, 1.0, 1.0, 1.0, 1., 1.0, 1., 1., 1., 1.};

  char fname[200];
  char fname2[200];
  char fname3[200]; //Y-axis name of histogramme
  
  char histname[200];
  char histname2[200]; //for isdouble rat
  char histname_covmatrix[200]; // if we want to calculate chiq2 using covarince matrix
  TMatrixD covmatrix;
  TMatrixD errormc_matrix;
  char hname[200];  //Name and title for original ratio histogramme
  char name[200];
  char htitle[200];
  char typname[200];
  TPad *dd1[24];
  int yyy;
  const int ndata=12;
  int ndata2 = ndata-8; //52; -2;

  switch(isame) 
    {
    case 0 :
      switch(ipteta) {
      case 2 : const char* ttl2[ndata]={"npv:1-3", "npv:4-6", "npv:7-9", "npv:10-15", "npv:16-30", "xx","xx","xx","xx","xx", "xx", "xx"}; ndata2=ndata-8; break;
      case 3 : const char* ttl2[ndata]={"# of jet=2", "# of jet=3", "# of jet>=4", "xx","xx","xx","xx","xx","xx","xx", "xx", "xx"}; ndata2=ndata-10; break;

      case 8 : const char* ttl2[ndata]={"bayes", "svd","xx", "xx","xx","xx","xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; idError=4; break; 

      case 10 : const char* Datasets[6]={"Pythia6", "Pythia8", "Herwig", "Madgraph", "Data", "Alpgen"};
	const char* ttl2[ndata]={"Gen", "it=1", "it=2", "it=3", "it=4", "it=5", "it=6", "it=7", "it=8","it=9","it=10", "it=11"}; ndata2 = ndata-1; break;
	//	const char* ttl2[ndata]={"Gen", "ir=4", "ir=5", "ir=6", "ir=7", "ir=8","ir=9", "ir=10", "ir=11","ir=12","ir=13","ir=14"}; ndata2 = ndata-4; break;


      case 11 : const char* ttl2[ndata]={"Madgraph", "Madgraph4%","xx","xx","xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-10; break;
	//      case 12 : const char* ttl2[ndata]={"Pythia6", "Pythia8", "Herwig", "Madgraph", "Data","xx","xx","xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-7; break;
      case 12 : const char* ttl2[ndata]={"Pythia6", "Pythia8", "Madgraph", "Data", "xx", "xx", "xx","xx","xx", "xx", }; ndata2 = ndata-9; idError=3;  break;

	//      case 13 : const char* ttl2[ndata]={"Data", "Pythia6", "D6T", "Z2*", "Pythia8", "Tune2C","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-6; break; //7TeV
	//      case 13 : const char* ttl2[ndata]={"Data", "Pythia6", "D6T", "Z2*", "Pythia8", "Tune1", "Tune2C","Tune2M", "Tune4Cx","Madgraph","xx", "xx"}; ndata2 = ndata-2; break; //7TeV
      case 13 : const char* ttl2[ndata]={"Data", "Pythia6", "D6T", "Z2*", "Pythia8", "Tune2C","Tune2M", "Tune4Cx","Herwig", "TuneEE3C", "Madgraph", "tun23"}; ndata2 = ndata-1; break; //7TeV	

	//      case 13 : const char* ttl2[ndata]={"Data", "Pythia6", "D6T", "Z2", "Tune2C", "Tune4C","Herwig++","Madgraph","xx","xx", "xx", "xx"}; ndata2 = ndata-4; break; //8TeV

      case 14 : const char* ttl2[ndata]={"Data", "Pythia8", "Madgraph", "Herwig++", "Data","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-8; break;

	//      case 14 : const char* ttl2[ndata]={"Herwig", "Data", "Herwig443", "Herwig_UE_EE_3C_7000GeV", "Hetwig8TeV", "xx", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-7; break;
	//      case 14 : const char* ttl2[ndata]={"Herwig", "Data", "443_7Tev", "Def8TeV", "553_8TeV", "xx", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-7; break;

      case 15 : const char* ttl2[ndata]={"Data", "Py8+CUETP8M1", "Py8+Monash", "MadGraph", "Herwig++","xx", "xx","xx","xx", "xx"}; ndata2 = ndata-7; break;
      //case 15 : const char* ttl2[ndata]={"Unf-Pythia8", "Gen-Pythia8", "Madgraph", "Herwigpp","xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; break;
      case 16 : const char* ttl2[ndata]={"Pythia8", "MadGraph", "Herwig++", "Pythia6_reco", "Pythia6_gen", "xx", "xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; break;
      //case 16 : const char* ttl2[ndata]={"Bayes", "SVD", "Pythia6_unf", "Pythia6_reco", "Pythia6_gen", "xx", "xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; break;
	//      case 16 : const char* ttl2[ndata]={"Data_unf", "Data_reco", "Pythia8_unf", "Pythia8_reco", "Pythia8_gen", "xx", "xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-9; break;
	//      case 16 : const char* ttl2[ndata]={"Data_unf", "Data_reco", "Madgraph_unf", "Madgraph_reco", "Madgraph_gen", "xx", "xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-9; break;

      case 17 : const char* ttl2[ndata]={"Data", "Pythia8", "Madgraph", "Herwig++", "Pythia6_reco", "xx", "xx", "xx","xx","xx", "xx", "xx"}; ndata2 = ndata-8; break;
      case 18 : const char* ttl2[ndata]={"DefSim", "Jetsmr", "JetsmrUp", "JetsmrDown", "JetExtSmr", "JetExtSmrUp", "xx", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-6; idError = 1; break;
      case 19 : const char* ttl2[ndata]={"PFjet", "JetESPlus", "JetESMinus", "xx", "xx", "xx", "xx", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-9; idError=2; break;	

//	case 20 :  const char* ttl2[ndata]={ "Data", "Py8+CUETP8M1", "mpi-off", "isr-off", "fsr-off","Pythia8-BeamRemnants-reconnectRange", "Madgraph", "Par-7"}; ndata2 = ndata-7; break;
	case 20 :  const char* ttl2[ndata]={ "Data_Unf_Pythia8", "Data_Unf_Madgraph", "Data_Unf_Herwigpp", "MPI-expPow-", "MPI-pT0Ref+", "MPI-pT0Ref-", "MPI-ecmPow+","MPI-ecmPow-"}; ndata2 = ndata-10; break;
	//      case 20 :  const char* ttl2[ndata]={"Data", "Pythia6", "Perugia", "D6T", "Pythia8", "Herwig", "Madgraph", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-5; break;
		//      case 20 :  const char* ttl2[ndata]={"Data", "Official", "Official_Gen", "xqcut_10", "xqcut_20", "xqcut_30", "xqcut_35", "xqcut_20_Final", "Set7", "Set8", "Set9", "xx"}; ndata2 = ndata-4; break;
//	 case 20 :  const char* ttl2[ndata]={"Data", "Pythia8-official", "Pythia8-new", "Monash", "Par-17", "Par-18", "Par-19","Par-20", "Tune4Cx","Herwig", "TuneEE3C", "Madgraph"}; ndata2 = ndata-4; break;
//	 case 21 :  const char* ttl2[ndata]={"Data", "Pythia8-official", "Pythia8-new", "Monash", "/alpha-set1", "/alpha-set2", "/alpha-set3","/alpha-set5", "/alpha-set2","Herwig", "TuneEE3C", "Madgraph","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"}; ndata2 = ndata; break;

//	case 20 :  const char* ttl2[ndata]={"Data", "Pythia8.185", "Monash", "set-16", "set-17", "set-18", "set-19","set-20", "xx","xx", "xx", "xx"}; ndata2 = ndata-4; break;
	case 21 :  const char* ttl2[ndata]={"Data","Monash", "Monash-20%", "Monash-6%","Monash+20%"}; ndata2 = ndata; break;  
    case 200 :  const char* ttl2[ndata]={"Data", "Pythia6", "Pythia8", "Herwig", "Madgraph", "xx", "xx", "xx", "xx", "xx", "xx", "xx"}; ndata2 = ndata-7; break;
      default : const char* ttl2[ndata]={"Pt>90", "Pt>150", "Pt>240", "Pt>300", "Pt>450","xx","xx","xx","xx","xx", "xx", "xx"}; ndata2=ndata-7; break;
	
      }
      break;
	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8", "MadGraph", "xx", "xx", "xx", "xx","xx","xx","xx" }; ndata2 = ndata-10; break;      
 
//	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8-official", "Pythia8-new", "Monash", "Par-9", "Par-10", "Par-11","Par-12","xx","xx", "xx","xx"}; ndata2 = ndata-4; break;
//	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8-official", "Pythia8-new", "Monash", "Par-13", "Par-14", "Par-15","Par-16","xx","xx", "xx","xx"}; ndata2 = ndata-4; break;
//	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8-official", "Pythia8-new", "Monash", "Par-17", "Par-18", "Par-19","Par-20","xx","xx", "xx","xx"}; ndata2 = ndata-4; break;
//	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8-official", "Pythia8-new", "Monash", "Par-17", "Par-18", "Par-19","Par-20","xx","xx", "xx","xx"}; ndata2 = ndata-4; break;
//	case -19 :  const char* ttl2[ndata]={ "Data", "Pythia8.185", "Monash", "Par-1", "Par-2", "Par-3", "Par-4","xx", "xx","xx", "xx", "xx"}; ndata2 = ndata-5; break;

   case -4 :  const char* ttl2[ndata]={"Data", "MC", "Official", "Private", "xx", "xx", "xx", "xx", "xx"}; ndata2=ndata-9; break;
    case -34 : const char* ttl2[ndata]={"Herwig",  "443_7Tev", "Def8TeV", "553_8TeV","xx", "xx","xx","xx","xx", "xx", "xx", "xx"}; ndata2 = ndata-8; break; // "NoPuNoTrg", "noTrg", "ReGen", "xx", "xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-6; break;
      //    case -34 : const char* ttl2[ndata]={"Def", "ReGen", "NoPU", "NoPuNoTrg", "noTrg", "xx", "xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; break;

    case -35 : const char* ttl2[ndata]={"Def", "NoPU", "NoPuNoTrg", "noTrg", "xx", "xx", "xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-8; break;
    case -36 :  const char* ttl2[ndata]={"Without bkg", "With bkg", "xx", "xx", "xx", "xx", "xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-10; break;
    case -37 : const char* ttl2[ndata]={"Pythia6", "D6T", "Z2*", "Pythia8", "Tune2C","Tune2M", "Tune4Cx","Herwig", "TuneEE3C", "Madgraph", "xx"}; ndata2 = ndata-2; break; //7TeV	
    case -38 : const char* ttl2[ndata]={"Data", "Pythia6", "D6T", "Z2*", "Pythia8", "Tune2C","Tune2M", "Tune4Cx","Herwig", "TuneEE3C", "Madgraph"}; ndata2 = ndata; break; //7TeV


    case 4 :  const char* ttl2[ndata]={"PFjet", "PFjetak7", "PFjetkt4","xx", "xx", "xx","xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-9; break;
    case 7 : 	  const char* ttl2[ndata]={ "Jetsmr", "JetsmrUp", "JetsmrDown", "JetExtSmr", "JetExtSmrUp", "PFjet","xx","xx","xx","xx", "xx", "xx"}; ndata2=ndata-6; break;
    case 61 : const char* ttl2[ndata] = {"Genjet", "PFjet", "PFjetak7", "PFjetkt4", "Calojet", "xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-7; break;
    case 62 : const char* ttl2[ndata] = {"PFjet", "PFjetak7", "PFjetkt4", "Calojet", "xx","xx","xx","xx","xx", "xx", "xx"}; ndata2 = ndata-8; break;

    default : break;
  }

  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat((istatoth>=0) ? 1110 : 0); //1100); //1110); //111110);  
  gStyle->SetStatFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");

  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineWidth(3);

  gStyle->SetStatY(0.99);
  
  gStyle->SetStatX(statxpos);
  if (idx<0) {
    if (ndata2<=3) { 
      gStyle->SetStatW( (nzone4>0) ? 0.34: 0.20); //pss
      gStyle->SetStatH(0.14); //6); //pss
    } else {
      gStyle->SetStatW( (nzone4>0) ? 0.30: 0.18); //pss
      gStyle->SetStatH(0.13); //6); //pss      
    }
  } else {
    gStyle->SetStatFontSize(0.015); //0.04); //pss
  }

  gStyle->SetPadTopMargin(0.01);


  //  if (gStyle->GetOptLogy()==1) {
    gStyle->SetPadLeftMargin((idx<0) ? 0.17 : 0.17); //0.14 : 0.17); // 0.07 : 0.10); // 0.18 for paper
    //  } else {
    //    gStyle->SetPadLeftMargin((idx<0) ? 0.14 : 0.19); // 0.17 : 0.19); // 0.10 : 0.12); 
    //  }
  gStyle->SetPadRightMargin(0.04);
  if (noYtitle) gStyle->SetPadLeftMargin(.09);
  if (mRatio==0) {
    gStyle->SetPadBottomMargin(0.17); //0.17
    gStyle->SetPadLeftMargin(0.20);
    gStyle->SetPadRightMargin(0.02);
  }
  if (isSlide && mRatio!=0) { gStyle->SetPadLeftMargin(.17);}

  TPad* c1_2;
  TPad* c1_3;
  TPad* c2;
  int xrng;
  int yrng;
  if (ipad<=1) {
    if (mRatio!=0) { // alltogether>0) {
     // xrng = 650;  yrng=800;
      ///      xrng = 650;  yrng=900; //900, 650; //1000;
            xrng = 800;  yrng=350; //for isSlide
      //      xrng = 800;  yrng=600;
      //      xrng = 900;  yrng=500;
      //      xrng = 900;  yrng=500; //600; //450;
      //      xrng = 950;  yrng=1100;
  setTDRStyle(gStyle);   
    } else {
      if (isArbitrary) {
	xrng=700; yrng=500; // xrng=600; yrng=300;
      } else {
        if (idx<0) {
	  if (nzone4>0) {
	    xrng = 700;  yrng=360;
	  } else {
	    xrng = 500; yrng = 500; 
	  }
        } else {
   	  xrng = 600; yrng = 750;
        }
      }
    }
  //}
    TCanvas *c1 = new TCanvas("c1", "runfile", xrng, yrng);
    //    c1->SetCanvasSize(1050,1400);
  int iPeriod=4;
  int iPos=0;
   // CMS_lumi(c1,iPeriod,iPos);
    c1->SetCanvasSize(600,600);
    //c1->SetCanvasSize(600,700);
    //CMS_lumi(c1,iPeriod,iPos);
   // c1->Update();
  }
  if (ipad==1) {
    if (isArbitrary) {
      if (nzone4>0) { c1->Divide(4,1); } else {c1->Divide(2,1);} 
    } else if (mRatio!=0) { // alltogether>0) { 
      if (nzone4>0) { c1->Divide(2,1,0.02,0.000005); } else {c1->Divide(22,1);} 
//      if (nzone4>0) { c1->Divide(3,2,0.003,0.000005); } else {c1->Divide(1,1);} 
    } else {
      //      if (nzone4>0) { c1->Divide(5,2); } else {c1->Divide(1,2);}
      if (nzone4>0) { c1->Divide(3,1); } else {c1->Divide(2,3);}
    }
  }

  if (ipad>=1) c1->cd(ipad);
  if (mRatio!=0) { // (alltogether>0) {
    gStyle->SetPadBottomMargin((isSepRatio) ? 0.362 : 0.27); //0.03);
    gStyle->SetPadTopMargin(1.53);
    //c1_3 = new TPad("c1_3", "newpad3",0., 0., 1., 0.988);
    c1_3 = new TPad("c1_3", "newpad3",0., 0., 1, 0.95);
    //    c1_3 = new TPad("c1_3", "newpad3",0., 0., 0.9, 0.9);
     c1_3->Draw();
    
    //CMS_lumi(c1,iPeriod,iPos);
    //    if (ipad==1){
//      cmsPrel2();
    //    }else if (ipad==2){
    //      cmsPrel3_2();
    //    }
    c1_3->Divide(1, (isSepRatio) ? ((mRatio>0) ? ndata2 : ndata2-1) : 2 , -1.e-40, -1.e-40,0);

    //CMS_lumi(c1,iPeriod,iPos);
    if (mRatio>0) { 
     //setTDRStyle(gStyle); 
      dd1[0] = (TPad*) c1_3->GetPad(1);
      dd1[0]->SetTopMargin(0.01);
      if (isUnitScale){
	dd1[0]->SetLeftMargin(0.23);
      }else{
	dd1[0]->SetLeftMargin(0.182);
      }
      dd1[0]->SetRightMargin(0.005);
      //      dd1[0]->SetLogy(0);
      if (isSepRatio) {
	double gapxx =  0.95/(ndata2+2); // 0.083; 0.93 ->0.94
	double yy2 = 0.995;
	double yy1= yy2- 3*gapxx;
	dd1[0]->SetPad(0.01,yy1,0.9985,yy2+0.0049);
      } else {
	dd1[0]->SetPad(0.0001,0.4,1,0.999);
	//dd1[0]->SetPad(0.0001,0.4,0.1,0.999);
      }
      dd1[0]->cd();
    }
  }

  double x1= (statxpos<0.6) ? 0.65 : 0.20; // 0.45; //0.60 : 0.20;
  double gap = (mRatio!=0) ? 0.35 : 0.35; // 0.25 : 0.35;
  double x2= x1 + gap;
  
  //  double y1 = 0.05, y2=y1+0.10+0.05*ndata2; //y1 = 0.28,
  ///  double y1 = 0.28, y2=y1+0.10+0.04*ndata2; //y2=y1+0.10+0.04*ndata2;
  double y1 = 0.28, y2=y1+0.10+0.02*ndata2;
  if (isSlide && mRatio!=0) y2+=.2;
  //  TLegend *tleg = new TLegend(x1,y1,x2,y2,"","brNDC");
  //TLegend *tleg = new TLegend(.3,0.3,0.5,0.5,"","brNDC");
  if (strstr(etapt,"typ_0_pt0_eta1_3"))  TLegend *tleg = new TLegend(0.63,0.2,0.80,0.55,"","brNDC");
  else if (strstr(etapt,"typ_0_pt1_eta1_3"))  TLegend *tleg = new TLegend(0.56,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt2_eta1_3"))  TLegend *tleg = new TLegend(0.56,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt3_eta1_3"))  TLegend *tleg = new TLegend(0.54,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt4_eta1_3"))  TLegend *tleg = new TLegend(0.55,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt5_eta1_3"))  TLegend *tleg = new TLegend(0.54,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt6_eta1_3"))  TLegend *tleg = new TLegend(0.54,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_0_pt7_eta1_3"))  TLegend *tleg = new TLegend(0.51,0.2,0.80,0.6,"","brNDC");
  else if (strstr(etapt,"typ_1_pt0_eta1_3"))  TLegend *tleg = new TLegend(0.60,0.2,0.80,0.55,"","brNDC");
  else if (strstr(etapt,"typ_2_pt0_eta1_3"))  TLegend *tleg = new TLegend(0.63,0.2,0.80,0.55,"","brNDC");
  else TLegend *tleg = new TLegend(0.50,0.2,0.67,0.6,"","brNDC");

  tleg->SetFillColor(10);
  tleg->SetTextSize((isSlide && mRatio!=0) ? 0.055 : (idx<0) ?0.06 : 0.058); //pss 0.07 : 0.058); //0.06); //pss
  tleg->SetBorderSize(0);
  tleg->SetTextFont(42);
  
  const int nmxbin=1000;

  TH1F* fdfz[ndata]; //copy from root file
  TH1F* fdfx[ndata]; //for all plots etc (normalised)
  TH1F* fdfy[ndata]; //for isExternalError==2, add systematic errors in data and MC
  TH1F* fdfy1[ndata]; //for isExternalError==2, add systematic errors in data and MC
  TH1F* fdfy11[ndata]; //for isExternalError==2, add systematic errors in data and MC
  TH1* h1 = new TH1I("h1", "h1 title", 100, 0.0, 4.0);
  TH1* h2 = new TH1I("h1", "h2 title", 100, 0.0, 4.0);
  TH1* h3 = new TH1I("h3", "h3 title", 100, 0.0, 4.0);

  TH1F* fdfzunf[ndata]; //To get error from maximum iteration //copy from root file

  TH1F* ratHist[ndata]={0}; //Histogramme for ratio plot;
  TH1F* ratHist1[ndata]={0}; //Histogramme for ratio plot;
  TH1F* ratHist2[ndata]={0}; //Histogramme for ratio plot;
  TH1F* ratHistOrg[ndata]={0}; //Histogramme for ratio plot with only statistical error;

  double scalex[ndata];
  TPaveStats* st1[ndata];
  TGraphAsymmErrors* gr4[ndata+5][ndata];

  double binrange[nmxbin+1]={0};	  
  double binwidth[nmxbin]={0};
  double bincenter[nmxbin]={0};
  double binentry[nmxbin][ndata]={0};
//  double staterrData[nmxbin][ndata]={0}; 
  double binerror[nmxbin][ndata]={0};
  double binerrorNor[nmxbin][ndata]={0};
  double minbincontent[ndata]={0};
  double efficiency[nmxbin]={0};
  double efficiencyforerr[nmxbin]={0};
  double errorinx[nmxbin]={0};
  double erroriny[nmxbin]={0};
  double erroriny2[nmxbin]={0}; //for isLargerDeno
  double errorinyExt[nmxbin]={0.}; //for normalised plots, error in ratio excluding error in data, for isExternalError=2,
  double errorinyDoub[nmxbin]={0}; //for double ratio, store maxmimum error of two
  double meanerr[ndata]={0};

  double errpl[nmxbin]; //Calculation
  double errmi[nmxbin];  

  //03/06/2010 store maximum deviation to combine total errors at the end
  double devEfficiency[nmxbin]={0.0};
  double efficiencyErr[nmxbin]={0.0};

  for (int ij=0; ij<nmxbin; ij++) { errpl[ij] = errmi[ij] = 0.0;}
  
  double error_pl[25]; //Total error
  double error_mi[25]; 
  
  double tterrpl[25]; //Total systematic error (from previous steps)
  double tterrmi[25];

  double chisqq=0;
  int ndff=0;
  float fact=1.;

  int icol = 0;
  int linestl=0;
  float ntot = 0.;
  float mean = 0.;
  float rms = 0.;
  float alow=100.;
  float ahg=-100.;
  double sysperr[600];
  double sysnerr[600];
  double sysperrMC[600];
  double sysnerrMC[600];
  double updnerrMC[600];
  double syserrData[600];
  double staterrData[600];
  double toterrData[600];
  double chi2[ndata] = {0};
  char namey[100];
  sprintf(namey, "outfile_%s.txt", title);
  ofstream file_out(namey);
  ofstream outfile; //write all output in same file
  ofstream outfile1; //write all output in same file
  ofstream file_stat_err; //write all stat err
  outfile1.open("chi.tex",ios::app);
  file_stat_err.open("stat_err.inc", ios::app);
  //if(typindex==0 && titleindex==3){ 
  TString s(etapt);
  //   TString s2( s(6,7) );
    if(strstr(etapt,"typ_0_pt0_eta1_3")){
  outfile1 << " \\documentclass{article}"<<endl;

 outfile1 << "\\usepackage{amstext}"<<endl;
 outfile1 << "\\usepackage{array}"<<endl;
 outfile1 << "\\begin{document}"<<endl;

outfile1 << "\\begin{table}"<<endl;
// outfile << "\\begin{sidewaystable} "<<endl;
 outfile1 << "\\begin{center} "<<endl;
 outfile1 << "\\begin{tabular}{ |c||c| } " <<endl;

 outfile1<< "\\hline "<<endl;


 outfile1<< "name " << "& " << "$\\chi^2$" <<"\\\\"<<endl;
 outfile1<< "\\hline"<< endl;
 }
  if(typindex==0 && titleindex==3) outfile.open("typ0_allht3_delta.log",ios::app);
  if(typindex==1 && titleindex==3) outfile.open("typ1_allht3_delta.log",ios::app);
  if(typindex==2 && titleindex==3) outfile.open("typ2_allht3_delta.log",ios::app);
  if(typindex==0 && titleindex==9) outfile.open("typ0_allht9_delta.log",ios::app);
  if(typindex==1 && titleindex==9) outfile.open("typ1_allht9_delta.log",ios::app);
  if(typindex==2 && titleindex==9) outfile.open("typ2_allht9_delta.log",ios::app);
  if(typindex==0 && titleindex==18) outfile.open("typ0_allht18_delta.log",ios::app);
  if(typindex==1 && titleindex==18) outfile.open("typ1_allht18_delta.log",ios::app);
  if(typindex==2 && titleindex==18) outfile.open("typ2_allht18_delta.log",ios::app);
  if(typindex==0 && titleindex==24) outfile.open("typ0_allht24_delta.log",ios::app);
  if(typindex==1 && titleindex==24) outfile.open("typ1_allht24_delta.log",ios::app);
  if(typindex==2 && titleindex==24) outfile.open("typ2_allht24_delta.log",ios::app);
  //outfile.open("thrust_typ0_delta.log",ios::app);
  ofstream file_figname("figname_jet.sh", ios::app);
  //outfile << " \\documentclass{article}"<<endl;

 //outfile << "\\usepackage{amstext}"<<endl;
 //outfile << "\\usepackage{array}"<<endl;
 //outfile << "\\begin{document}"<<endl;

//outfile << "\\begin{table}"<<endl;
// outfile << "\\begin{sidewaystable} "<<endl;
 //outfile << "\\begin{center} "<<endl;
 //outfile << "\\begin{tabular}{ |c||c|c|c|c||c|c|c|c|} " <<endl;

 //outfile<< "\\hline "<<endl;


 //outfile<<"Name" << " & " << "73<$H_{T,2}$>93" << " & " << "93<$H_{T,2}$>165" << " & " << "165<$H_{T,2}$>225" << " & " << "225<$H_{T,2}$>298" << "298<$H_{T,2}$>365" << " & " << "365<$H_{T,2}$>452" << " & " << "452<$H_{T,2}$>557" << " & " << "$H_{T,2}$ > 557" << "\\\\" <<endl;
 //outfile<< "\\hline"<< endl;

//  char[5]={Y23,Broad,ttmass,t3mass};
  //for matchin with unfolding distributions
  int nbinunf=0;
  double binunfs[600];
  // Roo_Data_thrustc_PFjet_c3_e1
  //    int nbinunf=18; //thrust
  //  double binunfs[19]={-10.87, -8.21, -7.23, -6.55, -6, -5.53, -5.1, -4.7, -4.32, -3.95, -3.59, -3.24, -2.9, -2.57, -2.25, -1.95, -1.67, -1.37, -0.98};
    //    int nbinunf=18; //minorc
  //  double binunfs[19]={-5.48, -4.32, -3.79, -3.42, -3.12, -2.86, -2.63, -2.42, -2.22, -2.02, -1.83, -1.65, -1.47, -1.3, -1.13, -0.97, -0.82, -0.69, -0.57};
  //  int nbinunf=24; //broadt
  //  double binunfs[25]={-3.63, -3.45, -3.27, -3.09, -2.91, -2.73, -2.56, -2.39, -2.23, -2.07, -1.92, -1.78, -1.64, -1.51, -1.39, -1.27, -1.16, -1.06, -0.96, -0.87, -0.78, -0.7, -0.62, -0.55, -0.48 };

  int nbinx = 0;
  char  pch1[300];

  for (int ij=0; ij<ndata2; ij++) {
    cout <<"ij "<< ij<<endl;
    icol++;
    if (icol%3==0 ||icol%5==0 || icol%7==0) icol++;
    if (icol>10) icol +=9;

    if (isame<0) {
      switch(ij) {
      case 0 : _file0->cd(); break;
      //case 1 : _file1->cd(); break;
      //case 2 : _file2->cd(); break;
/*      case 3 : _file3->cd(); break;
      case 4 : _file4->cd(); break;
      case 5 : _file5->cd(); break;	
      case 6 : _file6->cd(); break;
      case 7 : _file7->cd(); break;
      case 8 : _file8->cd(); break;
      case 9 : _file9->cd(); break;	
      case 10 : _file10->cd(); break;
      case 11 : _file11->cd(); break;	
*/
      default : _file0->cd(); break;
      }
      //      sprintf(histname, "%s", title);
      if (ij==0) {
//        sprintf(histname, "%s_Genjet_c0", title);
//          sprintf(histname, "qcdevt/%s_thrustce_0", title);
 //         sprintf(histname, "allcharged_thrustc_%s", title);
    //        sprintf(histname, "analyzeBasicPat/recojetallave_pt_0");
            //sprintf(histname, "analyzeBasicPat/%s", title);
      //       sprintf(histname, "analyzeBasicPat/hjetptbypl_0");
            sprintf(histname, "recojet1_pt_1");
   } else if(ij==1) {sprintf(histname, "%s", title);}
     else if(ij==2) {sprintf(histname, "%s", title);}  


else {
//	sprintf(histname, "%s_Genjet_c0", title);
//        sprintf(histname, "qcdevt/%s_thrustce_0", title);
        //sprintf(histname, "allcharged_thrustc_%s", title);
  //    sprintf(histname, "hjetptbypl_0");
       //sprintf(histname, "recojetallave_pt_0");
       sprintf(histname, "%s", title);
      }
      cout <<"histnamex "<<histname<<endl;
    } else if (isame==0) {
      
      switch (ipteta) {
      case 2 : sprintf(histname, "%s_c3_e1_vtx%i", title, ij);  break; //Primary vertex
      //case 3 : sprintf(histname, "%s_c3_e1_j%i", title, ij); break;  //Number of jets
      case 3 : sprintf(histname, "analyzeBasicPat/njets_0"); break;  //Number of jets


     case 8 :
       //Comparison of Bayesian and SVD
       //        root ../unfold/outfile_ff3alliter_Pythia6_1_c3_e1.root ../unfold/outfile_ff3alliter_Pythia8_1_c3_e1.root ../unfold/outfile_ff3alliter_Madgraph_1_c3_e1.root

       // _file0->cd();
       //       compare("Pythia6_thrustc_PFjet_c3_e1",1,"xx",-1, 0, 8);
       // _file1->cd();
       //       compare("Pythia8_thrustc_PFjet_c3_e1",2,"xx",-1, 0, 8);
       // _file2->cd();
       //       compare("Madgraph_thrustc_PFjet_c3_e1",3,"xx",-1, 0, 8);
       sprintf(histname,"%s_Data_PFjet_%i_unfold_Roo_%s", ttl2[ij], (ij==0)?istrt:9,title);

       break;

      case 9 :
	// unfold Set0 Comparison of Pythia two sets/iteration
	// root ../unfold/outfile_ff3alliter_Pythia6_1_c3_e1.root
	// _file0->cd();
	//      compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,9,"svd","Pythia6", 0);
	// _file1->cd();
	//      compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,9,"svd","Pythia8", 0);

	//	compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,9,"svd","Pythia6", 2);
	//      compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,9,"svd","Pythia6", 2);
	//	compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,9,"svd","Pythia6", 3);
 	//compare("minorc_PFjet_c3_e1",2,"xx",-1, 0,9,"bayes","Pythia6", 0);

	//	int bayessft = (strstr(unf, "svd")) ? 15 : 15;
	int bayessft=9;
	sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_%s_%s", unf, Datasets[istrt], bayessft-ij, etapt, title);
	cout <<"histname =================== "<< ij<< ""<< histname<<endl;
	break;
	

      case 10 :
	// unfold Set0 Comparison of Pythia two sets/iteration
	// root ../unfold/outfile_ff3alliter_Pythia6_1_c3_e1.root
	// svd_Madgraph_PFjet_7_unfold_Roo_Pythia6_thrustc_PFjet_c3_e1
	//	compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,10,"svd","Pythia6",0);
	//    compare("thrustc_PFjet_c3_e1",2,"xx",-1, 0,10,"bayes","Pythia6",0);

	//	_file0->cd(); 
	isame=0;

	if (ij==0) {
	  /*char  varnm[200];
	  char* stry = title;
	  char * pchy = strchr(stry,'_');
	  int len = pchy-stry;
	  strncpy (varnm,stry,len);
	  varnm[len]='\0';
	  char * pchz = strstr(stry,"_c");*/
	  //	  sprintf(histname, "Roo_%s_%s_Genjet%s",etapt, varnm,pchz);
	  //sprintf(histname, "Roo_%s_%s_Genjet%s",Datasets[istrt], varnm,pchz);
	  sprintf(histname, "Roo_Pythia8_gen_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24

	} else {
	  //int bayessft=(strstr(unf, "svd")) ? 3 : 0; //1 : 1; //4 : 5;
	  //sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_%s_%s", unf, Datasets[istrt], ij+bayessft, etapt, title);
          sprintf(histname, "bayes_Madgh_Jets_%i_unfold_Roo_Pythia8_reco_%s", ij+1, etapt); //bayes_Madgh_Jets_11_unfold_Roo_Pythia8_reco_typ_2_pt7_eta0_24
	}
	break;

       case 11 :
        isame=0;

        if (ij==0) {
          sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_typ_2_%s", istrt, etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
        } else {
          sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_typ_4_%s", istrt, etapt); //bayes_Madgh_Jets_11_unfold_Roo_Pythia8_reco_typ_2_pt7_eta0_24
        }
        break;

	
/*      case 11 :
	//xxxxxx unfold Set1 Same resposne matrix on different MC
	//root ../unfold/outfile_ff_Pythia6_c3_e1.root ../unfold/outfile_ff_Pythia8_c3_e1.root ../unfold/outfile_ff_Madgraph_c3_e1.root

	//compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,11, "svd", "cxx", 3, 1);
	//compare("thrustc_PFjet_c3_e1",2,"xx",-1, 0,11, "svd", "cxx", 3, 2);
	_file0->cd(); isame=0;
	if (isUnfolding==1) {
	  sprintf(histname, "svd_%s_PFjet_3_unfold_Roo_%s_%s", ttl2[ij], ttl2[0], title);
	} else {
	  sprintf(histname, "Roo_%s_%s", ttl2[ij], title);
	}
	  break;
*/	
      case 12 :
	//       unfold Set2 // Diferent response matrix on data

	// root ../unfold/outfile_bb_Pythia6_1_c0_e1.root ../unfold/outfile_bb_Pythia8_1_c0_e1.root ../unfold/outfile_bb_Madgraph_1_c0_e1.root
	//	compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,12, "bayes", "Data", 2, 1);
	//	compare("thrustc_PFjet_c3_e1",2,"xx",-1, 0,12, "bayes", "Pythia6", 2, 1); 
	//	compare("thrustc_PFjet_c3_e1",3,"xx",-1, 0,12, "bayes", "Pythia8", 2, 1);
	//	compare("thrustc_PFjet_c3_e1",4,"xx",-1, 0,12, "bayes", "Madgraph", 2, 1); 

	//	compare("thrustc_PFjet_c3_e1",1,"xx",-1, 0,12, "svd", "Data", 7, 1);
	//	compare("thrustc_PFjet_c3_e1",2,"xx",-1, 0,12, "svd", "Pythia6", 7, 1); 
	//	compare("thrustc_PFjet_c3_e1",3,"xx",-1, 0,12, "svd", "Pythia8", 7, 1);
	//	compare("thrustc_PFjet_c3_e1",4,"xx",-1, 0,12, "svd", "Madgraph", 7, 1); 

	switch(ij) {
	case 0 : _file0->cd(); break;
	case 1 : _file1->cd(); break;
	case 2 : _file2->cd(); break;
	default : _file0->cd(); break;
	}
	sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_%s_%s", unf, etapt, istrt, ttl2[ij], title);
	break;

      case 13 :
	//      Unfolded data and Genrator MC
	// root ../unfold/outfile_dd_Pythia6_1_c0_e1.root gen_qcdevt_pythia6_7tev_d6t_flat_all_cc.root gen_qcdevt_pythia6_7tev_z2star_flat_all_cc.root gen_qcdevt_pythia8_7tev_tune2c_flat_all_cc.root

	//       compare("thrustc_PFJet_c1_e1", 1, "xx", -1, 0,14)
	// compare("thrustc", 1, "xx", -1, 0,14) compare("minorc", 2, "xx", -1, 0,14)
	//   	compare("thrustc", 1, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)
	// compare("minorc", 2, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)
	//7TeV
	switch(ij) {
	case 0 : _file0->cd(); break;
	case 1 : _file0->cd(); break;
	case 2 : _file1->cd(); break;
	case 3 : _file2->cd(); break;
	case 4 : _file0->cd(); break;
	case 5 : _file3->cd(); break;
	case 6 : _file4->cd(); break;
	case 7 : _file5->cd(); break;	  
	case 8 : _file0->cd(); break;
	case 9 : _file6->cd(); break;	  
	case 10 : _file7->cd(); break;
	case 11 : _file8->cd(); break;	  

	default : _file0->cd(); break;
	}

	if (ij==0) { 
	  sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt);
	} else if (ij==1 || ij==4 || ij==8) {
	  sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
	} else {
	  sprintf(histname, "%s_Genjet_%s", title, etapt);
	}
	break;

      case 14 :
	//      Set4 : Unfolded data and Genrator MC (fro mRECO smaple)
	// root ../unfold/outfile_aa_Pythia6_1_c3_e1.root gen_qcdevt_herwigpp_fall11_443_all_redgen2.root list_rootfiles_sum12_herwigpp_ee3c_redbin.root  gen_qcdevt_herwigpp_8tev_flat_redgen.root 
	//       compare("thrustc_PFJet_c1_e1", 1, "xx", -1, 0,14)
	// compare("thrustc", 1, "xx", -1, 0,14) compare("minorc", 2, "xx", -1, 0,14)
	//   	compare("thrustc", 1, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)
	// compare("minorc", 2, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)

	if (ij==0) { 
	 // sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt);
	  sprintf(histname, "Roo_Data_reco_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
          //sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_%s", etapt);
	} 
         else if(ij==1) { // if (ij==0) {
	  //sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
          //sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt0_eta0_3"); //bayes_Madgh_Jets_11_unfold_Roo_Pythia8_reco_typ_2_pt7_eta0_24
	  sprintf(histname, "Roo_Pythia8_reco_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
	}
        else if(ij==2) {sprintf(histname, "Roo_Madgh_reco_%s", etapt);}	
        else {sprintf(histname, "Roo_Herwigpp_reco_%s", etapt);}	
	break;

 case 15 :
        //      Set4 : Unfolded data and Genrator MC (fro mRECO smaple)
        // root ../unfold/outfile_aa_Pythia6_1_c3_e1.root gen_qcdevt_herwigpp_fall11_443_all_redgen2.root list_rootfiles_sum12_herwigpp_ee3c_redbin.root  gen_qcdevt_herwigpp_8tev_flat_redgen.root 
        //       compare("thrustc_PFJet_c1_e1", 1, "xx", -1, 0,14)
        // compare("thrustc", 1, "xx", -1, 0,14) compare("minorc", 2, "xx", -1, 0,14)
        //      compare("thrustc", 1, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)
        // compare("minorc", 2, "xx", -1, 0, 14, "bayes", "c3_e1",8,2)

        if (ij==0) {
         // sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt);
          //sprintf(histname, "Roo_Data_reco_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
          //sprintf(histname, "bayes_Pythia8_Jets_%i_unfold_Roo_Madgh_reco_%s", istrt, etapt);
          sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt);
        }
         else if(ij==1) { // if (ij==0) {
          //sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
          //sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt0_eta0_3"); //bayes_Madgh_Jets_11_unfold_Roo_Pythia8_reco_typ_2_pt7_eta0_24
          sprintf(histname, "Roo_Pythia8_gen_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
          //sprintf(histname, "svd_Data_Jets_9_unfold_Roo_Pythia8_reco_%s", etapt);
          //sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt);
        }
        else if(ij==2) {
          sprintf(histname, "Roo_default_gen_%s", etapt);
         } 
        else if(ij==3) {
        sprintf(histname, "Roo_Madgh_gen_%s", etapt);
        }
        else {
         // sprintf(histname, "svd_Data_Jets_9_unfold_Roo_Pythia8_reco_%s", etapt);
         sprintf(histname, "Roo_Herwigpp_gen_%s", etapt);
         }
        break;
 
/*      case 15 : 
	// unfolded data with unfolded MC/ Reco Data with Reco MC
	// root ../unfold/outfile_aa_Pythia6_0_c3_e1.root ../unfold/outfile_aa_Pythia8_0_c3_e1.root ../unfold/outfile_aa_Madgraph_0_c3_e1.root
	// _file0->cd();
	// compare("thrustc_PFjet_c3_e1", 1, "xx", -1, 0, 15, "svd", "Pythia6",8,2)
	// _file0->cd();
	// compare("thrustc_PFjet_c3_e1", 1, "xx", -1, 0, 15, "bayes", "Pythia6",8,1)
	// _file1->cd();
	// compare("thrustc_PFjet_c3_e1", 2, "xx", -1, 0, 15, "bayes", "Pythia8",8,1)
	// _file2->cd();
	// compare("thrustc_PFjet_c3_e1", 3, "xx", -1, 0, 15, "bayes", "Madgraph",8,1)

	if (isUnfolding==1) { 
	  if (ij<5) {
	    sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_%s_%s", unf, ttl2[ij], istrt, etapt, title); 
          sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_%s", etapt);
	  } else {
	    //	    _file1->cd();
	  }
	} else {
	  if (ij<5) {
	    sprintf(histname, "Roo_%s_%s", ttl2[ij],title); 
	  } else {

	  }
	}
	break;
*/
 /*     case 16 : 
	//comparison of unfolded data and with unfolded MC and with reco data and MC
	//root ../unfold/outfile_aa_Pythia6_1_c3_e1.root
	//compare("thrustc", 1, "xx", -1, 0, 16, "bayes", "c3_e1",8,1)
	//compare("thrustc", 2, "xx", -1, 0, 16, "svd", "c3_e1",8,1)
	//compare("minorc", 1, "xx", -1, 0, 16, "bayes", "c3_e1",8,1)
	//compare("minorc", 2, "xx", -1, 0, 16, "svd", "c3_e1",8,1)
	//compare("y3c", 1, "xx", -1, 0, 16, "bayes", "c3_e1",8,1)
	//compare("y3c", 2, "xx", -1, 0, 16, "svd", "c3_e1",8,1)
	
	char  varnm[200];
	char* stry = ttl2[ij];
	char * pchy = strchr(stry,'_');
	int len = pchy-stry;
	strncpy (varnm,stry,len);
	varnm[len]='\0';
	_file0->cd();
	switch(ij) {
	case 0 : sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, varnm, istrt, title, etapt); break;
	case 1 : sprintf(histname, "Roo_%s_%s_PFjet_%s", varnm,title, etapt); break;
	case 2 : sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, varnm, istrt, title, etapt); break; 
	case 3 : sprintf(histname, "Roo_%s_%s_PFjet_%s", varnm,title, etapt); break;
	case 4 : sprintf(histname, "Roo_%s_%s_Genjet_%s", varnm,title, etapt); break;
	default : break;
	}
	cout <<"histname "<<histname<<endl;
	break;*/
case 16 :

        if (ij==0) {
          //sprintf(histname, "bayes_Pythia8_Jets_%i_unfold_Roo_Madgh_reco_%s", istrt, etapt); //clouser plot
          sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt); //Diff algo Bayes, SVD
        //sprintf(histname,"refold_bayes_effi_Pythia8_Jets_%i_unfold_Roo_Pythia8_reco_%s",istrt, etapt); //refold
        }
         else if(ij==1) { // if (ij==0) {
          //sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
          //sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt0_eta0_3"); //bayes_Madgh_Jets_11_unfold_Roo_Pythia8_reco_typ_2_pt7_eta0_24
          sprintf(histname, "Roo_Pythia8_reco_%s", etapt); // Reold //Roo_Madgh_gen_typ_2_pt7_eta0_24
          //sprintf(histname, "svd_Data_Jets_9_unfold_Roo_Pythia8_reco_%s", etapt);
        //  sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt);
        }
        //else {sprintf(histname, "Roo_Madgh_gen_%s", etapt);}
        break;
    /*  case 17 :
	char  varnm[200];
	char* stry = ttl2[ij];
	char * pchy = strchr(stry,'_');
	int len = pchy-stry;
	strncpy (varnm,stry,len);
	varnm[len]='\0';
	_file0->cd();
	switch(ij) {
	case 1 : sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, varnm, istrt, title, etapt); break;
	case 2 : sprintf(histname, "Roo_%s_%s_PFjet_%s", varnm,title, etapt); break;
	case 3 : sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, varnm, istrt, title, etapt); break; 
	case 4 : sprintf(histname, "Roo_%s_%s_PFjet_%s", varnm,title, etapt); break;
	case 0 : sprintf(histname, "Roo_%s_%s_Genjet_%s", varnm,title, etapt); break;
	default : break;
   	  }
	cout <<"histname "<<histname<<endl;
	break;*/
	case 17 :

        if (ij==0) {
          sprintf(histname, "gen_mean_unolded_bayes_%s", etapt);
        }
         else if(ij==1) { // if (ij==0) {
          sprintf(histname, "gen_mean_Pythia8_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24
        }
        else if(ij==2) {
          sprintf(histname, "gen_mean_Madgh_%s", etapt);
          }
        else {
         sprintf(histname, "gen_mean_Herwigpp_%s", etapt);
         }
        break;

      case 18 : 
	// reco/unfolded MC with different jet smearing
	//  root ../unfold/tmp130206/outfile_ff3alliter_Pythia8_1_c3_e1.root
	//compare("thrustc", 1, "xx", -1, 0, 18, "svd", "c3_e1",9,2)
	//compare("thrustc", 2, "xx", -1, 0, 18, "svd", "c3_e1",9,1)

	//	_file0->cd(); 
	if (isUnfolding==1) { 
	  // jet smearing effect
	  if (ij==0) {
	    sprintf(histname, "Roo_Pythia6_%s_Genjet_%s", title, etapt); 
	  } else {
	    sprintf(histname, "%s_Pythia6_%s_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt); 
	  }
	} else {
	  if (ij==0) {
	    sprintf(histname, "Roo_Pythia6_%s_PFjet_%s", title, etapt); 
	  } else {
	    sprintf(histname, "Roo_Pythia6_%s_%s_%s", title, ttl2[ij], etapt);
	  }
	}
	break;
	
      case 19 : 
	// reco/unfolded Data with different jet energy scale
	//  root ../unfold/tmp130206/outfile_ff3alliter_Pythia8_1_c3_e1.root
	//compare("thrustc", 1, "xx", -1, 0, 19, "svd", "c3_e1",9,2)
	//compare("thrustc", 2, "xx", -1, 0, 19, "svd", "c3_e1",9,1)
	//	_file0->cd(); 
	if (isUnfolding==1) { 
	  // jet smearing effect
	  sprintf(histname, "%s_Data_%s_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt); 
	} else {
	  sprintf(histname, "Roo_Data_%s_%s_%s", title, ttl2[ij], etapt);
	}
	break;

      case 200 :
	switch(ij) {
        case 0 : _file0->cd(); break;
        case 1 : _file0->cd(); break;
        case 2 : _file0->cd(); break;
        case 3 : _file0->cd(); break;
        case 4 : _file0->cd(); break;
        default : _file0->cd(); break;
        }
	if (ij==0) {
          sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt);
        } else {
          sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
        } 

        break;
	
      case 20 :  //    Unfolded data and Genrator MC	//7TeV
	switch(ij) {
	case 0 : _file0->cd(); 
                 cout << "FILE0=========="<<_file0->GetName() << endl;
                 //cout << _file0->IsOpen() << endl;
                 break;
	case 1 : _file1->cd(); 
                 cout <<"FILE1=========" <<_file1->GetName() << endl;
                 break;
//	case 2 : _file2->cd(); break;
//	case 3 : _file3->cd(); break;
//	case 4 : _file4->cd(); break;
/*	case 5 : _file5->cd(); break;
	case 6 : _file6->cd(); break;
	case 7 : _file7->cd(); break;	  
	case 8 : _file8->cd(); break;
	case 9 : _file9->cd(); break;	  
	case 10 : _file10->cd(); break;
	case 11 : _file11->cd(); break;	 
        case 12 : _file12->cd(); break;
        case 13 : _file13->cd(); break; */
	default : _file0->cd(); break;
	}

	if (ij==0) {
        //sprintf(histname,"bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt0_eta0_3"); 
        //sprintf(histname,"bayes_Data_Jets_%i_unfold_Roo_Pythia8_%s", istrt, etapt); 
       
   sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt);
//   sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Herwigpp_reco_%s", istrt, etapt); //Unfolded by Herwigpp
     
        //sprintf(histname, "Roo_Pythia8_gen_%s", etapt);
        // sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Pythia8_reco_%s", istrt, etapt);
            //sprintf(histname, "analyzeBasicPat/recojetallave_pt_0");
          // sprintf(histname, "bayes_Data_Jets_3_unfold_Roo_Pythia8_reco_typ_2_pt3_eta0_3");
          // sprintf(histname,"bayes_Data_Jets_5_unfold_Roo_Pythia8_reco_typ_2_pt4_eta0_15");
         cout<< "HistName 1=" <<histname<<endl;
	} 
       // else if(ij==1) sprintf(histname, "Roo_Pythia8_gen_%s", etapt);
          else if(ij==1) {
            sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Madgh_reco_%s", istrt, etapt);
         cout<< "HistName 2=" <<histname<<endl;
         cout << "What the f " << endl;
        }

        //else if(ij==2) sprintf(histname, "Roo_default_gen_%s", etapt);
 //       else if(ij==3) sprintf(histname, "Roo_MultipartonInteractions_pT0Ref_23_gen_%s", etapt);
        //else if(ij==4) sprintf(histname, "Roo_MultipartonInteractions_pT0Ref_2_gen_%s", etapt);
        //else if(ij==3) sprintf(histname, "Roo_MultipartonInteractions_ecmPow_23_gen_%s", etapt);
       // else if(ij==6) sprintf(histname, "Roo_MultipartonInteractions_ecmPow_19_gen_%s", etapt);
        //else if(ij==3) sprintf(histname, "Roo_MultipartonInteractions_expPow_2.2_gen_%s", etapt);
        //else if(ij==8) sprintf(histname, "Roo_MultipartonInteractions_expPow_17_gen_%s", etapt);
//        else if(ij==2) sprintf(histname, "Roo_mpioff_gen_%s", etapt);
//        else if(ij==3) sprintf(histname, "Roo_isr_off_gen_%s", etapt);
          

 /*      else {
//	  sprintf(histname, "bayes_Data_Jets_3_unfold_Roo_Pythia8_reco_typ_2_pt3_eta0_3");
         // sprintf(histname,"svd_Pythia8_Jets_9_unfold_Roo_Pythia8_reco_typ_2_pt4_eta0_4");   
         //sprintf(histname, "bayes_Data_Jets_4_unfold_Roo_Pythia8_reco_typ_0_pt4_eta0_3");
         // cout<< "HistName=" <<histname<<endl;
        //sprintf(histname,"refold_bayes_effi_Pythia8_Jets_5_unfold_Roo_Pythia8_reco_typ_0_pt4_eta0_3"); 
         // sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Madgh_reco_%s", istrt, etapt);
          sprintf(histname, "bayes_Data_Jets_%i_unfold_Roo_Herwigpp_reco_%s", istrt, etapt); //clouser plot
          
        //sprintf(histname, "Roo_ TimeShower_pTmin_6_gen_%s", etapt);
        //sprintf(histname, "Roo_fsr_off_gen_%s", etapt);
        //sprintf(histname, "Roo_MultipartonInteractions_expPow_17_gen_%s", etapt);
        //sprintf(histname, "Roo_MultipartonInteractions_pT0Ref_2_gen_%s", etapt);
        //sprintf(histname, "Roo_MultipartonInteractions_ecmPow_19_gen_%s", etapt);
        //sprintf(histname, "Roo_MultipartonInteractions_expPow_17_gen_%s", etapt);
      //    sprintf(histname, "Roo_Madgh_gen_%s", etapt); //Roo_Madgh_gen_typ_2_pt7_eta0_24

        //sprintf(histname, "Roo_fsr_off_gen_%s", etapt);
         cout<< "HistName 3=" <<histname<< " ; " << ij << endl;
        //sprintf(histname, "Roo_BeamRemnants_reconnectRange_22_gen_%s", etapt);
        //sprintf(histname, "Roo_StringZ_aLund_9_gen_%s", etapt);
       //sprintf(histname, "Roo_StringZ_bLund_12_gen_%s", etapt);
            //sprintf(histname, "recojetallave_pt_0");
	}
	break;*/
	/*		
      case 20 : // for madgrapg comp
	switch(ij) {
        case 0 : _file0->cd(); break;
        case 1 : _file1->cd(); break;
        case 2 : _file2->cd(); break;
        case 3 : _file3->cd(); break;
        case 4 : _file4->cd(); break;
        case 5 : _file5->cd(); break;
        case 6 : _file6->cd(); break;
	case 7 : _file7->cd(); break;
	case 8 : _file8->cd(); break;
        case 9 : _file9->cd(); break;
        case 10 : _file10->cd(); break;
        case 11 : _file11->cd(); break;

        default : _file0->cd(); break;
        }

        if (ij==0) {
          sprintf(histname, "%s_Data_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, istrt, title, etapt);
	  //        } else if (ij==1) {
	  //          sprintf(histname, "Roo_Madgraph_%s_Genjet_%s", title, etapt);
	  //	  sprintf(histname, "%s_PFjet_%s", title, etapt);
        } else {
          sprintf(histname, "%s_Genjet_%s", title, etapt);
        }
        break;
	*/

	case 21 :  //    Unfolded data and Genrator MC    //7TeV
        switch(ij) {
        case 0 : _file0->cd(); break;
        case 1 : _file1->cd(); break;
        case 2 : _file2->cd(); break;
        case 3 : _file3->cd(); break;
        case 4 : _file4->cd(); break;
        case 5 : _file5->cd(); break;
        case 6 : _file6->cd(); break;
        case 7 : _file7->cd(); break;
        case 8 : _file8->cd(); break;
        case 9 : _file9->cd(); break;
        case 10 : _file10->cd(); break;
        case 11 : _file11->cd(); break;
        case 12 : _file12->cd(); break;
        case 13 : _file13->cd(); break;
        case 14 : _file14->cd(); break;
        case 15 : _file15->cd(); break;
        case 16 : _file16->cd(); break;
        case 17 : _file17->cd(); break;
        case 18 : _file18->cd(); break;
        case 19 : _file19->cd(); break;
        case 20 : _file20->cd(); break;
        case 21 : _file21->cd(); break;
        case 22 : _file22->cd(); break;
        case 23 : _file23->cd(); break;
       	default : _file0->cd(); break;
        }

        if (ij==0) {
         sprintf(histname, "%s_Genjet_%s", title, etapt);
//          sprintf(histname, "%s_%s_PFjet_%i_unfold_Roo_Pythia6_%s_PFjet_%s", unf, ttl2[ij], istrt, title, etapt);
//          sprintf(histname, "%s_Genjet_%s", title, etapt);
//      } else if (ij==1 || ij==5 || ij==9) {
//        sprintf(histname, "Roo_%s_%s_Genjet_%s", ttl2[ij], title, etapt);
        } else {
          sprintf(histname, "%s_Genjet_%s", title, etapt);
        }
        break;	
      default :
	sprintf(histname, "%s_c%i_e1", title, ij); break;
      }
    } else { 
      if (ifile==0) {
	if (isUnfolding==1) {
	  //compare("7_unfold_Roo_PFJet_thrustc_PFJet_c1_e1", 1, "xx", -1, 33);
	  _file0->cd(); sprintf(histname, "svd_%s_%s", ttl2[ij],title);
	} else {
	  switch(ij) {
	  case 0 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  case 1 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  case 2 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  case 3 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  case 4 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  case 5 : _file0->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	  default : _file0->cd();sprintf(histname, "%s_%s_c1_e1", title, ttl2[0]); break;
	  }
	}
      } else {
	switch(ij) {
	case 0 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	case 1 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	case 2 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	case 3 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	case 4 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;
	case 5 : _file1->cd(); sprintf(histname, "%s_%s_c1_e1", title, ttl2[ij]); break;	  
	default : _file1->cd();sprintf(histname, "%s_%s_c1_e1", title, ttl2[0]); break;
	}
      }
    }
    if (isDoublerat && ij==0) { sprintf(histname2, "%s", histname);}

//TDirectory* dir = _file0->GetDirectory("qcdevt");
// dir->cd();
//fdfz[ij] = (TH1F*) _file1->Get(histname);
// fdfz[ij] = (TH1F*)dir->Get(histname);
//fdfz[ij] = (TH1F*) _file1->GetList()->FindObject("qcdevt");
    fdfz[ij] = (TH1F*)gDirectory->Get(histname); // title);
        cout <<"yxxxx "<<histname<< endl;
        cout<<"ijjjjjjjjjjjjjj"<<ij<<endl;
//TCanvas *c21 = new TCanvas("c21", "multipads1", 900, 700);
         
    //    fdfz[ij]->GetXaxis()->SetRangeUser(0,-7);
    //ordering of scale changes on 5th April 2013
    //if (irebin>1) fdfz[ij]->Rebin(irebin);
    if (ij==0 && nbinunf==0) {
     cout<<"Bins====================="<<endl;
      nbinunf = fdfz[ij]->GetNbinsX();
      for (int ix=0; ix<=nbinunf; ix++) {
  	binunfs[ix] = fdfz[ij]->GetBinLowEdge(ix+1);
      }
      fdfx[ij] = (TH1F*)fdfz[ij]->Clone();
    } else {
      int tmpnbin = fdfz[ij]->GetNbinsX();
      if (nbinunf!=tmpnbin) {
  	fdfx[ij] = rebin_hist( fdfz[ij], nbinunf, binunfs);
      } else {
  	fdfx[ij] = (TH1F*)fdfz[ij]->Clone();
      }
    }
    for (int ix=0; ix<=fdfx[ij]->GetNbinsX(); ix++) {
            cout << "ij= " << ij << " nbin= " << fdfx[ij]->GetNbinsX() << " ix= " << ix << " binledge " << fdfx[ij]->GetBinLowEdge(ix+1) << " bincontent "<< fdfx[ij]->GetBinContent(ix+1) << "Bin Error= " << fdfx[ij]->GetBinError(ix+1)/(fdfx[ij]->GetBinContent(ix+1)) <<endl;
  }
     
    /*
    for (int ix=0; ix<fdfx[ij]->GetNbinsX(); ix++) {
      //      cout << "11111111111111111111111  " << fdfx[ij]->GetBinContent(ix+1) << endl;
      binentry[ix][ij] = fdfx[ij]->GetBinContent(ix+1);
      binerror[ix][ij] = fdfx[ij]->GetBinError(ix+1);
      if (ij>0){
	double ratio = binentry[ix][0]/binentry[ix][ij];
	cout << "ratio " << ratio << " " << binentry[ix][ij] << " " << binentry[ix][0] << endl;
	errormc_matrix.ResizeTo(nbinx,nbinx);
	cout << "binerror " << binerror[ix][ij] << endl;
	//	errormc_matrix(ix,ix) = 0.0; // for wothout stat. error
	errormc_matrix(ix,ix) = pow(ratio*binerror[ix][ij],2);
	//	cout << "errormc_matrix(ix,ix)  " << errormc_matrix(ix,ix) << endl; 
      }
    }
    */
    //GMA to remove last bin
//    fdfx[ij] = rebin_histwithoutlast( fdfz[ij], nbinunf, binunfs);
//    fdfx[ij] = (TH1F*)fdfz[ij]->Clone();
//    cout << "model name " << ttl2[ij] << " Integral " << fdfx[ij]->Integral() << endl;
cout<<"yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy"<<endl;

    if (isUnitScale) { 
      scalex[ij] =  1./(fdfx[ij]->Integral());
    } else {
      if (ij==0) {
 	scalex[ij] = 1.;
 	fact = fdfz[ij]->Integral();
	
      } else {
	scalex[ij] =fact*1.0/(TMath::Max(fdfz[ij]->Integral(), 0.0001));
      }
    }
    double tmpwid=1;
    if (isUnfolding>-10) { 
      int tmpnbn = fdfx[ij]->GetNbinsX();
      for (int ix=0; ix<tmpnbn; ix++) {
	double awidth = fdfx[ij]->GetBinWidth(ix+1); // /tmpwid;
	fdfx[ij]->SetBinContent(ix+1, fdfx[ij]->GetBinContent(ix+1)/awidth);
	double error = fdfx[ij]->GetBinError(ix+1);
	fdfx[ij]->SetBinError(ix+1, error/awidth);
      }
   }
  

  meanerr[ij]=fdfx[ij]->GetMean(); 
  if (isExternalError==1 && ij==1) { 
   double mean1=meanerr[0];
   double mean2=meanerr[1];

  //cout << "Mean Error = " << " ; " << meanerr[0]  << endl;
    fileout<<"double trackmeanhist_"<<etapt<<"["<<8<<"] ={";
    fileout <<meanerr[0]-meanerr[1];
    fileout <<"};"<<endl;
 

   }


    //    fdfx[ij]->Scale(1./fdfx[ij]->Integral());
    //    for (int jk=0; jk<nbinx; jk++) {
    //      cout << "binentry " << fdfx[ij]->GetBinContent(jk+1) << " binerror " << fdfx[ij]->GetBinError(jk+1) << endl;
    //    }
    for (int jk=0; jk<nbinx; jk++) {
      //      cout <<  "binentry = " << fdfx[ij]->GetBinContent(jk+1) << "binerror = " << fdfx[ij]->GetBinError(jk+1) << endl;
    }
    fdfx[ij]->Scale(scalex[ij]);
    cout << "scale  " << scalex[ij] << endl;
    cout<<"yes"<<endl;
    nbinx = fdfx[ij]->GetNbinsX(); 
    minbincontent[ij]=int(1.e9);
    for (int jk=0; jk<nbinx; jk++) { 
      if (fdfx[ij]->GetBinContent(jk+1) >0 && fdfx[ij]->GetBinContent(jk+1)<minbincontent[ij]) { minbincontent[ij] = fdfx[ij]->GetBinContent(jk+1);} 

      binentry[jk][ij] = fdfx[ij]->GetBinContent(jk+1);
      binerror[jk][ij] = fdfx[ij]->GetBinError(jk+1);
      


            //fileout << "binentry we are inside " <<binentry[jk][ij] << " binerror " << binerror[jk][ij]/binentry[jk][ij] << endl; 
    }

    for (int ix=0; ix<fdfx[ij]->GetNbinsX(); ix++) {
      //      cout << "11111111111111111111111  " << fdfx[ij]->GetBinContent(ix+1) << endl;  
     // binentry[ix][ij] = fdfx[ij]->GetBinContent(ix+1);
     // binerror[ix][ij] = fdfx[ij]->GetBinError(ix+1);
      //cout<<"Tanmay"<<ij<<endl;

      if (ij>0){
        double ratio = binentry[ix][0]/binentry[ix][ij];
	cout << "ratio " << ratio << " " << binentry[ix][ij] << " " << binentry[ix][0] << endl;

//         cout<<"fine"<<endl;
        errormc_matrix.ResizeTo(nbinx,nbinx);
	//        cout << "binerror " << binerror[ix][ij] << endl;
        //      errormc_matrix(ix,ix) = 0.0; // for wothout stat. error  
        errormc_matrix(ix,ix) = pow(ratio*binerror[ix][ij],2);
        //      cout << "errormc_matrix(ix,ix)  " << errormc_matrix(ix,ix) << endl;   
      }
    }

    //    if (ij <2){
    //      file_out << " *************** THIS IS FOR " <<  ttl2[ij] << "  *******************" << endl;
    //	for (int ix=0; ix<fdfx[ij]->GetNbinsX(); ix++) {
    //	  file_out << " bin center " << fdfx[ij]->GetBinCenter(ix+1) << " bin content " << fdfx[ij]->GetBinContent(ix+1) << " bin error " << fdfx[ij]->	  GetBinError(ix+1) << endl;
    //	}
    //    }
    char  varnm[200];
    if (idx>=0) {
      char* strx = title;
      char * pchx = strchr(strx,'_c');
      int len = pchx-strx;
      strncpy (varnm,strx,len);
      varnm[len]='\0';
      
    }

    if (ij==0) {
      if (isitForPaper>0) { 
	//      if (isitForPaper==0) { 

	char* str = fdfx[ij]->GetName();
	int ix=0;
	for (ix=3; ix<nvar; ix++) {
	  if (strstr(str,varname[ix])) {
	    sprintf(pch1, "log_{e} (%s)", vartitle[ix]);  break;
	  }
	}
	
	if (ix==nvar) {
	  char* str1 = fdfx[ij]->GetTitle();
	  char * pchx = strchr(str1,'_');
	  int len = pchx-str1;
	  strncpy (pch1,str1,len);
	  pch1[len]='\0';
	  if (strstr(str1, "Pt") && (!strstr(str1, "sin"))) {  sprintf(pch1, "%s (GeV/c)", pch1);}
	  len=strlen(pch1);
	  pch1[len]='\0';
	}
      } else {
	char* str = fdfx[ij]->GetTitle();
	int len = strlen(str);
	
	if (isUnfolding==1) {
	  int lnx=0;
	  char* pch = strstr(str, "Roo_Pythia6_log"); 
	  if (pch<=0) pch = strstr(str, "Roo_Pythia8_log"); 
	  //	if (pch<=0) pch = strstr(str, "Roo_Herwig_log"); 
	  if (pch<=0) { lnx=1; pch = strstr(str, "Roo_Madgraph_log");}
	  if (pch<=0) pch = strstr(str, "Roo_Pythia6"); //only for cparamter
	  if (pch<=0) pch = strstr(str, "Roo_Pythia8");
	  if (pch<=0) pch = strstr(str, "Roo_Madgraph");
	  char* pxx = strstr(str, "_PFjet_c"); 
	  //	  cout <<pxx<<endl;
	  int len1 =pch-str+11;
	  if (istrt>9) len1++;
	  int len2=pch-pxx+1;
	  strncpy (pch1,str+len1, len-len1);
	  //	strncpy (pch1,str+len1, len2);
	  len1 = strlen(pch1);
	  pch1[len1]='\0';
	} else {
	  if (idx>=0) {
	  } else {
	    if (isame<=0) {
	      len = strstr(str, "_c");   
	    } else {
	      len = strstr(str, "__");   
	    }
	  }
	  pch1[len1]='\0';
	}
      }
    }
    
    if (icol==1) {
      nbinx = fdfx[ij]->GetNbinsX();
      fdfx[ij]->SetMarkerColor(1);
      mean = TMath::Max(1.,fdfx[ij]->GetMean());
      rms = TMath::Max(1.,fdfx[ij]->GetRMS());
      double mxall = fdfx[ij]->GetMaximum();
      for (int jk=0; jk<nbinx; jk++){
	////	cout <<"jk "<< jk<<" "<<fdfx[ij]->GetBinContent(jk+1)<<" "<<scalex[ij]<<" "<<mxall<<endl;
	if (fdfx[ij]->GetBinContent(jk+1) >= 0.000005*mxall) { // scalex[ij]) {
	  //	  alow = fdfx[ij]->GetBinCenter(jk+1)-0.5*fdfx[ij]->GetBinWidth(jk+1);
	  alow = fdfx[ij]->GetBinCenter(jk+1)-0.5*fdfx[ij]->GetBinWidth(jk+1)+1.e-10;
	  //	  cout <<"alowwwwwwwwwwwwwwwwwwwwww "<<alow<<" "<<fdfx[ij]->GetBinCenter(jk+1)<<" "<<fdfx[ij]->GetBinWidth(jk+1)<<" "<<fdfx[ij]->GetBinLowEdge(jk+1)<<endl;

	  break;
	}
      }
      for ( int jk=nbinx-1; jk>=0; jk--){
	if (fdfx[ij]->GetBinContent(jk+1) >= 0.000005*mxall) { // *scalex[ij]) {
	  //	  ahg = fdfx[ij]->GetBinCenter(jk+1)+0.5*fdfx[ij]->GetBinWidth(jk+1);
	  ahg = fdfx[ij]->GetBinCenter(jk+1)+0.5*fdfx[ij]->GetBinWidth(jk+1)-1.e-10;
	  //	  cout <<"ahg "<<ahg<<endl;
	  break;
	}
      } 
      
      //alow = fdfx[ij]->GetXaxis()->GetXmin()+fdfx[ij]->GetBinWidth(1);
      alow = fdfx[ij]->GetXaxis()->GetXmin();
      ahg = fdfx[ij]->GetXaxis()->GetXmax();
      cout <<"aloowww "<< alow<<" "<<ahg<< " binwidth " << fdfx[ij]->GetBinWidth(1)<<endl;
      //      alow = 250;
      //      ahg=550.;
       //CMS_lumi(c1,iPeriod,iPos);
      if (mRatio>=0) {
	
	fdfx[ij]->SetMarkerStyle(20);
	fdfx[ij]->SetMarkerSize(0.75);
	fdfx[ij]->GetXaxis()->SetTickLength(.075);
	fdfx[ij]->GetYaxis()->SetTickLength(.05);

	fdfx[ij]->GetXaxis()->SetRangeUser(alow, ahg);
	if (gStyle->GetOptLogy()==0) fdfx[ij]->SetMinimum(-0.02);
	

	fdfx[ij]->GetXaxis()->SetTitle(pch1);
	fdfx[ij]->GetXaxis()->SetTitleSize((idx<0) ? 0.08 : 0.40); //(.08); //pss
	fdfx[ij]->GetXaxis()->SetTitleOffset((idx<0) ? 0.92 : 0.96);
	fdfx[ij]->GetXaxis()->CenterTitle();
	
	fdfx[ij]->GetXaxis()->SetLabelFont(42);
	fdfx[ij]->GetXaxis()->SetTitleFont(42);
	fdfx[ij]->GetXaxis()->SetLabelSize(0.058);
	fdfx[ij]->GetXaxis()->SetLabelOffset(.01);

	if (mRatio==0) 	{
          fdfx[ij]->GetXaxis()->SetLabelOffset(-.005);
	  fdfx[ij]->GetXaxis()->SetTitleSize(0.055); //pss
	  fdfx[ij]->GetXaxis()->SetTitleOffset(0.93);
	}

	if (isUnitScale) {
         // char ytitle[200]= {"H_{T2}"};
         //char ytitle[200]=  {"Pt of leading jet (GeV/c)"};
        // char ytitle[200]=  {"Pt of second leading jet (GeV/c)"}; 
        // char ytitle[200]= {"#eta of leading jet"}; 
         //char ytitle[200]= {"#eta of second leading jet"}; 
         //char ytitle[200]= {"#phi of leading jet"}; 
         //char ytitle[200]= {"#phi of second leading jet"}; 
         //char ytitle[200]= {"#phi jet"}; 
	  // char ytitle[200]={"#Delta Pt of two leading jets (GeV/c)"};
         //char ytitle[200]=  {"#Delta#phi of Jets"};
         //char ytitle[200]= {"No. of jet"}; 
         //char ytitle[200]={"Pt2 x sin( #Delta #phi )/Pt1"};
         //char ytitle[200]=  {"#Delta#phi of Jets"};
         //sprintf(fname3, "1/N dN/d%s", pch1);
           //sprintf(fname3, "1/N dN/d%s", ytitle);
           if(basic_title) sprintf(fname3, "1/N dN/d%s", basic_vartitle[titleindex]);
           else sprintf(fname3, "1/N dN/d%s", vartitle[titleindex]);
	} else {
	  sprintf(fname3, "dN/d%s",  pch1);
	}

	if (!noYtitle) {
	  fdfx[ij]->GetYaxis()->SetTitle(fname3); 
	  double extrasize=(isSlide&& mRatio!=0) ? 0.03 : (gStyle->GetOptLogy()) ? 0.00 : 0.0;

	  //	  fdfx[ij]->GetYaxis()->SetTitleSize((mRatio!=0) ? extrasize+0.0366+0.0033*ndata2 :(idx<0) ? 0.07 : 0.06); //0.07 : 0.05; //pss
	  fdfx[ij]->GetYaxis()->SetTitleSize((mRatio!=0) ? extrasize+0.0514+0.0005*ndata2 :(idx<0) ? 0.07 : 0.05);
	  if (isUnitScale) {
	    fdfx[ij]->GetYaxis()->SetTitleOffset(1.0);
	  }else{
	    fdfx[ij]->GetYaxis()->SetTitleOffset((mRatio!=0) ? 1.951833-0.0433*ndata2-50*extrasize : 1.17);
	  ///	  fdfx[ij]->GetYaxis()->SetTitleOffset((mRatio!=0) ? 1.833-0.0633*ndata2-15*extrasize : 1.17); //0.0833*ndata2 0.97;
	  }
	  fdfx[ij]->GetYaxis()->CenterTitle();
	  ///	  fdfx[ij]->GetYaxis()->SetLabelSize((mRatio!=0) ? extrasize+0.03+0.003*ndata2: 0.06); 
	  fdfx[ij]->GetYaxis()->SetLabelSize((mRatio!=0) ? extrasize+0.040+0.001*ndata2: 0.05);
          

          
	  fdfx[ij]->GetYaxis()->SetLabelOffset(.001);
	  fdfx[ij]->GetYaxis()->SetLabelFont(42);
	  fdfx[ij]->GetYaxis()->SetTitleFont(42);
	  //	  fdfx[ij]->GetYaxis()->TGaxis::SetMaxDigits(4);
	}
      }
    } else { 
      if (istatoth>0) {
	gStyle->SetOptStat(1110);
	gStyle->SetStatX(statxpos);
	double yy = (idx>=0 ) ? 0.99-ij*0.10 : (ndata2<=3) ? 0.99-ij*0.17 : 0.99-ij*0.15;;
	gStyle->SetStatY(yy);
	gStyle->SetStatColor(icol);
      }

      chisqq=10;
      ndff = 0;
      int good = 0;
           // if (isrootchi>0 && mRatio>=0) fdfx[0]->Chi2TestX(fdfx[ij], chisqq, ndff, good, "UW"); //,res);
    }
   
    //setTDRStyle(gStyle);
       //CMS_lumi(c1,iPeriod,iPos);
    if (mRatio>=0) {
      fdfx[ij]->SetLineColor(icol);
      if(icol==8) fdfx[ij]->SetLineColor(icol+20);
      fdfx[ij]->SetLineWidth(3.0);
      //      if (isExternalError<=2) { fdfx[ij]->SetLineStyle(linestl++); }
            if (isExternalError<=2) { 
             if(ij==1) fdfx[ij]->SetLineStyle(1);
             else fdfx[ij]->SetLineStyle(icol); 
             }
      //if (isExternalError<=2) { fdfx[ij]->SetLineStyle(0); }
      if (ij==0) fdfx[ij]->SetMaximum(anorm*fdfx[ij]->GetMaximum());
    }
    if (isExternalError<2) {
      if (mRatio>=0){
       fdfx[ij]->SetLineColor(icol);
       file_out<< "Draw time  ij = " << ij << " ; icol = "<< icol << endl;
      //else fdfx[ij]->SetLineColor(5);
       //setTDRStyle(gStyle);
       //fdfx[ij]->GetXaxis()->SetRangeUser(-4.5, 0); //off previously
       //cout << "iPeriod ======" << iPeriod << "iPos ===  =" << iPos << endl;
       fdfx[ij]->Draw((icol==1)?"E": (istatoth>0) ? "sames:hist" : "same:hist");
       //fdfx[ij]->Draw((icol==1)?"E": (istatoth>0) ? "same:hist" : "same");
       //fdfx[ij]->Draw((ij==1 && isExternalError!=2) ? "hist" : "same:hist");
       //if(ij==3) fdfx[ij]->Draw("sames:hist");
       //if(ij==2)CMS_lumi(c1,iPeriod,iPos);
       //c1->Update();
  //c1->RedrawAxis();
  //c1->GetFrame()->Draw();
       }  
   } else {
      fdfy[ij] = (TH1F*)fdfx[ij]->Clone();
      //fdfy[ij]->GetXaxis()->SetRangeUser(-8.1., -1.7.);
      if (icol==1) {
	if (mRatio>=0) fdfy[ij]->Draw("E1xx");
//totsyserr_18_1_3_1
      if(typindex==0){	
	if (strstr(title,"thrustc")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_0_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_1_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_2_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_3_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_4_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_5_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_6_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_0_7_1[ix];}
	  }
	} else if (strstr(title,"rhotot")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_9")) {
            cout <<"RHOTOT  ============= " << endl;
            cout <<"RHOTOT  ============= " << endl;
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_0_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_1_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_2_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_3_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_4_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_5_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_6_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_0_7_1[ix];}
	  }
	} else if (strstr(title,"broadt")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_18")) {
             cout << "Broadt ======= " << endl;
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_0_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_1_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_2_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_3_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_4_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_5_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_6_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_0_7_1[ix];}
	  }
	} else if (strstr(title,"rhottot")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_24")) {
            cout<< " RHOTTOT ==" << endl;
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_0_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_1_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_2_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_3_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_4_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_5_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_6_1[ix];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_0_7_1[ix];}
	  }

	} else {
	  cout <<"Data : Variable name is not matching with the name "<<endl;
	  break;
	}
      }
    
      if(typindex==1){ 
      if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_1_pt0_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_0_1[ix];}
          } else if (strstr(etapt,"typ_1_pt1_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_1_1[ix];}
          } else if (strstr(etapt,"typ_1_pt2_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_2_1[ix];}
          } else if (strstr(etapt,"typ_1_pt3_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_3_1[ix];}
          } else if (strstr(etapt,"typ_1_pt4_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_4_1[ix];}
          } else if (strstr(etapt,"typ_1_pt5_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_5_1[ix];}
          } else if (strstr(etapt,"typ_1_pt6_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_6_1[ix];}
          } else if (strstr(etapt,"typ_1_pt7_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_1_7_1[ix];}
          } 
        } else if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_9")) {
            cout <<"RHOTOT  ============= " << endl;
            cout <<"RHOTOT  ============= " << endl;
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_0_1[ix];}
          } else if (strstr(etapt,"typ_1_pt1_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_1_1[ix];}
          } else if (strstr(etapt,"typ_1_pt2_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_2_1[ix];}
          } else if (strstr(etapt,"typ_1_pt3_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_3_1[ix];}
          } else if (strstr(etapt,"typ_1_pt4_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_4_1[ix];}
          } else if (strstr(etapt,"typ_1_pt5_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_5_1[ix];}
          } else if (strstr(etapt,"typ_1_pt6_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_6_1[ix];}
          } else if (strstr(etapt,"typ_1_pt7_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_1_7_1[ix];}
          }
          } else if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_1_pt0_eta1_18")) {
             cout << "Broadt ======= " << endl;
            for (int ix=0; ix<nbinx; ix++) {
                cout << "Broadt Nbinsx ====" << nbinx << endl;
syserrData[ix] = totsyserr_18_1_0_1[ix];
}
          } else if (strstr(etapt,"typ_1_pt1_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_1_1[ix];}
          } else if (strstr(etapt,"typ_1_pt2_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_2_1[ix];}
          } else if (strstr(etapt,"typ_1_pt3_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_3_1[ix];}
          } else if (strstr(etapt,"typ_1_pt4_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_4_1[ix];}
          } else if (strstr(etapt,"typ_1_pt5_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_5_1[ix];}
          } else if (strstr(etapt,"typ_1_pt6_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_6_1[ix];}
          } else if (strstr(etapt,"typ_1_pt7_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_1_7_1[ix];}
          }
         } else if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_24")) {
            cout<< " RHOTTOT TYP1==" << endl;
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_0_1[ix];}
          } else if (strstr(etapt,"typ_1_pt1_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_1_1[ix];}
          } else if (strstr(etapt,"typ_1_pt2_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_2_1[ix];}
          } else if (strstr(etapt,"typ_1_pt3_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_3_1[ix];}
          } else if (strstr(etapt,"typ_1_pt4_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_4_1[ix];}
          } else if (strstr(etapt,"typ_1_pt5_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_5_1[ix];}
          } else if (strstr(etapt,"typ_1_pt6_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_6_1[ix];}
          } else if (strstr(etapt,"typ_1_pt7_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_1_7_1[ix];}
          }
       

       }else {
          cout <<"Data : Variable name is not matching with the name "<<endl;
          break;
        }
      }

     if(typindex==2){ 
     if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_2_pt0_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_0_1[ix];}
          } else if (strstr(etapt,"typ_2_pt1_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_1_1[ix];}
          } else if (strstr(etapt,"typ_2_pt2_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_2_1[ix];}
          } else if (strstr(etapt,"typ_2_pt3_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_3_1[ix];}
          } else if (strstr(etapt,"typ_2_pt4_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_4_1[ix];}
          } else if (strstr(etapt,"typ_2_pt5_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_5_1[ix];}
          } else if (strstr(etapt,"typ_2_pt6_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_6_1[ix];}
          } else if (strstr(etapt,"typ_2_pt7_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_3_2_7_1[ix];}
          }
          } else if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_9")) {
            cout <<"RHOTOT  ============= " << endl;
            cout <<"RHOTOT  ============= " << endl;
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_0_1[ix];}
          } else if (strstr(etapt,"typ_2_pt1_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_1_1[ix];}
          } else if (strstr(etapt,"typ_2_pt2_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_2_1[ix];}
          } else if (strstr(etapt,"typ_2_pt3_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_3_1[ix];}
          } else if (strstr(etapt,"typ_2_pt4_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_4_1[ix];}
          } else if (strstr(etapt,"typ_2_pt5_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_5_1[ix];}
          } else if (strstr(etapt,"typ_2_pt6_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_6_1[ix];}
          } else if (strstr(etapt,"typ_2_pt7_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_9_2_7_1[ix];}
          }
         } else if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_2_pt0_eta1_18")) {
             cout << "Broadt ======= " << endl;
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_0_1[ix];}
          } else if (strstr(etapt,"typ_2_pt1_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_1_1[ix];}
          } else if (strstr(etapt,"typ_2_pt2_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_2_1[ix];}
          } else if (strstr(etapt,"typ_2_pt3_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_3_1[ix];}
          } else if (strstr(etapt,"typ_2_pt4_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_4_1[ix];}
          } else if (strstr(etapt,"typ_2_pt5_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_5_1[ix];}
          } else if (strstr(etapt,"typ_2_pt6_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_6_1[ix];}
          } else if (strstr(etapt,"typ_2_pt7_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_18_2_7_1[ix];}
          }
         } else if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_24")) {
            cout<< " RHOTTOT Typ2==" << endl;
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_0_1[ix];}
          } else if (strstr(etapt,"typ_2_pt1_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_1_1[ix];}
          } else if (strstr(etapt,"typ_2_pt2_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_2_1[ix];}
          } else if (strstr(etapt,"typ_2_pt3_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_3_1[ix];}
          } else if (strstr(etapt,"typ_2_pt4_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_4_1[ix];}
          } else if (strstr(etapt,"typ_2_pt5_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_5_1[ix];}
          } else if (strstr(etapt,"typ_2_pt6_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_6_1[ix];}
          } else if (strstr(etapt,"typ_2_pt7_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {syserrData[ix] = totsyserr_24_2_7_1[ix];}
          }

     }else {
          cout <<"Data : Variable name is not matching with the name "<<endl;
          break; 
        }
      }
	//	cout <<"nbinx "<<nbinx<<endl;
	//	file_out << " ********************* THE SYSTEMATIC ERROR OF DATA **************** " << endl;
        //file_stat_err<<"double stat_"<<etapt<<"["<<nbinx<<"] ={";
	for (int kl=0; kl<nbinx; kl++) {
	  //double yy = fdfy[ij]->GetBinError(kl+1);
         //if(kl<nbinx-1) file_stat_err <<yy<<", ";  
         //else file_stat_err <<yy<<"};"<< endl;  
	  //	  cout <<"errordata "<< kl<<" "<<syserrData[kl]<<endl;
	  //	  file_out << " " << syserrData[kl] << " "; 
	}
       //file_stat_err<<"};"<<endl;
/*
file_out<<"double jecvarpe_"<<ij<<"_"<<jk<<"_"<<mn<<"_"<<kl<<"["<<nbins+1<<"] ={";
        for (int ix=1; ix<nbins+1; ix++) {file_out <<poserr[ix]/max(1.e-6,mean[ix])<<", ";}
        file_out<<poserr[nbins]/max(1.e-6,mean[nbins])<<"};"<<endl;
        cout<<"jecposerr=" <<poserr[nbins]/max(1.e-6,mean[nbins]) << endl;
*/
	file_out << endl;
	for (int kl=0; kl<nbinx; kl++) {
	  double xx = fdfy[ij]->GetBinContent(kl+1);
	  double yy = fdfy[ij]->GetBinError(kl+1);
 
	  double errr = sqrt(yy*yy + pow(xx* syserrData[kl], 2.0));
	  fdfy[ij]->SetBinError(kl+1, errr); // fdfx[ij]->GetBinContent(kl+1)*syserrData[kl]);
          cout << "yes we solvrd it       ===============kl" << kl << yy << endl;
	}
    
	if (isAida) {
	  //	  for aida file
	  _filex->cd();
	  fdfy[ij]->Write();
	  _file0->cd();
	}

	//19th March 2013	if (mRatio) {
	  fdfy[ij]->SetFillColor(5); //4); //Green); //898); //5); 
	  fdfy[ij]->SetFillStyle(1111); 
	  fdfy[ij]->SetLineWidth(1);
	  fdfy[ij]->SetMarkerSize(0.25);
	  //	}
      } else { //MC samples
 cout << " yya" << endl;	
	if (strstr(title,"thrustc")) {
	  /*if (strstr(etapt,"typ_0_pt0_eta1_3")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
	    for (int kl=0; kl<nbinx; kl++) {
             sysperrMC[kl] = mcperr_3_0_1[kl];
            cout << " EROOOOOOOOOOOOOOOOOOOOOO" << mcperr_3_0_1[kl] << endl;
            }
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_0_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_1_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_2_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_3_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_4_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_5_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_6_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_3")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_3_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_3_7_1[kl];}
	  }*/
	} else if (strstr(title,"minorc")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_6_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_6_7_1[kl];}
	  }*/
	} else if (strstr(title,"rhotot")) {
	  /*if (strstr(etapt,"typ_0_pt0_eta1_18")) {
            cout <<"RHOTOT  ============= " << endl;
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_0_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_1_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_2_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_3_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_4_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_5_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_6_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_7_1[kl];}
	  }*/
	} else if (strstr(title,"y3c")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_15_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_15_7_1[kl];}
	  }*/
	} else if (strstr(title,"broadt")) {
	  /*if (strstr(etapt,"typ_0_pt0_eta1_18")) {
             cout << "Broadt ======= "
             cout << "Broadt ======= "
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_0_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_1_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_2_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_3_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_4_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_5_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_6_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_18")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_18_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_18_7_1[kl];}
	  }*/

	} else if (strstr(title,"broadw")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_21_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_21_7_1[kl];}
	  }*/

	} else if (strstr(title,"ttmass")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_7_1[kl];}
	  }*/

	} else if (strstr(title,"htmass")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_27_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_27_7_1[kl];}
	  }*/
	} else if (strstr(title,"sphericity")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_30_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_30_7_1[kl];}
	  }*/
	} else if (strstr(title,"cparameter")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_31_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_31_7_1[kl];}
	  }*/
	} else if (strstr(title,"rhottot")) {
	  /*if (strstr(etapt,"typ_0_pt0_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_0_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_1_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_2_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_3_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_4_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_5_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_6_1[kl];}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_24")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_24_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_24_7_1[kl];}
	  }*/
	} else if (strstr(title,"p3byp12c")) {
	  /*if (strstr(etapt,"c0_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_0_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_0_1[kl];}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_1_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_1_1[kl];}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_2_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_2_1[kl];}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_3_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_3_1[kl];}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_4_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_4_1[kl];}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_5_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_5_1[kl];}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_6_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_6_1[kl];}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int kl=0; kl<nbinx; kl++) {sysperrMC[kl] = mcperr_33_7_1[kl];}
	    for (int kl=0; kl<nbinx; kl++) {sysnerrMC[kl] = mcnerr_33_7_1[kl];}
	  }*/


	} else {
	  cout <<"MC : Variable name is not matching with the name "<<endl;
	  break;
	}
	///	cout << " before taking syserror errormc_matrix = " << endl;
	for (int ix=0; ix<errormc_matrix.GetNrows(); ix++){
	  for (int iy=0; iy<errormc_matrix.GetNcols(); iy++){
	    ///	    cout << errormc_matrix(ix,iy) << " " ;
	  }
	  ///	  cout << endl;
	}	
	//	if (ij < 2){
	//	  file_out << " ************** THE SYSTEMATIC ERROR OF MC " << ttl2[ij] << endl; 
	//	  for (int kl=0; kl<nbinx; kl++) {
	//	    file_out << " " << sysnerrMC[kl] << " " ;
	//	  }
	//	  file_out << endl;
	//	}
	///	cout << " sysnerrMC " ;

   // Pythia8 up down error 
if(mcsys){
 if(icol==2){
if(typindex==0){
   if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_0_pt0_eta1_3")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_0_1_3 >cuetp8m1varne_0_0_1_3) updnerrMC[kl] = cuetp8m1varpe_0_0_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_0_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_0_pt1_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_1_1_3 >cuetp8m1varne_0_1_1_3) updnerrMC[kl] = cuetp8m1varpe_0_1_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_1_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_0_pt2_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_2_1_3 >cuetp8m1varne_0_2_1_3) updnerrMC[kl] = cuetp8m1varpe_0_2_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_2_1_3[kl];
            }
          } else if (strstr(etapt,"typ_0_pt3_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_3_1_3 >cuetp8m1varne_0_3_1_3) updnerrMC[kl] = cuetp8m1varpe_0_3_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_3_1_3[kl];
               }
          } else if (strstr(etapt,"typ_0_pt4_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_4_1_3 >cuetp8m1varne_0_4_1_3) updnerrMC[kl] = cuetp8m1varpe_0_4_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_4_1_3[kl];
            }
          } else if (strstr(etapt,"typ_0_pt5_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_5_1_3 >cuetp8m1varne_0_5_1_3) updnerrMC[kl] = cuetp8m1varpe_0_5_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_5_1_3[kl];
             }
          } else if (strstr(etapt,"typ_0_pt6_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_6_1_3 >cuetp8m1varne_0_6_1_3) updnerrMC[kl] = cuetp8m1varpe_0_6_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_6_1_3[kl];
            }
          } else if (strstr(etapt,"typ_0_pt7_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_7_1_3 >cuetp8m1varne_0_7_1_3) updnerrMC[kl] = cuetp8m1varpe_0_7_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_7_1_3[kl];
            }
          }
        }
      
     if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_0_pt0_eta1_9")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_0_1_9 >cuetp8m1varne_0_0_1_9) updnerrMC[kl] = cuetp8m1varpe_0_0_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_0_1_9[kl];
            }
          } else if (strstr(etapt,"typ_0_pt1_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_1_1_9 >cuetp8m1varne_0_1_1_9) updnerrMC[kl] = cuetp8m1varpe_0_1_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_1_1_9[kl];
            }
          } else if (strstr(etapt,"typ_0_pt2_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_2_1_9 >cuetp8m1varne_0_2_1_9) updnerrMC[kl] = cuetp8m1varpe_0_2_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_2_1_9[kl];
            }
          } else if (strstr(etapt,"typ_0_pt3_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_3_1_9 >cuetp8m1varne_0_3_1_9) updnerrMC[kl] = cuetp8m1varpe_0_3_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_3_1_9[kl];
               }
          } else if (strstr(etapt,"typ_0_pt4_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_4_1_9 >cuetp8m1varne_0_4_1_9) updnerrMC[kl] = cuetp8m1varpe_0_4_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_4_1_9[kl];
            }
          } else if (strstr(etapt,"typ_0_pt5_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_5_1_9 >cuetp8m1varne_0_5_1_9) updnerrMC[kl] = cuetp8m1varpe_0_5_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_5_1_9[kl];
             }
          } else if (strstr(etapt,"typ_0_pt6_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_6_1_9 >cuetp8m1varne_0_6_1_9) updnerrMC[kl] = cuetp8m1varpe_0_6_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_6_1_9[kl];
            }
          } else if (strstr(etapt,"typ_0_pt7_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) { 
           if(cuetp8m1varpe_0_7_1_9 >cuetp8m1varne_0_7_1_9) updnerrMC[kl] = cuetp8m1varpe_0_7_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_7_1_9[kl];
            }
          }
        }
    
      if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_0_pt0_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_0_1_18 >cuetp8m1varne_0_0_1_18) updnerrMC[kl] = cuetp8m1varpe_0_0_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_0_1_18[kl];
            }
          } else if (strstr(etapt,"typ_0_pt1_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_1_1_18 >cuetp8m1varne_0_1_1_18) updnerrMC[kl] = cuetp8m1varpe_0_1_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_1_1_18[kl];
            }
          } else if (strstr(etapt,"typ_0_pt2_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_2_1_18 >cuetp8m1varne_0_2_1_18) updnerrMC[kl] = cuetp8m1varpe_0_2_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_2_1_18[kl];
            }
          } else if (strstr(etapt,"typ_0_pt3_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_3_1_18 >cuetp8m1varne_0_3_1_18) updnerrMC[kl] = cuetp8m1varpe_0_3_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_3_1_18[kl];
               }
          } else if (strstr(etapt,"typ_0_pt4_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_4_1_18 >cuetp8m1varne_0_4_1_18) updnerrMC[kl] = cuetp8m1varpe_0_4_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_4_1_18[kl];
            }
          } else if (strstr(etapt,"typ_0_pt5_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_5_1_18 >cuetp8m1varne_0_5_1_18) updnerrMC[kl] = cuetp8m1varpe_0_5_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_5_1_18[kl];
             }
          } else if (strstr(etapt,"typ_0_pt6_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_6_1_18 >cuetp8m1varne_0_6_1_18) updnerrMC[kl] = cuetp8m1varpe_0_6_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_6_1_18[kl];
            }
          } else if (strstr(etapt,"typ_0_pt7_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_7_1_18 >cuetp8m1varne_0_7_1_18) updnerrMC[kl] = cuetp8m1varpe_0_7_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_7_1_18[kl];
         }
       }
     }

   if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_0_pt0_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_0_1_24 >cuetp8m1varne_0_0_1_24) updnerrMC[kl] = cuetp8m1varpe_0_0_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_0_1_24[kl];
            }
          } else if (strstr(etapt,"typ_0_pt1_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_1_1_24 >cuetp8m1varne_0_1_1_24) updnerrMC[kl] = cuetp8m1varpe_0_1_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_1_1_24[kl];
            }
          } else if (strstr(etapt,"typ_0_pt2_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_0_2_1_24 >cuetp8m1varne_0_2_1_24) updnerrMC[kl] = cuetp8m1varpe_0_2_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_2_1_24[kl];
            }
          } else if (strstr(etapt,"typ_0_pt3_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_3_1_24 >cuetp8m1varne_0_3_1_24) updnerrMC[kl] = cuetp8m1varpe_0_3_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_3_1_24[kl];
               }
          } else if (strstr(etapt,"typ_0_pt4_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_0_4_1_24 >cuetp8m1varne_0_4_1_24) updnerrMC[kl] = cuetp8m1varpe_0_4_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_4_1_24[kl];
            }
          } else if (strstr(etapt,"typ_0_pt5_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_5_1_24 >cuetp8m1varne_0_5_1_24) updnerrMC[kl] = cuetp8m1varpe_0_5_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_5_1_24[kl];
             }
          } else if (strstr(etapt,"typ_0_pt6_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_0_6_1_24 >cuetp8m1varne_0_6_1_24) updnerrMC[kl] = cuetp8m1varpe_0_6_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_6_1_24[kl];
            }
          } else if (strstr(etapt,"typ_0_pt7_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_0_7_1_24 >cuetp8m1varne_0_7_1_24) updnerrMC[kl] = cuetp8m1varpe_0_7_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_0_7_1_24[kl];
         }
       }
     }

   } 

   if(typindex==1){
   if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_1_pt0_eta1_3")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_0_1_3 >cuetp8m1varne_1_0_1_3) updnerrMC[kl] = cuetp8m1varpe_1_0_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_0_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_1_pt1_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_1_1_3 >cuetp8m1varne_1_1_1_3) updnerrMC[kl] = cuetp8m1varpe_1_1_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_1_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_1_pt2_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_2_1_3 >cuetp8m1varne_1_2_1_3) updnerrMC[kl] = cuetp8m1varpe_1_2_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_2_1_3[kl];
            }
          } else if (strstr(etapt,"typ_1_pt3_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_3_1_3 >cuetp8m1varne_1_3_1_3) updnerrMC[kl] = cuetp8m1varpe_1_3_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_3_1_3[kl];
               }
          } else if (strstr(etapt,"typ_1_pt4_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_4_1_3 >cuetp8m1varne_1_4_1_3) updnerrMC[kl] = cuetp8m1varpe_1_4_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_4_1_3[kl];
            }
          } else if (strstr(etapt,"typ_1_pt5_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_5_1_3 >cuetp8m1varne_1_5_1_3) updnerrMC[kl] = cuetp8m1varpe_1_5_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_5_1_3[kl];
             }
          } else if (strstr(etapt,"typ_1_pt6_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_6_1_3 >cuetp8m1varne_1_6_1_3) updnerrMC[kl] = cuetp8m1varpe_1_6_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_6_1_3[kl];
            }
          } else if (strstr(etapt,"typ_1_pt7_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_7_1_3 >cuetp8m1varne_1_7_1_3) updnerrMC[kl] = cuetp8m1varpe_1_7_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_7_1_3[kl];
            }
          }
        }
      
     if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_9")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_0_1_9 >cuetp8m1varne_1_0_1_9) updnerrMC[kl] = cuetp8m1varpe_1_0_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_0_1_9[kl];
            }
          } else if (strstr(etapt,"typ_1_pt1_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_1_1_9 >cuetp8m1varne_1_1_1_9) updnerrMC[kl] = cuetp8m1varpe_1_1_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_1_1_9[kl];
            }
          } else if (strstr(etapt,"typ_1_pt2_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_2_1_9 >cuetp8m1varne_1_2_1_9) updnerrMC[kl] = cuetp8m1varpe_1_2_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_2_1_9[kl];
            }
          } else if (strstr(etapt,"typ_1_pt3_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_3_1_9 >cuetp8m1varne_1_3_1_9) updnerrMC[kl] = cuetp8m1varpe_1_3_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_3_1_9[kl];
               }
          } else if (strstr(etapt,"typ_1_pt4_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_4_1_9 >cuetp8m1varne_1_4_1_9) updnerrMC[kl] = cuetp8m1varpe_1_4_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_4_1_9[kl];
            }
          } else if (strstr(etapt,"typ_1_pt5_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_5_1_9 >cuetp8m1varne_1_5_1_9) updnerrMC[kl] = cuetp8m1varpe_1_5_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_5_1_9[kl];
             }
          } else if (strstr(etapt,"typ_1_pt6_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_6_1_9 >cuetp8m1varne_1_6_1_9) updnerrMC[kl] = cuetp8m1varpe_1_6_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_6_1_9[kl];
            }
          } else if (strstr(etapt,"typ_1_pt7_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) { 
           if(cuetp8m1varpe_1_7_1_9 >cuetp8m1varne_1_7_1_9) updnerrMC[kl] = cuetp8m1varpe_1_7_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_7_1_9[kl];
            }
          }
        }
    
      if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_1_pt0_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_0_1_18 >cuetp8m1varne_1_0_1_18) updnerrMC[kl] = cuetp8m1varpe_1_0_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_0_1_18[kl];
            }
          } else if (strstr(etapt,"typ_1_pt1_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_1_1_18 >cuetp8m1varne_1_1_1_18) updnerrMC[kl] = cuetp8m1varpe_1_1_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_1_1_18[kl];
            }
          } else if (strstr(etapt,"typ_1_pt2_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_2_1_18 >cuetp8m1varne_1_2_1_18) updnerrMC[kl] = cuetp8m1varpe_1_2_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_2_1_18[kl];
            }
          } else if (strstr(etapt,"typ_1_pt3_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_3_1_18 >cuetp8m1varne_1_3_1_18) updnerrMC[kl] = cuetp8m1varpe_1_3_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_3_1_18[kl];
               }
          } else if (strstr(etapt,"typ_1_pt4_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_4_1_18 >cuetp8m1varne_1_4_1_18) updnerrMC[kl] = cuetp8m1varpe_1_4_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_4_1_18[kl];
            }
          } else if (strstr(etapt,"typ_1_pt5_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_5_1_18 >cuetp8m1varne_1_5_1_18) updnerrMC[kl] = cuetp8m1varpe_1_5_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_5_1_18[kl];
             }
          } else if (strstr(etapt,"typ_1_pt6_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_6_1_18 >cuetp8m1varne_1_6_1_18) updnerrMC[kl] = cuetp8m1varpe_1_6_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_6_1_18[kl];
            }
          } else if (strstr(etapt,"typ_1_pt7_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_7_1_18 >cuetp8m1varne_1_7_1_18) updnerrMC[kl] = cuetp8m1varpe_1_7_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_7_1_18[kl];
         }
       }
     }

   if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_0_1_24 >cuetp8m1varne_1_0_1_24) updnerrMC[kl] = cuetp8m1varpe_1_0_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_0_1_24[kl];
            }
          } else if (strstr(etapt,"typ_1_pt1_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_1_1_24 >cuetp8m1varne_1_1_1_24) updnerrMC[kl] = cuetp8m1varpe_1_1_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_1_1_24[kl];
            }
          } else if (strstr(etapt,"typ_1_pt2_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_1_2_1_24 >cuetp8m1varne_1_2_1_24) updnerrMC[kl] = cuetp8m1varpe_1_2_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_2_1_24[kl];
            }
          } else if (strstr(etapt,"typ_1_pt3_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_3_1_24 >cuetp8m1varne_1_3_1_24) updnerrMC[kl] = cuetp8m1varpe_1_3_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_3_1_24[kl];
               }
          } else if (strstr(etapt,"typ_1_pt4_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_1_4_1_24 >cuetp8m1varne_1_4_1_24) updnerrMC[kl] = cuetp8m1varpe_1_4_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_4_1_24[kl];
            }
          } else if (strstr(etapt,"typ_1_pt5_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_5_1_24 >cuetp8m1varne_1_5_1_24) updnerrMC[kl] = cuetp8m1varpe_1_5_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_5_1_24[kl];
             }
          } else if (strstr(etapt,"typ_1_pt6_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_1_6_1_24 >cuetp8m1varne_1_6_1_24) updnerrMC[kl] = cuetp8m1varpe_1_6_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_6_1_24[kl];
            }
          } else if (strstr(etapt,"typ_1_pt7_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_1_7_1_24 >cuetp8m1varne_1_7_1_24) updnerrMC[kl] = cuetp8m1varpe_1_7_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_1_7_1_24[kl];
         }
       }
     }

   } 

if(typindex==2){
   if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_2_pt0_eta1_3")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_0_1_3 >cuetp8m1varne_2_0_1_3) updnerrMC[kl] = cuetp8m1varpe_2_0_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_0_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_2_pt1_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_1_1_3 >cuetp8m1varne_2_1_1_3) updnerrMC[kl] = cuetp8m1varpe_2_1_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_1_1_3[kl]; 
            }
          } else if (strstr(etapt,"typ_2_pt2_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_2_1_3 >cuetp8m1varne_2_2_1_3) updnerrMC[kl] = cuetp8m1varpe_2_2_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_2_1_3[kl];
            }
          } else if (strstr(etapt,"typ_2_pt3_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_3_1_3 >cuetp8m1varne_2_3_1_3) updnerrMC[kl] = cuetp8m1varpe_2_3_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_3_1_3[kl];
               }
          } else if (strstr(etapt,"typ_2_pt4_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_4_1_3 >cuetp8m1varne_2_4_1_3) updnerrMC[kl] = cuetp8m1varpe_2_4_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_4_1_3[kl];
            }
          } else if (strstr(etapt,"typ_2_pt5_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_5_1_3 >cuetp8m1varne_2_5_1_3) updnerrMC[kl] = cuetp8m1varpe_2_5_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_5_1_3[kl];
             }
          } else if (strstr(etapt,"typ_2_pt6_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_6_1_3 >cuetp8m1varne_2_6_1_3) updnerrMC[kl] = cuetp8m1varpe_2_6_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_6_1_3[kl];
            }
          } else if (strstr(etapt,"typ_2_pt7_eta1_3")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_7_1_3 >cuetp8m1varne_2_7_1_3) updnerrMC[kl] = cuetp8m1varpe_2_7_1_3[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_7_1_3[kl];
            }
          }
        }
      
     if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_9")) {
         cout << " EROOOOOOOOOOOOOOOOOOOOOO" << endl;
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_0_1_9 >cuetp8m1varne_2_0_1_9) updnerrMC[kl] = cuetp8m1varpe_2_0_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_0_1_9[kl];
            }
          } else if (strstr(etapt,"typ_2_pt1_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_1_1_9 >cuetp8m1varne_2_1_1_9) updnerrMC[kl] = cuetp8m1varpe_2_1_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_1_1_9[kl];
            }
          } else if (strstr(etapt,"typ_2_pt2_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_2_1_9 >cuetp8m1varne_2_2_1_9) updnerrMC[kl] = cuetp8m1varpe_2_2_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_2_1_9[kl];
            }
          } else if (strstr(etapt,"typ_2_pt3_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_3_1_9 >cuetp8m1varne_2_3_1_9) updnerrMC[kl] = cuetp8m1varpe_2_3_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_3_1_9[kl];
               }
          } else if (strstr(etapt,"typ_2_pt4_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_4_1_9 >cuetp8m1varne_2_4_1_9) updnerrMC[kl] = cuetp8m1varpe_2_4_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_4_1_9[kl];
            }
          } else if (strstr(etapt,"typ_2_pt5_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_5_1_9 >cuetp8m1varne_2_5_1_9) updnerrMC[kl] = cuetp8m1varpe_2_5_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_5_1_9[kl];
             }
          } else if (strstr(etapt,"typ_2_pt6_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_6_1_9 >cuetp8m1varne_2_6_1_9) updnerrMC[kl] = cuetp8m1varpe_2_6_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_6_1_9[kl];
            }
          } else if (strstr(etapt,"typ_2_pt7_eta1_9")) {
            for (int kl=0; kl<nbinx; kl++) { 
           if(cuetp8m1varpe_2_7_1_9 >cuetp8m1varne_2_7_1_9) updnerrMC[kl] = cuetp8m1varpe_2_7_1_9[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_7_1_9[kl];
            }
          }
        }
    
      if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_2_pt0_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_0_1_18 >cuetp8m1varne_2_0_1_18) updnerrMC[kl] = cuetp8m1varpe_2_0_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_0_1_18[kl];
            }
          } else if (strstr(etapt,"typ_2_pt1_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_1_1_18 >cuetp8m1varne_2_1_1_18) updnerrMC[kl] = cuetp8m1varpe_2_1_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_1_1_18[kl];
            }
          } else if (strstr(etapt,"typ_2_pt2_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_2_1_18 >cuetp8m1varne_2_2_1_18) updnerrMC[kl] = cuetp8m1varpe_2_2_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_2_1_18[kl];
            }
          } else if (strstr(etapt,"typ_2_pt3_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_3_1_18 >cuetp8m1varne_2_3_1_18) updnerrMC[kl] = cuetp8m1varpe_2_3_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_3_1_18[kl];
               }
          } else if (strstr(etapt,"typ_2_pt4_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_4_1_18 >cuetp8m1varne_2_4_1_18) updnerrMC[kl] = cuetp8m1varpe_2_4_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_4_1_18[kl];
            }
          } else if (strstr(etapt,"typ_2_pt5_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_5_1_18 >cuetp8m1varne_2_5_1_18) updnerrMC[kl] = cuetp8m1varpe_2_5_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_5_1_18[kl];
             }
          } else if (strstr(etapt,"typ_2_pt6_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_6_1_18 >cuetp8m1varne_2_6_1_18) updnerrMC[kl] = cuetp8m1varpe_2_6_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_6_1_18[kl];
            }
          } else if (strstr(etapt,"typ_2_pt7_eta1_18")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_7_1_18 >cuetp8m1varne_2_7_1_18) updnerrMC[kl] = cuetp8m1varpe_2_7_1_18[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_7_1_18[kl];
         }
       }
     }

   if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_0_1_24 >cuetp8m1varne_2_0_1_24) updnerrMC[kl] = cuetp8m1varpe_2_0_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_0_1_24[kl];
            }
          } else if (strstr(etapt,"typ_2_pt1_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_1_1_24 >cuetp8m1varne_2_1_1_24) updnerrMC[kl] = cuetp8m1varpe_2_1_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_1_1_24[kl];
            }
          } else if (strstr(etapt,"typ_2_pt2_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
             if(cuetp8m1varpe_2_2_1_24 >cuetp8m1varne_2_2_1_24) updnerrMC[kl] = cuetp8m1varpe_2_2_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_2_1_24[kl];
            }
          } else if (strstr(etapt,"typ_2_pt3_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_3_1_24 >cuetp8m1varne_2_3_1_24) updnerrMC[kl] = cuetp8m1varpe_2_3_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_3_1_24[kl];
               }
          } else if (strstr(etapt,"typ_2_pt4_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
              if(cuetp8m1varpe_2_4_1_24 >cuetp8m1varne_2_4_1_24) updnerrMC[kl] = cuetp8m1varpe_2_4_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_4_1_24[kl];
            }
          } else if (strstr(etapt,"typ_2_pt5_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_5_1_24 >cuetp8m1varne_2_5_1_24) updnerrMC[kl] = cuetp8m1varpe_2_5_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_5_1_24[kl];
             }
          } else if (strstr(etapt,"typ_2_pt6_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
            if(cuetp8m1varpe_2_6_1_24 >cuetp8m1varne_2_6_1_24) updnerrMC[kl] = cuetp8m1varpe_2_6_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_6_1_24[kl];
            }
          } else if (strstr(etapt,"typ_2_pt7_eta1_24")) {
            for (int kl=0; kl<nbinx; kl++) {
           if(cuetp8m1varpe_2_7_1_24 >cuetp8m1varne_2_7_1_24) updnerrMC[kl] = cuetp8m1varpe_2_7_1_24[kl];
             else updnerrMC[kl] = cuetp8m1varne_2_7_1_24[kl];
         }
       }
     }

   } 



   fdfy1[ij] = (TH1F*)fdfx[ij]->Clone(); 
   //fdfy11[ij] = (TH1F*)fdfx[ij]->Clone(); 
   for (int kl=0; kl<nbinx; kl++) {
          double xx = fdfy1[ij]->GetBinContent(kl+1);
          double yy = fdfy1[ij]->GetBinError(kl+1);
                   

 
          double errr = sqrt(yy*yy + pow(xx* updnerrMC[kl], 2.0));
          fdfy1[ij]->SetBinError(kl+1, errr); // fdfx[ij]->GetBinContent(kl+1)*syserrData[kl]);
          cout << "yes we solved it       ===============kl" << kl << errr << endl;
        } 
          fdfy1[ij]->SetFillColor(30); //4); //Green); //898); //5); 
          fdfy1[ij]->SetFillStyle(1111);
          fdfy1[ij]->SetLineWidth(1);
          fdfy1[ij]->SetMarkerSize(0.35);
        
}

}







   // Pythia8 up down error 



	for (int kl=0; kl<nbinx; kl++) {
	  double xx = fdfy[ij]->GetBinContent(kl+1);
	  double yy = fdfy[ij]->GetBinError(kl+1);
          //staterrData[kl][ij]=sqrt(yy*yy);
//	  double errr = sqrt(yy*yy + pow(xx* sysperrMC[kl], 2.0)); //GMA 29th MArch: changed positive to negative
	  //	  fdfy[ij]->SetBinError(kl+1, errr);
	  ///	  errormc_matrix(kl,kl) = sqrt(errormc_matrix(kl,kl)*errormc_matrix(kl,kl) + pow(xx*sysnerrMC[kl], 4.0));
	  //	  cout << "errormc_matrix " << errormc_matrix(kl,kl) << endl;
	  ///	  cout << sysnerrMC[kl] << " " ;
	}
	///	cout << endl;
	///	cout << " after taking syserror errormc_matrix = " << endl;
	for (int ix=0; ix<errormc_matrix.GetNrows(); ix++){
	  for (int iy=0; iy<errormc_matrix.GetNcols(); iy++){
	    ///	    cout << errormc_matrix(ix,iy) << " " ;
	  }
	  ///	  cout << endl;
	}
      }
    if(sysonmc && icol==2) fdfy1[ij]->Draw("E2 same");
      if (mRatio>=0) {
       //if(mcsys) fdfy1[1]->Draw("E2 same");
	if (icol==1) { 
	  gStyle->SetErrorX(0.5);
	  fdfy[ij]->Draw("E2 same");  //: (istatoth>0) ? "sames" : "same");
	  fdfx[ij]->Draw("same:E:x0");
	} else {
	  fdfx[ij]->Draw((istatoth>0) ? "sames:hist" : "same:hist");
	}
      }
    //if(mcsys && icol==2) fdfy1[ij]->Draw("E2 same");
    }
//covmatrix_3_1_2_21
    if (isCovariance_matrix){
      if (icol==1) {
	covmatrix.ResizeTo(nbinx,nbinx);
        if(typindex==0){
	if (strstr(title,"thrustc")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_3")) { 
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_0_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_1_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_2_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_3_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_4_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_5_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_6_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_3")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_0_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"minorc")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_6_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"rhotot")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_0_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_1_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_2_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_3_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_4_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_5_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_6_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_9")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_0_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"y3c")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_15_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"broadt")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_0_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_1_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_2_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_3_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_4_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_5_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_6_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_18")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_0_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"broadw")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_21_7_21[ix][iy];}}
	  }

	} else if (strstr(title,"ttmass")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_7_21[ix][iy];}}
	  }

	} else if (strstr(title,"htmass")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_27_7_21[ix][iy];}}
	  }
	} else if (strstr(title,"sphericity")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_30_7_21[ix][iy];}}
	  }

	} else if (strstr(title,"cparameter")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_31_7_21[ix][iy];}}
	  }

	} else if (strstr(title,"rhottot")) {
	  if (strstr(etapt,"typ_0_pt0_eta1_24")) {
            cout << "RHOTOT == " << endl;
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_0_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt1_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_1_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt2_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_2_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt3_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_3_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt4_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_4_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt5_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_5_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt6_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_6_21[ix][iy];}}
	  } else if (strstr(etapt,"typ_0_pt7_eta1_24")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_0_7_21[ix][iy];}}
	  }

	} else if (strstr(title,"p3byp12c")) {
	  if (strstr(etapt,"c0_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_0_21[ix][iy];}}
	  } else if (strstr(etapt,"c1_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_1_21[ix][iy];}}
	  } else if (strstr(etapt,"c2_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_2_21[ix][iy];}}
	  } else if (strstr(etapt,"c3_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_3_21[ix][iy];}}
	  } else if (strstr(etapt,"c4_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_4_21[ix][iy];}}
	  } else if (strstr(etapt,"c5_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_5_21[ix][iy];}}
	  } else if (strstr(etapt,"c6_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_6_21[ix][iy];}}
	  } else if (strstr(etapt,"c7_e")) {
	    for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_33_7_21[ix][iy];}}
	  }


	} else {
	  cout <<"Data : Variable name is not matching with the name "<<endl;
	  break;
	}
      }
      if(typindex==1){
       if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_1_pt0_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt1_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt2_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt3_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt4_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt5_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt6_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt7_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_1_7_21[ix][iy];}}
          }
         } else if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt1_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt2_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt3_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt4_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt5_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt6_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt7_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_1_7_21[ix][iy];}}
          }
         } else if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_1_pt0_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt1_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt2_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt3_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt4_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt5_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt6_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt7_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_1_7_21[ix][iy];}}
          }
          } else if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_1_pt0_eta1_24")) {
            cout << "RHOTOT == " << endl;
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt1_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt2_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt3_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt4_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt5_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt6_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_1_pt7_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_1_7_21[ix][iy];}}
          }
      } else {
          cout <<"Data : Variable name is not matching with the name "<<endl;
          break;
        }
      }
      if(typindex==2){
      if (strstr(title,"thrustc")) {
          if (strstr(etapt,"typ_2_pt0_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt1_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt2_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt3_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt4_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt5_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt6_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt7_eta1_3")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_3_2_7_21[ix][iy];}}
          }
          } else if (strstr(title,"rhotot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt1_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt2_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt3_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt4_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt5_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt6_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt7_eta1_9")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_9_2_7_21[ix][iy];}}
          }
          } else if (strstr(title,"broadt")) {
          if (strstr(etapt,"typ_2_pt0_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt1_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt2_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt3_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt4_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt5_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt6_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt7_eta1_18")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_18_2_7_21[ix][iy];}}
          }
         } else if (strstr(title,"rhottot")) {
          if (strstr(etapt,"typ_2_pt0_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_0_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt1_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_1_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt2_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_2_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt3_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_3_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt4_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_4_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt5_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_5_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt6_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_6_21[ix][iy];}}
          } else if (strstr(etapt,"typ_2_pt7_eta1_24")) {
            for (int ix=0; ix<nbinx; ix++) {for (int iy=0; iy<nbinx; iy++) {covmatrix(ix,iy) = covmatrix_24_2_7_21[ix][iy];}}
          }

      } else {
          cout <<"Data : Variable name is not matching with the name "<<endl;
          break;
        }
      } 


	///	cout <<  " before taking syserror covmatrix= " << endl;
	for (int ix=0; ix<nbinx; ix++) {
	  for (int iy=0; iy<nbinx; iy++) {
	    if (isUnitScale){
              covmatrix(ix,iy) = covmatrix(ix,iy)*scalex[0]*scalex[0];
             // covmatrix(ix,iy) = covmatrix(ix,iy);
            }
	    ///	    cout << covmatrix(ix,iy) << " " ;
	  }
	  ///	  cout << endl;
	}
         for (int ix=0; ix<nbinx; ix++) {
         if(isUnitScale) staterrData[ix]= binerror[ix][0]*scalex[0]; 
           else staterrData[ix]= binerror[ix][0];
        //syserrData[ix]=pow(syserrData[ix],2);  
         cout << "STAT ERROR =  " << staterrData[ix] << " ;SYS ERROR = " <<syserrData[ix] <<endl;
         }
     
         //if (ij>0 && mRatio>=0) fdfx[0]->Chi2TestX(fdfx[ij], chisqq, ndff, good, "UW");

	///	cout << "syserrData = " ;
         //double xx=0;
	for (int ix=0; ix<nbinx; ix++) {
	   double xx = fdfx[ij]->GetBinContent(ix+1);
                  cout << " STATTTTTTT " << fdfx[0]->GetBinError(ix+1) << endl;
	  //	  fileout << " covmatrix " << covmatrix(ix,ix) << " xx " << xx << " syserrData " << syserrData[ix] << endl;
	  covmatrix(ix,ix) = sqrt(covmatrix(ix,ix)*covmatrix(ix,ix) + pow(xx* syserrData[ix], 4.0) + (staterrData[ix]*staterrData[ix]));
	  //covmatrix(ix,ix) = sqrt((covmatrix(ix,ix)*covmatrix(ix,ix))+ (xx*xx*syserrData[ix])+ (staterrData[ix]*staterrData[ix]));
	  	  //fileout << "Cov Matrix= "<< covmatrix(ix,ix) <<" SysErr= " << syserrData[ix] << " StatErr= " << staterrData[ix]  << " xx= " << xx << endl;
	}
	///	cout << endl;
		cout <<  " after taking syserror covmatrix= " << endl;
        for (int ix=0; ix<nbinx; ix++) {
          for (int iy=0; iy<nbinx; iy++) {
	   //             fileout << " ; " <<covmatrix(ix,iy) << " " ;
          }
          cout << endl;
        }
      }//if (icol==1) 
      if (ij>0){
       double chisum=0;
	TVectorD diffdmc(nbinx);
	for (int kl=0; kl<fdfx[ij]->GetNbinsX(); kl++) {
	  diffdmc(kl) = (binentry[kl][0]-binentry[kl][ij]);
          //fileout << " OE =" << pow((binentry[kl][0]-binentry[kl][ij]),2) << " CHI2 = " << endl;
           chisum=chisum +(pow((binentry[kl][0]-binentry[kl][ij]),2))/binentry[kl][ij];
           //fileout <<" CHI2 = " << chisum <<endl; 
	  	  //fileout << "binentry[kl][0] " << binentry[kl][0] << " binentry[kl][ij] " << binentry[kl][ij] << " diffdmc " << diffdmc(kl) << endl;
	}
	//	TMatrixD comcovmatrix = covmatrix + errormc_matrix;
	TMatrixD comcovmatrix = covmatrix;                                                                                                 
	///	cout <<  " covmatrix= " << endl;
	for (int ix=0; ix<covmatrix.GetNrows(); ix++) {
	  for (int iy=0; iy<covmatrix.GetNcols(); iy++) {
	    	    //fileout << " covmatrix + errormc_matrix = " << covmatrix(ix,iy) << " " ;
	  }
	  	  cout << endl;
	}
		cout <<  " errormc_matrix= " << endl;
	for (int ix=0; ix<errormc_matrix.GetNrows(); ix++) {
	  for (int iy=0; iy<errormc_matrix.GetNcols(); iy++) {
	    ///	    cout << errormc_matrix(ix,iy) << " " ;
	  }
	  ///	  cout << endl;
	}
	///	cout <<  " com covmatrix= " << endl;
	for (int ix=0; ix<comcovmatrix.GetNrows(); ix++) {
	  for (int iy=0; iy<comcovmatrix.GetNcols(); iy++) {
	    	    cout << comcovmatrix(ix,iy) << " " ;
	  }
	  	  cout << endl;
	}
	comcovmatrix.Invert();
	TVectorD covdx = comcovmatrix*diffdmc;
	for (int ix=0; ix<comcovmatrix.GetNcols(); ix++) { 
	  chi2[ij] +=(diffdmc(ix)*covdx(ix));
	  	  //fileout<<"for ix = "<< ix << " Inverse of comcovmatrix " << comcovmatrix(ix,ix) << " covdx " << covdx(ix) << " diffmc " << 
                 //diffdmc(ix) <<" chi2 " << chi2[ij] << endl; 
	}
	//	chi2[ij] = chi2[ij]/(comcovmatrix.GetNcols()-1); 
	//	cout << chi2[ij] << endl;
	//	file_out << chi2[ij] << endl;
      }//(ij>0)
    } //(isCovariance_matrix)

    if (mRatio>=0) {
      if (icol==1) {
	if (idx<0) {
	  if (isttl2>0 && mRatio>=0) { 
	    //GMA 1302161	    if (mRatio>0 || isrootchi==0) {
	    tleg->AddEntry(fdfx[ij],ttl2[ij],"lpfe"); //pss
	    //	    } else {
	    //	      sprintf(fname, "%s   : #chi^{2}/ndf", ttl2[ij]);
	      //	      tleg->AddEntry(fdfx[ij],fname,"lpfe"); //pss
	      //	    }
	  }
	} else {
	  tleg->AddEntry(fdfx[ij],title,"lpfe");
	}
      } else {
	
	if (isrootchi>0) {
	  double xx=int(10*chisqq)/10.;
	  double prob = TMath::Prob(chisqq, ndff);
	  //GMA 130216	  if (mRatio!=0) {
	  tleg->AddEntry(fdfx[ij],ttl2[ij],"lf");
	} else {
	  if (isttl2>0) tleg->AddEntry(fdfx[ij],ttl2[ij],"lf"); //lpf
	}
	
	if (idx>=0 && isttl2>0) {
	  fileout<<ttl2[ij]<<"_"<<varnm<<"_ch["<<idx<<"]="<< chisqq <<"/;"<<endl;
	  //fileout<<ttl2[ij]<<"_"<<varnm<<"_nd["<<idx<<"]="<< ndff <<";s"<<endl;
	  fileout<<ttl2[ij]<<"_"<<varnm<<"_pr["<<idx<<"]="<< TMath::Prob(chisqq, ndff) <<";"<<endl;
	}
      }
      //    if (mRatio>=0) {
      if (istatoth>0) {
	gPad->Update();
	TPaveStats * stx = (TPaveStats*)fdfx[ij]->FindObject("stats");
	//      stx->SetTextColor(icol);
	//      stx->SetLineColor(icol);
	//      stx->SetFillColor(10);
      }
      //if (isttl2>0 && idx<0 && ipad==0) tleg->Draw();
 //        tleg->Draw();
	  //  tleg->AddEntry(fdfy[5],"Total Unc","lpfe"); //pss
      tleg->SetTextSize(0.050);
      TPaveText *ptstt = new TPaveText(.5,0.1,0.8,0.1,"brNDC");
      ptstt->SetFillColor(10);
    ptstt->SetBorderSize(0);
    ptstt->SetTextFont(42);
    //std::size_t found = etapt.find(typ0);
  //if (found!=std::string::npos) sprintf (fname, "%g", typ[0]);
   /* std::size_t found = etapt.find(typ1);
  if (found!=std::string::npos) sprintf (fname, "%g", typ[1]);
    std::size_t found = etapt.find(typ2);
  if (found!=std::string::npos) sprintf (fname, "%g", typ[2]);*/
    if(typnamerange && ij==ndata2-1){
     TString s(etapt);
     TString s2( s(6,7) );
    if(typindex==0) sprintf(typname,"%s",typ[0]);
    if(typindex==1) sprintf(typname,"%s",typ[1]);
    if(typindex==2) sprintf(typname,"%s",typ[2]);
/*    if(strstr(etapt,"pt0_eta1")) sprintf (fname, "%s: %s", typname,htrang[0]);
    if(strstr(etapt,"pt1_eta1")) sprintf (fname, "%s: %s", typname,htrang[1]);
    if(strstr(etapt,"pt2_eta1")) sprintf (fname, "%s: %s", typname,htrang[2]);
    if(strstr(etapt,"pt3_eta1")) sprintf (fname, "%s: %s", typname,htrang[3]);
    if(strstr(etapt,"pt4_eta1")) sprintf (fname, "%s: %s", typname,htrang[4]);
    if(strstr(etapt,"pt5_eta1")) sprintf (fname, "%s: %s", typname,htrang[5]);
    if(strstr(etapt,"pt6_eta1")) sprintf (fname, "%s: %s", typname,htrang[6]);
    if(strstr(etapt,"pt7_eta1")) sprintf (fname, "%s: %s", typname,htrang[7]);
*/
    if(strstr(etapt,"pt0_eta1")) sprintf (fname, "%s %s", htrang[0], "GeV/c");
    if(strstr(etapt,"pt1_eta1")) sprintf (fname, "%s %s", htrang[1], "GeV/c");
    if(strstr(etapt,"pt2_eta1")) sprintf (fname, "%s %s", htrang[2], "GeV/c");
    if(strstr(etapt,"pt3_eta1")) sprintf (fname, "%s %s", htrang[3], "GeV/c");
    if(strstr(etapt,"pt4_eta1")) sprintf (fname, "%s %s", htrang[4], "GeV/c");
    if(strstr(etapt,"pt5_eta1")) sprintf (fname, "%s %s", htrang[5], "GeV/c");
    if(strstr(etapt,"pt6_eta1")) sprintf (fname, "%s %s", htrang[6], "GeV/c");
    if(strstr(etapt,"pt7_eta1")) sprintf (fname, "%s %s", htrang[7], "GeV/c");
      TText* text = ptstt->AddText(fname);
      cout << "Jet Name catch=  " << fname << endl;
      text->SetTextSize(0.080);
      //if(ipad==0) ptstt->Draw();
    }
      if(ipad==0) ptstt->Draw();
     if(isExternalError==2){
     h1->SetFillColor(5);
     h2->SetFillColor(21);
     h3->SetFillColor(30);
     if(ij==4){  
      //  fdfy1[ij]->SetFillColor(30);
        tleg->AddEntry(h1,"Total Uncertainty","f");
       if(isstat) tleg->AddEntry(h2,"Statistical Uncertainty","f");
        if(sysonmc) tleg->AddEntry(h3,"Py8+CUETP8M1 Uncertainty","f");
         }
      }
         tleg->Draw();
      //tleg->SetTextSize(0.048);
       if(ij==ndata2-1){
       int iPeriod=4;
       int iPos =0; 
      CMS_lumi(c1,iPeriod,iPos);
      }
      //      tleg->SetTextColor(icol);
    }
  }    
  
  //  if (isitForPaper==2 && mRatio>=0) {
  if (isitForPaper==2 && mRatio>=0 && ipad==2) {
    //    if (mRatio==0) {cmsPrel(0.25, 0.58, .94, .2);} else { cmsPrel(0.38, 0.85, 0.065, 0.13);}
  }
  
  if (mRatio>=0 && ipad>=1 && isIndex>0) {
    //    if (ipad>2){
    //      TPaveText *ptst1 = new TPaveText(0.28,0.88,0.28,0.94,"brNDC"); //.18
    //    }else{
      TPaveText *ptst1 = new TPaveText(0.88,0.88,0.98,0.94,"brNDC");
      //    }
    ptst1->SetFillColor(10);
    ptst1->SetBorderSize(0);
    TText* text1 = ptst1->AddText(subttl[ipad-1]);
    text1->SetTextSize((isitForPaper>0)?0.07:0.08); //GMA 0.09
    text1->SetTextFont(42);
    //ptst1->Draw();
    TPaveText *ptst3 = new TPaveText(0.45,0.92,0.70,0.96,"brNDC");
    ptst3->SetFillColor(10);
    ptst3->SetBorderSize(0);
    //    TText* text3 = ptst3->AddText("320 < p_{T,1} < 390 GeV/c");
    //    TText* text3 = ptst3->AddText("p_{T,1} > 250 GeV/c"); 
    TText* text3 = ptst3->AddText(ptrange[ipad+3-1]); //for y23 plot
    text3->SetTextSize((isitForPaper>0)?0.06:0.08); //GMA 0.09 
    text3->SetTextFont(42);
    //    ptst3->Draw();
  } 
  
  int ipad2  = (nzone4>0) ? (TMath::Max(ipad,1)+4) : (TMath::Max(ipad,1)+1);
  if (isArbitrary) ipad2 = ipad+4;
  
  if (mRatio!=0) {
    TPaveText *ptst;
    icol = 0;
    for (int ij=0; ij<ndata2; ij++) {
      icol++;
      if (icol%3==0 || icol%5==0 || icol%7==0) icol++;
      if (icol>10) icol +=9;
      double xmx=-100.;
      double xmn=100.;
      double errsq = 0;
      double chisq =0;
      double avedev = 0;
      TVectorD errormc(nbinx);
      int nchn = 0;
      for (int kl=0; kl<nbinx; kl++) {
	if (isSepRatio==0) {
	  bincenter[kl] = fdfx[ij]->GetBinCenter(kl+1); //  + 0.2*fdfx[ij]->GetBinWidth(kl+1)*(ij-2)/(ndata2-1);  
	} else {
	  if (icol==1) bincenter[kl] = fdfx[ij]->GetBinCenter(kl+1);
	}
	
	binwidth[kl] = fdfx[ij]->GetBinWidth(kl+1);
	binrange[kl] = bincenter[kl] - 0.5*binwidth[kl];
	if (kl==nbinx-1) { binrange[kl+1] = bincenter[kl] + 0.5*binwidth[kl];}
	
	if (isExternalError==2) { 
	  binerrorNor[kl][ij] = fdfy[ij]->GetBinError(kl+1)/TMath::Max(minbincontent[ij], binentry[kl][ij]);
	}

	if (ij>0) {
	  double tmpp =0., tmpm=0.;
	  if (binentry[kl][ij] + binentry[kl][0]>0) {
	    
	    efficiency[kl] = binentry[kl][ij] / TMath::Max(minbincontent[0],binentry[kl][0]); 
	    
	    efficiencyforerr[kl] = TMath::Max(minbincontent[ij],binentry[kl][ij])/TMath::Max(minbincontent[0],binentry[kl][0]); 
	    erroriny[kl] = efficiency[kl]*sqrt(pow (binerror[kl][0]/(TMath::Max(minbincontent[0], binentry[kl][0])),2.0) +
					       pow (binerror[kl][ij]/(TMath::Max(minbincontent[ij], binentry[kl][ij])),2.0));
	    if (isLargerDeno) {
	      erroriny2[kl] = (1./(TMath::Max(1.e-6, efficiency[kl])))*
		sqrt(pow (binerror[kl][0]/TMath::Max(minbincontent[ij], binentry[kl][0]),2.0) +
		     pow (binerror[kl][ij]/TMath::Max(minbincontent[ij], binentry[kl][ij]),2.0));
	    }
	    tmpp = sqrt ( pow(TMath::Max(0.,efficiency[kl]-1.), 2.) + erroriny[kl]*erroriny[kl]);
	    tmpm = sqrt ( pow(TMath::Min(0.,efficiency[kl]-1.), 2.) + erroriny[kl]*erroriny[kl]); 
	    
	    if (binerror[kl][0]/(TMath::Max(minbincontent[0], binentry[kl][0]))> binerror[kl][ij]/(TMath::Max(minbincontent[ij], binentry[kl][ij]))) {
	      errorinyDoub[kl] = efficiency[kl]*binerror[kl][0]/(TMath::Max(minbincontent[0], binentry[kl][0]));
	    } else {
	      errorinyDoub[kl] = efficiency[kl]*binerror[kl][ij]/(TMath::Max(minbincontent[ij], binentry[kl][ij]));
	    }
	  } else { 
	    efficiency[kl] = 0; erroriny[kl]= errorinyExt[kl] = errorinyDoub[kl]=1.; //0.001; //1.e12;
	  }
	  
	  if (errpl[kl] < tmpp) errpl[kl] = tmpp;
	  if (errmi[kl] < tmpm) errmi[kl] = tmpm;
	  if (erroriny[kl]>0 && bincenter[kl]>=alow && bincenter[kl]<=ahg) {
	    
	    nchn++;
	    
	    if (isLargerDeno && efficiency[kl] >1.) {
	      chisq +=pow((1./efficiency[kl]-1.)/erroriny2[kl], 2.);
	      avedev += fabs(1./efficiency[kl]-1.)/pow(erroriny2[kl], 2.);
	      
	      if (fabs(1./efficiency[kl]-1.) > devEfficiency[kl]) {
		devEfficiency[kl] = TMath::Max(devEfficiency[kl], fabs(1./efficiency[kl]-1.));
		efficiencyErr[kl] = TMath::Max(efficiencyErr[kl], erroriny2[kl]);
               
	      }
	    } else {
	      chisq +=pow((efficiency[kl]-1.)/erroriny[kl], 2.);
	      avedev += fabs(efficiency[kl]-1.)/pow(erroriny[kl], 2.);
	      
	      if (fabs(efficiency[kl]-1.) > devEfficiency[kl]) {
		devEfficiency[kl] = TMath::Max(devEfficiency[kl], fabs(efficiency[kl]-1.));
		efficiencyErr[kl] = TMath::Max(efficiencyErr[kl],erroriny[kl]);
	      }
	    }
	    
	    if (isLargerDeno && efficiency[kl] >1.) { 
	      errsq +=pow(1./erroriny2[kl], 2.);
	    } else {
	      errsq +=pow(1./erroriny[kl], 2.);
	    }
	  }
	  
	  if (efficiency[kl] >0) {
	    if (xmx <efficiency[kl]+erroriny[kl]) xmx =efficiency[kl]+erroriny[kl];
	    if (xmn >efficiency[kl]-erroriny[kl]) xmn =efficiency[kl]-erroriny[kl]; 
	  }
	  
	  errorinx[kl] = 0.0001;
	  
	} // if (ij>0)
      } // for (int kl=0; kl<fdfx[ij]->GetNbinsX(); kl++)      

      if (isDoublerat && ij>0) {
	cout<<"  double "<<histname2<<"_"<<ttl2[ij]<<"_x"<<ifile<<"["<<nbinx<<"] ={";
        for (int kl=0; kl<nbinx-1; kl++) {
	  cout<<bincenter[kl]<<",";
	}
	cout<<bincenter[nbinx-1]<<"};"<<endl;

	cout<<"  double "<<histname2<<"_"<<ttl2[ij]<<"_y"<<ifile<<"["<<nbinx<<"] ={";
        for (int kl=0; kl<nbinx-1; kl++) {
	  cout<<TMath::Max(1.e-10, efficiency[kl])<<",";
	}
	cout<<TMath::Max(1.e-10,efficiency[nbinx-1])<<"};"<<endl;	

	cout<<"  double "<<histname2<<"_"<<ttl2[ij]<<"_dy"<<ifile<<"["<<nbinx<<"] ={";
        for (int kl=0; kl<nbinx-1; kl++) {
	  cout<<errorinyDoub[kl]<<",";
	}
	cout<<errorinyDoub[nbinx-1]<<"};"<<endl;
      }

      sprintf(hname, "xx_%s_%i", ttl2[ij], ij);
      sprintf(htitle, "xx_%s_%i",ttl2[ij], ij);
      ratHist[ij] = new TH1F(hname, htitle, nbinx,  binrange);
      sprintf(hname, "xxy_%s_%i", ttl2[ij], ij);
      sprintf(htitle, "xxy_%s_%i",ttl2[ij], ij);
      ratHist1[ij] = new TH1F(hname, htitle, nbinx,  binrange);
      sprintf(hname, "xxyy_%s_%i", ttl2[ij], ij);
      sprintf(htitle, "xxyy_%s_%i",ttl2[ij], ij);
      //ratHist2[ij] = (TH1F*)ratHist[ij]->Clone();
      //ratHist2[ij] = new TH1F(hname, htitle, nbinx,  binrange);

      //      if (ij>0) { 
      int ix = (mRatio>=0) ? ij : ij-1;
      //	if (isSepRatio || ij==1) {
      cout <<"xx "<< isSepRatio<<" "<< ij<<endl;
      if ((isSepRatio && ij>0) || (!isSepRatio && ij==0)) {
	cout <<"xx1 "<< isSepRatio<<" "<< ij<<endl;
	if (!isSepRatio) ix++;
	dd1[ix] = (TPad*)c1_3->GetPad(ix+1);
	dd1[ix]->SetRightMargin(0.005);
	if (isUnitScale){
	  dd1[ix]->SetLeftMargin(0.23);
	}else{
	  dd1[ix]->SetLeftMargin(0.182);
	}
	dd1[ix]->SetLogy(0);
	if (isSepRatio) {
	  double gapxx = 0.95/(ndata2+2);
	  double yy2 = (1. - gapxx*(ij+2))-0.005;
	  if (mRatio<0) {
	    double gapxx = 0.97/(ndata2-1);
	    double yy2 = (1. - gapxx*(ij-1));
	  }
	  double yy1= yy2- gapxx;
	  if (ij==ndata2-1) yy1 -=.1266 - 0.00666*ndata2; //0.1;
	  
	  dd1[ix]->SetPad(0.01,yy1,0.9985,yy2);
	} else if (ij==0){
	  dd1[ix]->SetPad(0.0001,0.0,0.999,0.4);
	}
	dd1[ix]->cd();
      }
      
      double width = 0.5*fdfx[ij]->GetBinWidth(10);
           //	sprintf(hname, "xx_%s_%i", ttl2[ij], ij);
      //	sprintf(htitle, "xx_%s_%i",ttl2[ij], ij);
      
      double shift = 0;
      if (isSepRatio==0) {
	double shift = (int(ij/2.))*0.2*width; // (ij + 1 - ndata2/2.)*0.2*width;
	if (ij%2==0) shift *=-1.;
      }
      
      //	ratHist[ij] = new TH1F(hname, htitle, nbinx,  binrange);
      
      ratHist[ij]->GetXaxis()->SetLabelFont(42);
      ratHist[ij]->GetXaxis()->SetTitleFont(42);
      ratHist[ij]->GetYaxis()->SetLabelFont(40);
      ratHist[ij]->GetYaxis()->SetTitleFont(42);
      ratHist[ij]->SetMaximum(1.4996); // 1.60); //(1.70); //1.09); //1.70); //1.07); // 1.85); //2.45); // 1.85); // 1.25); // 1.85); // 2.95); // 1.8\5);
      //      ratHist[ij]->SetMinimum(0.55); //(0.65); // 0.55); //0.91); //0.35); // 0.55); //0.85); // 0.55); //0.25); // 0.55);
      ratHist[ij]->SetMinimum(0.7); 

      ratHist[ij]->GetXaxis()->SetRangeUser(alow, ahg);
      
      if (isSepRatio==0) {
	ratHist[ij]->GetXaxis()->SetTickLength(.075);
	ratHist[ij]->GetYaxis()->SetTickLength(.05);
	
	ratHist[ij]->GetXaxis()->SetRangeUser(alow, ahg);
	ratHist[ij]->GetXaxis()->SetTitleSize(0.14); //(.08); //pss                                                                                   
	ratHist[ij]->GetXaxis()->SetTitleOffset(1.0);
	ratHist[ij]->GetXaxis()->CenterTitle();
	ratHist[ij]->GetXaxis()->SetLabelSize(0.12);
	ratHist[ij]->GetXaxis()->SetLabelOffset(-.000171);
	
	if (isExternalError==2) {
	  sprintf(fname3, "%s / xx", ttl2[0]);
	} else {
	  sprintf(fname3, "%s / %s", ttl2[1], ttl2[0]);
	}
sprintf(fname3, "%s", vartitle[titleindex]);
        ratHist[ij]->GetXaxis()->SetTitle(fname3);	
	ratHist[ij]->GetYaxis()->SetTitle("MC / Data");
	ratHist[ij]->GetYaxis()->SetTitleOffset(.70);
	ratHist[ij]->GetYaxis()->CenterTitle();
	ratHist[ij]->GetYaxis()->SetTitleSize((idx<0) ? 0.12 : 0.12);
	ratHist[ij]->GetYaxis()->SetLabelSize(0.11);
       
	if (ij==0 && isExternalError==2) {
        
           file_stat_err<<"double stat_"<<etapt<<"["<<nbinx<<"] ={";
	  for (int kl=0; kl<nbinx; kl++) {
	    //double error =sqrt(pow(sysnerrMC[kl],2.) + pow(binerrorNor[kl][0],2.));
	    //double error =sqrt(pow(syserrData[kl],2.) + pow(binerrorNor[kl][0],2.));
	    double error =sqrt(pow(syserrData[kl],2.) + pow((binerror[kl][0]/binentry[kl][0]),2.));
	    //double error1 =sqrt(pow(binerrorNor[kl][0],2.));
	    double error1 =sqrt(pow((binerror[kl][0]/binentry[kl][0]),2.));
            double error2 =sqrt(pow(updnerrMC[kl],2.));// + pow(binerrorNor[kl][1],2.)); 
	 //   double error1 =sqrt(pow(staterrData[kl],2.));
            if(isUnitScale) cout<< "Got or = "<< binerror[kl][0]*scalex[0] << endl;;
            cout << "ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR kl = " << kl << "BinError ==  " << error1 <<  "; A = " <<endl;
	    ratHist[ij]->SetBinContent(kl+1, 1.);
	    ratHist[ij]->SetBinError(kl+1, error);
            ratHist1[ij]->SetBinContent(kl+1, 1.);
            ratHist1[ij]->SetBinError(kl+1, error1);
           // ratHist2[ij]->SetBinContent(kl+1, 1.);
          //  ratHist2[ij]->SetBinError(kl+1, error2);

         //if(kl<nbinx-1) file_stat_err <<binerrorNor[kl][0]<<", ";
         //else file_stat_err <<binerrorNor[kl][0]<<"};"<< endl;
         if(kl<nbinx-1) file_stat_err <<binerror[kl][0]/binentry[kl][0]<<", ";
         else file_stat_err <<binerror[kl][0]/binentry[kl][0]<<"};"<< endl;

	  }
       
          //if (ij==1 && isExternalError==2){
          //double error2 =sqrt(pow(updnerrMC[kl],2.) + pow(binerrorNor[kl][1],2.)); 
	    //double error1 =sqrt(pow(binerrorNor[kl][0],2.));
            //ratHist2[ij]->SetBinContent(kl+1, 1.);
            //ratHist2[ij]->SetBinError(kl+1, error1);

          // }
         
	  gStyle->SetErrorX(0.5);
	  ratHist[0]->SetFillColor(5);
	  ratHist[0]->SetMarkerSize(0.1);
	  ratHist[0]->SetFillStyle(1111);

          ratHist1[0]->SetFillColor(21);
          ratHist1[0]->SetMarkerSize(0.1);
          ratHist1[0]->SetFillStyle(1111);
          //if(ij==1){ 
          //ratHist2[0]->SetFillColor(30);
          //ratHist2[0]->SetMarkerSize(0.1);
          //ratHist2[0]->SetFillStyle(1111);
         // }

	  for (int kl=0; kl<nbinx; kl++) { cout <<"kl "<< kl<<" "<<ratHist[0]->GetBinContent(kl+1)<<" "<<ratHist[0]->GetBinError(kl+1)<<endl;
	  }
	  ratHist[0]->Draw("E2");
	 if(isstat) ratHist1[0]->Draw("E2 same");
	 //if(mcsys) ratHist2[0]->Draw("E2 same");
	} // if (ij==0 && isExternalError==2)

         





      } //if (isSepRatio==0)
     
    




 
      if (ij>0) {
	for (int kl=0; kl< nbinx; kl++) {
	  ratHist[ij]->Fill(bincenter[kl], efficiency[kl]);
	}
	
	if ( isExternalError==2) {
	  for (int kl=0; kl<nbinx; kl++) {errorinyExt[kl] = efficiency[kl]*sqrt(pow(binerrorNor[kl][ij],2.) + pow(binerrorNor[kl][0],2.));
	  }
	    
	  for (int kl=0; kl< nbinx; kl++) { ratHist[ij]->SetBinError(kl+1, errorinyExt[kl]);}
	  sprintf(hname, "xx_%s_org_%i", ttl2[ij], ij);
	  sprintf(htitle, "xx_%s_org_%i",ttl2[ij], ij);
	  
	  ratHistOrg[ij] = new TH1F(hname, htitle, nbinx,  binrange); //bincenter[0]-width+shift,  bincenter[nbinx-1]+width+shift);
	  
	  for (int kl=0; kl< nbinx; kl++) {
	    ratHistOrg[ij]->Fill(bincenter[kl], efficiency[kl]);
	    ratHistOrg[ij]->SetBinError(kl+1, erroriny[kl]);
	  }
	  
	  ratHistOrg[ij]->SetLineWidth(1);
	  ratHistOrg[ij]->SetMarkerSize(0.65);
	  ratHistOrg[ij]->SetMarkerColor(1);
	  ratHistOrg[ij]->SetMarkerStyle(20);
	  ratHistOrg[ij]->SetMarkerSize(0.67); //1.0);
	  if (isSepRatio>0) {
	    ratHist[ij]->SetFillColor(icol); 
	    ratHist[ij]->SetFillStyle(1111);
	  } 
	  ratHist[ij]->SetLineWidth(1);
	  
	} else if (isDoublerat) {
	  for (int kl=0; kl< nbinx; kl++) {ratHist[ij]->SetBinError(kl+1, errorinyDoub[kl]);}
	} else {
	  for (int kl=0; kl< nbinx; kl++) { ratHist[ij]->SetBinError(kl+1, erroriny[kl]);}
	}
	
	if (isExternalError<=2) {
	  
	  ratHist[ij]->SetLineColor(1); 
	  ratHist[ij]->SetLineWidth(1);
	  ratHist[ij]->SetMarkerSize(0.65);
	  ratHist[ij]->SetMarkerColor(icol);
	  ratHist[ij]->SetMarkerStyle(20);
	  ratHist[ij]->SetMarkerSize(0.7); //1.0);
	}
	
	///	ratHist[ij]->SetMaximum(2.19); // 1.60); //(1.70); //1.09); //1.70); //1.07); // 1.85); //2.45); // 1.85); // 1.25); // 1.85); // 2.95); // 1.85);
	///	ratHist[ij]->SetMinimum(0.55); //(0.65); // 0.55); //0.91); //0.35); // 0.55); //0.85); // 0.55); //0.25); // 0.55);
	
	///	ratHist[ij]->GetXaxis()->SetRangeUser(alow, ahg);
	ratHist[ij]->GetXaxis()->SetTitle(pch1);
	ratHist[ij]->GetXaxis()->SetNdivisions(608); //608); //404);
	ratHist[ij]->GetYaxis()->SetNdivisions(404); //812); //505); //404);
	double xx = (isSepRatio>0) ? 0.01+0.065*(ndata2-1) : 0.1;
	//	double xx = (isSepRatio>0) ? 0.01+0.065*(ndata2-1) : .109;
	ratHist[ij]->GetYaxis()->SetLabelSize(xx);
	//	ratHist[ij]->Fit("pol1");
	if (isSepRatio==1) {
	  if (isExternalError==2) {
	      sprintf(fname3, "#frac{%s}{%s}", ttl2[ij], ttl2[0]);
	  } else {
	    sprintf(fname3, "#frac{%s}{%s}", ttl2[ij], ttl2[0]);
	  }
	  
	 // ratHist[ij]->GetYaxis()->SetTitle(fname3);
	 // ratHist[ij]->GetYaxis()->CenterTitle();
	  
	  float sclsize = 0.12 + 0.003*(ndata2-2); // 4; // 0.12;
	  float scloff = 0.45; // 0.45;
	    
	  sclsize *=0.95+0.1*(ndata2-3); scloff /=0.85+0.1*(ndata2-3);
	  
	  if (ij==ndata2-1) {
	    sclsize /=1.4+0.0333*ndata2; 
	    scloff *=1.65;
	  }
	  
	  if (isitForPaper>0 && mRatio!=0) {sclsize *=0.67; scloff *=1.65;} // sclsize *=0.8; scloff *=1.4;
	  // 3,12 sclsize *=0.67; scloff *=1.65;
	  //note sclsize *=0.73; scloff *=1.55;
	  // 4,12 
	  
	  //	  if (isSlide&& mRatio!=0) { sclsize *=1.07; scloff *=0.9;}
	  if (isNote) { sclsize *=1.25; scloff *=0.75;}
	  //	  if (isNote) { sclsize *=1.15; scloff *=0.85;}
	  
	  //	  ratHist[ij]->GetYaxis()->SetTitleSize(sclsize);
	  //	  ratHist[ij]->GetYaxis()->SetLabelSize(sclsize+0.001*ndata2);
	  scloff + = 0.04; //0.15; //0.10; //0.04;
	  
	  //	  ratHist[ij]->GetYaxis()->SetTitleOffset(xx);
	  
	  ratHist[ij]->GetXaxis()->SetTickLength( 0.02*(ndata2-1));
	  
	  if (ij==ndata2-1) {
	    double ttscl = (isitForPaper && mRatio!=0) ? 0.68 : 0.9; //0.85 for (3,12) 0.8; //1.1; //0.9; //1.0;
	    if (isSlide && mRatio!=0) ttscl=1;
	    if (isNote) ttscl=1.3;
	    ratHist[ij]->GetXaxis()->SetTitle(pch1);
          //ratHist[ij]->GetXaxis()->SetTitle("H_{t2}");
            //ratHist[ij]->GetXaxis()->SetTitle("#eta of leading jet");
            //ratHist[ij]->GetXaxis()->SetTitle("#phi of jet");
          // ratHist[ij]->GetXaxis()->SetTitle("#eta of second leading jet");
	   // ratHist[ij]->GetXaxis()->SetTitle("#Delta Pt of two leading jets (GeV/c)");
              //char ytitle[200]= {"No. of jet"}; 
          // ratHist[ij]->GetXaxis()->SetTitle("log_{e} (#tau_{#perp})");
           //ratHist[ij]->GetXaxis()->SetTitle("#rho_{Tot,C}"); 	
          //ratHist[ij]->GetXaxis()->SetTitle("#rho^{T}_{Tot,C}");  
        //ratHist[ij]->GetXaxis()->SetTitle("Y_{23,C}");    
         //ratHist[ij]->GetXaxis()->SetTitle("B_{ T,C}");
          //ratHist[ij]->GetXaxis()->SetTitle("#Delta#phi of Jets"); 
   //ratHist[ij]->GetXaxis()->SetTitle("Pt2 x sin( #Delta #phi )/Pt1");
          // ratHist[ij]->GetXaxis()->SetTitle("Pt of second leading jet (GeV/c)");
          // ratHist[ij]->GetXaxis()->SetTitle("Pt of leading jet (GeV/c)");
           if(basic_title)sprintf(fname3, "%s", basic_vartitle[titleindex]);
           else sprintf(fname3, "%s", vartitle[titleindex]);
          ratHist[ij]->GetXaxis()->SetTitle(fname3);
           ratHist[ij]->GetXaxis()->SetLabelOffset(0.02); // 1.01); // .02); //.01);
	    
	    ratHist[ij]->GetXaxis()->SetLabelSize(ttscl*((idx>=0) ? 0.03+0.02*ndata2 :0.08+0.02*ndata2));
	    
	    if (isSlide&& mRatio!=0) {
	      ratHist[ij]->GetXaxis()->SetTitleOffset((1.475 - 0.076*ndata2)/ttscl); //0.055*ndata2 //1.275
	    } else {
	      ratHist[ij]->GetXaxis()->SetTitleOffset((1.275 - 0.02*ndata2)/ttscl); //0.030*ndata2 0.055*ndata2 //1.275
	    }
	    double ttsiz=(idx>=0) ? 0.055+0.003*ndata2 : 0.060+0.0125*ndata2;
	    if (isitForPaper) ttsiz *=.8;
	    ratHist[ij]->GetXaxis()->SetTitleSize(ttscl*((idx>=0) ? 0.055+0.004*ndata2 : 0.065+0.0125*ndata2));
	    ratHist[ij]->GetXaxis()->CenterTitle();
	    scloff +=0.04; //Change Y-axis title/labale size for final ratio plot, which has different pad size                                       
	    sclsize *=1.05;
	  }
	  ratHist[ij]->GetYaxis()->SetTitleSize(2.20*sclsize);
	  ratHist[ij]->GetYaxis()->SetLabelSize(8.1*sclsize); //+0.0001*ndata2);
	  ratHist[ij]->GetYaxis()->SetTitleOffset(0.35*scloff);
          
 
	  if (isExternalError==2) {
           setTDRStyle(gStyle); 
	    gStyle->SetErrorX(0.5);
	    ratHist[ij]->Draw("E2");
            //ratHist[ij]->GetXaxis()->SetRangeUser(-4.5., 0);
	    ratHistOrg[ij]->Draw("same:E1:x0");
	  } else {
            setTDRStyle(gStyle);
	    //ratHist[ij]->GetXaxis()->SetRangeUser(-4.5, 0);
	    ratHist[ij]->Draw("E1:x0");
	  }
	} else {
	  
           setTDRStyle(gStyle); 
	   ratHist[ij]->SetLineColor(icol);
//	   if(icol==8) ratHist[ij]->SetLineColor(icol+20);
//	   if(ij==1) ratHist[ij]->SetLineStyle(1);
//	    else ratHist[ij]->SetLineStyle(icol);
           if(ij==1 && sysonmc) {
		ratHist2[ij] = (TH1F*)ratHist[ij]->Clone();
                ratHist2[1]->SetFillColor(30);
                ratHist2[1]->SetMarkerSize(0.1);
                ratHist2[1]->SetFillStyle(1111);        
             for (int kl=0; kl<nbinx; kl++) {
               double error2 =sqrt(pow(updnerrMC[kl],2.));
               ratHist2[1]->SetBinError(kl+1, error2);   
               ratHist2[1]->Draw("E2 same");
               }
               
             }
	  ratHist[ij]->SetLineStyle(0);
	  ratHist[ij]->SetLineWidth(1.0);
           sprintf(fname3, "%s", vartitle[titleindex]);
          ratHist[ij]->GetXaxis()->SetTitle(fname3);
           //ratHist[ij]->GetXaxis()->SetLabelOffset(0.02); // 1.01); // .02); //.01);
	  ratHist[ij]->Draw((ij==1 && isExternalError!=2) ? "hist" : "same:hist");
	}
	
	double xamx=0.60; // 0.70;
	if (nzone4>0 || ndata2<=3) xamx = 0.75; //0.85 ;
	///	double yamn = 0.9 - 0.02*(ndata2-1);
	double yamn = 0.9 - 0.042*(ndata2-1); 
	
	if (isSepRatio==1 || ij==1) {
	  	  ptst = new TPaveText(0.35,yamn,xamx,0.95,"brNDC");
	 // ptst = new TPaveText(0.45,0.50,0.50,0.68,"brNDC");
	  //	  ptst = new TPaveText(0.645,0.91*yamn,0.55,0.99,"brNDC");
	  ptst->SetFillColor(10);
	  ptst->SetBorderSize(0);
	  ptst->SetTextFont(42);
	}
	avedev = int(1000*avedev/errsq)/1000.;
	chi2[ij] = int(100*chi2[ij])/100.;
	if (isSepRatio==0) {
	  if (isCovariance_matrix) {
            sprintf (fname, "%s : #chi^{2}=%g", ttl2[ij], chi2[ij]);
            cout << chi2[ij] << endl;
            file_out << chi2[ij] << endl;
	   outfile <<  chi2[ij] <<endl; 
             outfile1<< "$"<<etapt<<"$" <<  " & " <<chi2[ij] <<"\\\\"<<endl;
                    outfile1<< "\\hline"<< endl;
          }
          else {	  
	    sprintf (fname, "%s  :  <#Delta> = %g", ttl2[ij], avedev);
	    cout << ttl2[ij] <<" ;AVEDEV =" << avedev << endl;
	    file_out << avedev << endl;
          }
	} else {
	  if (isCovariance_matrix) {
	    sprintf (fname, "#chi^{2}=%g", chi2[ij]);
	    cout << chi2[ij] << endl;
	    file_out << chi2[ij] << endl;
	   outfile <<  chi2[ij] <<endl; 
           outfile1<< "$"<<etapt<<"$" <<  " & " <<chi2[ij] <<"\\\\"<<endl;
                    outfile1<< "\\hline"<< endl;
             
	}
//	  outfile << "\n";
	  else {
	    sprintf (fname, "<#Delta>=%g", avedev);
	    cout << "Delta Value ===============================" <<avedev << endl;
	    file_out << avedev << endl;
	    //outfile<< " Bayes -"<< " & " << avedev <<" &  " <<endl;
                    //outfile<< "\\hline"<< endl; 
            outfile<<avedev<<endl;
         //if(titleindex==3 || titleindex==9 || titleindex==18 || titleindex==24){
           // outfile <<"\n" <<"\n";
            //}
	  }
	}
	TText* text = ptst->AddText(fname);
	if (isSepRatio==0) text->SetTextColor(icol);
	
	double txtsize = .06; // (isSepRatio>0) ? (alltogether>0) ? ((ij==ndata2-1) ? 0.05+0.01*ndata2 : 0.045+0.135*ndata2) : 0.051;
	
	if (isSepRatio>0) { 
	  if (ij==ndata2-1) {
	    txtsize = 0.05+0.023*ndata2;
	  } else {
	    txtsize = 0.045+0.043*ndata2;
	  }
	}
	if (isitForPaper>0 && mRatio!=0) txtsize *=0.73;
	if (isSlide&& mRatio!=0) txtsize *=1.5;
	if (isNote) txtsize *=1.4; 
	text->SetTextSize(1.1*txtsize); //GMAA
	
	ptst->SetBorderSize(0);
	if ((isSepRatio==0 || ij==ndata2-1)&& isExternalError<=2) {
         
        cout << "FNAME == " << fname << endl;
        ptst->Draw();
        }	
	if (isExternalError<=2) {
	  float yyv=1.0;
	  if (isAssym) yyv = 0.0;
	  
	  TLine lnn(alow, yyv, ahg+.5, yyv);
	  ratHist[ij]->GetXaxis()->GetXmin();
	  
	  TLine lnn();
	  lnn.SetLineColor(1);
	  lnn.SetLineWidth(1);
	  lnn.SetLineStyle(0);
	  
	  lnn.DrawLine(alow, yyv, ahg+.5, yyv);
//	  	  lnn.DrawLine(ratHist[ij]->GetXaxis()->GetXmin()+1.8, yyv, ratHist[ij]->GetXaxis()->GetXmax(), yyv);
	          cout << "RATTT HIST X MIN ======" << ratHist[ij]->GetXaxis()->GetXmin() << endl; 
       //CMS_lumi(c1,iPeriod,iPos);
	}
      } // if (ij>0) 
    } // for (int ij=0; ij<ndata; ij++) 
    
    if (ipad2>=1 && isIndex>0 && (ndata2<=2 || isSepRatio==0)) {
      TPaveText *ptst2 = new TPaveText(0.28,0.88,0.28,0.94,"brNDC");
      //	TPaveText *ptst2 = new TPaveText(0.78,0.88,0.88,0.94,"brNDC");
      ptst2->SetFillColor(10);
      ptst2->SetBorderSize(0);
      TText* text2 = ptst2->AddText(subttl[ipad2-1]);
      //      TText* text2 = ptst2->AddText(subttl[ipad]);
      text2->SetTextSize(.07);
      text2->SetTextFont(42);
      //ptst2->Draw();

      //ptst = new TPaveText(0.35,yamn,xamx,0.95,"brNDC");
	 // ptst = new TPaveText(0.45,0.50,0.50,0.68,"brNDC");
      TPaveText *ptst3 = new TPaveText(0.75,0.84,0.80,0.87,"brNDC");
      ptst3->SetFillColor(10);
      ptst3->SetBorderSize(0);
      sprintf(fname, "ndf = %g", nbinx-1);
      TText* text3 = ptst3->AddText(fname);
      text3->SetTextSize(.065);
      text3->SetTextFont(42);
      //ptst3->Draw();
    }
  }
  
  if (isExternalError==1) {

   cout<< "Test 1" << endl;
    char pcha[100];
    sprintf(pcha,"%s",histname);
    int len = strlen(pcha);
    pcha[len]='\0';
    len = strstr(pcha, "_"); 
    int len1 = len-int(pcha);
    char  histname2[50];
    strncpy (histname2, pcha, len1);
    histname2[len1]='\0';

    cout <<pcha<<" "<<endl;

    char* str = strstr(pcha, "_c"); 
 //   strcat(histname2, str);
   len = strlen(histname2);
    histname2[len]='\0';
    
    fileout<<"  double "<<histname<<"_xx"<<idError<<"["<<nbinx<<"] ={";
    cout<<"  double "<<histname<<"_xx"<<idError<<"["<<nbinx<<"] ={";

    for (int kl=0; kl<nbinx-1; kl++) {
      fileout<<bincenter[kl]<<",";
      cout<<bincenter[kl]<<",";
    }
    fileout<<bincenter[nbinx-1]<<"};"<<endl;
    cout<<bincenter[nbinx-1]<<"};"<<endl;

    fileout<<"  double " <<histname<< "_yy"<<idError<<"["<<nbinx<<"] ={";
    cout<<"  double " <<histname<< "_yy"<<idError<<"["<<nbinx<<"] ={";

    for (int kl=0; kl<nbinx-1; kl++) {
      fileout<<TMath::Max(1.e-10, devEfficiency[kl])<<",";
      cout<<TMath::Max(1.e-10, devEfficiency[kl])<<",";
      //      fileout<<TMath::Max(1.e-10, -devEfficiency[kl])<<",";
    }
    fileout<<TMath::Max(1.e-10,devEfficiency[nbinx-1])<<"};"<<endl;	
    cout<<TMath::Max(1.e-10,devEfficiency[nbinx-1])<<"};"<<endl;
    //    fileout<<TMath::Max(1.e-10,-devEfficiency[nbinx-1])<<"};"<<endl;

    fileout<<"  double "<<histname<<"_dyy"<<idError<<"["<<nbinx<<"] ={";
    cout<<"  double "<<histname<<"_dyy"<<idError<<"["<<nbinx<<"] ={";
    
    for (int kl=0; kl<nbinx-1; kl++) {
      fileout<<efficiencyErr[kl]<<",";
      cout<<efficiencyErr[kl]<<",";
    }
    fileout<<efficiencyErr[nbinx-1]<<"};"<<endl;
    cout<<efficiencyErr[nbinx-1]<<"};"<<endl;
  }
  if (isame2>0) { 
    sprintf(fname, "p%i", isame2, abs(ipteta));
  } else {
    sprintf(fname, "m%i", -isame2, abs(ipteta));
  }
  if (isame==0) {
    sprintf(fname, "%s_%s_%s_%i", fname, unf, etapt, abs(ipteta));
  }

  if (isame>0) { 
    sprintf(fname2, "com_%i_%s_%s_%s_ratio_16feb_8tevsys_%s.eps", isame, title, ttl2[0], ttl2[1], fname);
  } else {
    //    sprintf(fname2, "%s_12thjuly_v4_%s_sepratio.eps", title, fname);
    sprintf(fname2, "%s_%s_sepratio.eps", title, etapt);
  }
  cout <<"fname2 "<<fname2<<endl;
    //if (ipad==4)  c1->SaveAs(fname2);
     if (strstr(etapt,"typ_0_pt0_eta1_3")) {
      file_figname << "#!/bin/bash" << endl;
      file_figname << "pdftk " ;
     }
      sprintf(fname2, "%s.pdf",etapt);
      c1->SaveAs(fname2);
      file_figname << fname2 << " "; 
       if (strstr(etapt,"typ_0_pt7_eta1_24")) file_figname << "cat output all_var_jet.pdf"; 
    if ((ipad==0 || ipad==2) && idx>=0 ) {  c1->Update();  delete c1;}

  //  for (int ij=0; ij<ndata; ij++) {
  //    if (ratHist[ij]) { delete ratHist[ij]; ratHist[ij]=0;}
  //    if (ratHistOrg[ij]) { delete ratHistOrg[ij]; ratHistOrg[ij]=0;}
  //  }
//outfile << " \\end{tabular}" <<endl;
 //outfile << " \\end{center} "<<endl;

 //outfile << "\\end{table}"<<endl;
 //outfile << "\\end{document}"<<endl;
  if(typindex==2 && titleindex==24){
  outfile1 << " \\end{tabular}" <<endl;
  outfile1 << " \\end{center} "<<endl;

  outfile1 << "\\end{table}"<<endl;
  outfile1 << "\\end{document}"<<endl;
 }
file_figname.close();
outfile.close();
outfile1.close();
//cin>>yyy;
} //void compare



void convert_error(TH1F* thin) {
  int nbins = thin->GetNbinsX();
  for (int ix=1; ix<=nbins; ix++) {
    //    double val= thin->GetBinContent(ix);
    double err = thin->GetBinError(ix);
    thin->SetBinContent(ix, err);
    thin->SetBinError(ix, 0.0);
  }
}

void convert_relerror(TH1F* thin) {
  int nbins = thin->GetNbinsX();
  for (int ix=1; ix<=nbins; ix++) {
    double val= thin->GetBinContent(ix);
    double err = thin->GetBinError(ix);
    cout <<"== "<< ix<<" "<<thin->GetName()<<"; Val = "<<val <<  " ; Mean = " << thin->GetBinCenter(ix) <<" ; Error = "<< err << " ; Rel Error = " << err/TMath::Max(1.e-6,val)<<endl;
    if (val>1.e-6) {

      thin->SetBinContent(ix, err/val);
    } else {
      thin->SetBinContent(ix, 0.0);
    }
    thin->SetBinError(ix, 0.0);
  }

  //  thin->SetLineColor(ifil);
  //  thain->SetLineStyle(2*ior+1);
  //  thin->SetLineWidth(2*ior+1);
}




void unfolding_error(const char* unf="bayes", const char* varname="reco_typ_0_pt3_eta0_3", int var=5, int ipd=1) {
  // unfolding_error()
  // unfolding_error("svd", "Pythia6_thrustc_PFjet_c3_e1", 2);
  //  gStyle->SetPadGridX(3);
  //  gStyle->SetPadGridY(3);
  //  gStyle->SetGridStyle(0);
  /*gStyle->SetOptLogy(1);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetOptStat(0);
*/
 /* TStyle *gStyle = new TStyle("gStyle","Style for P-TDR");
  SetStyle st;
  st.SetPars(gStyle);*/
  setTDRStyle(gStyle);

  TH1F* histx[20];
  TH1F* histy[20];
  TH1F* hister[20];
  TH1F* histreler[20];
  char name[100];
  char  pch1[200];

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.055);
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // align right

  TLegend *tleg = new TLegend(0.25, 0.55, 0.55, 0.85,"","brNDC");
  tleg->SetFillColor(10);
  tleg->SetTextSize(0.09); //pss 0.07 : 0.058); //0.06); //pss
  tleg->SetBorderSize(0);
  tleg->SetTextFont(42);
  if (ipd==1) {
    TCanvas *c1 = new TCanvas("c1", "runfile", 700, 350);
    c1->Divide(3,1);
  }
  c1->cd(ipd);
  int icol=0;

  for (int ij=0; ij<8; ij++) {
    icol++; if (icol%5==0) icol++;
    // if (ij==0 || ij==3 ) continue;
    if (ij==0) {
      //      sprintf(name, "Roo_%s", varname);
      if (ipd==1) {
	//  printf(name, "Roo_Data_%s", varname);      
	sprintf(name, "Roo_Data_%s", varname);
	// sprintf(name, "reco_typ_0_pt0_eta0_3"); 
      } else if (ipd==2) {
	sprintf(name, "Roo_Data_%s", varname);
      }
      else {
       sprintf(name, "Roo_Data_%s", varname);
      }
      // cout<< "Hist name1=" << name << endl;
    }
    else if (ij==1) {
      // sprintf(name, "Roo_Data_reco_typ_0_pt4_eta0_3");
      sprintf(name, "bayes_Data_Jets_1_unfold_Roo_Pythia8_%s",varname);
      }
      else if (ij==2) {
    //  sprintf(name, "Roo_Madgh_%s",varname);
     sprintf(name, "bayes_Data_Jets_2_unfold_Roo_Pythia8_%s",varname);
     }   
    
       else if (ij==3) {
       sprintf(name, "bayes_Data_Jets_3_unfold_Roo_Pythia8_%s",varname);
       }   

       else if (ij==4) {
       sprintf(name, "bayes_Data_Jets_4_unfold_Roo_Pythia8_%s",varname);
       }   

       else if (ij==5) {
       sprintf(name, "bayes_Data_Jets_5_unfold_Roo_Pythia8_%s",varname);
       }   

       else if (ij==6) {
       sprintf(name, "bayes_Data_Jets_6_unfold_Roo_Pythia8_%s",varname);
       }   
    
    else {
      //sprintf(name, "%s_Data_Jets_%i_unfold_Roo_Pythia8_%s", unf, ij, varname);
      sprintf(name, "bayes_Data_Jets_7_unfold_Roo_Pythia8_%s",varname);
      //  sprintf(name, "Roo_Madgh_reco_typ_0_pt4_eta0_3");
      //      sprintf(name, "Roo_Madgh_reco_typ_0_pt3_eta0_3");

      cout<< "Hist name2=" << name << endl;
    }
    //else {
    // sprintf(name, "Roo_Pythia8_reco_typ_0_pt4_eta0_3", varname);
    //}
    // cout<< "Hist name1=" << name << endl;
//  }

  histx[ij] = (TH1F*)gDirectory->Get(name);
  histx[ij]->SetLineColor(icol);
  histy[ij] = (TH1F*)histx[ij]->Clone();

  if (ij==0) {
    histy[ij]->GetYaxis()->CenterTitle();
    histy[ij]->GetYaxis()->SetTitleSize(0.055);
    histy[ij]->GetYaxis()->SetTitleOffset(1.55);
    histy[ij]->GetYaxis()->SetLabelSize(0.055);
    histy[ij]->GetXaxis()->SetTickLength(0.05);
    histy[ij]->GetXaxis()->SetLabelSize(0.05);
    if (ipd==1) {
      if(var==3) histy[ij]->GetXaxis()->SetTitle("log #tau_{_{#perp} _{   ,C}} "); //histy[ij]->GetTitle());
      if(var==9) histy[ij]->GetXaxis()->SetTitle("#rho_{Tot,C}"); //histy[ij]->GetTitle());
      if(var==18) histy[ij]->GetXaxis()->SetTitle("B_{ T,C}"); //histy[ij]->GetTitle());
      if(var==24) histy[ij]->GetXaxis()->SetTitle("#rho^{T}_{Tot,C}"); //histy[ij]->GetTitle());
    } else if (ipd==2) {
      if(var==3) histy[ij]->GetXaxis()->SetTitle("log #tau_{_{#perp} _{   ,C}} "); //histy[ij]->GetTitle());
      if(var==9) histy[ij]->GetXaxis()->SetTitle("#rho_{Tot,C}"); //histy[ij]->GetTitle());
      if(var==18) histy[ij]->GetXaxis()->SetTitle("B_{ T,C}"); //histy[ij]->GetTitle());
      if(var==24) histy[ij]->GetXaxis()->SetTitle("#rho^{T}_{Tot,C}"); //histy[ij]->GetTitle());     
    } else {
      if(var==3) histy[ij]->GetXaxis()->SetTitle("log #tau_{_{#perp} _{   ,C}} "); //histy[ij]->GetTitle());
      if(var==9) histy[ij]->GetXaxis()->SetTitle("#rho_{Tot,C}"); //histy[ij]->GetTitle());
      if(var==18) histy[ij]->GetXaxis()->SetTitle("B_{ T,C}"); //histy[ij]->GetTitle());
      if(var==24) histy[ij]->GetXaxis()->SetTitle("#rho^{T}_{Tot,C}"); //histy[ij]->GetTitle());
    }
    histy[ij]->GetXaxis()->SetTitleColor(1);
    histy[ij]->GetXaxis()->CenterTitle();
    histy[ij]->GetXaxis()->SetTitleSize(0.055);
    histy[ij]->GetXaxis()->SetTitleOffset(1.10);
    histy[ij]->GetXaxis()->SetLabelFont(42);
    histy[ij]->GetXaxis()->SetTitleFont(42);
  }

  hister[ij] = (TH1F*)histy[ij]->Clone();
  histreler[ij] = (TH1F*)histy[ij]->Clone();
  convert_error(hister[ij]);
  convert_relerror(histreler[ij]);
  histreler[ij]->SetMaximum(0.21); //1.032);
  histreler[ij]->SetTitle(" ");
  setTDRStyle();
  histreler[ij]->GetYaxis()->SetTitle("Rel error (%)");
  histreler[ij]->Draw((ij==0)?"":"same");


  //hister[ij]->GetYaxis()->SetTitle("Rel error (%)");
  // hister[ij]->Draw((ij==0)?"":"same");
  if (ij==0) {
    tleg->AddEntry(hister[ij], "Data","lf"); //pss
  }

 /* else if (ij==1) {
    tleg->AddEntry(hister[ij], "Madgh","lf"); //pss
  }
  else if (ij==2){
    tleg->AddEntry(hister[ij], "Madgh","lf");
  }*/
  /* 
     if (ij==3) {
     tleg->AddEntry(hister[ij], "iter 3","lf"); //pss
     } 

     if (ij==4) {
     tleg->AddEntry(hister[ij], "iter 4","lf"); //pss
     } 

     if (ij==5) {
     tleg->AddEntry(hister[ij], "iter 5","lf"); //pss
     }

     if (ij==6) {
     tleg->AddEntry(hister[ij], "iter 6","lf"); //pss
     }*/
   else {
  tleg->AddEntry(hister[ij], Form("iter=%i", ij),"lf"); //pss
  }
}
//  latex.DrawLatex(.26, .91, subttlx[ipd-1]);
tleg->Draw();
}


void tdrGrid(bool gridOn) {
  gStyle->SetPadGridX(gridOn);
  gStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle(TStyle *gStyle) {
  TStyle *gStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

// For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

// For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  
// For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  // gStyle->SetErrorMarker(20);
  //gStyle->SetErrorX(0.);
  
  gStyle->SetMarkerStyle(20);
  
//For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

//For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

// Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02);

// For the Global title:

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // gStyle->SetTitleH(0); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

// For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

// Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(1);
  gStyle->SetOptLogz(0);

// Postscript options:
  gStyle->SetPaperSize(20.,20.);
  // gStyle->SetLineScalePS(Float_t scale = 3);
  // gStyle->SetLineStyleString(Int_t i, const char* text);
  // gStyle->SetHeaderPS(const char* header);
  // gStyle->SetTitlePS(const char* pstitle);

  // gStyle->SetBarOffset(Float_t baroff = 0.5);
  // gStyle->SetBarWidth(Float_t barwidth = 0.5);
  // gStyle->SetPaintTextFormat(const char* format = "g");
  // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // gStyle->SetTimeOffset(Double_t toffset);
  // gStyle->SetHistMinimumZero(kTRUE);

  gStyle->SetHatchesLineWidth(5);
  gStyle->SetHatchesSpacing(0.05);

  gStyle->cd();

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
      latex.SetTextFont(9cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextFont(42);
  latex.SetTextAlign(31);
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(0.495-l,1-t+lumiTextOffset*t,cmsText); //all other plot
      //latex.DrawLatex(0.50-l,1-t+lumiTextOffset*t,cmsText); //only for delta
      //latex.DrawLatex(1-r,1-t+lumiTextOffset*t,cmsText);
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
      latex.DrawLatex(0.535-posX_, posY_, extraText);  //on for all plots  
      //latex.DrawLatex(0.24-posX_, posY_, extraText);   // only for delta    
    }
  return;
}
