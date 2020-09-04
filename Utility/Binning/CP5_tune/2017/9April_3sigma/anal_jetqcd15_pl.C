const char* varname[32]={"thrust", "minor", "y3", "thrustc", "thrustce", "thrustcr",
                           "minorc", "minorce", "minorcr", "tmass", "tmasse", "tmassr",
                           "hmass", "hmasse", "hmassr", "y3c", "y3ce", "y3cr",
                           "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
                           "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
                           "sphericity", "cparameter"};
Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kTRUE));
}

const int nfilemx=16;
const int colcode[nfilemx]={1,2,3,4,6,7,11,kCyan+2,kCyan-10, kBlue-7, kMagenta-9, kRed-5, kYellow-6, kGreen-6, kRed-7, kGreen+3};
void def_setting() {

gStyle->SetCanvasBorderMode(0);
gStyle->SetPadBorderMode(0);
gStyle->SetStatBorderSize(1);
gStyle->SetStatStyle(1001);
gStyle->SetCanvasColor(10);
gStyle->SetPadColor(10);
gStyle->SetStatColor(10);

gStyle->SetTitleFontSize(0.04); //Global histogram title
gStyle->SetStatFont(22);        // Times New Roman
gStyle->SetTextFont(22);        // Times New Roman
gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
gStyle->SetLabelFont(22,"XYZ"); // Times New Roman

gStyle->SetPadBottomMargin(0.12);
gStyle->SetPadTopMargin(0.05);

gStyle->SetPadLeftMargin(0.10);
gStyle->SetPadRightMargin(0.08);
gStyle->SetStatX(0.95);
gStyle->SetStatY(0.95);
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.1);

gStyle->SetTitlePS("error parametrization");
gStyle->SetOptStat(0011100);   // print fit results in the stat box
gStyle->SetOptFit(0110);
gStyle->SetOptDate(0);
 gStyle->SetOptStat(1); //to have error, remove this, also change gStyle->SetStatH(0.4); 
gStyle->SetOptTitle(0);
gStyle->SetMarkerColor(4444);
gStyle->SetMarkerStyle(4);
gStyle->SetMarkerSize(0.2);
gStyle->SetHistLineColor(1);
gStyle->SetHistLineWidth(2);


  gStyle->SetPalette(1,0); 
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.02); 
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetTitleFontSize(0.09);
  gStyle->SetTitleOffset(-0.09);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetLabelSize(0.095,"XY");  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");

  gStyle->SetNdivisions(404,"XYZ");

  gStyle->SetStatTextColor(1);
  gStyle->SetStatX(.99);
  gStyle->SetStatY(.99);  
  gStyle->SetStatW(.3);
  gStyle->SetStatH(.1);

  gStyle->SetOptLogy(0);
  gStyle->SetOptLogx(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //  gStyle->SetStripDecimals (kTRUE);
  // h->GetXaxis()->SetTimeFormat("%d/%m/%Y");

  gStyle->SetHistLineWidth(3);
}




void plot2d(int ipt=3, int ieta=0, int ivar=3) {
  // ipt : pt bins
  // ieta : eta bin
  // ivar : variable
  
  def_setting();
  
  TH2F* histxx[20];
  TH2F* histyy[20];
  const int ntype=5;
  const char* typname[ntype]={"Jets", "All particle",   "Charged particle", "All particle: P_{T}>0.5", "Charged particle: P_{T}>0.5"};

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 500);
  c1->Divide(3,2, 1.e-6, 1.e-6);
  char name[100];

  for (int ij=0; ij<ntype; ij++) {
    sprintf(name,"analyzebasicpat/corr_typ_%i_pt%i_eta%i_%i", ij, ipt, ieta, ivar);
    cout<<"name "<<name<<endl;
    histxx[ij] = (TH2F*)gDirectory->Get(name);
    
    c1->cd(ij+1);
    histxx[ij]->Draw("colz");
  }
}


void resolution(int ityp=0, int ivar=3, int ipt=0, int ieta=0, double alow=-10){ 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitle(0);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(.12);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetTitleSize(0.095,"XYZ");
  gStyle->SetLabelSize(0.075,"XYZ");
  gStyle->SetLabelOffset(-0.02, "X");
  gStyle->SetTitleOffset(0.5, "X");
  gStyle->SetTitleOffset(0.75, "Y");

  gStyle->SetNdivisions(406,"x");
  gStyle->SetNdivisions(1212,"Y");

 

  char varname2[100]=" log #tau_{_{#perp} _{   ,C}}";
  TH2F* histx[10];
  TH1F* hist1x[10][200];
  TF1*  fitx[10][200];

  TGraphErrors* gry[10];
  TGraphErrors* gr1[10];
  TGraphErrors* regr[10];

  double peak[200]={0};
  double xval[200]={0};
  double xerr[200]={0};
  double yval[200]={0};
  double yerr[200]={0};

  double ywid[200]={0};
  double ywer[200]={0};

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // align right

  char name[200];
  char title[200];
  
 

  //  if (ipt==0) {
    TCanvas *c0 = new TCanvas("c0", "runfile0", 500, 350);
    gStyle->SetPadRightMargin(0.06);
    TCanvas *c1 = new TCanvas("c1", "runfile", 800, 650);
    c1->Divide(6,5);
    //  }

  for (int ij=ipt; ij<ipt+1; ij++) {
    //    c1->cd(ij+1);
    // sprintf(name,"corr_typ_%i_pt%i_eta%i_%i", ityp, ipt, ieta, ivar);
    sprintf(name,"analyzeBasicPat/corr_typ_%i_pt%i_eta%i_%i", ityp, ipt, ieta, ivar);
//corr_typ_0_pt7_eta0_9
     cout << "Corr name =  "<< name << endl;
    histx[ij] =(TH2F*)gDirectory->Get(name);
    //    histx[ij]->GetXaxis()->SetRangeUser(-9., -2.);
    //    histx[ij]->GetYaxis()->SetRangeUser(-9.0, -2.0);
    c0->cd();
    //    histx[ij]->SetMinimum(0.01*histx[ij]->GetMaximum());
    histx[ij]->Draw("colz");

    int nhstbins = histx[ij]->GetNbinsX();
    cout <<"name "<<name<<" "<<nhstbins<<endl;
    double apeakmx=0;
    int nbins = 0;
    int igap=2;
    for (int jk=1; jk<=nhstbins; jk +=igap) {
      apeakmx=0;
      c1->cd(nbins+1);
      sprintf(title, "%s_%i", varname[ivar], jk);

      hist1x[ij][jk] = (TH1F*)histx[ij]->ProjectionY(title, jk, jk+igap-1, "e");
      hist1x[ij][jk]->GetXaxis()->SetLabelSize(0.10);
      hist1x[ij][jk]->GetYaxis()->SetLabelSize(0.07);

      double bincon = hist1x[ij][jk]->GetMaximum();
      if ( bincon<1.e-10) continue;
      double amn =hist1x[ij][jk]->GetXaxis()->GetXmin();
      double amx =hist1x[ij][jk]->GetXaxis()->GetXmax();

      sprintf(name, "fitx_%i_%i", ij, jk);
      fitx[ij][jk] = new TF1(name, gausX, amn, amx, 3);  
      double tmpbinwid=hist1x[ij][jk]->GetBinWidth(nhstbins/2.);
      double par[3]={hist1x[ij][jk]->GetMaximum(), hist1x[ij][jk]->GetBinCenter(hist1x[ij][jk]->GetMaximumBin()), 4*tmpbinwid}; 
      fitx[ij][jk]->SetParameters(par);
      fitx[ij][jk]->SetParLimits(2, tmpbinwid, 10*tmpbinwid);
      fitx[ij][jk]->SetLineColor(2);
      fitx[ij][jk]->SetLineWidth(1);

      hist1x[ij][jk]->Fit(name, "BQ"); //"gaus","SQ");
      double hght=fitx[ij][jk]->GetParameter(0);
      double mean=fitx[ij][jk]->GetParameter(1);
      double rms=fitx[ij][jk]->GetParameter(2);
      fitx[ij][jk]->SetParLimits(0, 0.5*hght, 2.0*hght);
      fitx[ij][jk]->SetParLimits(1, 0.5*mean, 2.0*mean);
      fitx[ij][jk]->SetParLimits(2, 0.5*rms, 5.0*rms);

      hist1x[ij][jk]->Fit(name, "BQ");
      hght=fitx[ij][jk]->GetParameter(0);
      mean=fitx[ij][jk]->GetParameter(1);
      rms=fitx[ij][jk]->GetParameter(2);
      fitx[ij][jk]->SetParLimits(0, 0.5*hght, 2.0*hght);
      fitx[ij][jk]->SetParLimits(1, 0.5*mean, 2.0*mean);
      fitx[ij][jk]->SetParLimits(2, 0.5*rms, 2.0*rms);

      hist1x[ij][jk]->Fit(name, "BQ");
      //      TFitResultPtr ptr = hist1x[ij][jk]->Fit(name, "SQ"); //"gaus","SQ");
      
      peak[nbins] = fitx[ij][jk]->GetParameter(0);
      if (peak[nbins]>apeakmx) apeakmx = peak[nbins];
      if (igap%2==0) { 
	xval[nbins] = histx[ij]->GetXaxis()->GetBinLowEdge(jk+igap/2); //Center(jk);
      } else {
	xval[nbins] = histx[ij]->GetXaxis()->GetBinCenter(jk+igap/2);
      }

      yval[nbins] = fitx[ij][jk]->GetParameter(1);
      yerr[nbins] = fitx[ij][jk]->GetParError(1);
      ywid[nbins] = fitx[ij][jk]->GetParameter(2);
      ywer[nbins] = fitx[ij][jk]->GetParError(2);
      hist1x[ij][jk]->GetXaxis()->SetNdivisions(204);
      hist1x[ij][jk]->GetYaxis()->SetNdivisions(1212);
      hist1x[ij][jk]->GetXaxis()->SetRangeUser(yval[nbins]-6*ywid[nbins], yval[nbins]+6*ywid[nbins]);
      //      hist1x[ij][jk]->SetMaximum(1.1*hist1x[ij][jk]->GetMaximum());
      hist1x[ij][jk]->Draw();
      
      latex.DrawLatex(.6, .9, Form("[%g : %g]",histx[ij]->GetXaxis()->GetBinLowEdge(jk),histx[ij]->GetXaxis()->GetBinLowEdge(jk+igap/2))); 


      yval[nbins] -=xval[nbins];


      nbins++;
    }

    if (nbins<4) continue; 
    //    if (ipt==0) {
      TCanvas *c2 = new TCanvas("c2", "runfile2", 700, 350);
      c2->Divide(3,1);
      //    }
    gry[ij] = new TGraphErrors(nbins, xval, yval, xerr, yerr);
    gr1[ij] = new TGraphErrors(nbins, xval, ywid, xerr, ywer);
    
    gry[ij]->SetMarkerSize(0.8); gry[ij]->SetMarkerColor(2.0);
    gry[ij]->SetMarkerStyle(4); 
//    gry[ij]->SetMarkerSize(1.5);
    gry[ij]->GetXaxis()->SetTitle(varname2);
    gry[ij]->GetXaxis()->CenterTitle();
    gry[ij]->GetYaxis()->SetTitle("#Delta#mu");
    gry[ij]->GetYaxis()->CenterTitle();

    gr1[ij]->SetMarkerSize(0.8); gr1[ij]->SetMarkerColor(2.0); 
    gr1[ij]->SetMarkerStyle(4);  
    gr1[ij]->GetXaxis()->SetTitle(varname2);
    gr1[ij]->GetXaxis()->CenterTitle();
    gr1[ij]->GetYaxis()->SetTitle("#sigma");
    gr1[ij]->GetYaxis()->CenterTitle();
    
    c2->cd(1); gry[ij]->Draw("ALP");// latex.DrawLatex(0.8, 0.9, "(a)");
    c2->cd(2); gr1[ij]->Draw("ALP");// latex.DrawLatex(0.8, 0.9, "(b)");
    
    double hbinwid=0.5*histx[ij]->GetXaxis()->GetBinWidth(1);
    int nrebin=0;
    double xrewid[200];
    double yrewid[200];
    double yreerr[200];
    for (int jk=10; jk<nbins; jk++) {
      if (fabs(yval[jk])<0.2 && yerr[jk]<.1 && peak[jk]>0.0001*apeakmx) {
	xrewid[nrebin] = xval[jk];
	yrewid[nrebin] = ywid[jk];
	yreerr[nrebin] = ywer[jk];

	nrebin++;
      }
    }
    
    for (int jk=0; jk<nbins; jk++) {
      cout<<xval[jk]<<", ";
    }
    cout<<endl;
    
    for (int jk=0; jk<nbins; jk++) {
      cout<<ywid[jk]<<", ";
    }
    cout<<endl;
    

    c2->cd(3);
    regr[ij] = new TGraphErrors(nrebin, xrewid, yrewid, xerr, yreerr);
    regr[ij]->SetMarkerSize(1.1); regr[ij]->SetMarkerColor(2.0);
    regr[ij]->GetXaxis()->SetTitle(varname2);
    regr[ij]->GetXaxis()->CenterTitle();
    regr[ij]->GetYaxis()->SetTitle("#sigma");
    regr[ij]->GetYaxis()->CenterTitle();

    TFitResultPtr ptrre = regr[ij]->Fit("pol3","SQ");
    regr[ij]->Draw("ALP"); // latex.DrawLatex(0.8, 0.9, "(c)");
    double p0 = ptrre->Parameter(0);
    double p1 = ptrre->Parameter(1);
    double p2 = ptrre->Parameter(2);
    double p3 = ptrre->Parameter(3);
    
    double width = xrewid[nrebin-1] - xrewid[0] + 2*hbinwid;
    int nfibins=0;
    double xfibins[200];
    
    double gap = 0;

    double xpoint = xrewid[0] + 2*hbinwid;
    xpoint = int(100*xpoint+0.05)/100.;

    while(gap<width) {
      xfibins[nfibins] = xpoint;
      double tmpwd = 3*(p0 + p1*xpoint + p2*xpoint*xpoint+p3*xpoint*xpoint*xpoint);
      if (tmpwd<hbinwid) tmpwd =hbinwid; 
      tmpwd = int(100*tmpwd+0.5)/100.0;
      gap +=tmpwd;
      xpoint +=tmpwd;

      nfibins++;
      if (nfibins>=100) break;
    }
    xfibins[nfibins] = xpoint;

    cout <<"double "<<varname[ivar]<<"_"<<ij<<"["<<nfibins+1<<"] ={";
    for (int jk=0; jk<nfibins; jk++) {
      cout <<xfibins[jk]<<", ";
    }
    cout <<xfibins[nfibins]<<" }; "<<endl;
  }



}

void resolution2(int ivar=3, int ipt=0, int ieta=0, double alow=-10){ 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitle(0);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(.12);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetTitleSize(0.095,"XYZ");
  gStyle->SetLabelSize(0.075,"XYZ");
  gStyle->SetLabelOffset(-0.02, "X");
  gStyle->SetTitleOffset(0.5, "X");
  gStyle->SetTitleOffset(0.75, "Y");

  gStyle->SetNdivisions(406,"x");
  gStyle->SetNdivisions(1212,"Y");

 

  char varname2[100]=" log #tau_{_{#perp} _{   ,C}}";
  TH2F* histx[40];
  TH1F* hist1x[40][200];
  TF1*  fitx[40][200];
  TF1*  fitxy[40][200]={0};

  TGraphErrors* gry[10];
  TGraphErrors* gr1[10];
  TGraphErrors* regr[10];

  double peak[200]={0};
  double xval[200]={0};
  double xerr[200]={0};
  double yval[200]={0};
  double yerr[200]={0};

  double ywid[200]={0};
  double ywer[200]={0};

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.SetTextFont(42);
  latex.SetTextAlign(31); // align right

  char name[200];
  char title[200];
  
 

  TCanvas *c0 = new TCanvas("c0", "runfile0", 500, 350);
  gStyle->SetPadRightMargin(0.06);
  TCanvas *c1 = new TCanvas("c1", "runfile", 800, 650);
  c1->Divide(6,5);
  TCanvas *c2 = new TCanvas("c2", "runfile2", 700, 350);

  int icolx=-1;
  for (int ityp=0; ityp<3; ityp++) {
    icolx++;
  for (int ij=ipt; ij<ipt+1; ij++) {
    //    c1->cd(ij+1);
     sprintf(name,"analyzeBasicPat/corr_typ_%i_pt%i_eta%i_%i", ityp, ipt, ieta, ivar);

    histx[ij] =(TH2F*)gDirectory->Get(name);
    //    histx[ij]->GetXaxis()->SetRangeUser(-9., -2.);
    //    histx[ij]->GetYaxis()->SetRangeUser(-9.0, -2.0);
    c0->cd();
    //    histx[ij]->SetMinimum(0.01*histx[ij]->GetMaximum());
    histx[ij]->Draw("colz");

    int nhstbins = histx[ij]->GetNbinsX();
    cout <<"name "<<name<<" "<<nhstbins<<endl;
    double apeakmx=0;
    int nbins = 0;
    int igap=2;
    for (int jk=1; jk<=nhstbins; jk +=igap) {
      apeakmx=0;
      c1->cd(nbins+1);
      sprintf(title, "%s_%i", varname[ivar], jk);

      hist1x[ij][jk] = (TH1F*)histx[ij]->ProjectionY(title, jk, jk+igap-1, "e");
      hist1x[ij][jk]->GetXaxis()->SetLabelSize(0.12);
      hist1x[ij][jk]->GetYaxis()->SetLabelSize(0.12);

      double bincon = hist1x[ij][jk]->GetMaximum();
      if ( bincon<1.e-10) continue;
      double amn =hist1x[ij][jk]->GetXaxis()->GetXmin();
      double amx =hist1x[ij][jk]->GetXaxis()->GetXmax();

      sprintf(name, "fitx_%i_%i", ij, jk);
      fitx[ij][jk] = new TF1(name, gausX, amn, amx, 3);  
      double tmpbinwid=hist1x[ij][jk]->GetBinWidth(nhstbins/2.);
      double par[3]={hist1x[ij][jk]->GetMaximum(), hist1x[ij][jk]->GetBinCenter(hist1x[ij][jk]->GetMaximumBin()), 4*tmpbinwid}; 
      fitx[ij][jk]->SetParameters(par);
      fitx[ij][jk]->SetParLimits(2, tmpbinwid, 10*tmpbinwid);
      fitx[ij][jk]->SetLineColor(2);
      fitx[ij][jk]->SetLineWidth(1);

      hist1x[ij][jk]->Fit(name, "BQ"); //"gaus","SQ");
      double hght=fitx[ij][jk]->GetParameter(0);
      double mean=fitx[ij][jk]->GetParameter(1);
      double rms=fitx[ij][jk]->GetParameter(2);
      fitx[ij][jk]->SetParLimits(0, 0.5*hght, 2.0*hght);
      fitx[ij][jk]->SetParLimits(1, 0.5*mean, 2.0*mean);
      fitx[ij][jk]->SetParLimits(2, 0.5*rms, 5.0*rms);

 /*     hist1x[ij][jk]->Fit(name, "BQ");
      hght=fitx[ij][jk]->GetParameter(0);
      mean=fitx[ij][jk]->GetParameter(1);
      rms=fitx[ij][jk]->GetParameter(2);
      fitx[ij][jk]->SetParLimits(0, 0.5*hght, 2.0*hght);
      fitx[ij][jk]->SetParLimits(1, 0.5*mean, 2.0*mean);
      fitx[ij][jk]->SetParLimits(2, 0.5*rms, 2.0*rms);

      hist1x[ij][jk]->Fit(name, "BQ");*/
      //      TFitResultPtr ptr = hist1x[ij][jk]->Fit(name, "SQ"); //"gaus","SQ");
      
      peak[nbins] = fitx[ij][jk]->GetParameter(0);
      if (peak[nbins]>apeakmx) apeakmx = peak[nbins];
      if (igap%2==0) { 
	xval[nbins] = histx[ij]->GetXaxis()->GetBinLowEdge(jk+igap/2); //Center(jk);
      } else {
	xval[nbins] = histx[ij]->GetXaxis()->GetBinCenter(jk+igap/2);
      }

      yval[nbins] = fitx[ij][jk]->GetParameter(1);
      yerr[nbins] = fitx[ij][jk]->GetParError(1);
      ywid[nbins] = fitx[ij][jk]->GetParameter(2);
      ywer[nbins] = fitx[ij][jk]->GetParError(2);
      hist1x[ij][jk]->GetXaxis()->SetNdivisions(204);
      hist1x[ij][jk]->GetYaxis()->SetNdivisions(1212);
      hist1x[ij][jk]->GetXaxis()->SetRangeUser(yval[nbins]-6*ywid[nbins], yval[nbins]+6*ywid[nbins]);
      //      hist1x[ij][jk]->SetMaximum(1.1*hist1x[ij][jk]->GetMaximum());
      hist1x[ij][jk]->Draw();
      
      latex.DrawLatex(.6, .9, Form("[%g : %g]",histx[ij]->GetXaxis()->GetBinLowEdge(jk),histx[ij]->GetXaxis()->GetBinLowEdge(jk+igap/2))); 


      yval[nbins] -=xval[nbins];


      nbins++;
    }

    if (nbins<4) continue; 

    gry[ij] = new TGraphErrors(nbins, xval, yval, xerr, yerr);
    gr1[ij] = new TGraphErrors(nbins, xval, ywid, xerr, ywer);
    
    gry[ij]->SetMarkerSize(0.8); gry[ij]->SetMarkerColor(2.0); 
    gry[ij]->GetXaxis()->SetTitle(varname2);
    gry[ij]->GetXaxis()->CenterTitle();
    gry[ij]->GetYaxis()->SetTitle("#Delta#mu");
    gry[ij]->GetYaxis()->CenterTitle();

    gr1[ij]->SetMarkerSize(0.8); gr1[ij]->SetMarkerColor(2.0); 
    gr1[ij]->GetXaxis()->SetTitle(varname2);
    gr1[ij]->GetXaxis()->CenterTitle();
    gr1[ij]->GetYaxis()->SetTitle("#sigma");
    gr1[ij]->GetYaxis()->CenterTitle();
    
    //    c2->cd(1); gry[ij]->Draw("ALP");// latex.DrawLatex(0.8, 0.9, "(a)");
    //    c2->cd(2); gr1[ij]->Draw("ALP");// latex.DrawLatex(0.8, 0.9, "(b)");
    
    double hbinwid=0.5*histx[ij]->GetXaxis()->GetBinWidth(1);
    int nrebin=0;
    double xrewid[200];
    double yrewid[200];
    double yreerr[200];
    for (int jk=60; jk<nbins; jk++) {
      if (fabs(yval[jk])<0.2 && yerr[jk]<.1 && peak[jk]>0.0001*apeakmx) {
	xrewid[nrebin] = xval[jk];
	yrewid[nrebin] = ywid[jk];
	yreerr[nrebin] = ywer[jk];
      cout << "jk =" << jk << " nrebin = " << nrebin << "; xrewid " << xrewid[nrebin] << " yrewid " << yrewid[nrebin]  << endl;

	nrebin++;
      }
    }
    
    for (int jk=0; jk<nbins; jk++) {
      cout<<xval[jk]<<", ";
    }
    cout<<endl;
    
    for (int jk=0; jk<nbins; jk++) {
      cout<<ywid[jk]<<", ";
    }
    cout<<endl;
    

    c2->cd();
    regr[ij] = new TGraphErrors(nrebin, xrewid, yrewid, xerr, yreerr);
    regr[ij]->SetMarkerSize(1.1);// regr[ij]->SetMarkerColor(2.0);
    regr[ij]->GetXaxis()->SetTitle(varname2);
    regr[ij]->GetXaxis()->CenterTitle();
    regr[ij]->GetYaxis()->SetTitle("#sigma");
    regr[ij]->GetYaxis()->CenterTitle();
    regr[ij]->SetMarkerColor(colcode[icolx]);
    regr[ij]->SetMarkerStyle(24);

    //    gStyle->SetFuncColor(colcode[icolx]);
    //    gStyle->SetLineColor(colcode[icolx]);
    regr[ij]->SetMaximum(1.0);
    TFitResultPtr ptrre = regr[ij]->Fit("pol3","SQ0");

    if (ityp==0) {
      regr[ij]->Draw("AP"); // latex.DrawLatex(0.8, 0.9, "(c)");
    } else {
      regr[ij]->Draw("P:same");
    }
    double p0 = ptrre->Parameter(0);
    double p1 = ptrre->Parameter(1);
    double p2 = ptrre->Parameter(2);
    double p3 = ptrre->Parameter(3);
   int jk=0; 
    sprintf(name, "tmp_%i", ityp);
    fitxy[ityp][jk] = new TF1(name, "[0]+[1]*x+[2]*x*x+[3]*x*x*x", -30.0, 0.0);
    fitxy[ityp][jk]->SetParameter(0, p0);
    fitxy[ityp][jk]->SetParameter(1, p1);
    fitxy[ityp][jk]->SetParameter(2, p2);
    fitxy[ityp][jk]->SetParameter(3, p3);
    fitxy[ityp][jk]->SetLineColor(colcode[icolx]);
    
    fitxy[ityp][jk]->Draw("same");

    /*
    double width = xrewid[nrebin-1] - xrewid[0] + 2*hbinwid;
    int nfibins=0;
    double xfibins[200];
    
    double gap = 0;

    double xpoint = xrewid[0] + 2*hbinwid;
    xpoint = int(100*xpoint+0.05)/100.;

    while(gap<width) {
      xfibins[nfibins] = xpoint;
      double tmpwd =1.50*(p0 + p1*xpoint + p2*xpoint*xpoint+p3*xpoint*xpoint*xpoint);
      if (tmpwd<hbinwid) tmpwd =hbinwid; 
      tmpwd = int(100*tmpwd+0.5)/100.0;
      gap +=tmpwd;
      xpoint +=tmpwd;

      nfibins++;
      if (nfibins>=100) break;
    }
    xfibins[nfibins] = xpoint;

    cout <<"double "<<varname[ivar]<<"_"<<ij<<"["<<nfibins+1<<"] ={";
    for (int jk=0; jk<nfibins; jk++) {
      cout <<xfibins[jk]<<", ";
    }
    cout <<xfibins[nfibins]<<" }; "<<endl;
    */
  }
  }


}

