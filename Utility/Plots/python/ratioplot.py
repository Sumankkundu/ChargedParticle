from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT, gSystem
from ROOT import kBlack, kBlue, kRed
import array as arr
gROOT.SetBatch()

#from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
#from matplotlib.figure import Figure
#import numpy as np
#from python
#from PyPDF2 import PdfFileMerger

#from sys 
import PyPDF2


pdfWriter = PyPDF2.PdfFileWriter()
pdfOutputFile = open('MergedFiles.pdf', 'wb')


nHLTmx =8 #HT2 Range
nusedvar = 5   #Event Shape variables used
ntype = 2 
#histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100]
color = arr.array('i',[2,4,6,30,6,46,3,28,38,42])  #define the color for different histograms
var = [3,9,15,18,24]
HT2range = [84, 111, 172, 241, 309, 376, 462, 568, 3000] # 2017 1June21 UL17
#  int  HT2range[nHLTmx+1]={83, 109, 172, 241, 309, 377, 462, 570, 3000}; //2017
#  int HT2range[nHLTmx+1]={83, 109, 176, 247, 318, 387, 477, 573, 3000}; //2018 25May21 UL18



unfolddata = TFile('Unfolded_Result_Data_PY8.root',"READ") # 

Data = TFile( 'Data_19UL18_JECV5_25May21.root' )
PYCP5 = TFile( 'Pythia8_Flat_CP5_Tune.root',"READ" )
PYCP52018 = TFile( 'PY8_19UL18_Flat_JRV2_25May21.root',"READ" )
PYMPIoff = TFile('Pythia8_Flat_CP5_mpioff_20June21.root',"READ" )
PYISRoff = TFile('Pythia8_Flat_CP5_ISRoff_20June21.root', "READ" )
PYFSRoff = TFile('Pythia8_Flat_CP5_FSRoff_20June21.root', "READ" )
inputroot = [PYCP52018, PYMPIoff, PYISRoff, PYFSRoff]
modelname = ["Py8 CP5 tune", "Py8 MPI off", "Py8 ISR off", "Py8 FSR off"]  #define the model name

Esvsym = ["#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"]
Esvname = ["Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"]
htrang=["83 < H_{T,2} < 109", "109 < H_{T,2} < 172", "172 < H_{T,2} < 241", "241 < H_{T,2} < 309","309 < H_{T,2} < 377", "377 < H_{T,2} < 462", "462 < H_{T,2} <570","H_{T,2}"]

Esvlogx = ["ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"]
Esvlogy = ["1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"]
Unit = ["GeV","fb^{-1}","Pb","#%"]



if not Data :
   print ("Failed to get data")

def ReadH1(froot,histname):
    h1 = froot.Get(histname)
    #h1.SetLineColor(kBlue+1)
    h1.SetLineWidth(2)
    #h1.FillRandom("gaus")
    h1.GetYaxis().SetTitleSize(20)
    h1.GetYaxis().SetTitleFont(43)
    h1.GetYaxis().SetTitleOffset(1.55)
    h1.SetStats(0)
    return h1
 
 
def createH2():
    h2 = TH1F("h2", "h2", 100, -5, 5)
    h2.FillRandom("gaus")
    h2.SetLineColor(kRed)
    h2.SetLineWidth(2)
    return h2
 
 
def createRatio(h1, h2, ymax, ymin, h1name, h2name ):
    #h1 -> Data
    #h2->MC we use MC/ Data 
    h3 = h2.Clone("h3")
    h3.SetLineColor(h2.GetLineColor())
    h3.SetMarkerStyle(h2.GetMarkerStyle())
    h3.SetTitle("")
    #h3.Sumw2()
    #h3.SetStats(0)
    h3.Divide(h1)
    h3.SetMinimum(ymin)
    h3.SetMaximum(ymax)
 
    return h3

def sethist(h1, xtitle, ytitle, xtitlesize, ytitlesize, xlablesize, ylablesize):
    h1.SetTitle("")
    y = h1.GetYaxis()
    y.SetTitle(ytitle)
    y.CenterTitle()
    y.SetNdivisions(505)
    y.SetTitleSize(ytitlesize)
    y.SetTitleFont(43)
    y.SetTitleOffset(2)
    y.SetLabelFont(43)
    y.SetTickLength(0.06)
    y.SetLabelSize(ylablesize)

    x = h1.GetXaxis()
    x.CenterTitle()
    x.SetTitle(xtitle)
    x.SetTitleSize(xtitlesize)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetTickLength(0.06)
    x.SetLabelSize(xlablesize)
    
    return 0

def divBybinWidth(h1):
    tmpnbn = h1.GetNbinsX()
    for ix in range(tmpnbn):
        awidth = h1.GetBinWidth(ix+1)
        h1.SetBinContent(ix+1, h1.GetBinContent(ix+1)/awidth)
        error = h1.GetBinError(ix+1)
        h1.SetBinError(ix+1, error/awidth)
    return 0




def createCanvasPads():
    c = TCanvas("c", "canvas", 680, 750)
    # Upper histogram plot is pad1
    c.SetBorderSize(0)
    #c.SetLeftMargin(0.1);
    c.SetRightMargin(0.02);
    c.SetTopMargin(0.04);
    c.SetBottomMargin(0.0);

    pad1 = TPad("pad1", "pad1", 0.05, 0.35, 1, 1.0)
    pad1.SetTopMargin(0.04)  # joins upper and lower plot
    pad1.SetBottomMargin(0.0)
    pad1.SetLeftMargin(0.17)
    pad1.SetRightMargin(0.02)
    pad1.SetTickx(1)
    pad1.SetTicky(1)

    #pad1.SetGridy()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0.05, 0.02, 1.0, 0.35)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.5)
    pad2.SetLeftMargin(0.17)
    pad2.SetRightMargin(0.02)
    pad2.SetGridy()
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()
 
    return c, pad1, pad2
 

name = "test.pdf"
canvas = TCanvas("can", "can", 800, 800)
canvas.Print(name+"[") # Add a dummy canvas

def ratioplot():
   for ity in range(0,2):
      for ivar in range(0,5):
         for ipt in range(8):
            #h1 = ReadH1(unfolddata,"Unfold2D"+"/"+"Edd_TUnfold_NoReg_typ_"+str(ity)+"_eta0_"+str(var[ivar])+"_pt"+str(ipt))
            h1 = ReadH1(unfolddata,"Unfold2D"+"/"+"Edd_TUnfold_NoReg_typ_"+str(ity)+"_eta0_"+str(var[ivar])+"_pt"+str(ipt))
            h1.Scale(1/h1.Integral())
            divBybinWidth(h1)
            h1.SetLineColor(kBlack)                
            h1.SetMarkerStyle(9)                
            h1.SetMarkerSize(.9)                
            h1.SetLineWidth(3)                
            c, pad1, pad2 = createCanvasPads()
            pad1.cd()
            sethist(h1,Esvlogx[ivar],Esvlogy[ivar], 22,30, 25,25)
            h1.Draw()
            hc1 = []
            legend = TLegend(0.3 ,.15 ,0.7 ,0.4)
            legend.SetFillStyle(0);
            legend.SetBorderSize(0);
            legend.AddEntry (h1,"Data")
            for index in range(len(inputroot)):
               h2 = ReadH1(inputroot[index],"analyzeBasicPat"+"/"+"gen_typ_"+str(ity)+"_pt"+str(ipt)+"_eta0_"+str(var[ivar]))
               h2.Scale(1/h2.Integral())
               divBybinWidth(h2)
               h2.SetLineColor(color[index])
               h2.SetLineWidth(3)
               legend.AddEntry (h2, modelname[index])
               legend.SetLineWidth (0)
               legend.Draw("same")
               rh2 = createRatio(h1, h2,2.0, 0.1,"MC", "Data")
               sethist(rh2,Esvlogx[ivar],"MC/Data", 28,28,25,25)
               hc1.append(rh2)
               pad2.cd()
               rh2.Draw("same e hist")
               #cantest.cd()
               #rh2.Draw("same e")
               #if index==2:
               #cantest.Print("hello.pdf")
               pad1.cd()
               pad1.SetLogy ( True )
               h2.Draw("same hist")
               
            c.Print(name+"[")
        
   c.Print(name+"]")

            
            
if __name__ == "__main__":
   ratioplot()
   
