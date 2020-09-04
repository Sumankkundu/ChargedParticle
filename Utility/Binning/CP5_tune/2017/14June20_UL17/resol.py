import ROOT
from ROOT import TCanvas, TPaveLabel, TPaveText, TPavesText, TText
from ROOT import TArrow, TLine
from ROOT import gROOT, gBenchmark

nHT = 8
nvar =5
ntype = 2
var=[3,9,15,18,24]

f = ROOT.TFile.Open('PY8_UL2017_120Bin_14june20.root', 'read')
ROOT.gROOT.LoadMacro("anal_jetqcd15_pl.C")
for j in range(ntype):
    for k in range(nvar):
        for l in range(nHT):
            #if (j!=1 or l!=1 or k!=1):
                if (j!=1 or l!=7 or k!=2 ):
                 ROOT.resolution(j,var[k],l,0) 
                 print "   "
#dest.write(ROOT.resolution(0,3,0,0))
f.Close()


