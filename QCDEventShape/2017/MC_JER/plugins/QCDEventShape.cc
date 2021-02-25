// -*- C++ -*-
//
// Package:    Test/QCDEventShape
// Class:      QCDEventShape
// 
/**\class QCDEventShape QCDEventShape.cc Test/QCDEventShape/plugins/QCDEventShape.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tanmay Sarkar
//         Created:  Wed, 01 Jul 2015 10:24:21 GMT
//
//


// system include files

#define DIJETAVE 

////for data
//#define JETENERGY
//#define TRIGGER

// //for Madgraph
//#define LHAPDF
//#define JETRESO
//#define TRACKSYS
//#define TRIGGER

////for Pythia8
#define JETRESO
#define TRIGGER

//For Flat
#define FLAT


////For GenParticle only
//#define GENPART


#include <memory>
#include <map>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TFormula.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <cmath>
#include "TMath.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"
#include "TUnfoldBinningXML.h"
#include "TUnfold.h"
#include "TUnfoldSys.h"

#include "TH2F.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <time.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "CondFormats/DataRecord/interface/L1GtStableParametersRcd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "Test/QCDEventShape/plugins/EventShape_vector.h" 


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <CondFormats/DataRecord/interface/JetResolutionRcd.h>
#include <CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
using namespace edm;
using namespace reco;
using namespace std;
using namespace CLHEP;
using namespace trigger;
using namespace math;
static const int nvar=32;
static const int nhist=10;
static const int typen=2;

static const int nHLTmx=8; 
const char* varname[nvar]={"y3anti", "y3ceanti", "y3cranti", "thrustc", "thrustce", "thrustcr",
                           "minorc", "minorce", "minorcr", "tmass", "tmasse", "tmassr",
                           "hmass", "hmasse", "hmassr", "y3c", "y3ce", "y3cr",
                           "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
                           "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
                           "sphericity", "cparameter"};

const char* vartitle[nvar]={"Anti-Y_{23,C} ", "Anti-Y_{23,E} ", "Anti-Y_{23,R} ",
                            "#tau_{_{#perp} _{   ,C}} ", "#tau_{_{#perp} _{   ,E}} ", "#tau_{_{#perp} _{   ,R}} ",
                            "T_{ m,C} ", "T_{ m,E} ", "T_{ m,R} ",
                            "#rho_{Tot,C} ", "#rho_{Tot,E} ", "#rho_{Tot,R} ",
                            "#rho_{H,C} ", "#rho_{H,E} ", "#rho_{H,R} ",
                            "Y_{23,C} ", "Y_{23,E} ", "Y_{23,R} ",
                            "B_{ T,C} ", "B_{ T,E} ", "B_{ T,R} ",
                            "B_{ W,C} ", "B_{ W,E} ", "B_{ W,R} ",
                            "#rho^{T}_{Tot,C} ", "#rho^{T}_{Tot,E} ", "#rho^{T}_{Tot,R} ",
                            "#rho^{T}_{H,C} ", "#rho^{T}_{H,E} ", "#rho^{T}_{H,R} ",
                            "S_{_{#perp} _{   ,C}}", "C-parameter_{C}"};

//--------------------------For fixed Binning
int nbinsx[nvar]={120, 120, 120, 120, 120, 120,  //Number of Bins 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120};

double endxt[nvar]={8.0, 8.0, 6.0, 8.0, 6.0,  5.0, //Lower Edges
                    4.0, 4.0,  2.0, 7.0, 7.0, 4.0,
                    7.0, 7.0, 4.0, 8.0, 8.0, 6.0,
                    5.0, 5.0, 4.0, 5.0, 5.0, 4.0,
                    10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                    0.0, 0.0};
double startx[nvar]={2., 0.5, 1., 1.0, 0.0, 2.,   //Lower Edges
                     0., 0., 1., 0.2, -1., -1.,
                     0., 0., 1., 1.0, 2., 0.5,
                     0.2, -0.5, -1., 0., -1., -1.,
                     1.0, 0., -0.5, 0., 0., -0.5,
                    -1., -1.};
double endx[nvar]={10.0, 10.0, 10.0, 8.0, 12.0, 12.0,
                   19.0, 19.0, 10.0, 7.0, 10.0, 10.0,
                   13.0, 13.0,  10.0, 8.0, 10.0, 6.0,
                   6.0, 5.0, 5.0, 8.0, 8.0, 8.0,
                   7.0, 8.0, 8.0, 12.0, 12.0, 8.0,
                   0.0, 0.0};
//--------------------------------------------------------------------
/*
//-------------------------No. of Bins For Reco Level for Fixed Bin
//For Jet 
const int rnmxbins=26;
int rnbinsx0[nvar]={0,0,0,26,0,0,
                   0,0,0,20,0,0,
                   0,0,0,12,0,0,
                   20,0,0,0,0,0,
                   22,0,0,0,0,0,0,0};

//For charge Particles
int rnbinsx1[nvar]={0,0,0,20,0,0,
                   0,0,0,16,0,0,
                   0,0,0,8,0,0,
                   22,0,0,0,0,0,
                   14,0,0,0,0,0,0,0};
*/
//-----------------------------------For Reco Level
const int rnmxbins=32;  //Maximum Bins in bellow array
int rnbinsx0[nvar]={0,0,0,32,0,0,0,0,0,
                  32,0,0,0,0,0,14,0,
                  0,24,0,0,0,0,0,26,
                  0,0,0,0,0,0,0};

//For charge Particles
int rnbinsx1[nvar]= {0,0,0,18,0,0,0,0,
                  0,18,0,0,0,0,0,8,
                  0,0,16,0,0,0,0,0,
                  14,0,0,0,0,0,0,0};

//-----------------------------------For Gen Level
const int nmxbins=16;  //Maximum Bins in bellow array
int nbinsx0[nvar]={0,0,0,16,0,0,0,0,
                  0,16,0,0,0,0,0,7,
                  0,0,12,0,0,0,0,0,
                  13,0,0,0,0,0,0,0};

 int  nbinsx1[nvar]={0,0,0,9,0,0,0,0,
                  0,9,0,0,0,0,0,4,
                  0,0,8,0,0,0,0,0,
                  7,0,0,0,0,0,0,0};


 //Reco level  /Jet
double rbinrngs0[nvar][rnmxbins+1] ={{},{},{},
                                      {-6.71, -6.4, -6.11, -5.83, -5.56, -5.3, -5.05, -4.8, -4.56, -4.33, -4.11, -3.89, -3.68, -3.48, -3.28, -3.09, -2.91, -2.74, -2.58, -2.43, -2.29, -2.15, -2.02, -1.9, -1.79, -1.69, -1.6, -1.51, -1.43, -1.36, -1.3, -1.24, -1.19},//32
                                      {},{},{},{},{},
                                      {-5.75, -5.45, -5.18, -4.92, -4.68, -4.45, -4.22, -4, -3.79, -3.58, -3.37, -3.17, -2.97, -2.77, -2.57, -2.38, -2.19, -2.01, -1.84, -1.67, -1.51, -1.36, -1.22, -1.09, -0.97, -0.86, -0.76, -0.67, -0.59, -0.52, -0.45, -0.39, -0.34}, //32
                                      {},{},{},{},{},
                                      {-6.71, -6.18, -5.67, -5.18, -4.73, -4.31, -3.92, -3.56, -3.22, -2.9, -2.59, -2.3, -2.01, -1.73, -1.44}, //14
                                      {},{},
                                      {-3.58, -3.34, -3.1, -2.87, -2.65, -2.44, -2.24, -2.05, -1.88, -1.72, -1.57, -1.43, -1.3, -1.18, -1.06, -0.95, -0.85, -0.76, -0.67, -0.59, -0.51, -0.43, -0.36, -0.29, -0.23}, //24
                                      {},{},{},{},{},
                                      {-5.89, -5.5, -5.17, -4.88, -4.61, -4.36, -4.12, -3.89, -3.67, -3.45, -3.24, -3.03, -2.83, -2.64, -2.45, -2.27, -2.1, -1.94, -1.79, -1.66, -1.54, -1.43, -1.34, -1.26, -1.19, -1.13, -1.08},//26
                                     {},{},{},{},{},{},{}};

double rbinrngs1[nvar][rnmxbins+1] = {{},{},{},
                                        {-5.66, -5.36, -5.06, -4.76, -4.46, -4.16, -3.86, -3.56, -3.27, -2.98, -2.7, -2.44, -2.19, -1.96, -1.75, -1.56, -1.4, -1.26, -1.14}, //18
                                        {},{},{},{},{},
                                        {-5.75, -5.3, -4.87, -4.45, -4.04, -3.65, -3.27, -2.9, -2.54, -2.2, -1.87, -1.56, -1.28, -1.02, -0.79, -0.59, -0.41, -0.26, -0.13 }, //18
                                        {},{},{},{},{},
                                        {-6.71, -6.04, -5.36, -4.69, -4.04, -3.41, -2.79, -2.15, -1.58 }, //8
                                        {},{},
                                        {-3.77, -3.49, -3.21, -2.93, -2.66, -2.39, -2.13, -1.87, -1.62, -1.39, -1.17, -0.97, -0.79, -0.62, -0.47, -0.34, -0.23 }, //16
                                        {},{},{},{},{},
                                        {-5.89, -5.4, -4.95, -4.52, -4.11, -3.72, -3.35, -2.99, -2.65, -2.33, -2.03, -1.75, -1.49, -1.25, -1.04},
                                        {},{},{},{},{},{},{}};//14

//Gen Level
double binrngs0[nvar][nmxbins+1] = {{},{},{},
                                      {-6.71,-6.11,-5.56,-5.05,-4.56,-4.11,-3.68,-3.28,-2.91,-2.58,-2.29,-2.02,-1.79,-1.6,-1.43,-1.3,-1.19}, //16
                                      {},{},{},{},{},
                                      {-5.75,-5.18,-4.68,-4.22,-3.79,-3.37,-2.97,-2.57,-2.19,-1.84,-1.51,-1.22,-0.97,-0.76,-0.59,-0.45,-0.34}, //16
                                      {},{},{},{},{},
                                      {-6.71,-5.67,-4.73,-3.92,-3.22,-2.59,-2.01,-1.44},
                                      {},{},
                                      {-3.58,-3.1,-2.65,-2.24,-1.88,-1.57,-1.3,-1.06,-0.85,-0.67,-0.51,-0.36,-0.23},//12
                                      {},{},{},{},{},
                                      {-5.89,-5.17,-4.61,-4.12,-3.67,-3.24,-2.83,-2.45,-2.1,-1.79,-1.54,-1.34,-1.19,-1.08},
                                      {},{},{},{},{},{},{}};

double binrngs1[nvar][nmxbins+1]={{},{},{},
                                    {-5.66,-5.06,-4.46,-3.86,-3.27,-2.7,-2.19,-1.75,-1.4,-1.14},
                                    {},{},{},{},{},
                                    {-5.75,-4.87,-4.04,-3.27,-2.54,-1.87,-1.28,-0.79,-0.41,-0.13},
                                    {},{},{},{},{},
                                    {-6.71,-5.36,-4.04,-2.79,-1.58},
                                    {},{},
                                    {-3.77,-3.21,-2.66,-2.13,-1.62,-1.17,-0.79,-0.47,-0.23},
                                    {},{},{},{},{},
                                    {-5.89,-4.95,-4.11,-3.35,-2.65,-2.03,-1.49,-1.04},
                                    {},{},{},{},{},{},{}};



//----------------------------------------HT2 Binning For 2D unfold 
double recohtbins[nHLTmx+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0};

//---------------------------------------------------------------------------------------------------------
const int nusedvar=5;
double usedvars[nusedvar]={3, 9, 15, 18, 24};

int isItUsed(int ival) {
	for (int ij=0; ij<nusedvar; ij++) {
		if (ival==usedvars[ij]) {return 1;}
	}
	return 0;
} 


const int npileupmx=99; //49;
double rat_pileup[nHLTmx][npileupmx]={{0}};
//clock_t t1,t2;
//UL PU

double mcpileup[npileupmx] ={ 1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
    0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
    0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
    0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
    0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
    0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
    0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
    0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
    0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
    0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
    0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
    0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
    0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
    3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
    1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
//HLT Path PileUP

//20APril Ultra Legacy true
double datpileup[nHLTmx][npileupmx] ={{1.53675e-05, 4.23286e-05, 0.000116697, 0.000208638, 0.000142335, 0.000182092, 0.000180895, 0.000155309, 0.000289895, 0.00125488, 0.00271058, 0.00492075, 0.00633856, 0.00618827, 0.00658225, 0.00866802, 0.0132994, 0.0207452, 0.0306966, 0.042707, 0.054913, 0.0633769, 0.064582, 0.0587072, 0.0494767, 0.0410557, 0.0353318, 0.0320279, 0.0301134, 0.0287383, 0.0275101, 0.0263219, 0.025095, 0.0237351, 0.0222264, 0.020644, 0.0190731, 0.017558, 0.016128, 0.0148406, 0.0137844, 0.0130433, 0.0126617, 0.0126278, 0.0128688, 0.0132535, 0.013604, 0.0137256, 0.0134542, 0.012701, 0.0114786, 0.00990479, 0.00815462, 0.00641259, 0.00482923, 0.00349585, 0.00244328, 0.00165653, 0.00109469, 0.000708293, 0.000450588, 0.000282908, 0.000175936, 0.000108744, 6.70361e-05, 4.13631e-05, 2.56357e-05, 1.60105e-05, 1.01022e-05, 6.45021e-06, 4.16938e-06, 2.72647e-06, 1.80078e-06, 1.19867e-06, 8.02187e-07, 5.38496e-07, 3.61838e-07, 2.42943e-07, 1.62748e-07, 1.0865e-07, 7.22147e-08, 4.77475e-08, 3.1384e-08, 2.04947e-08, 1.329e-08, 8.55381e-09, 5.46224e-09, 3.4594e-09, 2.17225e-09, 1.352e-09, 8.3386e-10, 5.09523e-10, 3.08394e-10, 1.84863e-10, 1.09732e-10, 6.44909e-11, 3.75233e-11, 2.16122e-11, 1.23213e-11},
{2.09998e-05, 5.81646e-05, 0.000160839, 0.000286491, 0.00018581, 0.000229088, 0.000226163, 0.000190714, 0.000376732, 0.00172036, 0.0037652, 0.00683484, 0.00847828, 0.00754201, 0.00760162, 0.00998974, 0.0149877, 0.0217875, 0.0287254, 0.0352424, 0.0413189, 0.0459001, 0.0475604, 0.0461527, 0.0429618, 0.0395778, 0.036789, 0.0346623, 0.0329776, 0.031479, 0.0300664, 0.0287391, 0.0274254, 0.0260041, 0.0244393, 0.0228039, 0.0211878, 0.0196366, 0.0181764, 0.0168634, 0.015791, 0.0150545, 0.0147118, 0.0147593, 0.015122, 0.0156514, 0.0161383, 0.0163486, 0.0160821, 0.0152277, 0.0137974, 0.0119314, 0.0098409, 0.00775013, 0.00584327, 0.00423333, 0.00295992, 0.00200666, 0.00132526, 0.000856429, 0.000543802, 0.000340559, 0.000211103, 0.000129975, 7.97684e-05, 4.89777e-05, 3.01962e-05, 1.87576e-05, 1.1773e-05, 7.47962e-06, 4.81303e-06, 3.13509e-06, 2.06388e-06, 1.37012e-06, 9.14963e-07, 6.13151e-07, 4.1144e-07, 2.7594e-07, 1.84684e-07, 1.232e-07, 8.18311e-08, 5.40743e-08, 3.55241e-08, 2.31874e-08, 1.50297e-08, 9.66978e-09, 6.17265e-09, 3.90803e-09, 2.45321e-09, 1.52643e-09, 9.41188e-10, 5.74961e-10, 3.4792e-10, 2.0851e-10, 1.23742e-10, 7.27107e-11, 4.2298e-11, 2.43579e-11, 1.38842e-11},
{2.35835e-05, 6.54134e-05, 0.000180301, 0.000320928, 0.000205544, 0.000250447, 0.000250517, 0.000210823, 0.000417088, 0.00192635, 0.00422647, 0.0076718, 0.00945682, 0.00827455, 0.00821479, 0.0106622, 0.0156455, 0.0218821, 0.027118, 0.0306929, 0.0334707, 0.035936, 0.0377992, 0.0387448, 0.0387444, 0.0380018, 0.0367266, 0.0351854, 0.0336028, 0.0320463, 0.0305733, 0.0292453, 0.0280075, 0.0267245, 0.0253176, 0.0238007, 0.0222186, 0.0206137, 0.019053, 0.017653, 0.0165595, 0.0158982, 0.0157354, 0.0160584, 0.0167693, 0.0176839, 0.0185443, 0.0190606, 0.01898, 0.0181552, 0.0165886, 0.0144445, 0.0119806, 0.00947694, 0.00716879, 0.00520516, 0.00364357, 0.00247026, 0.00162965, 0.00105071, 0.000664744, 0.000414189, 0.000255033, 0.000155707, 9.45924e-05, 5.73942e-05, 3.49191e-05, 2.13875e-05, 1.32343e-05, 8.29563e-06, 5.27501e-06, 3.40264e-06, 2.22358e-06, 1.46875e-06, 9.77936e-07, 6.54513e-07, 4.39182e-07, 2.94795e-07, 1.97582e-07, 1.32033e-07, 8.78643e-08, 5.81736e-08, 3.82901e-08, 2.50388e-08, 1.62581e-08, 1.04773e-08, 6.6984e-09, 4.247e-09, 2.6696e-09, 1.66318e-09, 1.02674e-09, 6.27936e-10, 3.80387e-10, 2.28202e-10, 1.35563e-10, 7.97317e-11, 4.64246e-11, 2.67578e-11, 1.52651e-11},
{2.94733e-05, 8.13293e-05, 0.000220365, 0.000392672, 0.000236212, 0.000279115, 0.000270768, 0.000209453, 0.000465055, 0.00234836, 0.00520614, 0.00944458, 0.0115626, 0.00990584, 0.0093488, 0.011373, 0.0158949, 0.0216814, 0.0265205, 0.0296053, 0.0318658, 0.0340989, 0.0360623, 0.0372195, 0.0373273, 0.0365257, 0.0351071, 0.0334582, 0.0318641, 0.0303735, 0.0289862, 0.0277199, 0.0265125, 0.0252664, 0.0239896, 0.0228013, 0.0218002, 0.0209744, 0.0202384, 0.0195305, 0.0188757, 0.0183729, 0.0181361, 0.0182298, 0.0186249, 0.019185, 0.0196831, 0.019852, 0.0194589, 0.018375, 0.0166156, 0.0143478, 0.0118222, 0.00930409, 0.00701102, 0.00507632, 0.00354639, 0.0024012, 0.00158275, 0.00101992, 0.000645035, 0.000401783, 0.000247304, 0.000150912, 9.16116e-05, 5.5528e-05, 3.37368e-05, 2.06267e-05, 1.27361e-05, 7.96353e-06, 5.05005e-06, 3.24823e-06, 2.11657e-06, 1.39414e-06, 9.25781e-07, 6.18057e-07, 4.13758e-07, 2.77137e-07, 1.85382e-07, 1.23657e-07, 8.21547e-08, 5.4311e-08, 3.5698e-08, 2.3314e-08, 1.51203e-08, 9.73349e-09, 6.21665e-09, 3.9379e-09, 2.47315e-09, 1.53955e-09, 9.49692e-10, 5.80397e-10, 3.51349e-10, 2.10645e-10, 1.25055e-10, 7.35079e-11, 4.27762e-11, 2.46413e-11, 1.40502e-11},
{2.27285e-05, 6.28966e-05, 0.000173563, 0.000309334, 0.00019925, 0.000245037, 0.000241188, 0.000201572, 0.000402568, 0.00185697, 0.00406832, 0.00738526, 0.00916867, 0.00816646, 0.00816752, 0.0105166, 0.0153225, 0.0213648, 0.0264376, 0.0298462, 0.0324307, 0.0347949, 0.0367969, 0.0381212, 0.0386142, 0.0383485, 0.0374593, 0.0361929, 0.034786, 0.0333258, 0.0318747, 0.0304856, 0.0290943, 0.0275767, 0.025903, 0.0241592, 0.0224459, 0.0208117, 0.019284, 0.0179258, 0.0168454, 0.0161577, 0.0159359, 0.0161807, 0.0168032, 0.0176209, 0.0183747, 0.0187773, 0.0185859, 0.0176723, 0.0160574, 0.0139144, 0.011498, 0.00907326, 0.00685589, 0.00497797, 0.00348701, 0.00236626, 0.00156193, 0.00100688, 0.000636315, 0.000395703, 0.000243057, 0.000148046, 8.97933e-05, 5.44685e-05, 3.31896e-05, 2.03973e-05, 1.26845e-05, 7.9984e-06, 5.11723e-06, 3.31919e-06, 2.17862e-06, 1.44338e-06, 9.62578e-07, 6.44463e-07, 4.3217e-07, 2.89705e-07, 1.93825e-07, 1.29258e-07, 8.58331e-08, 5.67061e-08, 3.72455e-08, 2.43066e-08, 1.57525e-08, 1.01333e-08, 6.46759e-09, 4.09421e-09, 2.56976e-09, 1.59876e-09, 9.85681e-10, 6.02078e-10, 3.64293e-10, 2.18301e-10, 1.29541e-10, 7.61114e-11, 4.42725e-11, 2.54928e-11, 1.453e-11},
{2.21934e-05, 6.14397e-05, 0.000170057, 0.000303246, 0.000197536, 0.0002444, 0.000241052, 0.000203391, 0.000400335, 0.00182211, 0.00398326, 0.00723306, 0.00903071, 0.00816199, 0.00825969, 0.0106703, 0.015557, 0.021685, 0.0268075, 0.0302318, 0.0328348, 0.0352344, 0.0372748, 0.0386301, 0.0391408, 0.0388729, 0.0379606, 0.0366589, 0.0352149, 0.0337199, 0.0322385, 0.0308263, 0.0294206, 0.0278944, 0.0262118, 0.0244535, 0.0227168, 0.0210497, 0.0194795, 0.0180684, 0.0169217, 0.0161491, 0.0158185, 0.0159281, 0.0163937, 0.0170467, 0.0176497, 0.0179386, 0.0176884, 0.0167767, 0.0152187, 0.0131722, 0.0108731, 0.00857026, 0.00646747, 0.00468983, 0.00328168, 0.00222587, 0.00147, 0.000949319, 0.000601914, 0.00037612, 0.000232471, 0.000142641, 8.72146e-05, 5.33446e-05, 3.27671e-05, 2.02861e-05, 1.26953e-05, 8.04623e-06, 5.16788e-06, 3.3614e-06, 2.21049e-06, 1.46627e-06, 9.78561e-07, 6.55446e-07, 4.39639e-07, 2.94747e-07, 1.97208e-07, 1.31517e-07, 8.7332e-08, 5.76951e-08, 3.78942e-08, 2.47292e-08, 1.6026e-08, 1.03089e-08, 6.57951e-09, 4.16496e-09, 2.61411e-09, 1.62632e-09, 1.00265e-09, 6.12428e-10, 3.70548e-10, 2.22046e-10, 1.31761e-10, 7.74142e-11, 4.50295e-11, 2.59283e-11, 1.4778e-11},
{2.39377e-05, 6.62419e-05, 0.000182187, 0.000324618, 0.000206085, 0.000250018, 0.000248141, 0.000206365, 0.000415288, 0.00194489, 0.00427255, 0.00776125, 0.00959895, 0.00845094, 0.00833784, 0.0106032, 0.0152607, 0.021134, 0.0261417, 0.0295838, 0.0322083, 0.0346006, 0.0366529, 0.0380563, 0.0386329, 0.0384361, 0.0375983, 0.0363664, 0.0349795, 0.0335317, 0.0320955, 0.0307324, 0.0293802, 0.027909, 0.026279, 0.0245667, 0.0228674, 0.0212302, 0.0196865, 0.0183016, 0.0171808, 0.0164312, 0.016118, 0.0162384, 0.0167085, 0.0173608, 0.0179586, 0.0182383, 0.0179748, 0.0170444, 0.0154616, 0.0133844, 0.0110505, 0.00871168, 0.00657494, 0.00476799, 0.00333644, 0.00226314, 0.00149486, 0.000965675, 0.000612562, 0.000382973, 0.000236812, 0.000145332, 8.88344e-05, 5.42849e-05, 3.32886e-05, 2.05587e-05, 1.28264e-05, 8.10091e-06, 5.18406e-06, 3.36003e-06, 2.20249e-06, 1.45689e-06, 9.70042e-07, 6.48505e-07, 4.34317e-07, 2.9082e-07, 1.94386e-07, 1.29527e-07, 8.59505e-08, 5.67485e-08, 3.72529e-08, 2.42994e-08, 1.57408e-08, 1.01216e-08, 6.45774e-09, 4.08656e-09, 2.56413e-09, 1.59478e-09, 9.82948e-10, 6.00249e-10, 3.63094e-10, 2.17531e-10, 1.29054e-10, 7.58084e-11, 4.40869e-11, 2.53807e-11, 1.44632e-11},
{1.80294e-05, 6.04463e-05, 0.000172336, 0.000310608, 0.000202004, 0.000249998, 0.000245965, 0.000206972, 0.000410215, 0.0018777, 0.00411048, 0.00746252, 0.00925659, 0.00822804, 0.00826162, 0.0107372, 0.0157553, 0.0220336, 0.0272888, 0.030851, 0.0336004, 0.0361119, 0.038222, 0.0396231, 0.0401661, 0.0399152, 0.0389974, 0.0376669, 0.0361757, 0.0346178, 0.0330587, 0.0315524, 0.0300297, 0.0283585, 0.0265074, 0.0245749, 0.0226808, 0.0208881, 0.0192247, 0.017741, 0.0165249, 0.0156695, 0.015231, 0.0152047, 0.0155138, 0.0160077, 0.0164738, 0.0166726, 0.016397, 0.0155296, 0.014078, 0.0121816, 0.0100536, 0.00792233, 0.00597609, 0.00433119, 0.00302905, 0.00205369, 0.00135617, 0.000876136, 0.000556013, 0.00034792, 0.000215417, 0.000132427, 8.11144e-05, 4.96854e-05, 3.05477e-05, 1.8918e-05, 1.18359e-05, 7.49598e-06, 4.80957e-06, 3.12491e-06, 2.05289e-06, 1.36062e-06, 9.0753e-07, 6.07657e-07, 4.07523e-07, 2.73217e-07, 1.82823e-07, 1.21946e-07, 8.09947e-08, 5.35216e-08, 3.51618e-08, 2.29518e-08, 1.48776e-08, 9.5724e-09, 6.11078e-09, 3.86903e-09, 2.42883e-09, 1.51133e-09, 9.31914e-10, 5.69318e-10, 3.44518e-10, 2.06478e-10, 1.22541e-10, 7.2007e-11, 4.18899e-11, 2.41236e-11, 1.37511e-11}};



static const int nsrc = 27;   // Change form 26 as for 2015 data .  See JEC for 2017 94X
const char* srcnames[nsrc] =
  {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF","RelativePtBB", "RelativePtEC1", "RelativePtEC2","RelativePtHF","RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"};


double intlumi[nHLTmx]={1., 1, 1, 1, 1, 1,1,1};
double lumiwt[nHLTmx]={1., 1, 1, 1, 1, 1,1,1};
//unsigned int l1trg[4], hlttr[8], tetrg[2];
unsigned int mypow_2[32];

//std::ofstream myfile;
//  myfile.open("txt.log");


//const bool m_trigeff = true;
const int njetptmn=nHLTmx; // 8;
const int njetptbin=120;

#ifdef DIJETAVE
const char* jethlt_name[nHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};
//double leadingPtThreshold[njetptmn+1] ={90.0, 120.0, 180.0, 250.0, 320.0, 400.0, 480.0, 600.0, 2000.0};

double leadingPtThreshold[njetptmn+1] ={83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger

//double compres[njetptmn] = {1630, 5320, 62.1, 38.9, 27.0, 4.33, 1.23, 1.0};
//double compres[njetptmn] = {1630, 5320, 62.1, 38.9, 27.0, 4.33, 1.23, 1.0};

const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};

//#else

#endif

#ifdef DIJETAVE
double jethlt_thr[nHLTmx]={60,80,140,200,260,320,400,500};
//#else

#endif
double prescl[nHLTmx];

#ifdef TRACKSYS
const int ntype=3;
#else
const int ntype=2;
#endif

const int njetetamn=1; // GMA 4;
#ifdef  LHAPDF
       const int nnnmx=101;
        double pdfwt[nnnmx];
  TH1F* h_genevtvarpdf[ntype][njetptmn][njetetamn][nvar][nnnmx];
  TH1* h_genevtvarpdf_2D[ntype][njetetamn][nvar][nnnmx];
#endif

#ifdef  JETENERGY
        //const int nsrc = 26;
        const int njecmx=2*nsrc+1;
  TH1F* h_recoevtvarjec[ntype][njetptmn][njetetamn][nvar][njecmx];
  TH1* h_recoevtvarjec_2D[ntype][njetetamn][nvar][njecmx]; //For 2D
#elif defined(JETRESO)
        const int njecmx = 3;
  TH1F* h_recoevtvarres[ntype][njetptmn][njetetamn][nvar][njecmx];
  TH1* h_recoevtvarres_2D[ntype][njetetamn][nvar][njecmx]; //For 2D //No use
#else
  const int njecmx=1;
#endif


//#ifdef  JETRESO
//  const int nGenReso = 3;
//  TH1F* h_genevtvarres[ntype][njetptmn][njetetamn][nvar][nGenReso];
//#else
  const int nGenReso=1;
  //const int nGenReso=njecmx;
//#endif

//int trgbit[nHLTmx]={10,11,12,13,14,16};
//double trgpas[nHLTmx+1]={0,0,0,0,0,0,0,0,0};

//const int njetetamn=3;
double etarange[njetetamn] ={2.4}; //{3.0, 2.4, 1.8, 1.3};
double resetarange[njetetamn+4] ={0, 0.5, 1.0, 1.5}; //, 2.0, 2.5, 3.0, 3.5};
double par0[njetetamn+4]={1.02, 1.02, 1.022, 1.017, 0.98}; //, 0.9327};
double par1[njetetamn+4]={7.3e-6, -7.3e-6, -5.66e-6, -9.9e-6, 1.41e-4}; //, 4.6e-4};
double par2[njetetamn+4]={-8.2e-9, -8.2e-9, -3.58e-9, -4.18e-9, -6.104e-8}; //, -4.041e-7};
double particlept[4]={0.0, 0.25, 0.50, 1.00};

//const int ntype=4;
//const char* typname[ntype]={"Jets", "All particle", "All particle: P_{T}>0.25", "All particle: P_{T}>0.50", "All particle: P_{T}>1"};

#ifdef TRACKSYS
const char* typname[ntype]={"Jets", "Charged Particles"};
#else
const char* typname[ntype]={"Jets", "Charged Particles"};
#endif
static const int njetmx =30;

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double Phi_0_2pi(double x) {
  while (x >= 2*M_PI) x -= 2*M_PI;
  while (x <     0.)  x += 2*M_PI;
  return x;
}

double Phi_mpi_pi(double x) {
  while (x >= M_PI) x -= 2*M_PI;
  while (x < -M_PI) x += 2*M_PI;
  return x;
}



double dPhi(double phi1,double phi2){
  phi1=Phi_0_2pi(phi1);
  phi2=Phi_0_2pi(phi2);
  return Phi_mpi_pi(phi1-phi2);
}

 int sbitx(unsigned ival, int ibit) {
      unsigned den = mypow_2[ibit]; // unsigned(pow(2., double(ibit)));
      int isel = unsigned(ival/den)%2;
 //  int isel = unsigned(ival/den);
       //cout <<"iv "<< ival<<" "<<ibit<<" "<<den<<" "<<ival/den<<" "<<unsigned(ival/den)<<" "<<isel<<endl;

      return isel;
    }

double respfun(double a, double b, double c, double x){
  double func=a+b*x+c*x*x;
  return func;
}


struct triggervar{
  HepLorentzVector trg4v;
  bool		  both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
};

//
// class declaration
//

class QCDEventShape : public edm::EDAnalyzer {
   public:
      explicit QCDEventShape(const edm::ParameterSet&);
      ~QCDEventShape();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 // int sbitx(unsigned ival, int ibit);

  bool isHistFill;
  bool isTrigger;
  bool isRECO[ntype][njetetamn];
  bool isRECO_JER[njecmx][ntype][njetetamn];//for JER
  bool isMC;
  //  bool isParticle; //Do we want particle level informations, other than jets ?
  //  bool isGenParticle; //Do we want Simulated particle level informations, other than jets ?
  bool isReconstruct; // otherwise Only generator level informations  
  //  bool isPartQCD; //For tracker variables, recosntruct QCD EVT variables
  bool isJetQCD;  //For Jet variables, recosntruct QCD EVT variables
  bool isGenJET; // Genjet information or note (for herwig/alpgen, donot store this ?)
  //  double trackPtThreshold; //Threshold of track Pt to store it in root file, -ve implies don't store

  //  double etarange; //Eta range of all jets
  double ptthreshold; //Pt threshold of JEC jets
  double leadingPtthreshold; //Pt threshold of JEC leading jet
  bool   isOtherAlgo; // store Kt4 and ak7 variables or not
  double weight=1; //weight for histogramme fit
  double weight2=1;

  std::string m_resolutions_file;
  std::string scalefile;

  std::string theHLTTag;
//  unsigned int mypow_2[32];
  int nevt;

  std::string theRootFileName;
  //TFile* //theFile;
 // TTree* //T1;

  //ifstream myfile ("example.txt");
  //std::ofstream myfile;
  //myfile.open("txt.log");
  TDirectoryFile *TUnfoldBinng2D =new TDirectoryFile("analyzeBasicPat2D","TUnfoldBinning 1D Historgams");

  //TH1F* h_recoevtvar[10][ntype][njetptmn][njetetamn][nvar];
  TH1F* h_recoevtvar[ntype][njetptmn][njetetamn][nvar];
  TH1F* h_recoevtfake[ntype][njetptmn][njetetamn][nvar];  //For fake
  TH1F* h_genevtmiss[ntype][njetptmn][njetetamn][nvar];  //For miss

  TH1F* h_genevtvar[ntype][njetptmn][njetetamn][nvar];  //For Gen
  TH1F* h_genevtvar2[ntype][njetptmn][njetetamn][nvar];
  
  TH2F* h_2devtvar[ntype][njetptmn][njetetamn][nvar];  //For RM
  TH2F* h_2ht;

  TH1F* vec_anglex[nhist];

  //static const int njetmx =30;
  int npfjets;
  int pfjetmul[njetmx];
  float pfjetpx[njetmx], pfjetpy[njetmx], pfjetpz[njetmx], pfjeten[njetmx],  pfjetenuc[njetmx], neuemf[njetmx], neuhad[njetmx];
  float pfjetenscl[njetmx], pfjetensmr[njetmx];

  float jetpt, jeteta, jetphi; 
  int nallpf, ncharged;
  float thphi[nhist], thrust[nhist], anglex[nhist];
  float jtthan;
  int irunhlt, l1pres[nHLTmx],  hltpres[nHLTmx], compres[nHLTmx]; 
  static const int nprimx=150;
  int nprim, ntkpm[nprimx];
  //  float  primdx[nprimx], primdy[nprimx], primdz[nprimx], 
  float primpr[nprimx];
  int irun, ilumi, ibrnc;
  unsigned int ievt;
  float inslumi;
  int nsicls, ntottrk;
//#ifdef FLAT 
//  bool isFlat=1;
//#else 
    bool isFlat=0;
//#endif

    float defweight=1.0, weighttrg=1., qlow=-10., qhigh=100000.;
  //=============****=========================

//----------------------------------------------------------------TunfoldBinning -----------------------------
  //----------------------------------2D Bining using TUnfoldBinning
   TUnfoldBinning *binsRec2D[ntype][njetetamn][nvar];
   TUnfoldBinning *RecoBinning2D[ntype][njetetamn][nvar];
   TUnfoldBinning *binsGen2D[ntype][njetetamn][nvar];
   TUnfoldBinning *GenBinning2D[ntype][njetetamn][nvar];

   TH1* h_recovar_2D[ntype][njetetamn][nvar]; //Reco
   TH1* h_recofake_2D[ntype][njetetamn][nvar];//For fake
  // TH1* h_recofakeOutE_2D[type][njetetamn][nvar];//For fake
  // TH1* h_recofakeOutHT_2D[type][njetetamn][nvar];//For fake


   TH1* h_genvar_2D[ntype][njetetamn][nvar]; //Gen
   TH1* h_genmiss_2D[ntype][njetetamn][nvar]; //For Miss
 //  TH1* h_genmissOutE_2D[type][njetetamn][nvar]; //For Miss
 //  TH1* h_genmissOutHT_2D[type][njetetamn][nvar]; //For Miss

   TH2* RM_2D[ntype][njetetamn][nvar];
   
   TH2* RM_JER_2D[ntype][njetetamn][nvar][njecmx];
   TH1* h_gen_JER_miss_2D[ntype][njetetamn][nvar][njecmx]; //For Miss
   TH1* h_reco_JER_fake_2D[ntype][njetetamn][nvar][njecmx];//For fake



  //=============****=========================
  //TH1F* recojt_hist;
//  TH1F* recojt_pt[njetetamn][nHLTmx];
  TH1F* recojt_pt[njetetamn];
  TH1F* recojt_eta;
  TH1F* recojt_phi;

  TH1F* recojtallave_pt[njetetamn];
  TH1F* recojtallavewt1_pt[njetetamn];

  TH1F* recojtave_pt[njetetamn][nHLTmx];
  TH1F* recojtavewt1_pt[njetetamn][nHLTmx];
  TH1F* recojt1_pt[njetetamn];
  TH1F* recojt1_eta;
  TH1F* recojt1_phi;

  TH1F* recojt2_pt[njetetamn];
  TH1F* recojt2_eta;
  TH1F* recojt2_phi;

  TH1F* recojt3_pt[njetetamn];
  TH1F* recojt3_eta;
  TH1F* recojt3_phi;

  TH1F* recoht2_pt[njetetamn];




  TH1F* hjetdpt[njetetamn];
  TH1F* hjetdphi[njetetamn];
  TH1F* hjetptbypl[njetetamn];
  TH1F* hjetpt2bypt1[njetetamn];
  TH1F* hjetpt3bypt2[njetetamn];
  // TH1F* recochg_hist;
  TH1F* recochg_pt;
  TH1F* recochg_eta;
  TH1F* recochg_phi;

  TH1F* recochg1_pt;
  TH1F* recochg1_eta;
  TH1F* recochg1_phi;

  TH1F* recochg2_pt;
  TH1F* recochg2_eta;
  TH1F* recochg2_phi;

  TH1F* recochg3_pt;
  TH1F* recochg3_eta;
  TH1F* recochg3_phi;

  //===============****==============================
  //  TH1F* genjt_hist;
  TH1F* genjt_pt[njetetamn];
  TH1F* genjt_eta;
  TH1F* genjt_phi;
  TH1F* genjtallave_pt[njetetamn];


  TH1F* genjt1_pt[njetetamn];
  TH1F* genjt1_eta;
  TH1F* genjt1_phi;

   TH1F* genjt2_pt[njetetamn];
  TH1F* genjt2_eta;
  TH1F* genjt2_phi;

  TH1F* genjt3_pt[njetetamn];
  TH1F* genjt3_eta;
  TH1F* genjt3_phi;

  TH1F* genjetdpt[njetetamn];
  TH1F* genjetdphi[njetetamn];
  TH1F* genjetptbypl[njetetamn];
  TH1F* genjetpt2bypt1[njetetamn];
  TH1F* genjetpt3bypt2[njetetamn];

  // TH1F* genchg_hist;
  TH1F* genchg_pt;
  TH1F* genchg_eta;
  TH1F* genchg_phi;

  TH1F* genchg1_pt;
  TH1F* genchg1_eta;
  TH1F* genchg1_phi;

  TH1F* genchg2_pt;
  TH1F* genchg2_eta;
  TH1F* genchg2_phi;

  TH1F* genchg3_pt;
  TH1F* genchg3_eta;
  TH1F* genchg3_phi;
/*
  // TH1F* genneu_hist;
  TH1F* genneu_pt;
  TH1F* genneu_eta;
  TH1F* genneu_phi;

  TH1F* genjt_oth_pt[njetetamn];
  TH1F* genjt_oth_eta;
  TH1F* genjt_oth_phi;

  //  TH1F* genchg_oth_hist;
  TH1F* genchg_oth_pt;
  TH1F* genchg_oth_eta;
  TH1F* genchg_oth_phi;

  //  TH1F* genneu_oth_hist;
  TH1F* genneu_oth_pt;
  TH1F* genneu_oth_eta;
  TH1F* genneu_oth_phi;
*/
  //Response hist
  TH2F* resp_jet[njetetamn+2];
  TH1F* resp_jet1[njetetamn+2];

  TH1F* prim_hist[nHLTmx+1];
  TH1F* prim_sel[nHLTmx+1];

  TH1F* prim_hist_rewt[nHLTmx+1];
  TH1F* prim_sel_rewt[nHLTmx+1];

  TH2F* prim_correl;

  TH1F* prim_alltrk[2];
  TH1F* prim_seltrk[2];
  TH1F* prim_goodtrk[2];
  TH1F* prim_dx[2];
  TH1F* prim_dy[2];
  TH2F* prim_dxy[2];
  TH1F* prim_dz[2];  
  TH1F* prim_prob[2];

  TH1F* h_jetpt[nHLTmx][njetetamn];
  TH1F* h_jeteta[nHLTmx];
  TH1F* h_jetphi[nHLTmx][njetetamn];
  TH1F* h_njets[njetetamn];
  TH1F* h_nchg[njetetamn];

  TH1F* gen_njets[njetetamn];




  TH1F* trgjet_angle[nHLTmx][2];
  TH2F* trgjet_2dangle[nHLTmx][2];
  TH1F* trgjet_pt[nHLTmx][2];
  TH1F* trgjet_eta[nHLTmx][2];
  TH1F* trgjet_phi[nHLTmx][2];
  TH1F* prbjet_pt[nHLTmx][2];
  TH1F* prbjet_eta[nHLTmx][2];
  TH1F* prbjet_phi[nHLTmx][2];


  //Dijet trigger efficiency
  TH1F* hlt_dijettag[nHLTmx][njetetamn];
  TH1F* hlt_dijetprob[nHLTmx][njetetamn];

  //Trigger Normal case

  TH1F* counthist;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GenEventInfoProduct> generator1_;
  edm::EDGetTokenT<pat::JetCollection> jetSrcToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genSrcToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> PFSrcToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<reco::GenJetCollection> genjetToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  edm::EDGetTokenT<reco::PFJetCollection> ak5PFjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak5GenJetToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfCTEQWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfMMTHWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfNNPDFWeightsInputToken_;
  const edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<double> m_rho_token;

  float qscale;
  float wtfact; //MC : eventinfo->weight(); data : hltpres[ihltfill]*l1pres[ihltfill];
  int procid, npilup1, npilup2; //1:-5 to -1, 2:0 to 3

  int idall;
  float xfrac1, xfrac2, xpdf1, xpdf2;  

  //HLTConfigProvider hltConfig_;
   HLTPrescaleProvider hltPrescaleProvider_;
  int nreco, naa, nbb, ncc;

	std::vector<JetCorrectionUncertainty*> vsrc; // (nsrc);
reweight::PoissonMeanShifter PShiftUp_;
reweight::PoissonMeanShifter PShiftDown_;
edm::LumiReWeighting *LumiWeights_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
QCDEventShape::QCDEventShape(const edm::ParameterSet& iConfig):
  generator1_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("evtinfo"))),
  jetSrcToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  genSrcToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genSrc"))),
  PFSrcToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsSrc"))),
  genjetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjetSrc"))),
  pileup_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSrc"))),
  ak5PFjetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak5pfJetSrc"))),
  ak5GenJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak5genJetSrc"))),
  pdfCTEQWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFCTEQWeightsInputTag"))),
  pdfMMTHWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFMMTHWeightsInputTag"))),
  pdfNNPDFWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFNNPDFWeightsInputTag"))),
  LHERunInfoToken_(consumes<LHERunInfoProduct, edm::InRun >(iConfig.getParameter<edm::InputTag>("LHERunInfoProductInputTag"))),
  lheEventProductToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInputTag"))),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  m_rho_token = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  //m_resolutions_file = iConfig.getParameter<edm::FileEEInPath>("resolutionsFile").fullPath();
 // scalefile = iConfig.getParameter<edm::FileInPath>("scaleFactorsFile").fullPath();
  isHistFill = iConfig.getUntrackedParameter<bool>("HistFill", true);
  //  isHistFill2 = pset.getUntrackedParameter<bool>("HistFill2", false);                                            
  isTrigger = iConfig.getUntrackedParameter<bool>("Trigger", true);
	//  isRECO = iConfig.getUntrackedParameter<bool>("RECO", false);
  isMC = iConfig.getUntrackedParameter<bool>("MonteCarlo", false);
  isReconstruct = iConfig.getUntrackedParameter<bool>("Reconstruct", true);
  isJetQCD = iConfig.getUntrackedParameter<bool>("JetQCD", false);
  isGenJET = iConfig.getUntrackedParameter<bool>("GenJET", false);
  //  etarange = iConfig.getUntrackedParameter<double>("EtaRange", 5.0);
  ptthreshold = iConfig.getUntrackedParameter<double>("PtThreshold", 10.0);
  //leadingPtthreshold = iConfig.getUntrackedParameter<double>("LeadingPtThreshold", 40.0);
  isOtherAlgo = iConfig.getUntrackedParameter<bool>("OtherAlgo", false);
  weight2 = iConfig.getUntrackedParameter<double>("HistWeight", 1.0);
  weight = weight2;
  theHLTTag = iConfig.getUntrackedParameter<string>("HLTTag", "HLT");
  theRootFileName = iConfig.getUntrackedParameter<string>("RootFileName");
  //theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  //theFile->cd();
  //T1 = new TTree("T1", "QCDEvt");

  //T1->Branch("irun", &irun, "irun/I");  
  //T1->Branch("ilumi", &ilumi, "ilumi/I");  
  //T1->Branch("ievt", &ievt, "ievt/i");
  //T1->Branch("ibrnc", &ibrnc, "ibrnc/I");  
  //T1->Branch("nsicls", &nsicls, "nsicls/I");  //to pf neutral paritcle (excluding HF)
  //T1->Branch("ntottrk", &ntottrk, "ntottrk/I");  //total pfcharged particle (HF)

  //  //T1->Branch("jetpt", &jetpt, "jetpt/F");
  //  //T1->Branch("jeteta", &jeteta, "jeteta/F");
  //  //T1->Branch("jetphi", &jetphi, "jetphi/F");

  //T1->Branch("nallpf", &nallpf, "nallpf/I");
  //T1->Branch("ncharged", &ncharged, "ncharged/I");
  //  //T1->Branch("jtthan",&jtthan,"jtthan/F");
  //  //T1->Branch("thphi",thphi,"thphi[10]/F");

  //T1->Branch("thrust",thrust,"thrust[10]/F");

  //  //T1->Branch("anglex",anglex,"anglex[10]/F");

  //cout << "Testing 1 ==== " <<njecmx<< endl;


//===============================
//--------------------------------------------Define TUnfoldBinning--------------------------
  char RecoBinName[100], GenBinName[100], Axisname[100]; 
   for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
             if (isItUsed(ij)) {
//-----------------------------------2D Binning--------------------------------------------------------
         sprintf(RecoBinName, "Detector2d_typ_%i_eta%i_%i", ityp, iet, ij);
         binsRec2D[ityp][iet][ij] = new TUnfoldBinning(RecoBinName);

         sprintf(RecoBinName, "Recobin2d_typ_%i_eta%i_%i",ityp, iet, ij);
         RecoBinning2D[ityp][iet][ij]= binsRec2D[ityp][iet][ij]->AddBinning(RecoBinName);

         sprintf(Axisname, "var_%i", ij);
         RecoBinning2D[ityp][iet][ij]->AddAxis(Axisname,(ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij],(ityp==0) ? rbinrngs0[ij]: rbinrngs1[ij],false,false);
         sprintf(Axisname, "ht");
         RecoBinning2D[ityp][iet][ij]->AddAxis(Axisname, nHLTmx, recohtbins, false,false);

          //--------------------------------------------------------------------
         sprintf(GenBinName, "Generator2d_typ_%i_eta%i_%i",ityp, iet, ij);
         binsGen2D[ityp][iet][ij] = new TUnfoldBinning(GenBinName);

         sprintf(GenBinName, "Genbin2d_typ_%i_eta%i_%i", ityp, iet, ij);
         GenBinning2D[ityp][iet][ij]= binsGen2D[ityp][iet][ij]->AddBinning(GenBinName);

         sprintf(Axisname, "var_%i", ij);
         GenBinning2D[ityp][iet][ij]->AddAxis(Axisname,(ityp==0) ? nbinsx0[ij] : nbinsx1[ij],(ityp==0) ? binrngs0[ij]: binrngs1[ij],false,false);
         sprintf(Axisname, "ht");
         GenBinning2D[ityp][iet][ij]->AddAxis(Axisname, nHLTmx, recohtbins,false,false);
	   }
	 }
      }
    }


  char name[200];
  char title[200];

  for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
             if (isItUsed(ij)) {
	     if (isReconstruct) { 
              sprintf(name, "dd_reco_typ_%i_eta%i_%i", ityp, iet, ij);
	      sprintf(title, "2D Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
              h_recovar_2D[ityp][iet][ij] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title); //false : global bin ID
              h_recovar_2D[ityp][iet][ij]->Sumw2();              

              sprintf(name, "dd_fake_reco_typ_%i_eta%i_%i", ityp, iet, ij);
              sprintf(title, "2D Fake Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
              h_recofake_2D[ityp][iet][ij] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);//false : global bin ID
              h_recofake_2D[ityp][iet][ij]->Sumw2();              
	     }

#ifdef  LHAPDF
            for (int ix=1; ix<nnnmx; ix++) {
              sprintf(name, "dd_genpdf_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Genpdf %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
             h_genevtvarpdf_2D[ityp][iet][ij][ix] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
             h_genevtvarpdf_2D[ityp][iet][ij][ix]->Sumw2();
            }
#endif
#ifdef  JETENERGY
            for (int ix=1; ix<njecmx; ix++) {
              sprintf(name, "dd_recojec_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Recojec %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_recoevtvarjec_2D[ityp][iet][ij][ix] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_recoevtvarjec_2D[ityp][iet][ij][ix]->Sumw2();
            }
#elif defined(JETRESO)
            for (int ix=1; ix<njecmx; ix++ ) {
              sprintf(name, "dd_recoreso_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Recoreso %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_recoevtvarres_2D[ityp][iet][ij][ix] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_recoevtvarres_2D[ityp][iet][ij][ix]->Sumw2();
            
             
              sprintf(name, "dd_corr_reso_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Corr reso %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              RM_JER_2D[ityp][iet][ij][ix] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D[ityp][iet][ij], binsGen2D[ityp][iet][ij], name ,0,0, title);
              RM_JER_2D[ityp][iet][ij][ix]->Sumw2();

              sprintf(name, "dd_fake_reso_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D fake reso %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_reco_JER_fake_2D[ityp][iet][ij][ix] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_reco_JER_fake_2D[ityp][iet][ij][ix]->Sumw2();

              sprintf(name, "dd_miss_reso_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D miss reso %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_gen_JER_miss_2D[ityp][iet][ij][ix] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_gen_JER_miss_2D[ityp][iet][ij][ix]->Sumw2();
             }
#endif
	      sprintf(name, "dd_gen_typ_%i_eta%i_%i", ityp, iet, ij);
              sprintf(title, "2D Gen %s %g %s", typname[ityp] , etarange[iet], vartitle[ij]);
              h_genvar_2D[ityp][iet][ij] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);//false : global bin ID
              h_genvar_2D[ityp][iet][ij]->Sumw2();              

	      sprintf(name, "dd_miss_gen_typ_%i_eta%i_%i", ityp, iet, ij);
              sprintf(title, "Miss Gen %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
	      h_genmiss_2D[ityp][iet][ij] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title); //false : global bin ID
              h_genmiss_2D[ityp][iet][ij]->Sumw2();              

	      if (isReconstruct) {
              sprintf(name, "dd_corr_typ_%i_eta%i_%i", ityp , iet, ij);
              sprintf(title, "Gen_Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
	      RM_2D[ityp][iet][ij] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D[ityp][iet][ij], binsGen2D[ityp][iet][ij], name ,0,0, title);
              RM_2D[ityp][iet][ij]->Sumw2();              
            
	      } //if (isReconstruct)
	    }//if (isItUsed(ij))
          }
        }
      }
//-----------------------------------------Root Histogram
  for (int ityp=0; ityp<ntype; ityp++) {
    for (int ipt=0; ipt<njetptmn; ipt++) {
      for (int iet=0; iet<njetetamn; iet++) {
	for (int ij=0; ij<nvar; ij++) {
	  if (isItUsed(ij)) { 
	    if (isReconstruct) { 
	      sprintf(name, "reco_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	      sprintf(title, "Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
	      h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij] );  
	      //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );  
              //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );
              //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsx[ij], -endx[ij], -startx[ij]);
	      h_recoevtvar[ityp][ipt][iet][ij]->Sumw2();
	     //For Fake
	     sprintf(name, "fake_reco_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
             sprintf(title, "Fake Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
             h_recoevtfake[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);  //Reco : only Var
	     // h_recoevtfake[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );//Reco each HT,var
              h_recoevtfake[ityp][ipt][iet][ij]->Sumw2();
	    }
#ifdef  LHAPDF
	    for (int ix=1; ix<nnnmx; ix++) {
	      sprintf(name, "genpdf_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	      sprintf(title, "Genpdf %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
	     h_genevtvarpdf[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
	     // h_genevtvarpdf[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);//Reco each HT, var
	      h_genevtvarpdf[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#endif
	    
#ifdef  JETENERGY
	    for (int ix=1; ix<njecmx; ix++) {
	      sprintf(name, "recojec_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	      sprintf(title, "Recojec %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
	      h_recoevtvarjec[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
	      //h_recoevtvarjec[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
	      h_recoevtvarjec[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#elif defined(JETRESO)
	    for (int ix=1; ix<njecmx; ix++ ) {
	      sprintf(name, "recoreso_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	      sprintf(title, "Recoreso %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
              h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
              //h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij] , (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
              //h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
	      h_recoevtvarres[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#endif
	    
	    sprintf(name, "gen_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	    sprintf(title, "Gen %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );
              h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsxgen[ij], -endx[ij], -startx[ij]);
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsx[ij], -endx[ij], -startx[ij]); //For Equal Binning
             h_genevtvar[ityp][ipt][iet][ij]->Sumw2();
             //for miss
             sprintf(name, "miss_gen_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
             sprintf(title, "Miss Gen %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
             h_genevtmiss[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );
             //h_genevtmiss[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
             h_genevtmiss[ityp][ipt][iet][ij]->Sumw2();
	    
	    
	    sprintf(name, "gen2_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	    sprintf(title, "Gen2 %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            h_genevtvar2[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij] , (ityp==0) ? binrngs0[ij] : binrngs1[ij]);	   
            //h_genevtvar2[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);	   
	    h_genevtvar2[ityp][ipt][iet][ij]->Sumw2();

	    if (isReconstruct) { 
	      sprintf(name, "corr_typ_%i_pt%i_eta%i_%i", ityp , ipt, iet, ij);
	      sprintf(title, "Gen_Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij], (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
            //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij], (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
           //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? rnbinsx0[ij][ipt] : rnbinsx1[ij][ipt], (ityp==0) ? rbinrngs0[ij][ipt] : rbinrngs1[ij][ipt], (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
            //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt], (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]); //Reco =Gen bin Each HT, Var
          //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, nbinsxgen[ij], -endx[ij], -startx[ij], nbinsx[ij], -endx[ij], -startx[ij]);
	  //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, nbinsx[ij], -endx[ij], -startx[ij], nbinsx[ij], -endx[ij], -startx[ij]);
	    h_2devtvar[ityp][ipt][iet][ij]->Sumw2();
	    }
	  }
	  //cout <<"ijx "<< ityp<<" "<< ipt<<" "<< iet<<" "<<ij<<endl;
	  
	}
      }
    }
  }

 sprintf(name, "corr_jet");
 sprintf(title, "Gen_Reco_HT2");
 h_2ht=fs->make<TH2F>(name, title, 8, leadingPtThreshold, 8, leadingPtThreshold);


  //==============****================================  
#ifndef GENPART                     
  //  recojt_hist = fs->make<TH1F>("recojt_hist","# of recojets",20,-0.5, 19.5);
  // recojt_hist->Sumw2();
  //recojt_pt = fs->make<TH1F>("recojt_pt","Et_{recojets}",100,20., 2020.);
  //recojt_pt->Sumw2();
  recojt_eta = fs->make<TH1F>("recojt_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt_eta->Sumw2();
  recojt_phi = fs->make<TH1F>("recojt_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt_phi->Sumw2();

  //recojt1_pt = fs->make<TH1F>("recojet1_pt","Et_{recojets}",100,20., 2020.);
  //recojt1_pt->Sumw2();
  recojt1_eta = fs->make<TH1F>("recojet1_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt1_eta->Sumw2();
  recojt1_phi = fs->make<TH1F>("recojet1_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt1_phi->Sumw2();

  //recojt2_pt = fs->make<TH1F>("recojet2_pt","Et_{recojets}",100,20., 2020.);
  //recojt2_pt->Sumw2();
  recojt2_eta = fs->make<TH1F>("recojet2_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt2_eta->Sumw2();
  recojt2_phi = fs->make<TH1F>("recojet2_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt2_phi->Sumw2();

  // recojt3_pt = fs->make<TH1F>("recojet2_pt","Et_{recojets}",100,20., 2020.);
  //recojt3_pt->Sumw2();
  recojt3_eta = fs->make<TH1F>("recojet2_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt3_eta->Sumw2();
  recojt3_phi = fs->make<TH1F>("recojet2_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt3_phi->Sumw2();


  for(int jk=0; jk<njetetamn; jk++){
    sprintf(name, "recojetallave_pt_%i",jk);
    sprintf(title, "Et_{recojetsallave}_%g", etarange[jk]);
    recojtallave_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    recojtallave_pt[jk]->Sumw2();

    sprintf(name, "recojetallavewt1_pt_%i",jk);
    sprintf(title, "Et_{recojetsallavewt1}_%g", etarange[jk]);
    recojtallavewt1_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    recojtallavewt1_pt[jk]->Sumw2();

    sprintf(name, "recojt_pt_%i",jk);
    sprintf(title, "Et_{recojets}_%g", etarange[jk]);
    recojt_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt_pt[jk]->Sumw2();

    sprintf(name, "recojet1_pt_%i",jk);
    sprintf(title, "Et_{recojets1}_%g", etarange[jk]);
    recojt1_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1_pt[jk]->Sumw2();


    sprintf(name, "recojet2_pt_%i",jk);
    sprintf(title, "Et_{recojets2}_%g", etarange[jk]);
    recojt2_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2_pt[jk]->Sumw2();


    sprintf(name, "recojet3_pt_%i",jk);
    sprintf(title, "Et_{recojets3}_%g", etarange[jk]);
    recojt3_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt3_pt[jk]->Sumw2();


    for (int kl=0; kl<nHLTmx; kl++) { 
      //  sprintf(name, "recojt_pt_%i_%i",jk, kl);
      //  sprintf(title, "Et_{recojets}_%g_%i", etarange[jk], kl);
      //  recojt_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      //  recojt_pt[jk][kl]->Sumw2();

      sprintf(name, "recojetave_pt_%i_%i",jk, kl);
      sprintf(title, "Et_{recojetsave}_%g_%i", etarange[jk], kl);
      recojtave_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      recojtave_pt[jk][kl]->Sumw2();

      sprintf(name, "recojetavewt1_pt_%i_%i",jk, kl);
      sprintf(title, "Et_{recojetsavewt1}_%g_%i", etarange[jk], kl);
      recojtavewt1_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      recojtavewt1_pt[jk][kl]->Sumw2();
    }

    sprintf(name, "recojetHT2_%i",jk);
    sprintf(title, "recojetsHT2_%g", etarange[jk]);

    recoht2_pt[jk] = fs->make<TH1F>(name, title, 400,20., 1500.);
    recoht2_pt[jk]->Sumw2();


    sprintf(name, "hjetdpt_%i",jk);
    sprintf(title, "dpt_{recojets12}_%g", etarange[jk]);

    hjetdpt[jk] = fs->make<TH1F>(name, title, 100,20., 500.);
    hjetdpt[jk]->Sumw2();

    sprintf(name, "hjetpt2bypt1_%i",jk);
    sprintf(title, "hjetpt2bypt1 reco jet_%g", etarange[jk]);


    hjetpt2bypt1[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetpt2bypt1[jk]->Sumw2();

    sprintf(name, "hjetpt3bypt2_%i",jk);
    sprintf(title, "hjetpt3bypt2 reco jet_%g", etarange[jk]);
    hjetpt3bypt2[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetpt3bypt2[jk]->Sumw2();


    sprintf(name, "hjetdphi_%i",jk);
    sprintf(title, "#phi_{recojets}_%g", etarange[jk]);
    hjetdphi[jk] = fs->make<TH1F>(name,title,100,-M_PI, M_PI);
    hjetdphi[jk]->Sumw2();
    sprintf(name, "hjetptbypl_%i",jk);
    sprintf(title, "1st recojet Pt*sin/1st Recojet_%g", etarange[jk]);
    hjetptbypl[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetptbypl[jk]->Sumw2();

    //hjetpt2bypt1 = fs->make<TH1F>("hjetpt2bypt1", "hjetpt2bypt1 reco jet", 60, 0., 1.0);
    //hjetpt2bypt1->Sumw2();
    //hjetpt3bypt2 = fs->make<TH1F>("hjetpt2bypt1", "hjetpt2bypt1 reco jet", 60, 0., 1.0);
    //hjetpt3bypt2->Sumw2();

  }

  recochg_pt = fs->make<TH1F>("recochg_pt","Et_{recocharge_alljet}",100, 1., 101.);
  recochg_pt->Sumw2();
  recochg_eta = fs->make<TH1F>("recochg_eta","#eta_{recocharge_alljet}",100,-3., 3.);
  recochg_eta->Sumw2();
  recochg_phi = fs->make<TH1F>("recochg_phi","#phi_{recocharge_alljet}",100,-M_PI, M_PI);
  recochg_phi->Sumw2();


  recochg1_pt = fs->make<TH1F>("recochg1_pt","Et_{recocharge_jet1}",100, 1., 101.);
  recochg1_pt->Sumw2();
  recochg1_eta = fs->make<TH1F>("recochg1_eta","#eta_{recocharge_jet1}",100,-3., 3.);
  recochg1_eta->Sumw2();
  recochg1_phi = fs->make<TH1F>("recochg1_phi","#phi_{recocharge_jet1}",100,-M_PI, M_PI);
  recochg1_phi->Sumw2();

  recochg2_pt = fs->make<TH1F>("recochg2_pt","Et_{recocharge_jet2}",100, 1., 101.);
  recochg2_pt->Sumw2();
  recochg2_eta = fs->make<TH1F>("recochg2_eta","#eta_{recocharge_jet2}",100,-3., 3.);
  recochg2_eta->Sumw2();
  recochg2_phi = fs->make<TH1F>("recochg2_phi","#phi_{recocharge_jet2}",100,-M_PI, M_PI);
  recochg2_phi->Sumw2();

  recochg3_pt = fs->make<TH1F>("recochg3_pt","Et_{recocharge_jet3}",100, 1., 101.);
  recochg3_pt->Sumw2();
  recochg3_eta = fs->make<TH1F>("recochg3_eta","#eta_{recocharge_jet3}",100,-3., 3.);
  recochg3_eta->Sumw2();
  recochg3_phi = fs->make<TH1F>("recochg3_phi","#phi_{recocharge_jet3}",100,-M_PI, M_PI);
  recochg3_phi->Sumw2();



#endif

 //================*****===================================
  for (int ij=0; ij<nhist; ij++) {
    sprintf(name, "anglex_%i", ij);
    vec_anglex[ij] = fs->make<TH1F>(name, name, 240, 0.7, 1.0);
  }
  // genjt_hist = fs->make<TH1F>("genjt_hist","# of genjets",20,-0.5, 19.5);
  // genjt_hist->Sumw2();
  for(int jk=0; jk<njetetamn; jk++){
    sprintf(name, "genjetallave_pt_%i",jk);
    sprintf(title, "Et_{genjetsallave}_%g", etarange[jk]);
    genjtallave_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    genjtallave_pt[jk]->Sumw2();

    sprintf(name, "genjt_pt_%i",jk);
    sprintf(title, "Et_{genjets}_%g", etarange[jk]);
    genjt_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt_pt[jk]->Sumw2();

    sprintf(name, "genjet1_pt_%i",jk);
    sprintf(title, "Et_{genjets1}_%g", etarange[jk]);
    genjt1_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt1_pt[jk]->Sumw2();

    sprintf(name, "genjet2_pt_%i",jk);
    sprintf(title, "Et_{genjets2}_%g", etarange[jk]);
    genjt2_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt2_pt[jk]->Sumw2();

    sprintf(name, "genjet3_pt_%i",jk);
    sprintf(title, "Et_{genjets3}_%g", etarange[jk]);
    genjt3_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt3_pt[jk]->Sumw2();

    /*sprintf(name, "genjt_oth_pt_%i",jk);
    sprintf(title, "#Et_{genjets_oth}_%g", etarange[jk]);

    genjt_oth_pt[jk] = fs->make<TH1F>(name,title,100, 20., 2020.);
    genjt_oth_pt[jk]->Sumw2();
   */
    sprintf(name, "genjetdpt_%i",jk);
    sprintf(title, "dpt_{genjets12}_%g", etarange[jk]);

    genjetdpt[jk] = fs->make<TH1F>(name, title, 100,20., 500.);
    genjetdpt[jk]->Sumw2();

    sprintf(name, "genjetpt2bypt1_%i",jk);
    sprintf(title, "jetpt2bypt1 gen jet_%g", etarange[jk]);

    genjetpt2bypt1[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetpt2bypt1[jk]->Sumw2();

    sprintf(name, "genjetpt3bypt2_%i",jk);
    sprintf(title, "hjetpt3bypt2 gen jet_%g", etarange[jk]);
    genjetpt3bypt2[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetpt3bypt2[jk]->Sumw2();

    sprintf(name, "genjetdphi_%i",jk);
    sprintf(title, "#phi_{genjets}_%g", etarange[jk]);
    genjetdphi[jk] = fs->make<TH1F>(name,title,100,-M_PI, M_PI);
    genjetdphi[jk]->Sumw2();
    sprintf(name, "genjetptbypl_%i",jk);
    sprintf(title, "1st genjet Pt*sin/1st genjet_%g", etarange[jk]);
    genjetptbypl[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetptbypl[jk]->Sumw2();
  }

  for(int jk=0; jk<njetetamn+2; jk++){
    sprintf(name, "response_jet_distribution_%i",jk);
    sprintf(title, "response_2d_ratio_eta_%g_%g", resetarange[jk], resetarange[jk+1]);
    resp_jet[jk] = fs->make<TH2F>(name,title,400, 20., 2020., 50 , 0.0, 5.0);
    resp_jet[jk]->Sumw2();
    sprintf(name, "response_jet_distribution_resolution%i",jk);
    sprintf(title, "response_1d_resulation_eta_%g_%g", resetarange[jk], resetarange[jk+1]);
    resp_jet1[jk] = fs->make<TH1F>(name,title, 100 , 0.0, 1);
    resp_jet1[jk]->Sumw2();
  }

  //genjt_pt = fs->make<TH1F>("genjt_pt","Et_{genjets}",100,20., 2020.);
  // genjt_pt->Sumw2();
  genjt_eta = fs->make<TH1F>("genjt_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt_eta->Sumw2();
  genjt_phi = fs->make<TH1F>("genjt_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt_phi->Sumw2();

  genjt1_eta = fs->make<TH1F>("genjet1_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt1_eta->Sumw2();
  genjt1_phi = fs->make<TH1F>("genjet1_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt1_phi->Sumw2();

  genjt2_eta = fs->make<TH1F>("genjet2_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt2_eta->Sumw2();
  genjt2_phi = fs->make<TH1F>("genjet2_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt2_phi->Sumw2();

  genjt3_eta = fs->make<TH1F>("genjet2_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt3_eta->Sumw2();
  genjt3_phi = fs->make<TH1F>("genjet2_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt3_phi->Sumw2();
  //  genjt_oth_pt = fs->make<TH1F>("genjt_oth_pt","Et_{genjets_oth}",100, 20., 2020.);
  //  genjt_oth_pt->Sumw2();
/*  genjt_oth_eta = fs->make<TH1F>("genjt_oth_eta","#eta_{genjets_oth}",100,-5., 5.);
  genjt_oth_eta->Sumw2();
  genjt_oth_phi = fs->make<TH1F>("genjt_oth_phi","#phi_{genjets_oth}",100,-M_PI, M_PI);
  genjt_oth_phi->Sumw2();
*/
  // genchg_hist = fs->make<TH1F>("genchg_hist","# of genchargeds",120,-0.5, 239.5);
  // genchg_hist->Sumw2();
  genchg_pt = fs->make<TH1F>("genchg_pt","Et_{gencharge_alljet}",100, 1., 101.);
  genchg_pt->Sumw2();
  genchg_eta = fs->make<TH1F>("genchg_eta","#eta_{gencharge_alljet)",100,-3., 3.);
  genchg_eta->Sumw2();
  genchg_phi = fs->make<TH1F>("genchg_phi","#phi_{gencharge_alljet}",100,-M_PI, M_PI);
  genchg_phi->Sumw2();

  genchg1_pt = fs->make<TH1F>("genchg1_pt","Et_{gencharge_jet1}",100, 1., 101.);
  genchg1_pt->Sumw2();
  genchg1_eta = fs->make<TH1F>("genchg1_eta","#eta_{gencharge_jet1}",100,-3., 3.);
  genchg1_eta->Sumw2();
  genchg1_phi = fs->make<TH1F>("genchg1_phi","#phi_{gencharge_jet1}",100,-M_PI, M_PI);
  genchg1_phi->Sumw2(); 

  genchg2_pt = fs->make<TH1F>("genchg2_pt","Et_{gencharge_jet2}",100, 1., 101.);
  genchg2_pt->Sumw2();
  genchg2_eta = fs->make<TH1F>("genchg2_eta","#eta_{gencharge_jet2}",100,-3., 3.);
  genchg2_eta->Sumw2();
  genchg2_phi = fs->make<TH1F>("genchg2_phi","#phi_{gencharge_jet2}",100,-M_PI, M_PI);
  genchg2_phi->Sumw2();

  genchg3_pt = fs->make<TH1F>("genchg3_pt","Et_{gencharge_jet3}",100, 1., 101.);
  genchg3_pt->Sumw2();
  genchg3_eta = fs->make<TH1F>("genchg3_eta","#eta_{gencharge_jet3}",100,-3., 3.);
  genchg3_eta->Sumw2();
  genchg3_phi = fs->make<TH1F>("genchg3_phi","#phi_{gencharge_jet3}",100,-M_PI, M_PI);
  genchg3_phi->Sumw2();

 

  // genchg_oth_hist = fs->make<TH1F>("genchg_oth_hist","# of genchargeds (others)",120,-0.5, 239.5);
  // genchg_oth_hist->Sumw2();
 /* genchg_oth_pt = fs->make<TH1F>("genchg_oth_pt","Et_{genchargeds_oth}",100,1., 101.);
  genchg_oth_pt->Sumw2();
  genchg_oth_eta = fs->make<TH1F>("genchg_oth_eta","#eta_{genchargeds_oth}",100,-5., 5.);
  genchg_oth_eta->Sumw2();
  genchg_oth_phi = fs->make<TH1F>("genchg_oth_phi","#phi_{genchargeds_oth}",100,-M_PI, M_PI);
  genchg_oth_phi->Sumw2();
  // genneu_hist = fs->make<TH1F>("genneu_hist","# of genneutrals",120,-0.5, 239.5);
  // genneu_hist->Sumw2();
  genneu_pt = fs->make<TH1F>("genneu_pt","Et_{genneutrals}",100,1., 101.);
  genneu_pt->Sumw2();
  genneu_eta = fs->make<TH1F>("genneu_eta","#eta_{genneutrals}",100,-3., 3.);
  genneu_eta->Sumw2();
  genneu_phi = fs->make<TH1F>("genneu_phi","#phi_{genneutrals}",100,-M_PI, M_PI);
  genneu_phi->Sumw2();

  // genneu_oth_hist = fs->make<TH1F>("genneu_oth_hist","# of genneutrals (others)",120,-0.5, 239.5);
  // genneu_oth_hist->Sumw2();
  genneu_oth_pt = fs->make<TH1F>("genneu_oth_pt","Et_{genneutrals_oth}",100, 1., 101.);
  genneu_oth_pt->Sumw2();
  genneu_oth_eta = fs->make<TH1F>("genneu_oth_eta","#eta_{genneutrals_oth}",100,-5., 5.);
  genneu_oth_eta->Sumw2();
  genneu_oth_phi = fs->make<TH1F>("genneu_oth_phi","#phi_{genneutrals_oth}",100,-M_PI, M_PI);
  genneu_oth_phi->Sumw2();
*/
  for (int ij=0; ij<nHLTmx; ij++) { 
    sprintf(name, "nprimall_%i", ij);
    sprintf(title, "# of primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_hist[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_hist[ij]->Sumw2();

    sprintf(name, "nprimsel_%i", ij);
    sprintf(title, "Selected # of primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_sel[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_sel[ij]->Sumw2();

    sprintf(name, "nprimall_rewt_%i", ij);
    sprintf(title, "# of rewighted primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_hist_rewt[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_hist_rewt[ij]->Sumw2();

    sprintf(name, "nprimsel_rewt_%i", ij);
    sprintf(title, "Selected # of reweighted primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_sel_rewt[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_sel_rewt[ij]->Sumw2();
  }

  prim_correl = fs->make<TH2F>("correl", "Correlation of all and Selected # of primary vtx", 60, -0.5, 59.5, 60, -0.5, 59.5);
  const char* namex[2]={"Selected", "Rejected"};
  for (int ij=0; ij<2; ij++) {
    sprintf(name, "primalltrk_%i", ij);
    sprintf(title, "All tracks in primary vtx (%s)", namex[ij]);
    prim_alltrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primgoodtrk_%i", ij);
    sprintf(title, "Good tracks in primary vtx (%s)", namex[ij]);
    prim_goodtrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primseltrk_%i", ij);
    sprintf(title, "Selected tracks in primary vtx (%s)", namex[ij]);
    prim_seltrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primdx_%i", ij);
    sprintf(title, "#Delta x of prim wrt beam spot (%s)", namex[ij]);
    prim_dx[ij] = fs->make<TH1F>(name, title, 120, -2.4, 2.4);

    sprintf(name, "primdy_%i", ij);
    sprintf(title, "#Delta y of prim wrt beam spot (%s)", namex[ij]);
    prim_dy[ij] = fs->make<TH1F>(name, title, 120, -2.4, 2.4);

    sprintf(name, "primdxy_%i", ij);
    sprintf(title, "#Delta y vs #Delta x of prim (%s)", namex[ij]);
    prim_dxy[ij] = fs->make<TH2F>(name, title, 60, -0.15, 0.15, 60, -0.15, 0.15);

    sprintf(name, "primdz_%i", ij);
    sprintf(title, "#Delta z of prim wrt beam spo (%s)", namex[ij]);
    prim_dz[ij] = fs->make<TH1F>(name, title, 120, -30.0, 30.0); 

    sprintf(name, "primprob_%i", ij);
    sprintf(title, "log10(vertex fit prob) (%s)", namex[ij]);
    prim_prob[ij] = fs->make<TH1F>(name, title, 120, -20.0, 0.0);   
  }

  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "njets_%i",ij);
    sprintf(title, "No of Jets_eta range_%gs", etarange[ij]);
    h_njets[ij] = fs->make<TH1F>(name, title, 60, 0, 30);
    h_njets[ij]->Sumw2();
  }

  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "ncharges_%i",ij);
    sprintf(title, "No of charge particles_eta range_%gs", etarange[ij]);
    h_nchg[ij] = fs->make<TH1F>(name, title, 800, 0, 400);
    h_nchg[ij]->Sumw2();
  }


  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "gennjets_%i",ij);
    sprintf(title, "No of GenJets_eta range_%gs", etarange[ij]);
    gen_njets[ij] = fs->make<TH1F>(name, title, 60, 0, 30);
    gen_njets[ij]->Sumw2();
  }


#ifdef TRIGGER
  const char* trigvar[2]={"L1", "HLT"};
  for(int ij=0; ij<nHLTmx; ij++){
    for(int jk=0; jk<2; jk++){
      sprintf(name, "trgjet_pt_%i_%i", ij, jk);
      sprintf(title, "trgjet_pt_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_pt[ij][jk] = fs->make<TH1F>(name, title, njetptbin, 20,1500);
      trgjet_pt[ij][jk]->Sumw2();

      sprintf(name, "trgjet_eta_%i_%i", ij, jk);
      sprintf(title, "trgjet_eta_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_eta[ij][jk] = fs->make<TH1F>(name, title, njetptbin, -5., 5.);
      trgjet_eta[ij][jk]->Sumw2();

      sprintf(name, "trgjet_phi_%i_%i", ij, jk);
      sprintf(title, "trgjet_phi_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_phi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
      trgjet_phi[ij][jk]->Sumw2();

      sprintf(name, "prbjet_pt_%i_%i", ij, jk);
      sprintf(title, "prbjet_pt_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_pt[ij][jk] = fs->make<TH1F>(name, title, njetptbin, 20,1500);
      prbjet_pt[ij][jk]->Sumw2();

      sprintf(name, "prbjet_eta_%i_%i", ij, jk);
      sprintf(title, "prbjet_eta_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_eta[ij][jk] = fs->make<TH1F>(name, title, 100,-5., 5.);
      prbjet_eta[ij][jk]->Sumw2();

      sprintf(name, "prbjet_phi_%i_%i", ij, jk);
      sprintf(title, "prbjet_phi_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_phi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
      prbjet_phi[ij][jk]->Sumw2();


    }
  } 
#endif
//Trigger special

//=================================

	if (isReconstruct) { 
		for(int ij=0; ij<nHLTmx; ij++){
			for(int jk=0; jk<njetetamn; jk++){
				sprintf(name, "jetpt_%i_%i",jk,ij);
				sprintf(title, "jetpt_%s_%g", jethlt_name[ij], etarange[jk]);
				h_jetpt[ij][jk] = fs->make<TH1F>(name, title, 300, 50, 1550);
				h_jetpt[ij][jk]->Sumw2();
				
				sprintf(name, "jetphi_%i_%i",jk, ij);
				sprintf(title, "jetphi_%s_%g", jethlt_name[ij],etarange[jk]);
				h_jetphi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
				h_jetphi[ij][jk]->Sumw2();
				
			}
		}
	}
#ifdef TRIGGER
  for(int ij=0; ij<nHLTmx; ij++){

    sprintf(name, "jeteta_%i", ij);
    sprintf(title, "jetphi_%s", jethlt_name[ij]);//, jetvar[ij]);
    h_jeteta[ij] = fs->make<TH1F>(name, title, 100, -5, 5);
    h_jeteta[ij]->Sumw2();

    for (int jk=0; jk<2; jk++){ 
      sprintf(name, "angle1d_%s_%i", jethlt_name[ij], jk);
      sprintf(title, "Angle%s_%i", jethlt_name[ij], jk);
      trgjet_angle[ij][jk] = fs->make<TH1F>(name, title, 90 , 0.1, 2.5);

      sprintf(name, "angle2d_%s_%i", jethlt_name[ij], jk);
      sprintf(title, "Angle_2d_hist%s_%i", jethlt_name[ij], jk);
      trgjet_2dangle[ij][jk] = fs->make<TH2F>(name, title, njetptbin, 20, 1500, 30 , 0.1, 2.5);
    }
  }

  for (int ij=0; ij<nHLTmx; ij++) {
    for (int jk=0; jk<njetetamn; jk++) {
      sprintf(name, "hlt_dijettag_%i_%i", ij, jk);
      sprintf(title, "dijet tagged P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
      hlt_dijettag[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettag[ij][jk]->Sumw2();

      sprintf(name, "hlt_dijetprob_%i_%i", ij, jk);
      sprintf(title, "dijet probed P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
      hlt_dijetprob[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijetprob[ij][jk]->Sumw2();
    }
  }
#endif
  counthist = fs->make<TH1F>("count","No of events",2,0,2); 



  for (int ix=0; ix<32; ix++) { mypow_2[ix] = pow(2,ix);}
  nevt = 0;
  // irun_old=-1;
  //trig_init=0;

  nreco=naa= nbb= ncc=0;

}




QCDEventShape::~QCDEventShape()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void QCDEventShape::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
 // t1=clock();
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  //gRandom->SetSeed(19919925);
  //float rn=gRandom->Uniform();
  //cout << " Random Number ini = " << rn << endl;
  //if (rn >0.90) return;
 // cout << " Random Number = " << rn << endl;
  //cout << "Time = " << t1 << "; " << t2 << endl;
  nevt++;
  int ievt = iEvent.id().event();
  counthist->Fill(1); 
  if (nevt%10000==1)   std::cout<<"QCDEventShape::analyze "<< nevt<<" IRUN= "<<iEvent.id().run()<<" ievt= "<< iEvent.id().event()<<" "<<ievt<<endl;
  //" ilumi" <<
 //iEvent.luminosityBlock() << " ibunch " << iEvent.bunchCrossing() <<std::endl;
 // cout << "NEvent = " <<  nevt << endl;
 // if(iEvent.luminosityBlock()==9881 || iEvent.luminosityBlock()==23185 || iEvent.luminosityBlock()==25334 || iEvent.luminosityBlock()== 26584 ||iEvent.luminosityBlock()== 35674 || iEvent.luminosityBlock()==32764 || iEvent.luminosityBlock()== 35675 || iEvent.luminosityBlock()==53681) return ;
 //if(iEvent.luminosityBlock()==2 || iEvent.luminosityBlock()==7175 || iEvent.luminosityBlock()==41151 || iEvent.luminosityBlock()==7389697 || iEvent.luminosityBlock()==60334 || iEvent.luminosityBlock()==51317 || iEvent.luminosityBlock()==53654 || iEvent.luminosityBlock()==10333 || iEvent.luminosityBlock()==54778 || iEvent.luminosityBlock()==10082 || iEvent.luminosityBlock()==54322 || iEvent.luminosityBlock()==64667 || iEvent.luminosityBlock()==65977 || iEvent.luminosityBlock()==55534 || iEvent.luminosityBlock()==55781 || iEvent.luminosityBlock()==55782 || iEvent.luminosityBlock()==55783 || iEvent.luminosityBlock()==61360 || iEvent.luminosityBlock()==61370 ||iEvent.luminosityBlock()==68258 || iEvent.luminosityBlock()==62147 || iEvent.luminosityBlock()==67194 || iEvent.luminosityBlock()==43070 || iEvent.luminosityBlock()==49429 || iEvent.luminosityBlock()==15102 || iEvent.luminosityBlock()==23306 || iEvent.luminosityBlock()==14242|| iEvent.luminosityBlock()==19080 || iEvent.luminosityBlock()==9312025) return;
  npfjets = 0;
//  if(iEvent.luminosityBlock()<4401) return; 
 //   if(nevt<3442) return;
//  if(nevt!=3080) return;
// cout << "Write test 1 = ok " << endl;
 //=======================*****======================================                                               
  vector<double> recovar;
  vector<double> recovar1;
  vector<double> recovar2[njecmx];   //for JER
  std::vector<HepLorentzVector> recomom[njecmx][ntype][njetetamn];
  std::vector<HepLorentzVector> tmpjt4v; 
  std::vector<HepLorentzVector> tmpcand4v;                
  std::vector<HepLorentzVector> tmpgen4v;  
                                                                                             
  std::vector<HepLorentzVector> genmom[nGenReso][ntype][njetetamn];
  vector<double> genvar;

  //====================*****===========================================        

  wtfact=1.0;
//  double px=0;
//  double py=0;
//	double ptxy=0;

 // int ncount=0;
  unsigned ncount=0;
//  double recterm=0;
//  int ithird=-1;
  int irecoht=-1;
	//#ifdef JETENERGY
	int irecohtjec[njecmx];
	for (int ij=0; ij<njecmx; ij++) { irecohtjec[ij]=-1;}
	//#endif	
  double aveleadingptjec[njecmx] ={0};//14Sep20

  int igenht=-1;
	//#ifdef  JETRESO
	int igenhtres[nGenReso];
	for (int ij=0; ij<nGenReso; ij++) { igenhtres[ij]=-1;}
	//#endif
      double avegenptres[nGenReso]={0};//14Sep20

#ifdef TRIGGER
  const char* variab1;
#endif
#ifndef DIJETAVE
  const char* variab2; 
#endif

  if (isMC) {
#ifdef LHAPDF
    edm::Handle<LHEEventProduct> EvtHandle ;
    iEvent.getByToken( lheEventProductToken_ , EvtHandle ) ;
		
		for ( unsigned int weightIndex = 0; weightIndex < EvtHandle->weights().size(); ++weightIndex ) {
			//      cout<< EvtHandle->weights()[weightIndex].wgt <<endl;
      //systematicWeightIDs->push_back( atoi(EvtHandle->weights()[weightIndex].id.c_str()) );
			if (weightIndex>=9 && weightIndex<=109) {
				pdfwt[weightIndex-9] = EvtHandle->weights()[weightIndex].wgt/EvtHandle->originalXWGTUP(); 
				//				std::cout << weightIndex << " " << EvtHandle->weights()[weightIndex].id << " " << EvtHandle->weights()[weightIndex].wgt <<" "<<pdfwt[weightIndex-9]<< std::endl;
			}
    }
#endif
		//		cout<<"AAAAAAAAAAAAAA"<<endl;
    edm::Handle<GenEventInfoProduct> eventinfo;
    iEvent.getByToken(generator1_, eventinfo);
    if (eventinfo.isValid()) { 
      qscale = eventinfo->qScale(); 
      wtfact = eventinfo->weight();
      // weight = weight2*wtfact;
      procid = eventinfo->signalProcessID();
       //cout << " qscale = " <<setw(14)<< qscale << " ; wtfact = " << wtfact << " ; procid = " << procid  << endl;

      if (eventinfo->hasPDF()) {
	const gen::PdfInfo* xpdf = eventinfo->pdf();
	
	int id1 = xpdf->id.first;
	int id2 = xpdf->id.second;
	
	idall = 100*(id1+50)+ (id2+50); 
	
	qscale = xpdf->scalePDF;
	
	xfrac1 = xpdf->x.first;
	xfrac2 = xpdf->x.second;
	
	xpdf1 = xfrac1*xpdf->xPDF.first;
	xpdf2 = xfrac2*xpdf->xPDF.second; 
      }
    }
  }
  
#ifdef TRIGGER
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerBits_, trigRes);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  //---------------------------------------------------------------------Trigger
  const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  //	int ihltfill = -1;
#endif
  
  tmpjt4v.clear();
  tmpcand4v.clear();
  tmpgen4v.clear();
  double aveleadingpt =0;
  bool isInEtaRange[njetetamn]={0}; //GMA{0,0,0,0};
#ifndef GENPART
  edm::Handle<pat::JetCollection> ak4PFJets;
  if (isReconstruct) { 
    iEvent.getByToken(jetSrcToken_, ak4PFJets);
  }
  //  cout<<"1 aveleadingpt"<<endl;
  if (isReconstruct && ((!ak4PFJets.isValid()) ||  ak4PFJets->size() <2)) return; //GMA, do we use this
  
  if (ak4PFJets.isValid() &&  ak4PFJets->size()>=2) {
#ifdef DIJETAVE
    //    aveleadingpt = 0.5*((*ak4PFJets)[0].pt() + (*ak4PFJets)[1].pt());
    //cout<<"1 aveleadingpt"<<aveleadingpt<<endl;
    for (int iet=0; iet<njetetamn; iet++) {
      isInEtaRange[iet] = true;
    }
    
    for (int ij=0; ij<2; ij++) { 
      for (int iet=0; iet<njetetamn; iet++) {
	if (abs((*ak4PFJets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
      }
      
      //Jet ID ================= 2017 jetID recomendation 
      double NHF = (*ak4PFJets)[ij].neutralHadronEnergyFraction();
      double NEMF = (*ak4PFJets)[ij].neutralEmEnergyFraction();
      double CHF = (*ak4PFJets)[ij].chargedHadronEnergyFraction();
      //double MUF = (*ak4PFJets)[ij].muonEnergyFraction();
      //double CEMF = (*ak4PFJets)[ij].chargedEmEnergyFraction();
      int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
      //int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
      int CHM = (*ak4PFJets)[ij].chargedMultiplicity();

      bool TightJetID =false;
      //if(abs((*ak4PFJets)[ij].eta())<=2.7){  //Updated for ReReco 17
      // if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && (abs((*ak4PFJets)[ij].eta())<=2.4 && CHF>0 && CHM>0) ) TightJetID =true;
      //      } else {
      //                           TightJetID =false;}
      //Updated for UL17 : 27Aug20
      if(abs((*ak4PFJets)[ij].eta())<=2.6){
      if(NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0)  TightJetID =true;
                      } else {
                                 TightJetID =false;
                                 }
      if (abs((*ak4PFJets)[ij].eta())>2.6) {TightJetID = false;}
      if ((*ak4PFJets)[ij].pt()<30.0) {TightJetID = false;}

      if (TightJetID) {
				aveleadingpt +=(*ak4PFJets)[ij].pt();
      } else {
				aveleadingpt -=100000;
      }
    }
    aveleadingpt /=2.0;
#else

#endif
  }
#endif
  //  cout<<"2 aveleadingpt"<<endl;
  if (isReconstruct && isMC && aveleadingpt>3*qscale) return;

  irecoht = getbinid(aveleadingpt, nHLTmx, leadingPtThreshold);

#ifdef TRIGGER
  bool trgpas[nHLTmx]={0,0,0,0,0,0,0,0};
  //  if (isMC && ak4PFJets.isValid() &&  ak4PFJets->size()>=2) {
  //    aveleadingpt = (*ak4PFJets)[0].pt();
  //    irecoht = getbinid(aveleadingpt, nHLTmx, leadingPtThreshold);
  //  }

  // cout<<"ave"<<aveleadingpt<<" ; Jet1 Pt= " <<   (*ak4PFJets)[0].pt() <<" ; Jet2 Pt= " <<  (*ak4PFJets)[1].pt()<<endl; 
  //std::pair<std::vector<std::pair<std::string,int> >,int> prescaleValuesInDetail(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& trigger) const;
/* 
  //Calcualte Trigger Efficiency for dijet events Manas Sir
 bool trg_fired1=false;
  for (int jk=0; jk<nHLTmx-1; jk++) {
      bool trg_fired=false;
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str();
      //bool trg_fired=false;
      for (int iet=0; iet<njetetamn; iet++) {
	if(aveleadingpt*2>jethlt_thr[jk] && (strstr(variab1,jethlt_name[jk]) && (strlen(variab1)-strlen(jethlt_name[jk])<5)) ){ 
	//if(aveleadingpt*2>jethlt_thr[jk] && aveleadingpt*2<jethlt_thr[jk] && (strstr(variab1,jethlt_name[jk]) && (strlen(variab1)-strlen(jethlt_name[jk])<5)) ){ 
	  //for (int iet=0; iet<njetetamn; iet++) {
	  if (isInEtaRange[iet]) {
	    if (trigRes->accept(ij)) {
                hlt_dijettag[jk][iet]->Fill(aveleadingpt);
	        trg_fired=true;
               } 
	    if(trg_fired) trg_fired1=true;
            else trg_fired1=false;
          if(iet==0)  cout << jk << " ; ij" <<ij <<" ; Prb name "  << variab1 << endl;
	 // }
	//}
	
	if(trg_fired1 && (strstr(variab1,jethlt_name[jk+1]) && (strlen(variab1)-strlen(jethlt_name[jk+1])<5))){
	  if (trigRes->accept(ij) && isInEtaRange[iet]) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt);}
         if(iet==0) cout << jk <<" Tag name "  << variab1 << endl;
	}
       }
      } //for (int iet=0; iet<njetetamn; iet++) {
    }
 }
}
*/
//Calcualte Trigger Efficiency for dijet events Test
/*      
 bool trg_fired1=false; 
  for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
    std::string name = names.triggerName(ij);
    variab1 = name.c_str();
//   std::cout << "Trigger " << names.triggerName(ij) << 
  //              ", prescale " << triggerPrescales->getPrescaleForIndex(ij) <<
    //            ": " << (trigRes->accept(ij) ? "PASS" : "fail (or not run)") 
      //          << std::endl;
    bool trg_fired=false;
    for (int iet=0; iet<njetetamn-2; iet++) {
      if (isInEtaRange[iet]) {
	if((strstr(variab1,jethlt_name[2]) && (strlen(variab1)-strlen(jethlt_name[2])<5))){
	  if (trigRes->accept(ij)) {
	  //if (isL3) {
	    hlt_dijettag[2][iet]->Fill(aveleadingpt);
	    trg_fired=true;
	  if(iet==0) cout <<" Tag name "  << variab1 << endl;
	  }
	  if(trg_fired) trg_fired1=true;
	  else trg_fired1=false;
	}
	if(trg_fired1 && (strstr(variab1,jethlt_name[3]) && (strlen(variab1)-strlen(jethlt_name[3])<5))){
	  //if (trigRes->accept(ij) && isInEtaRange[iet]) {hlt_dijetprob[2][iet]->Fill(aveleadingpt);}
	   if (trigRes->accept(ij)) {
          hlt_dijetprob[2][iet]->Fill(aveleadingpt);
	  if(iet==0) cout <<"Prob name "  << variab1 << endl;
         }
	}   
      }    
    }
  }
*/

  //Calcualte Trigger Efficiency for dijet events
  bool trg_prev=false;

  //   if (!isMC) {
  for (int jk=-1; jk<nHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str(); 
      if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) || 
	  (jk>=0 && strstr(variab1,jethlt_name[jk]) && strlen(variab1)-strlen(jethlt_name[jk])<5)) {
	
	//const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltConfig_.prescaleValuesInDetail(iEvent,iSetup, variab1));
	const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltPrescaleProvider_.prescaleValuesInDetail(iEvent,iSetup,variab1));
	if (jk>=0) { 
          //cout<<variab1<<endl;
	  //==============================================================================
	  // double tmpp1= prescalesInDetail.first[0].second;
	  // double tmpp2 = prescalesInDetail.first[1].second;
	  
	  // l1pres[jk] =min(tmpp1, tmpp2);
	  //=====================================================================================
	  l1pres[jk] = prescalesInDetail.first[0].second;
	  
	 // if (jk>=3 && l1pres[jk]>1) { l1pres[jk]=1.0;}
	 if(l1pres[jk]<=0){l1pres[jk]=1.0;}

	  
          hltpres[jk] = prescalesInDetail.second;	  
	  
	  //compres[jk] = (l1pres[jk])*(triggerPrescales->getPrescaleForIndex(ij)); 
	  //compres[jk] = triggerPrescales->getPrescaleForIndex(ij);
	  compres[jk] = (l1pres[jk])*(hltpres[jk]);

	  //
	  //cout<<"Run NO= "<< iEvent.id().run()<<" ; Event No = "<< iEvent.id().event()<< " ; ilumi = " << iEvent.luminosityBlock() << 
	  //	" ; ibunch = " << iEvent.bunchCrossing()<<" ; L1 Pres0 = " << l1pres[jk] <<" "<<
	  //            " ; HLT Path= "<<name <<" ; HLT Pres = " <<hltpres[jk]<<" ; compres ="<<compres[jk] <<"; irecoht = "<< irecoht <<"; Pt=" <<aveleadingpt<<endl;
	  if (trigRes->accept(ij)) {trgpas[jk] = true;} // ihltfill = jk;}
	  
	  //if (trg_prev && compres[jk]>0.99) {
	  if (trg_prev){
	    for (int iet=0; iet<njetetamn; iet++) {
	      if (isInEtaRange[iet]) { 
		hlt_dijettag[jk][iet]->Fill(aveleadingpt,compres[jk]);
		if (trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt, compres[jk]);} //{, (isMC) ? 1.0 : compres[jk]);}
	      }
	    }
	  }
	  /*for (int iet=0; iet<njetetamn; iet++) {
	    if (isInEtaRange[iet]) { 
	    if(trg_prev) hlt_dijettag[jk][iet]->Fill(aveleadingpt);
	    if (trg_prev && trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt);} 
	    }    
	    }*/
	  // if (trg_prev) cout << "Accept =" << " name = " <<name <<endl;
	  trg_prev = trigRes->accept(ij);
	  //	  trg_prev = trg_prev|trigRes->accept(ij);
	  //	  if (!trg_prev) { trg_prev = trigRes->accept(ij);}
	  break;
	} else {
	  trg_prev = trigRes->accept(ij);
	  break;
	}
      }
    }
  }
#endif
  // cout<<"ihltfill "<<ihltfill<<endl;
  
  //	cout<<"3 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
  
  //	if ((irecoht <0 || irecoht >=nHLTmx) || ((!isMC) && (!trgpas[irecoht]))) return; //GMA remopve this condition
  //  cout <<"irecoht = "<<irecoht<<endl;
  //  if (irecoht==-3) return;
#ifdef TRIGGER
  if (irecoht>=0 && ((!isMC) && (!trgpas[irecoht]))) return;
  if (irecoht==-2 && ((!isMC) && (!trgpas[0]))) return;
#endif
  
  if (!isMC) {
    if (irecoht>=0) {
      wtfact = compres[irecoht];
    } else if (irecoht==-2) {
      wtfact = compres[0];
    } else {
      return ;
    }
  }
  //  for (int ij=0; ij<nHLTmx; ij++) {lumiwt[ij]=intlumi[nHLTmx-1]/intlumi[ij];}// cout<<"nt "<<datpileup[ij][0]<<endl;}
  if (isMC) {
#ifndef GENPART
    //Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    //iEvent.getByLabel("addPileupInfo", PupInfo);
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileup_, PupInfo);
    int npu = -1;
//    int tnpv  = -1;
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    if (PupInfo.isValid()) {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	if (PVI->getBunchCrossing()==0) {
//	  npu = PVI->getPU_NumInteractions();
	  npu = PVI->getTrueNumInteractions();
//	  tnpv  = PVI->getTrueNumInteractions();
	  break;
	}
      }
    }
    // double MyWeight = LumiWeights_->weight(npu);
    
    
    // cout << "Main weight = " <<MyWeight << endl;
    //double TotalWeight_plus = MyWeight*PShiftUp_.ShiftWeight( npu );
//double TotalWeight_plus = PShiftUp_.ShiftWeight( npu );
//double TotalWeight_minus = PShiftDown_.ShiftWeight( npu ); 

//cout << "Plus " << wtfact*TotalWeight_plus << " Mi = " << endl;
//cout << "wt= " <<  wtfact << " : weightmi" <<wtfact*TotalWeight_minus << " Mi = " << endl;
 //wtfact=wtfact*TotalWeight_plus; 
// wtfact=wtfact-TotalWeight_minus; 
  //  cout << "npu Number of interactions : " << npu << endl; 
//    cout << "tnpv Number of true interactions : " << tnpv << endl; 
    if (npu<0) return; //GMA  
    if (isFlat) {
      weight =weight2*wtfact; // for flat MC sample
    } else {
      weight =weight2;
    }
#endif
    defweight = weight;

#ifndef GENPART
        int tmprecht = (irecoht>=0) ? irecoht : 0; //GMA
    
    if (npu<npileupmx) {
          weight *=rat_pileup[tmprecht][npu]; //GMA
    } else {
            weight *=rat_pileup[tmprecht][npileupmx-1]; //GMA
    }
#endif
    
    weighttrg = weight;
    //    cout <<"weight  "<<weight<<" "<< weight2<<endl;
    //sar 3D PU reweighting 111028
  } else {
    weight = weight2;
    defweight = weight2;
    weighttrg = weight*wtfact; // *lumiwt[irecoht];
    //    weighttrg = weight*lumiwt[3];
    // cout <<"TEST2  weighttrg "<< weighttrg<<" ; weight "<<weight<<" ; "<< wtfact<<endl;
  }
  //=====================================
#ifndef GENPART
  if(!isMC){
    reco::TrackBase::Point beamPoint(0,0, 0);
    // math::XYZPoint beamPoint(0,0, 0); 
    
    edm::Handle<reco::BeamSpot> beamSpotH;
    iEvent.getByToken(beamSpot_,beamSpotH);
    if (beamSpotH.isValid()){
      beamPoint = beamSpotH->position();
    }
    
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vtxToken_, primaryVertices);
    
    int tmpvert=0;
    nprim=0;
    if (primaryVertices.isValid()) {
      tmpvert = primaryVertices->size();
      //cout<<"temp"<<tmpvert<<endl;
      for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
	int isel = (vert->isValid() && !vert->isFake()) ? 1 : 0;
	int ngoodtrk = 0;
	int nseltrk = 0;
	double prob = ChiSquaredProbability(vert->chi2(),vert->ndof());
	for (reco::Vertex::trackRef_iterator reftrk =vert->tracks_begin(); reftrk<vert->tracks_end(); reftrk++) {
	  if ((*reftrk)->quality(TrackBase::highPurity) && vert->trackWeight(*reftrk)>0) {
	    ngoodtrk++; 
	    if ((*reftrk)->normalizedChi2()<100000 && 
		abs((*reftrk)->dxy()) < 10000 && 
		(*reftrk)->pt() >0.50) {nseltrk++; } 
	  }
	}
	prim_alltrk[isel]->Fill(vert->tracksSize());
	prim_goodtrk[isel]->Fill(ngoodtrk);
	prim_seltrk[isel]->Fill(nseltrk);
	prim_dx[isel]->Fill(vert->position().x() - beamPoint.x());
	prim_dy[isel]->Fill(vert->position().y() - beamPoint.y());
	prim_dxy[isel]->Fill(vert->position().x() - beamPoint.x(), vert->position().y() - beamPoint.y());
	prim_dz[isel]->Fill(vert->position().z() - beamPoint.z());
	prim_prob[isel]->Fill(max(-20.0, log10(prob)));
	
	if (isel==1 && nprim < nprimx-1) {
	  primpr[nprim] = prob;
	  ntkpm[nprim] = 1000*(1000*min(int(vert->tracksSize()),999) + min(ngoodtrk,999)) + min(999, nseltrk);
	  nprim++;
	}
      }
    }
    
    prim_hist[0]->Fill(tmpvert);
    prim_sel[0]->Fill(nprim);
    
    prim_hist_rewt[0]->Fill(tmpvert, weighttrg);
    prim_sel_rewt[0]->Fill(nprim, weighttrg);

    if (irecoht>=0 && irecoht<nHLTmx) { 
      prim_hist[irecoht]->Fill(tmpvert);
      prim_sel[irecoht]->Fill(nprim);
      
      prim_hist_rewt[irecoht]->Fill(tmpvert, weighttrg);
      prim_sel_rewt[irecoht]->Fill(nprim, weighttrg);   
    }
    prim_correl->Fill(tmpvert, nprim);
  } 
#endif 
  //	cout<<"2 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl; 
  
  
  
  vector<double> jetptx[njecmx];
  vector<double> jetscl[njecmx];
  vector<int> jetindx[njecmx];

#ifndef GENPART
  if (ak4PFJets.isValid()) { 
    for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
      double pt = (*ak4PFJets)[ijet].pt();
      
      //#ifndef JETENERGY
      //#ifdef JETRESO
      
#if defined(JETRESO)&&(!defined(JETENERGY))
      // resolution file 
      JME::JetResolution resolution;
    //resolution = JME::JetResolution("Fall17_V2_MC_PtResolution_AK4PFchs.txt");
      resolution = JME::JetResolution("Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt");
    //resolution = JME::JetResolution("Fall17_V3b_MC_PtResolution_AK4PFchs.txt");
    //resolution = JME::JetResolution("Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt");
 
      // Scalefactor file
      JME::JetResolutionScaleFactor res_sf;
      //		cout<<"Filename="<<scalefile<<endl;
    res_sf = JME::JetResolutionScaleFactor("Summer19UL17_JRV2_MC_SF_AK4PFchs.txt");
    //res_sf = JME::JetResolutionScaleFactor("Fall17_V3b_MC_SF_AK4PFchs.txt");
    //res_sf = JME::JetResolutionScaleFactor("Fall17_V2_MC_SF_AK4PFchs.txt");
    //res_sf = JME::JetResolutionScaleFactor("Fall15_25nsV2_MC_SF_AK4PFchs.txt");
      
      edm::Handle<double> rho;
      iEvent.getByToken(m_rho_token, rho);
      //		cout<< "  rho=" << *rho << endl;
      
      //      cout << "Write test 3 = ok " << endl;
      double eta = (*ak4PFJets)[ijet].eta();
      double reso = 1;
      JME::JetParameters parameters_5 = {{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, *rho}};
      float rp = resolution.getResolution(parameters_5);
      float sf = res_sf.getScaleFactor({{JME::Binning::JetEta, eta}});
      float sf_up= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::UP);
      float sf_dn= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::DOWN);
      //#endif
      //#endif
#endif
      
      for (int isrc = 0; isrc < njecmx; isrc++) {
	double sup = 1;
#ifdef JETENERGY
	double eta = (*ak4PFJets)[ijet].eta();
	if (isrc>0 && isrc<=nsrc) {
	  JetCorrectionUncertainty *jecUnc = vsrc[isrc-1];
	  jecUnc->setJetEta(eta);
	  jecUnc->setJetPt(pt);
	  
	  sup += jecUnc->getUncertainty(true);
	} else if (isrc>nsrc) {
	  JetCorrectionUncertainty *jecUnc = vsrc[isrc-nsrc-1];
	  jecUnc->setJetEta(eta);
	  jecUnc->setJetPt(pt);
	  sup -= jecUnc->getUncertainty(false);
	}
#elif defined(JETRESO)
	if (isrc==0) {  
	  reso = sqrt(sf*sf - 1)*rp;
	} else if (isrc==1) {
	  reso = sqrt(sf_up*sf_up - 1)*rp;
	} else if (isrc==2) {
	  reso = sqrt(sf_dn*sf_dn - 1)*rp;
	}
	sup = gRandom->Gaus(1.0, reso);			
#endif
	jetptx[isrc].push_back(sup*pt);
	jetscl[isrc].push_back(sup);
	jetindx[isrc].push_back(ijet);
      }
    }
    //#if defined(JETENERGY)||defined(JETRESO)
    
    for (int isrc = 0; isrc < njecmx; isrc++) {
      for (unsigned int ij=0; ij<jetptx[isrc].size()-1; ij++) {
	for (unsigned int jk=ij+1; jk<jetptx[isrc].size(); jk++) {
	  if (jetptx[isrc][jk]>jetptx[isrc][ij]) {
	    double tmppt = jetptx[isrc][ij];
	    double tmpscl = jetscl[isrc][ij];
	    int tmpindx = jetindx[isrc][ij];
	    
	    jetptx[isrc][ij] = jetptx[isrc][jk];
	    jetscl[isrc][ij] = jetscl[isrc][jk];
	    jetindx[isrc][ij] = jetindx[isrc][jk];					
	    
	    jetptx[isrc][jk] = tmppt;
	    jetscl[isrc][jk] = tmpscl;
	    jetindx[isrc][jk] = tmpindx;
	  }
	}
      }
    }
    //#endif
    
    //    cout<<"1 aveleadingpt "<<endl; //aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
 //   double aveleadingptjec[njecmx] ={0};
    for (int isrc = 0; isrc < njecmx; isrc++) {
      if (jetptx[isrc].size()>=2) {
	aveleadingptjec[isrc] = 0.5*(jetptx[isrc][0] + jetptx[isrc][1]);
	irecohtjec[isrc] = getbinid(aveleadingptjec[isrc], nHLTmx, leadingPtThreshold);
      } else {
	irecohtjec[isrc] = -1;
      }
    }
   
 
    //GMA Need the corection on aveleadingpt
    if (ak4PFJets.isValid() && ak4PFJets->size() >=2) { //  && aveleadingpt >leadingPtThreshold[0]) { //GMA look on this
      
      for (int iet=0; iet<njetetamn; iet++) {
	for (int isrc = 0; isrc < njecmx; isrc++) {
	  if (aveleadingptjec[isrc] >leadingPtThreshold[0]) { 
	    int njets=0;
	    ncount=0;
	    //recterm=0;
	   // ithird=-1;
	  //  double sup = 1;	
	   // px=0;
	   // py=0;
	  //  ptxy=0;
	    tmpjt4v.clear();
	    tmpcand4v.clear();
	    tmpgen4v.clear();
	    
	    //				if (abs((*ak4PFJets)[0].eta())<etarange[iet] && abs((*ak4PFJets)[1].eta())<etarange[iet]) {
	    //					for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
	    
	    for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
	      if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet]) {
		int ireorjt = jetindx[isrc][ijet];
		
		double pt = jetptx[isrc][ijet];
		double sup = jetscl[isrc][ijet];
		double abseta = abs((*ak4PFJets)[ireorjt].eta());
	        if (pt<30.0 || abseta >etarange[iet]) continue;	
		//							if (iet==0 && isrc==0) cout <<"pteta "<<pt<<" "<<abseta<<endl;
		bool isEta = (abseta<2.4) ? true : false;
		
		if (isEta && pt>30.0) { njets++;}
		if (abseta>5.0) continue;
		bool isPt = (pt>30.0) ? true : false;
		if (isEta && isPt) {ncount++;}
		
		//cout<< "ncount = " << ncount << endl;
		//Jet ID ================= Tight ID 2017 Recomendation  check for 2018
		double NHF = (*ak4PFJets)[ireorjt].neutralHadronEnergyFraction();
		double NEMF = (*ak4PFJets)[ireorjt].neutralEmEnergyFraction();
		double CHF = (*ak4PFJets)[ireorjt].chargedHadronEnergyFraction();
//		double MUF = (*ak4PFJets)[ireorjt].muonEnergyFraction();
//		double CEMF = (*ak4PFJets)[ireorjt].chargedEmEnergyFraction();
		int NumConst = (*ak4PFJets)[ireorjt].chargedMultiplicity()+(*ak4PFJets)[ireorjt].neutralMultiplicity();
		//int NumNeutralParticles =(*ak4PFJets)[ireorjt].neutralMultiplicity();
		int CHM = (*ak4PFJets)[ireorjt].chargedMultiplicity();
                //cout<<"NHF== "<< NHF << "; NEF== " << NEMF <<" ; CHF==" <<CHF <<" ;cef==" << CEMF <<"; no= " << NumConst <<" ; nch==" << CHM <<" ; NO of part==" << NumNeutralParticles <<endl;
		bool TightJetID =false;
		//if(abs((*ak4PFJets)[ireorjt].eta())<=2.7){       //Update for ReReco 2017 
                //if( (NHF<0.90 && NEMF<0.90 && NumConst>1 ) && (abs((*ak4PFJets)[ireorjt].eta())<=2.4 && CHF>0 && CHM>0 ) ) TightJetID =true;
                //       } else {
                //                           TightJetID =false; }
                //Update for UL17 : 27Aug20
                if(abs((*ak4PFJets)[ireorjt].eta())<=2.6){
                if(NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0)  TightJetID =true;
                      } else {
                                 TightJetID =false;
                                 }
                
                 if (abs((*ak4PFJets)[ireorjt].eta())>2.6) {TightJetID = false;}
                 if ((*ak4PFJets)[ireorjt].pt()<30.0) {TightJetID = false;}
		
		if( ireorjt<=1 && !TightJetID) break;
		if (!TightJetID) continue;
		//JetID ================
	/*	
		if (ncount <=2 && ncount !=ijet+1) {
		  for (int ix=0; ix<ntype; ix++) { 
		    recomom[isrc][ix][iet].clear(); 
		  }
		  break;
		}*/
		
		//if (isrc==0 && iet==0) {
		//cout <<"recomom[isrc][0][iet].size() "<<iet<<" "<<isrc<<" "<<ijet<<" "<<ncount<<" "<<recomom[isrc][0][iet].size()<<" "<<recomom[isrc][1][iet].size()<<" "<<ncount<<" "<<ijet<<endl;
		//						}
		
		HepLorentzVector tmp4v((*ak4PFJets)[ireorjt].px(), (*ak4PFJets)[ireorjt].py(), (*ak4PFJets)[ireorjt].pz(), (*ak4PFJets)[ireorjt].energy());
		
		tmp4v *=sup;
		/*double respfact=respfun(1.02, 0.000004799, 0.000000007044,tmp4v.perp()); 
		  cout <<"Response factor = " <<respfact << endl;
		  tmp4v /=respfact;
		*/
		//cout << "Pt before correction = "<< tmp4v.perp()<< endl;
		/*double  respfact=0.;
		  bool isCorrect=false;   					
		  for (int iresp=0; iresp<7; iresp++){
		  if(abs(tmp4v.eta())> resetarange[iresp] && abs(tmp4v.eta())<resetarange[iresp+1]){
		  respfact=respfun(par0[iresp], par1[iresp], par2[iresp], tmp4v.perp());
		  isCorrect =true;
		  }
		  //		cout << "iresp = "<< iresp << " Eta = " <<tmp4v.eta() <<endl;
		  if(isCorrect) break;
		  }
		  //       cout <<"Response factor = " <<respfact << endl;
		  double invrespfact=0;
		  if (respfact!=0) invrespfact =1/respfact;
		  tmp4v*=invrespfact;*/
		//	cout <<"Response factor = " <<respfact << "Inv Response factor =" << invrespfact << " Corrected Pt = "<< tmp4v.perp() <<endl;
		
		//						cout <<"perp "<<iet<<" "<<isrc<<" "<<ijet<<" "<<sup<<" "<<tmp4v<<" "<<tmp4v.eta()<<" "<<tmp4v.perp()<<" "<<pt<<endl;
		
		if (isEta && isPt) { tmpjt4v.push_back(tmp4v);}
		//tmpjt4v.push_back(tmp4v);	  
		//	 if (isEta && isPt) {allrecojetmom.push_back(tmp4v);}
	//	if (ncount<=2) {  //change for all jet 26th June
		  if (isEta && isPt) {
		    recomom[isrc][0][iet].push_back(tmp4v);
		  }
		  //}
		  //								cout <<"ncount filled "<<ncount<<" "<<isrc<<" "<<iet<<" "<<recomom[isrc][0][iet].size()<<endl;
	//	  px +=tmp4v.px();
	//	  py +=tmp4v.py();
	//	  ptxy +=tmp4v.perp();
		  if (isrc==0) { 
		    if ((isInEtaRange[iet])) {recojt_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt_phi->Fill(tmp4v.phi(), weighttrg);}
		    if (isEta && ncount==1) {recoht2_pt[iet]->Fill(aveleadingpt,weighttrg);}
		  }
		//} else {
		/*  if (isrc==0) { 
		    if ((isInEtaRange[iet])) {recojt_oth_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt_oth_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt_oth_phi->Fill(tmp4v.phi(), weighttrg);}
		  }*/
		/*  if (isEta && isPt) {
		    double tmppx = px + tmp4v.px();
		    double tmppy = py + tmp4v.py();
		    double tmppt = ptxy + tmp4v.perp();
		    double tmprec = sqrt(pow(tmppx, 2)+pow(tmppy, 2))/tmppt;
		    
		    if (tmprec>recterm) {
		      recterm = tmprec;
		      ithird = ireorjt;
		      //cout <<"ithird Data : "<< ijet<<endl;
		    }
		  }*/
		//}
		if (isrc==0) { 
		  if(ijet==0) { 
		    if (isInEtaRange[iet]) {recojt1_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt1_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt1_phi->Fill(tmp4v.phi(), weighttrg);}
		  } else if(ijet==1){
		    if (isInEtaRange[iet]) {recojt2_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt2_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {recojt2_phi->Fill(tmp4v.phi(), weighttrg);}
		    if (isInEtaRange[iet] && ncount==2) { 
		      if (irecoht>=0 && irecoht<nHLTmx) { 
			recojtave_pt[iet][irecoht]->Fill(aveleadingpt, weighttrg);
			recojtavewt1_pt[iet][irecoht]->Fill(aveleadingpt);
		      }
		      
		      recojtallavewt1_pt[iet]->Fill(aveleadingpt);
		      recojtallave_pt[iet]->Fill(aveleadingpt, weighttrg);
		    }
		    
		  } else if(ijet==2) {
		    if (isInEtaRange[iet]) {recojt3_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0 ) {recojt3_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {recojt3_phi->Fill(tmp4v.phi(), weighttrg);}
		  }
		  
		  if (tmpjt4v.size()==2 && isInEtaRange[iet]) { 
		    double dphi = dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi());
		    double dpt = tmpjt4v[0].perp() - tmpjt4v[1].perp();
		    double dperp = fabs(tmpjt4v[1].perp()*sin(dphi))/tmpjt4v[0].perp();
		    hjetdphi[iet]->Fill(dphi, weighttrg);
		    hjetdpt[iet]->Fill(dpt, weighttrg);
		    hjetptbypl[iet]->Fill(dperp, weighttrg);
		    hjetpt2bypt1[iet]->Fill(tmpjt4v[1].perp()/tmpjt4v[0].perp(), weighttrg);
		  }
		  
		  if (tmpjt4v.size()==3) {hjetpt3bypt2[iet]->Fill(tmpjt4v[2].perp()/tmpjt4v[1].perp(), weighttrg);}
		 } //if (isrc==0) {
		
	        int nchg=0;	
		std::vector<reco::CandidatePtr> daus((*ak4PFJets)[ireorjt].daughterPtrVector());           
		std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });                                                                                                  
		for (unsigned int i2 = 0; i2< daus.size(); ++i2) {   
		  const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
		  int charge = pfcand.charge();
		  HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
		  tmpcand4v.push_back(cand4v);	
                  nchg++;
                  h_nchg[iet]->Fill(nchg, weighttrg);
                  
		  //	   if (cand4v.perp()<0.5) continue;
		//  if (ncount<=2 && isEta && isPt) { 
		    //recomom[isrc][1][iet].push_back(cand4v);
		    
		    if (charge !=0) {
		      recomom[isrc][1][iet].push_back(cand4v);
#ifdef TRACKSYS
		      if (gRandom->Uniform() < 0.96) {recomom[isrc][2][iet].push_back(cand4v); }
#endif
		    }
		   /* if (charge==0) { //other option if need open
		      if (cand4v.perp()>1.0) {
			recomom[isrc][3][iet].push_back(cand4v);
		      }  
		    } else {
		      if (cand4v.perp()>0.5) {
			recomom[isrc][3][iet].push_back(cand4v);
		      }
		    }*/
		    
		    
		    
		    
		    //   double dphi = dPhi(recomom[0][0][0].phi(), recomom[0][0][1].phi());
		    // double dpt = recomom[0][0][0].perp() - recomom[0][0][1].perp();
		    //    double dperp = fabs(tmpcand4v[1].perp()*sin(dphi))/tmpjt4v[0].perp();
		    
		    //	if (dpt<0) cout <<" "<< jk<<" "<<ij<<" "<<mn<<" "<<seljtvar4v[0]<<" "<<seljtvar4v[1]<<" "<<seljtvar4v[0].perp()<<" "<<seljtvar4v[1].perp()<<" "<<endl;
		    
		    //    hjet1dphi->Fill(dphi, weighttrg);
		    //    hjet1dpt->Fill(dpt, weighttrg);
		 // }
		  if (isrc==0) { 
		    //if (isEta && isPt) {
		      if (charge !=0) {
			recochg_phi->Fill(cand4v.phi(), weighttrg);
			recochg_pt->Fill(cand4v.perp(), weighttrg);
			recochg_eta->Fill(cand4v.eta(), weighttrg);
		      }
		    
		      if (ijet==0 && charge !=0) {
                        recochg1_phi->Fill(cand4v.phi(), weighttrg);
                        recochg1_pt->Fill(cand4v.perp(), weighttrg);
                        recochg1_eta->Fill(cand4v.eta(), weighttrg);
		      }
                     else if (ijet==1 && charge !=0) {
                        recochg2_phi->Fill(cand4v.phi(), weighttrg);
                        recochg2_pt->Fill(cand4v.perp(), weighttrg);
                        recochg2_eta->Fill(cand4v.eta(), weighttrg);
                      }
                     else if (ijet==2 && charge !=0) {
                        recochg3_phi->Fill(cand4v.phi(), weighttrg);
                        recochg3_pt->Fill(cand4v.perp(), weighttrg);
                        recochg3_eta->Fill(cand4v.eta(), weighttrg);
                      }                     

		    }//if (isrc==0) {
		//  }
		} //for (unsigned int i2 = 0; i2< daus.size(); ++i2)
		//  if (isEta && isPt) {ncount++;}
	   //   } //if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet])
	  //  } // for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++)
/*	    if (ithird>=0) {
	      
	      recomom[isrc][0][iet].push_back(tmp4v);
	      //					cout <<"recomom[isrc][0][iet] "<< isrc<<" "<<iet<<" "<<recomom[isrc][0][iet].size()<<endl;
	      // tmpjt4v.push_back(tmp4v);   
	      
	      std::vector<reco::CandidatePtr> daus((*ak4PFJets)[ithird].daughterPtrVector());
	      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2 ->pt(); });
	      for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
		int charge = pfcand.charge();
		HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
		//      if (cand4v.perp()<0.5) continue;                                                             
		recomom[isrc][1][iet].push_back(cand4v);
		
		if (charge !=0) {
		  recomom[isrc][2][iet].push_back(cand4v);
#ifdef TRACKSYS
		  if (gRandom->Uniform() < 0.96) {recomom[isrc][4][iet].push_back(cand4v); }
#endif
		  
		}
		if (charge==0){
		  if (cand4v.perp()>1.0) {
		    recomom[isrc][3][iet].push_back(cand4v);
		  }
		} else{
		  if (cand4v.perp()>0.5) {
		    recomom[isrc][3][iet].push_back(cand4v);
		  }
		}
	      }
	    }*/ //if (ithird>=0) 
	    //if (isrc==0) {h_njets[iet]->Fill(ncount, weighttrg);}
	    h_njets[iet]->Fill(ncount, weighttrg);
              } //if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet])
            } // for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++)
	  } //if (aveleadingptjec[isrc] >leadingPtThreshold[0])
	} // 	for (int isrc = 0; isrc < njecmx; isrc++)
      } //for (int iet=0; iet<njetetamn; iet++)	   
    } // if (ak4PFJets.isValid() && ak4PFJets->size()>=2 && (*ak4PFJets)[0].pt()>leadingPtThreshold[0])
  } // if (ak4PFJets.isValid())
#endif
  
  //  cout << "Write test 31 = ok " << endl;
  //===================********Trigger****============================================================
  
  //t2=clock();
  //float diff ((float)t2-(float)t1);
  //if(diff>30000) return;
  //cout << "Time T2 = " << t2 << " ;Time Diff to Run ="<< diff << endl;
/*  
#ifdef TRIGGER
  if(!isMC){
#ifndef DIJETAVE
    //  vector<triggervar> alltrgobj;
    if (trigRes.isValid() && isReconstruct  &&
	(tmpjt4v.size() ==2 || (tmpjt4v.size()>=3 && tmpjt4v[2].perp()<30.0)) &&
	abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0){
      
      //  if (trigRes.isValid() && isReconstruct  &&
      //   (tmpjt4v.size() ==2) && abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0){
      
      int ijet = int(2*gRandom->Uniform())%2;
      int ijet2 = (ijet==0) ? 1 : 0;
      //cout <<"gRandom= "<<gRandom->Uniform() << " ijet" <<ijet<<endl; 
      HepLorentzVector tagjet4v = tmpjt4v[ijet];
      HepLorentzVector probjet4v = tmpjt4v[ijet2];

      vector<triggervar> alltrgobj;
      alltrgobj.clear(); 
      const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	obj.unpackPathNames(names);
	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	for (unsigned ih = 0, n = pathNamesAll.size(); ih < n; ++ih) {
	  variab2 = pathNamesAll[ih].c_str(); 
	  for (int jk=0; jk<nHLTmx; jk++) {
	    if (strstr(variab2,jethlt_name[jk]) && strlen(variab2)-strlen(jethlt_name[jk])<5){
	      triggervar tmpvec;
	      if( obj.pt()<jethlt_thr[jk] ) continue;
	      tmpvec.both = obj.hasPathName( pathNamesAll[ih], true, true );
	      if(obj.pt()>10){
		tmpvec.both = obj.hasPathName( pathNamesAll[ih], true, true );
		tmpvec.highl  = obj.hasPathName( pathNamesAll[ih], false, true );
		tmpvec.level1 = obj.hasPathName( pathNamesAll[ih], true, false );
		tmpvec.trg4v = HepLorentzVector(obj.px(), obj.py(), obj.pz(), obj.energy());
		tmpvec.prescl = 1;
		tmpvec.ihlt = jk;
		alltrgobj.push_back(tmpvec);
	      }
	    }
	  }
	}
      }
      
      for (unsigned ij=0; ij<alltrgobj.size(); ij++) {
	HepLorentzVector trigger4v = alltrgobj[ij].trg4v;
	int ihlt = -1;
	int tmphlt = alltrgobj[ij].ihlt;
	if( trigger4v.perp()<jethlt_thr[tmphlt] ) continue;
	//      bool isBoth=alltrgobj[ij].both;
	bool isLF =alltrgobj[ij].level1;
	bool isL3 =alltrgobj[ij].highl;
	double angle = deltaR(tagjet4v, trigger4v);
	if (isLF) { 
	  trgjet_angle[tmphlt][0]->Fill(angle);
	  trgjet_2dangle[tmphlt][0]->Fill(trigger4v.perp(), angle);	
	}
	if (isL3) { 
	  trgjet_angle[tmphlt][1]->Fill(angle);
	  trgjet_2dangle[tmphlt][1]->Fill(trigger4v.perp(), angle);		
	}	
	// bool tag=false;
	if (deltaR(tagjet4v, trigger4v)<0.2) {
	  // tag=true;
	  ihlt = alltrgobj[ij].ihlt;
	  if (isLF)  {
	    //        if (isLF && !isBoth)  {
	    trgjet_pt[ihlt][0]->Fill(probjet4v.perp());
	    trgjet_eta[ihlt][0]->Fill(probjet4v.eta());
	    trgjet_phi[ihlt][0]->Fill(probjet4v.phi());
	  }
	  if (isL3) {
	    // if (isLF && !isBoth)  {
	    trgjet_pt[ihlt][1]->Fill(probjet4v.perp());
	    trgjet_eta[ihlt][1]->Fill(probjet4v.eta());
	    trgjet_phi[ihlt][1]->Fill(probjet4v.phi());
	  }
	  
	  for (unsigned jk=0; jk<alltrgobj.size(); jk++) {
	    if (ij==jk || alltrgobj[jk].ihlt !=ihlt) continue;
	    //        if( trigprbjet4v.perp()<jethlt_thr[ihlt] ) continue;
	    HepLorentzVector trigprbjet4v = alltrgobj[jk].trg4v;
	    if( trigprbjet4v.perp()<jethlt_thr[ihlt] ) continue;
	    double angle1 = deltaR(probjet4v, trigprbjet4v);
	    if (isLF && angle1<0.5 ) {
	      //          if (isLF && !isBoth && angle1<0.5 ) {
	      prbjet_pt[ihlt][0]->Fill(probjet4v.perp());
	      prbjet_eta[ihlt][0]->Fill(probjet4v.eta());
	      prbjet_phi[ihlt][0]->Fill(probjet4v.phi());
	      isLF = false;
	    }
	    
	    if (isL3 && angle1<0.2) {
	      //          if (isL3 && !isBoth && angle1<0.2)
	      prbjet_pt[ihlt][1]->Fill(probjet4v.perp());
	      prbjet_eta[ihlt][1]->Fill(probjet4v.eta());
	      prbjet_phi[ihlt][1]->Fill(probjet4v.phi());
	      isL3 = false;
	    }
	    
	    if ((!isL3) && (!isLF)) continue;
	    
	  } //for (unsigned jk=0; jk<alltrgobj.size(); jk++) 
	} //if (deltaR(tagjet4v, trigger4v)<0.2)
	if (ihlt>=0) continue;
      } //for (int ij=0; ij<alltrgobj.size(); ij++)
    } // if (trigRes.isValid() && m_trigeff && isReconstruct  && tmpjt4v.size() ==2) &&
    //      abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0
#endif
  }
#endif*/
  //======******Trigger Efficiency Normal=======================
  
  //cout << "Write test 1 = ok " << endl;
  //==================================***GenJets*****=================================
  //	cout<<"0 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
  if(isMC) {
    
    edm::Handle<reco::GenJetCollection> genjets;
    iEvent.getByToken(genjetToken_,genjets);
    
    double avegenpt =0;
    //    cout <<"HGebjet "<<endl;
    if (genjets.isValid() &&  genjets->size()>=2) {
#ifdef DIJETAVE
      for (int iet=0; iet<njetetamn; iet++) {
	isInEtaRange[iet] = true;
      }
      
      for (int ij=0; ij<2; ij++) {
	for (int iet=0; iet<njetetamn; iet++) {
	  if (abs((*genjets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
	}

	
	if (abs((*genjets)[ij].eta())<2.4 && (*genjets)[ij].pt()>30.0 ) { 
	  avegenpt +=(*genjets)[ij].pt();
	} else {
	  avegenpt -=100000;
	}
      }
      avegenpt /=2.0;
      
#else 

#endif
    }
    
    igenht = getbinid(avegenpt, njetptmn, leadingPtThreshold);

    
    //    cout << "Write test 2 = ok " << endl;
    //cout << "Write test 321 = ok " << endl;
    vector<double> genjetptx[nGenReso];
    vector<double> genjetscl[nGenReso];
    vector<int> genjetindx[nGenReso];
    
    for(unsigned ijet = 0; ijet != genjets->size(); ijet++) {
      double pt = (*genjets)[ijet].pt();
      //#ifdef JETRESO		
      //			double eta = (*genjets)[ijet].eta();
      //			double reso = 1;
      //			JME::JetParameters parameters_5 = {{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, *rho}};
      //			float rp = resolution.getResolution(parameters_5);
      //			float sf = res_sf.getScaleFactor({{JME::Binning::JetEta, eta}});
      //			float sf_up= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::UP);
      //			float sf_dn= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::DOWN);
      
      //#endif		
      for (int isrc = 0; isrc < nGenReso; isrc++) {
	double sup = 1.0;
	//#ifdef JETRESO
	//				if (isrc==0) {  
	//					reso = sqrt(sf*sf - 1)*rp;
	//				} else if (isrc==1) {
	//					reso = sqrt(sf_up*sf_up - 1)*rp;
	//				} else if (isrc==2) {
	//					reso = sqrt(sf_dn*sf_dn - 1)*rp;
	//				}
	
	//				sup = gRandom->Gaus(1.0, reso);
	//				//				cout<<"isrc "<< ijet<<" "<< pt<<" "<<eta<<" "<<isrc<<" rp "<<rp<<" "<<sf<<" "<<sf_dn<<" "<<sf_up<<" "<<reso<<" "<<sup<<endl;
	//#endif
	
	genjetptx[isrc].push_back(sup*pt);
	genjetscl[isrc].push_back(sup);
	genjetindx[isrc].push_back(ijet);
      }
    }
    
    //    cout << "Write test 3 = ok " << endl;
    //    cout << "Write test 322 = ok "<<nGenReso << endl;
    //////#ifdef JETRESO
    for (int isrc = 0; isrc < nGenReso; isrc++) {
      //      cout << "Write test 31 = ok "<<isrc << " ; " << genjetptx[isrc].size() <<endl;
      if(genjetptx[isrc].size()==0) break;
      for (unsigned int ij=0; ij<genjetptx[isrc].size()-1; ij++) {
	//cout << "Write test 32 = ok "<<nGenReso << endl;
	for (unsigned int jk=ij+1; jk<genjetptx[isrc].size(); jk++) {
	  
	  //    if(jk<genjetptx[isrc].size()) return;
	  //cout << "Write test 33 = ok "<<nGenReso << endl;
	  if (genjetptx[isrc][jk]>genjetptx[isrc][ij]) {
	    //cout << "Write test 34 = ok "<<nGenReso << endl;
	    double tmppt = genjetptx[isrc][ij];
	    double tmpscl = genjetscl[isrc][ij];
	    int tmpindx = genjetindx[isrc][ij];
	    
	    genjetptx[isrc][ij] = genjetptx[isrc][jk];
	    genjetscl[isrc][ij] = genjetscl[isrc][jk];
	    genjetindx[isrc][ij] = genjetindx[isrc][jk];					
	    
	    genjetptx[isrc][jk] = tmppt;
	    genjetscl[isrc][jk] = tmpscl;
	    genjetindx[isrc][jk] = tmpindx;
	    //	    cout << "Write test 35 = ok "<<nGenReso << endl;
	  }
	}
      }
    }
    //////#endif
    //    cout << "Write test 4 = ok " << endl;
   // double avegenptres[nGenReso]={0};
    
    for (int isrc = 0; isrc < nGenReso; isrc++) {
      if (genjetptx[isrc].size()>=2) {
	avegenptres[isrc] = 0.5*(genjetptx[isrc][0] + genjetptx[isrc][1]);
	igenhtres[isrc] = getbinid(avegenptres[isrc], njetptmn, leadingPtThreshold);
      } else {
	igenhtres[isrc] = -1;
      }
    }
    
    if(genjets.isValid() && genjets->size() >=2) { //  && avegenpt>leadingPtThreshold[0]) { 
      for (int iet=0; iet<njetetamn; iet++) {
	for (int isrc=0; isrc<nGenReso; isrc++) { 
	  if (avegenptres[isrc] > leadingPtThreshold[0]) {
	    //double px =0;
	    //double py =0;
	    //double ptxy =0;
	    
	    ncount=0;
	   //int recterm=0;
	   // int ithird=-1;
	    
	    for(unsigned ijet = 0; ijet < genjets->size(); ijet++) {
	      int igenjt = genjetindx[isrc][ijet];
	     /* if ((*genjets)[igenjt].pt()>25.0) {
		cout<<"ievt "<<ievt<<" "<<ijet<<" "<<igenjt<<" "<<genjetptx[isrc][ijet]<<" "<<(*genjets)[igenjt].pt()<<" "<<(*genjets)[igenjt].eta()<<" "<<(*genjets)[igenjt].phi()<<endl;
	      }*/

	      if (abs((*genjets)[genjetindx[isrc][0]].eta())<etarange[iet] && 
		  abs((*genjets)[genjetindx[isrc][1]].eta())<etarange[iet]) {
		
		
		double pt = genjetptx[isrc][ijet];
		double sup = genjetscl[isrc][ijet];
		double abseta = abs((*genjets)[igenjt].eta());
		if (pt<30.0 || abseta >etarange[iet]) continue;
		
		//								if (iet==0 && isrc==0) 
		//		cout <<"MC:pteta "<<ijet<<" "<<pt<<" "<<abseta<<endl;
		if (abseta>5.0) continue;
		bool isEta = (abseta<2.4) ? true : false;
		
		HepLorentzVector tmp4v((*genjets)[igenjt].px(), (*genjets)[igenjt].py(), (*genjets)[igenjt].pz(), (*genjets)[igenjt].energy());
		tmp4v *=sup;
		bool isPt = (pt>30.0) ? true : false;
		//Response 
		/*if(isPt && isReconstruct) {
		  for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
		    HepLorentzVector tmp4vreco((*ak4PFJets)[ijet].px(), (*ak4PFJets)[ijet].py(), (*ak4PFJets)[ijet].pz(), (*ak4PFJets)[ijet].energy());
		    bool isResp=false;
		    for (int iresp=0; iresp<7; iresp++){
		      bool isEtaMatch=false; 
		      if(abs(tmp4v.eta())> resetarange[iresp] && abs(tmp4v.eta())<resetarange[iresp+1]){
			isEtaMatch=true;
			// cout << "Eta Match = " << tmp4v.eta() << " iresp = "<< iresp << endl; 
			double respangle=deltaR(tmp4v,tmp4vreco);
			if(respangle <0.2) {
			  resp_jet[iresp]->Fill(tmp4v.perp(), tmp4vreco.perp()/tmp4v.perp(), weighttrg);
			  resp_jet1[iresp]->Fill(abs((tmp4v.perp()-tmp4vreco.perp())/tmp4v.perp()), weighttrg);
			  //		cout << "Resolution = " << abs((tmp4v.perp()-tmp4vreco.perp())/tmp4v.perp()) << endl;
			  isResp=true;
			}
		      }		
		      if(isEtaMatch) break;
		    }
		    if(isResp) break;
		  }
		}*/
		//Response
		
		//								cout <<"isrc "<<iet<<" "<< isrc <<" "<<igenjt<<" "<<pt<<" " <<tmp4v.perp()<<" "<<tmp4v.eta()<<endl;
		//							pt = tmp4v.perp();
		
		//bool isPt = (pt>30.0) ? true : false;
		if (isEta && isPt) {ncount++;}
		
		/*if (ncount <=2 && ncount !=ijet+1) {
		  for (int ix=0; ix<ntype; ix++) { 
		    genmom[isrc][ix][iet].clear(); 
		  }
		  break;
		}*/
		//							if (isrc==0 && iet==0) cout <<"	genmom[isrc][0][iet].size() "<<genmom[isrc][0][iet].size()<<" "<<ncount<<" "<<igenjt<<" "<<endl;
		
		if (isEta && isPt) { tmpgen4v.push_back(tmp4v);} 
		//if (ncount<=2) {
		  if (isEta && isPt) {
		    genmom[isrc][0][iet].push_back(tmp4v);
		   // cout <<"twoijx "<<isrc<<" "<<iet<<" "<< genmom[isrc][0][iet].size()<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<" "<<tmp4v.phi()<<endl;
		  }
		  //px +=tmp4v.px();
		 // py +=tmp4v.py();
		 // ptxy +=tmp4v.perp();
		  if (isrc==0) { 
		    if (isInEtaRange[iet]) {genjt_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt_phi->Fill(tmp4v.phi(), weighttrg);}
		    
		  }
		//} else {
		  /*if (isrc==0) { 
		    if (isInEtaRange[iet]) {genjt_oth_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt_oth_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt_oth_phi->Fill(tmp4v.phi(), weighttrg);}
		  }*/
		 /* if (isEta && isPt) {
		    double tmppx = px + tmp4v.px();
		    double tmppy = py + tmp4v.py();
		    double tmppt = ptxy + tmp4v.perp();
		    double tmprec = sqrt(pow(tmppx, 2)+pow(tmppy, 2))/tmppt;
		    
		    if (tmprec>recterm) {
		      recterm = tmprec;
		      ithird = igenjt;
		      // 		      cout <<"ithird MC : "<< igenjt<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<endl;
		    }
		  }*/
		//}
		if (isrc==0) { 
		  if(ijet==0) {
		    //		    cout<<"Gen Pt= " << avegenpt <<endl;
		    if (isInEtaRange[iet]) {genjt1_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt1_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt1_phi->Fill(tmp4v.phi(), weighttrg);}
		  } else if(ijet==1){
		    //		    cout<<"okkkkkkkk" <<endl;
		    if (isInEtaRange[iet]) {genjt2_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt2_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {genjt2_phi->Fill(tmp4v.phi(), weighttrg);}
		    if (isInEtaRange[iet] && ncount==2) {
		      //cout<<"Gen Pt 1= " << avegenpt <<endl;
		      genjtallave_pt[iet]->Fill(avegenpt, weighttrg);
		    }
		  } else if(ijet==2) {
		    if (isInEtaRange[iet]) {genjt3_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0 ) {genjt3_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {genjt3_phi->Fill(tmp4v.phi(), weighttrg);}
		  }
		  if (tmpgen4v.size()==2 && isInEtaRange[iet]) {
		    double dphi = dPhi(tmpgen4v[0].phi(), tmpgen4v[1].phi());
		    double dpt = tmpgen4v[0].perp() - tmpgen4v[1].perp();
		    double dperp = fabs(tmpgen4v[1].perp()*sin(dphi))/tmpgen4v[0].perp();
		    genjetdphi[iet]->Fill(dphi, weighttrg);
		    genjetdpt[iet]->Fill(dpt, weighttrg);
		    genjetptbypl[iet]->Fill(dperp, weight);
		    genjetpt2bypt1[iet]->Fill(tmpgen4v[1].perp()/tmpgen4v[0].perp(), weight);
		  }
		  
		  if (tmpgen4v.size()==3) {genjetpt3bypt2[iet]->Fill(tmpgen4v[2].perp()/tmpgen4v[1].perp(), weight);}
		}
#ifdef GENPART
		std::vector <const GenParticle*> daus ((*genjets)[igenjt].getGenConstituents ());
		//								std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); 
		
		for (unsigned int i2 =0; i2< daus.size(); ++i2) {
		  const GenParticle* pfcand = daus[i2];
		  int charge = pfcand->charge();
		  HepLorentzVector cand4v(pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy());
		  //									int pdgid = pfcand->pdgId();
		  
#else								
		  std::vector<reco::CandidatePtr> daus((*genjets)[igenjt].daughterPtrVector());
		  std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });                               
		  
		  for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		    const pat::PackedCandidate &pfcand = static_cast<const pat::PackedCandidate &>(*daus[i2]);
		    int charge = pfcand.charge();
		    
		    HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
#endif
		    //	    if (cand4v.perp()<0.5) continue;
		    
		    //if (ncount<=2 && isEta && isPt) {
		      //genmom[isrc][1][iet].push_back(cand4v);
		      if (charge !=0) {
			genmom[isrc][1][iet].push_back(cand4v);
#ifdef TRACKSYS
			if (gRandom->Uniform() < 0.96) {genmom[isrc][2][iet].push_back(cand4v); }
#endif
		      }
		      
		      //   if (charge ==0) {genmom[isrc][2][iet].push_back(cand4v);}
		      
		     /* if(charge ==0) {
			if (cand4v.perp()>1.0) {
			  genmom[isrc][3][iet].push_back(cand4v);
			}
		      } else {
			if (cand4v.perp()>0.5) {
			  genmom[isrc][3][iet].push_back(cand4v);
			}
		      }*/
		    //}
		    if (isrc==0) { 
			if (charge !=0) {
		      if (isEta && isPt) {
			  genchg_phi->Fill(cand4v.phi(), weighttrg);
			} 
                      if (isEta) {
                         genchg_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg_eta->Fill(cand4v.eta(), weighttrg);
                        }

                       if(ijet==0) {
                      if (isEta && isPt) {
                          genchg1_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg1_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg1_eta->Fill(cand4v.eta(), weighttrg);
                        }                      
                      }

                      if(ijet==1) {
                      if (isEta && isPt) {
                          genchg2_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg2_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg2_eta->Fill(cand4v.eta(), weighttrg);
                        }
                      }
                  
                      if(ijet==2) {
                      if (isEta && isPt) {
                          genchg3_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg3_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg3_eta->Fill(cand4v.eta(), weighttrg);
                        }
                      }

		      }//if (charge !=0) { 
                      
                    /*if (isEta) {
                        if (charge !=0) {
                          genchg_pt->Fill(tmp4v.perp(), weighttrg);
                        }
                    if (isPt) {
                        if (charge !=0) {
                          genchg_eta->Fill(tmp4v.eta(), weighttrg);
                        }                      

else {
			if (charge !=0) {
			  genchg_oth_phi->Fill(tmp4v.phi(), weighttrg);
			} else {
			  genneu_oth_phi->Fill(tmp4v.phi(), weighttrg);
			}
		      }
		      
		      if (isEta) {
			if (charge !=0) {
			  genchg_pt->Fill(tmp4v.perp(), weighttrg);
			} else {
			  genneu_pt->Fill(tmp4v.perp(), weighttrg);
			}
		      } else {
			if (charge !=0) {
			  genchg_oth_pt->Fill(tmp4v.perp(), weighttrg);
			} else {
			  genneu_oth_pt->Fill(tmp4v.perp(), weighttrg);
			}
		      }
		      if (isPt) {
			if (charge !=0) {
			  genchg_eta->Fill(tmp4v.eta(), weighttrg);
			} else {
			  genneu_eta->Fill(tmp4v.eta(), weighttrg);
			}
		      } else {
			if (charge !=0) {
			  genchg_oth_eta->Fill(tmp4v.eta(), weighttrg);
			} else {
			  genneu_oth_eta->Fill(tmp4v.eta(), weighttrg);
			}
		      }*/
		    } //if (isrc==0)
		  } //for (unsigned int i2 = 0; i2< daus.size(); ++i2)
		  //  if (isEta && isPt) {ncount++;}
		} // if (abs((*genjets)[genjetindx[isrc][0]].eta())<etarange[iet] && 
		//								abs((*genjets)[genjetindx[isrc][1]].eta())<etarange[iet])
	      } //	for(unsigned ijet = 0; ijet != genjets->size(); ijet++) 
	      //cout << "Write test 324 = ok " << endl;
	     /* if (ithird>=0) {
		//							cout <<"ithird "<<isrc<<" "<< iet<<" "<< ithird<<endl;
		
		HepLorentzVector tmp4v((*genjets)[ithird].px(), (*genjets)[ithird].py(), (*genjets)[ithird].pz(), (*genjets)[ithird].energy());
		genmom[isrc][0][iet].push_back(tmp4v);
		//cout <<"thirdijxxx "<<isrc<<" "<<iet<<" "<< genmom[isrc][0][iet].size()<<" "<<genjets->size()<<" "<<ithird<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<" "<<tmp4v.phi()<<" "<<setprecision(14)<<weighttrg<<endl;
#ifdef GENPART
		std::vector <const GenParticle*> daus ((*genjets)[ithird].getGenConstituents ());
		//								std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); 
		
		for (unsigned int i2 =0; i2< daus.size(); ++i2) {
		  const GenParticle* pfcand = daus[i2];
		  int charge = pfcand->charge();
		  HepLorentzVector cand4v(pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy());
		  //								int pdgid = pfcand->pdgId();
		  
#else
		  std::vector<reco::CandidatePtr> daus((*genjets)[ithird].daughterPtrVector());
		  std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });    
		  
		  for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		    const pat::PackedCandidate &pfcand = static_cast<const pat::PackedCandidate &>(*daus[i2]);
		    
		    int charge = pfcand.charge();
		    HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
#endif
		    
		    
		    //      if (cand4v.perp()<0.5) continue;                                                                                                                     
		    genmom[isrc][1][iet].push_back(cand4v);
		    if (charge !=0) {
		      genmom[isrc][2][iet].push_back(cand4v);
#ifdef TRACKSYS
		      if (gRandom->Uniform() < 0.96) {genmom[isrc][4][iet].push_back(cand4v); }
#endif
		    }
		    
		    if(charge ==0) {
		      if (cand4v.perp()>1.0) {
			genmom[isrc][3][iet].push_back(cand4v);
		      }
		    } else {
		      if (cand4v.perp()>0.5) {
			genmom[isrc][3][iet].push_back(cand4v);
		      }
		    }
		  } //for (unsigned int i2 = 0; i2< daus.size(); ++i2) 
		}*/// if (ithird>=0)
		gen_njets[iet]->Fill(ncount,weighttrg); 
	      } // if (avegenptres[isrc] > leadingPtThreshold[0])
	    } //	for (int isrc=0; isrc<nGenReso; isrc++)
	  } //for (int iet=0; iet<njetetamn; iet++)
	} // if(genjets.isValid() && genjets->size()>=2 && (*genjets)[0].pt()>leadingPtThreshold[0])
	// } //if (genjets.isValid() &&  genjets->size()>=2) 
	h_2ht->Fill(aveleadingpt,avegenpt, weighttrg);
	///////Response
      } //isMC
      //	cout<<"22 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
      // if(isMC) h_2ht->Fill(aveleadingpt,avegenpt, weighttrg);
      //cout << "Write test 325 = ok " << endl;
      //for(int rnum=0; rnum<10; rnum++) {
      /*double rand=gRandom->Uniform();
      int k = rand/0.1;
//      cout << "Rand Number " << k << endl;*/
  
//-----------------------------------------------Calculate And Fill the  EventShape Variables--------------------------------
  for (int itp=0; itp<ntype; itp++) {
	for (int iet=0; iet<njetetamn; iet++) {
	  if (isReconstruct) {
              recovar2[0].clear(); 
              recovar2[1].clear(); 
              recovar2[2].clear(); 
	      recovar1.clear();
	   for (int isrc=0; isrc<njecmx; isrc++) { 
	   // for (int isrc=0; isrc<1; isrc++) { 
	      recovar.clear();
	      //recovar1.clear();
	      if (isrc==0) {isRECO[itp][iet]=false;}
	      isRECO_JER[isrc][itp][iet]=false;
	      
                 if (irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) {
		EventShape_vector  recoevtshape(recomom[isrc][itp][iet], 2.4, 0, 2, 1);
		recovar =  recoevtshape.getEventShapes();
                if(isrc==0){recovar1 =  recoevtshape.getEventShapes();}
                recovar2[isrc] = recoevtshape.getEventShapes();
	        isRECO_JER[isrc][itp][iet]=true;
		if (recovar[nvar]>=2) {
		  if (isrc==0) {isRECO[itp][iet] = true;}
		  for (int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 
		      if (isrc==0) { 
			if (int(recovar[nvar])>=2) {
			nreco++;
                      //if(ij==3 && itp==0 ){cout<<"reco: "<<ievt<<" "<<"Ty:" << itp  << " Nvar : "<<recovar[nvar]<<" "<<recomom[isrc][itp][iet].size() << " Ht2 Bins :" <<irecohtjec[isrc];}
                       //if(itp==0){ cout <<" Var :  " << ij <<" : "<< recovar[ij];}
		       h_recoevtvar[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar[ij], weighttrg); 
                       int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
                       h_recovar_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
			}
			/*for (int irand=0; irand<10; irand++) {
			  if(irand !=k ) h_recoevtvar[irand][itp][irecohtjec[isrc]][iet][ij]->Fill(recovar[ij], weighttrg); 
//#ifdef LHAPDF
			  //for (int ix=1; ix<nnnmx; ix++) {
			  //h_recoevtvarpdf[itp][irecohtjec[isrc]][iet][ij][ix]->Fill(recovar[ij], weighttrg*pdfwt[ix]); 
			  //			 	}
//#endif
			  }*/	
		      } else {
#ifdef JETENERGY
			if (int(recovar[nvar])>=2) {h_recoevtvarjec[itp][irecohtjec[isrc]][iet][ij][isrc]->Fill(recovar[ij], weighttrg);
                           int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
                           h_recoevtvarjec_2D[itp][iet][ij][isrc]->Fill(irecbin, weighttrg); }
#elif defined(JETRESO)
			if (int(recovar[nvar])>=2) {h_recoevtvarres[itp][irecohtjec[isrc]][iet][ij][isrc]->Fill(recovar[ij], weighttrg);
                           int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
                           h_recoevtvarres_2D[itp][iet][ij][isrc]->Fill(irecbin, weighttrg);}
		         
#endif
		      }
		    }
		  }
		}
	      }
	    }
	  } // if (isReconstruct)
//	  cout << endl;
	  if(isMC) {
	    for (int isrc=0; isrc<nGenReso; isrc++) {
	      //for (int isrc=0; isrc<1; isrc++) { 
	      genvar.clear();
	      bool isGEN=false;
	      if (isMC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1) { 
		EventShape_vector  genevtshape(genmom[isrc][itp][iet], 2.4, 0, 2, 1);
		
		genvar =  genevtshape.getEventShapes();
		if (genvar[nvar]>=2) {
		  isGEN = true;
		  for (int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 
		      if (isrc==0) { 
			if (int(genvar[nvar])>=2) {
			h_genevtvar[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);
			int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                        h_genvar_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
                        //if(ij==3 && itp==0 ){cout<<"Gen: "<<ievt<<" "<<"Ty:" << itp  << " Nvar : "<<genvar[nvar]<<" "<<genmom[isrc][itp][iet].size() << " Ht2 Bins :" <<igenhtres[isrc];}
                        //if(itp==0){ cout <<" Var :  " << ij <<" : "<< genvar[ij];}
			} //else {
			  //h_genevtvar2[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);
			//}
#ifdef JETRESO
			//	} else {
			//    	h_genevtvarres[itp][igenhtres[isrc]][iet][ij][isrc]->Fill(genvar[ij], weighttrg);	
#endif
#ifdef LHAPDF
			for (int ix=1; ix<nnnmx; ix++) {
			if (int(genvar[nvar])>=2) {h_genevtvarpdf[itp][igenhtres[isrc]][iet][ij][ix]->Fill(genvar[ij], weighttrg*pdfwt[ix]);
                        int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                        h_genevtvarpdf_2D[itp][iet][ij][ix]->Fill(igenbin, weighttrg*pdfwt[ix]);   
                                                                    }
			}
#endif
		      }
		    }
		  }
		}
	      }

  //                 cout << " JER ok : " <<isrc <<endl; 
///cout <<endl;
//cout << " is reco: " << isRECO_JER[0][itp][iet] << " " << isRECO_JER[1][itp][iet] << " "<< isRECO_JER[2][itp][iet] ;	
  	      if(isrc==0 && isReconstruct){ 
int ijer =2;
	  for(int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 	
        //if(isRECO[itp][iet] && isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) { 
//if(isRECO_JER[ijer][itp][iet]){	
	      if(isRECO_JER[ijer][itp][iet] && isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) { 
			naa++;
	                if(recovar2[ijer][nvar]>=2 &&  genvar[nvar]>=2){
          //               cout << " isrc : " <<isrc << " Reco var : " << recovar1[ij]<<" Gen Var " << genvar[ij]  <<endl; 
   //                      cout << " ijer: " <<ijer <<" HT Nominal "<< aveleadingptjec[isrc] <<endl;
     //                   cout << irecohtjec[ijer]<< " "<< recomom[ijer][itp][iet].size() << " " ; 
      //                   cout<< "HT JER "<< " Reco var : " << recovar2[0][ij]<< " " << recovar2[ijer][ij] << " " << recovar2[ijer][ij] <<" Gen Var " << genvar[ij]  <<endl; 
			 h_2devtvar[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], genvar[ij], weighttrg);
			 int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                         //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
                        int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
      //                   cout << "XId  "<< irecbin << " YID " << igenbin << " weg :" << weighttrg <<endl; 
                         RM_2D[itp][iet][ij]->Fill(irecbin, igenbin, weighttrg);
                        }else if (recovar2[isrc][nvar]>=2) {
			 
                         //h_2devtvar[itp][igenht][iet][ij]->Fill(recovar[ij],-10.0, weighttrg);				
			  h_recoevtfake[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], weighttrg);
                          //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
                         int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
                         h_recofake_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
                        }else if (genvar[nvar]>=2) {
			//h_2devtvar[itp][igenht][iet][ij]->Fill(-10.0, genvar[ij], weighttrg);	//Fill in Reco Underflow
			  h_genevtmiss[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);	
                          int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                          h_genmiss_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
                            }
			//  h_2devtvar[itp][0][iet][ij]->Fill(recovar[ij], genvar[ij], weighttrg);
		      } else {
			//if (isRECO[itp][iet] && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1 && recovar1[nvar]>=2) {
			if (isRECO_JER[ijer][itp][iet] && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1 && recovar2[ijer][nvar]>=2) {
			  nbb++;
			  //h_2devtvar[itp][igenht][iet][ij]->Fill(recovar[ij],-10.0, weighttrg); //Fill Fake in Gen Underflow
			    h_recoevtfake[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], weighttrg);
                            //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
                            int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
                            h_recofake_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
			}
			if (isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && genvar[nvar]>=2) {
			  ncc++;
                             
			   h_2devtvar[itp][igenht][iet][ij]->Fill(-10.0, genvar[ij], weighttrg);	//Fill Miss in Reco Underflow
			   h_genevtmiss[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);	
			   int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                           h_genmiss_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
			}
		      }
		    } //if (isItUsed(ij)) 
		  } // for(int ij=0; ij<nvar; ij++)	
		} // if (isrc==0 && isReconstruct)
/*
//-----------------------------------------------------------------for JER
#ifdef JETRESO
if(isrc==0 && isReconstruct){
   for(int ijer=0 ; ijer < njecmx ; ijer++){
                  for(int ij=0; ij<nvar; ij++) {
                    if (isItUsed(ij)) {
                      if(isRECO_JER[ijer][itp][iet] && isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[ijer]>=0 && irecohtjec[ijer]<njetptmn && recomom[ijer][itp][iet].size()>1) {                  
                        
                        if(recovar2[ijer][nvar]>=2 &&  genvar[nvar]>=2){
                         cout << " isrc : " <<isrc << " ijer : " << ijer << " HT : " << aveleadingptjec[ijer] << " " <<genvar[ij]<<" " <<recovar2[ijer][ij] <<"   ";
                         int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                         int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
                         cout << " OK " << igenbin<<" " << irecbin<<endl;
 //                        RM_JER_2D[itp][iet][ij][ijer]->Fill(irecbin, igenbin, weighttrg);
                        }else if (recovar2[ijer][nvar]>=2) {
                          int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
  //                        h_reco_JER_fake_2D[itp][iet][ij][ijer]->Fill(irecbin, weighttrg);
                        }else if (genvar[nvar]>=2) {
                          int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
  //                        h_gen_JER_miss_2D[itp][iet][ij][ijer]->Fill(igenbin, weighttrg);      
  
                          }
                      } else {
                        if (isRECO_JER[ijer][itp][iet] && irecohtjec[ijer]>=0 && irecohtjec[ijer]<njetptmn && recomom[ijer][itp][iet].size()>1 && recovar2[ijer][nvar]>=2) {
                         
                            int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar2[ijer][ij],aveleadingptjec[ijer]);
//                            h_reco_JER_fake_2D[itp][iet][ij][ijer]->Fill(irecbin, weighttrg);
                        }
                        if (isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && genvar[nvar]>=2) {
                         

                           int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
 //                          h_gen_JER_miss_2D[itp][iet][ij][ijer]->Fill(igenbin, weighttrg);
                        }

                      }
                    } //if (isItUsed(ij)) 
                  } // for(int ij=0; ij<nvar; ij++)     
                } //ijer <3
                } // if (isrc==0 && isReconstruct)
#endif
*/
//-----------------------------------------------------------


	      //} // if (igenht>=0 && igenht<njetptmn && genmom[isrc][itp][iet].size()>1)
	    } // for (int isrc=0; isrc<nGenReso; isrc++)
	  }//isMC
	} // for (int iet=0; iet<njetetamn; iet++)
      } //for (int itp=0; itp<ntype; itp++) 
      
      //if (nevt%1000==1) { std::cout <<"nevt "<<nevt<<" naa "<<naa<<" nbb "<<nbb<<" ncc "<<ncc<< std::endl;}
if(nevt==100){   cout <<igenht <<endl;}

    cout <<"end event" << endl;
    }

// ------------ method called once each job just before starting event loop  ------------
void 
QCDEventShape::beginJob() {
//  t1=clock();
  nevt = 0;
  if (isMC) { 
    double dattot[nHLTmx]={0};
    double mctot=0;
    for (int ij=0; ij<npileupmx; ij++) {
      for (int jk=0; jk<nHLTmx; jk++) {
				dattot[jk] +=datpileup[jk][ij];
      }
      mctot +=mcpileup[ij];
    }
    
    for (int ij=0; ij<npileupmx; ij++) {
      mcpileup[ij] /=max(1.e-6,mctot);
      for (int jk=0; jk<nHLTmx; jk++) {
				datpileup[jk][ij] /=max(1.e-6,dattot[jk]);
				
				rat_pileup[jk][ij] =  datpileup[jk][ij]/mcpileup[ij];
      }
    }
  }

#ifdef JETENERGY
  for (int isrc = 0; isrc < nsrc; isrc++) {
    const char *name = srcnames[isrc];
//    JetCorrectorParameters *p = new JetCorrectorParameters("Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt", name);
 // JetCorrectorParameters *p = new JetCorrectorParameters("Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt", name);    
  JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt", name);    
//  JetCorrectorParameters *p = new JetCorrectorParameters("Fall17_17Nov2017F_V6_DATA_UncertaintySources_AK4PFchs.txt", name);    
JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
		//    vsrc[isrc] = unc;
		vsrc.push_back(unc);
  }
#endif  


//cout << "Write test 34 = ok " << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDEventShape::endJob() 
{

         TUnfoldBinng2D->cd();
     for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
             if (isItUsed(ij)) {

        h_recovar_2D[ityp][iet][ij]->Write();
        h_recofake_2D[ityp][iet][ij]->Write();
        h_genvar_2D[ityp][iet][ij]->Write();
        h_genmiss_2D[ityp][iet][ij]->Write();
        RM_2D[ityp][iet][ij]->Write();

#ifdef  LHAPDF
            for (int ix=1; ix<nnnmx; ix++) {h_genevtvarpdf_2D[ityp][iet][ij][ix]->Write(); }
#endif
#ifdef  JETENERGY
            for (int ix=1; ix<njecmx; ix++) {h_recoevtvarjec_2D[ityp][iet][ij][ix]->Write();   }
#elif defined(JETRESO)
            for (int ix=1; ix<njecmx; ix++ ) {   
             h_recoevtvarres_2D[ityp][iet][ij][ix]->Write();  
//             RM_JER_2D[ityp][iet][ij][ix]->Write();
//             h_reco_JER_fake_2D[ityp][iet][ij][ix]->Write();
//             h_gen_JER_miss_2D[ityp][iet][ij][ix]->Write();

 }
#endif

           }
         }       
       }
     }
  //theFile->cd();
  //theFile->Write();
  //theFile->Close();
  //myfile1->Close();
  //fs->Write();
  //fs->Close();
}

// ------------ method called when starting to processes a run  ------------

void 
QCDEventShape::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
// Initialize hltConfig

#ifdef TRIGGER

// cout << "Write test 4 = ok " << endl;
	bool changed(true);
  if (hltPrescaleProvider_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
    hltConfig.dump("Triggers");
    hltConfig.dump("PrescaleTable");

    for (unsigned int ij=0; ij<nHLTmx; ij++) {
      l1pres[ij] = hltpres[ij]=-7;
    }

       } else {
         }



/*   bool changedConfig;
   if (!hltConfig_.init(iRun, iSetup, theHLTTag.c_str(), changedConfig)) {
     LogError("HLTMuonVal") << "Initialization of HLTConfigProvider failed!!"; 
     return;
 
    //for (unsigned int ij=0; ij<nHLTmx; ij++) {
      //l1pres[ij] = hltpres[ij]=-7;
    //}

  }*/
 
/* 
bool changed(true);
  if (hltConfig_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
//    cout <<"Trigger tables "<<endl;
    hltConfig_.dump("Triggers");
//    cout <<"Prescale tables "<<endl;
    hltConfig_.dump("PrescaleTable");

    for (unsigned int ij=0; ij<nHLTmx; ij++) {
      l1pres[ij] = hltpres[ij]=-7;
    }

    //trig_init=0; //GMA
    //    // ..
       } else {
              // ..
         }

*/




/* bool changed(true);
   if (hltConfig_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
     if (changed) {
      // check if trigger name in (new) config
       if (triggerName_!="@") { // "@" means: analyze all triggers in config
     const unsigned int n(hltConfig_.size());
     const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
     if (triggerIndex>=n) {
       LogVerbatim("HLTEventAnalyzerAOD") << "HLTEventAnalyzerAOD::analyze:"
            << " TriggerName " << triggerName_ 
            << " not available in (new) config!" << endl;
       LogVerbatim("HLTEventAnalyzerAOD") << "Available TriggerNames are: " << endl;
       hltConfig_.dump("Triggers");
     }
       }
       hltConfig_.dump("ProcessName");
       hltConfig_.dump("GlobalTag");
       hltConfig_.dump("TableName");
       hltConfig_.dump("Streams");
       hltConfig_.dump("Datasets");
       hltConfig_.dump("PrescaleTable");
       hltConfig_.dump("ProcessPSet");
     }
   } else {
     LogVerbatim("HLTEventAnalyzerAOD") << "HLTEventAnalyzerAOD::analyze:"
      << " config extraction failure with process name "
      << processName_ << endl;
   }
 
*/
#endif
 
  std::cout<<" End of QCDEventShape::beginRun"<<std::endl; //"nevt "<<nevt<<" naa "<<naa<<" nbb "<<nbb<<" ncc "<<ncc<< std::endl;
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
QCDEventShape::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
std::cout<<" End of QCDEventShape::beginRun"<<std::endl;
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
QCDEventShape::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{


}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
QCDEventShape::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QCDEventShape::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double PhiInRange(const double& phi) {
      double phiout = phi;
      
      if( phiout > 2*M_PI || phiout < -2*M_PI) {
	phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;
      
      return phiout;
}

template <class T, class U>
double deltaR(const T& t, const U& u) {
  return sqrt(pow(t.eta()-u.eta(),2) +pow(PhiInRange(t.phi()-u.phi()),2));
}



//define this as a plug-in
DEFINE_FWK_MODULE(QCDEventShape);

/*
L1_ZeroBias 	29989	5327	5327	5327	5327	1601	801	801	801	801	801	801	801
17	L1_SingleJet52 	3000	10000	6000	4000	3000	1500	800	500	400	300	150	100	262139
18	L1_SingleJet68 	1500	1500	1000	750	500	300	150	100	75	50	30	15	262139
19	L1_SingleJet92 	3000	3000	2000	2000	1500	800	400	300	200	150	80	40	262139
20	L1_SingleJet128 1000	1000	1	1	1	1	1	1	1	1	1	1	262139
21	L1_SingleJet176 300	300	1	1	1	1	1	1	1	1	1	1	262139


66	HLT_DiPFJetAve40_v2 (2013430) 	25 	25 	16 	12 	8 	5 	3 	2 	1 	1 	1 	1 	0 	L1_ZeroBias
69	HLT_DiPFJetAve60_v2 (2013431) 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	0 	L1_ZeroBias
71	HLT_DiPFJetAve80_v2 (2013432) 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	L1_SingleJet52
58	HLT_DiPFJetAve140_v2 (2013433) 	2 	2 	2 	2 	1 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet92
60	HLT_DiPFJetAve200_v2 (2013434) 	250 	250 	125 	85 	60 	35 	16 	12 	9 	6 	4 	2 	1 	L1_SingleJet128
62	HLT_DiPFJetAve260_v2 (2013435) 	85 	85 	85 	60 	42 	24 	12 	8 	6 	4 	2 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
64	HLT_DiPFJetAve320_v2 (2013436) 	15 	15 	15 	10 	6 	4 	2 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
65	HLT_DiPFJetAve400_v2 (2013437) 	5 	5 	5 	3 	2 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
67	HLT_DiPFJetAve500_v2 (2013438) 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
*/
