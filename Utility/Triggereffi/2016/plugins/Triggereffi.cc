// -*- C++ -*-
//
// Package:    Test/Triggereffi
// Class:      Triggereffi
//
/**\class Triggereffi Triggereffi.cc Test/Triggereffi/plugins/Triggereffi.cc
   
   Description: [one line class summary]
   
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Suman Kumar Kundu
//         Created:  Sun, 19 Jan 2020 10:50:07 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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

#include "TH2F.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>

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

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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

#include "FWCore/Framework/interface/EDAnalyzer.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace reco;
using namespace std;
using namespace CLHEP;
using namespace trigger;
using namespace math;

const int nHLTmx = 553;               //Total Number of Trigger ** previous 448
static const  int nMuHLTmx=11;      // muon trigger number
static const int nDiJetHLTmx=9;     //DiJet Trigger no.
static const int nJetHLTmx=8;       //DiJet Trigger no.

const int njetetamn=1;                     // one eta space is choosen 

const char* jethlt_name[nJetHLTmx]={"HLT_PFJet80_v",
                                    "HLT_PFJet140_v",
                                    "HLT_PFJet200_v",
                                    "HLT_PFJet260_v",
                                    "HLT_PFJet320_v",
                                    "HLT_PFJet400_v",
                                    "HLT_PFJet450_v",
                                    "HLT_PFJet500_v"};

//const char* dijethlt_name[nDiJetHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};                           // Di-Jet Trigger Name
const char* dijethlt_name[nDiJetHLTmx]={"HLT_DiPFJetAve40_v","HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};                           // Di-Jet Trigger Name


const char* dijethlt_label[nDiJetHLTmx]={"hltDiPFJetAve40","hltDiPFJetAve60","hltDiPFJetAve80", "hltDiPFJetAve140", "hltDiPFJetAve200", "hltDiPFJetAve260", "hltDiPFJetAve320", "hltDiPFJetAve400", "hltDiPFJetAve500"};                           // Di-Jet Trigger Name

const char* muhlt_name[nMuHLTmx]={"HLT_IsoMu17_eta2p1_", //63
                                  "HLT_IsoMu20_", //72
                                  "HLT_IsoMu30_", // 71
                                  "HLT_IsoMu24_eta2p1_", //77
                                  "HLT_IsoMu27_", //78
                                  "HLT_IsoTkMu20_", //79
                                  "HLT_IsoTkMu20_eta2p1_", //80
                                  "HLT_IsoTkMu24_eta2p1_", //81
                                  "HLT_IsoTkMu27_", //82
                                  "HLT_IsoMu18_",   //449
                                  "HLT_IsoMu22_"};   //453     // Muon trigger name

double l1Pt[nDiJetHLTmx] = {0,35,60,90,120,170,170,170};
//double jethlt_thr[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
double jethlt_thr[nDiJetHLTmx]={40,60,80,140,200,260,320,400,500};
double leadingPtThreshold[nDiJetHLTmx]={40,60,80,140,200,260,320,400,500};
//double leadingPtThreshold[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
double etarange[njetetamn] ={2.5};
bool trgpas[nDiJetHLTmx];//={0,0,0,0,0,0,0,0};


int l1pres[nDiJetHLTmx], hltpres[nDiJetHLTmx], compres[nDiJetHLTmx];    // For jet
int l1mupres[nMuHLTmx], hltmupres[nMuHLTmx];                           // For Muons

double dR(double eta1, double phi1, double eta2, double phi2);
static double PhiInRange(const double& phi);

int hlt_list[nHLTmx];
int hltmu_list[nHLTmx]; //[nMuHLTmx];
int hltjet_list[nHLTmx]; //[nDiJetHLTmx];
int hltdijet_list[nHLTmx]; //[nDiJetHLTmx];
int hltdijet_mimic[nDiJetHLTmx]; //[nDiJetHLTmx];

std::map<std::string, bool> fired;
int nevt;
int ntrig;
//int nSingleMutig =0;

// DeltaR function
// //template <class T, class U> static double deltaR(const T& t, const U& u);
// //double deltaR(double eta1, double phi1, double eta2, double phi2);
// int sbitx(unsigned ival, int ibit);
//

struct triggervar{
  HepLorentzVector trg4v;
  bool            both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
};



//using reco::TrackCollection;

class Triggereffi : public edm::EDAnalyzer {
public:
  explicit Triggereffi(const edm::ParameterSet&);
  ~Triggereffi();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
 
 virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
 
  //Dijet trigger efficiency
//  TH1F* hlt_dijettag[nDiJetHLTmx][njetetamn];
//  TH1F* hlt_dijetprob[nDiJetHLTmx][njetetamn];
  
  TH1F* hlt_dijettrg_fired[nDiJetHLTmx][njetetamn];
  TH1F* hlt_dijettrg_all_evt[nDiJetHLTmx][njetetamn];
  
  std::string theRootFileName;
  std::string theHLTTag;
  // ----------member data ---------------------------
  edm::EDGetTokenT<GenEventInfoProduct> generator1_;
  edm::EDGetTokenT<pat::MuonCollection> MuonToken_;     // Edited
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
  HLTPrescaleProvider hltPrescaleProvider_;
  
  
  //    edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
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
Triggereffi::Triggereffi(const edm::ParameterSet& iConfig):
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
  generator1_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("evtinfo"))),
  MuonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),  //Edited
  jetSrcToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  genSrcToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genSrc"))),
  PFSrcToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),                           // This is related to Trigger result and trigger fired decision
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),  // Trigger object 
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),         // for prescale 
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
  theHLTTag = iConfig.getUntrackedParameter<string>("HLTTag", "HLT");
  char name[200];
  char title[200];
  
  //define histogram
  //These are for Tag and probe
/*  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
    for (int unsigned jk=0; jk<njetetamn; jk++) {
      sprintf(name, "hlt_dijettag_%i_%i", ij, jk);
      sprintf(title, "dijet tagged P_T : (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      hlt_dijettag[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettag[ij][jk]->Sumw2();
      
      sprintf(name, "hlt_dijetprob_%i_%i", ij, jk);
      sprintf(title, "dijet probed P_T : (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      hlt_dijetprob[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijetprob[ij][jk]->Sumw2();
    }
  }*/
  
  
  //Define Histogram for efficiency
  // These are for New 
  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
    for (int unsigned jk=0; jk<njetetamn; jk++) {
      
      sprintf(name, "hlt_dijet_effi_%i_%i", ij, jk);
      sprintf(title, "dijet trigger fired: (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      hlt_dijettrg_fired[ij][jk] = fs->make<TH1F>(name, title, 1200, 40, 1240);
      //hlt_dijettrg_fired[ij][jk] = fs->make<TH1F>(name, title, 100, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettrg_fired[ij][jk]->Sumw2();
    }
  }
  
  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
    for (int unsigned jk=0; jk<njetetamn; jk++) {
      
      sprintf(name, "hlt_dijet_all_evt_%i_%i", ij, jk);
      sprintf(title, "dijet trigger All event: (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      hlt_dijettrg_all_evt[ij][jk] = fs->make<TH1F>(name, title, 1200, 40, 1240);
      hlt_dijettrg_all_evt[ij][jk]->Sumw2();
    }
  }
  
  nevt=0;  //event counter
  ntrig=0;  //event counter
  
  
}       //End constructor    //Triggereffi::Triggereffi(const edm::ParameterSet& iConfig): 


Triggereffi::~Triggereffi()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Triggereffi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  
  /*
    #ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
    #endif
    
    #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    #endif
  */
  //  cout<< "New Event  <<    **********************************************New Event"<< nevt <<endl;
  //   cout <<endl;
  
  if (nevt%10000==1) cout <<"TriggerEfficiency::analyze run number = "<<nevt<<endl;
//  cout << " Ok1 " << endl;  
  const char* variab1;
  const char* variab2;
  double aveleadingpt =0;
  
  
  
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerBits_, trigRes);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(MuonToken_,muons);
  
  
  const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  bool isInEtaRange[njetetamn]={0};
  
  edm::Handle<pat::JetCollection> ak4PFJets;
  iEvent.getByToken(jetSrcToken_, ak4PFJets);
  
  //select events where one of the SingleMu triggers fires 
  //  unsigned int musz=muons->size();
  //if(musz <= 0) return;
  // cout << "Event with muon Numbers : " <<musz<< endl;
  
  
  
  
  
  if (!trigRes.isValid()) return;
  bool mutrig=false;
  bool dijtrig=false;
  edm::TriggerNames triggerNames = iEvent.triggerNames(*trigRes);   //Get The trigger name for TriggerResult Input
  
  if (trigRes.isValid()) {
    unsigned int size = trigRes->size();
    for (unsigned jk=0; jk<nHLTmx; jk++) {
      hlt_list[jk] = hltmu_list[jk] =hltdijet_list[jk] =-1;
    }         
    for(unsigned ij = 0; ij != size; ij++) {
      std::string name = triggerNames.triggerName(ij);
      if(nevt==1)  {cout <<"inclu:hltobject "<<" "<<ij<<" "<<name<<endl;}
      //  if(trigRes->accept(ij)){cout << "trigger accepted  for : " << name << endl;}  
      variab1 = name.c_str();
      //      if(trigRes->accept(ij)){cout << "trigger accepted  for (cross checked) : " << variab1 << endl;} 
      
      //Check for muon trigger
      for (unsigned kl=0; kl<nMuHLTmx; kl++) {
	if ((strstr(variab1,muhlt_name[kl])) && strlen(variab1)-strlen(muhlt_name[kl])<5) {
	  //            cout << "singleMu Trigger Name : " << variab1 << endl;
	  hltmu_list[ij] = kl;
	  mutrig=true;
	  break;
	}
      }
      
      //check for DijetHLT trigger
      for (unsigned kl=0; kl<nDiJetHLTmx; kl++) {
	if ((strstr(variab1,dijethlt_name[kl])) && strlen(variab1)-strlen(dijethlt_name[kl])<5) {
	  hltdijet_list[ij] = kl; 
    //     cout << variab1 <<  endl; 
	  dijtrig=true; 
	  break;
	}
      }
      
    } //Trigger size
    
  }//Triger result valid
	  

  mutrig=true; //no function of this here
  
  if (!mutrig) return; // No need of this condition 
  if (!dijtrig) return;
  
  //*********************************************
  
  
  
  //**************************************************
  //Calculate average Pt 
  if (ak4PFJets.isValid() &&  ak4PFJets->size()>1) {
    for (int iet=0; iet<njetetamn; iet++) {
      isInEtaRange[iet] = true;
    }
    
    for (int ij=0; ij<2; ij++) {
      for (int iet=0; iet<njetetamn; iet++) {
        if (abs((*ak4PFJets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
      }
      //Jet ID ================= Tight ID 2017 Recomendation added as https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2017
      double NHF = (*ak4PFJets)[ij].neutralHadronEnergyFraction();
      double NEMF = (*ak4PFJets)[ij].neutralEmEnergyFraction();
      double CHF = (*ak4PFJets)[ij].chargedHadronEnergyFraction();
  //  double  MUF  = (*ak4PFJets)[ij].muonEnergyFraction();
  //  double CEMF = (*ak4PFJets)[ij].chargedEmEnergyFraction();
      int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
     // int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
      int CHM = (*ak4PFJets)[ij].chargedMultiplicity();
      bool TightJetID =false;
//-------------------------------2017 UL Update
//      if(abs((*ak4PFJets)[ij].eta())<=2.7){
//                     if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && (abs((*ak4PFJets)[ij].eta())<=2.4 && CHF>0 && CHM>0 )) TightJetID =true;
//      } else {
//                         TightJetID =false;
//      }
//      
//      if (abs((*ak4PFJets)[ij].eta())>2.7) {TightJetID = false;}
//      if ((*ak4PFJets)[ij].pt()<30.0) {TightJetID = false;}
//----------------------------------2018 UL ID update 
      if(abs((*ak4PFJets)[ij].eta())<=2.7){
      if(NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0  && abs((*ak4PFJets)[ij].eta())<=2.4 )  TightJetID =true;
      if(NHF<0.90 && NEMF<0.99 && NumConst>1 && CHF>0 && CHM>0  && abs((*ak4PFJets)[ij].eta())>2.4 )  TightJetID =true;
                                   }else{ 
                                 TightJetID =false;
                                 }
      if (abs((*ak4PFJets)[ij].eta())>2.7) {TightJetID = false;}
      if ((*ak4PFJets)[ij].pt()<30.0) {TightJetID = false;}

      
      if (TightJetID) {
        aveleadingpt +=(*ak4PFJets)[ij].pt();
      } else {
        aveleadingpt -=100000;
      }
    }
    aveleadingpt /=2.0;
    
  } // calculation of average pt
  //***********************************************************
  if( aveleadingpt < 0) return;
  nevt++;
  

  //Prescale Factor Calculation
  int preL1(-1), preHLT(-1), prescale(-1);

  // trgpas[nDiJetHLTmx]={false};
  for (int jk=-1; jk<nDiJetHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str();
      if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) ||
          (jk>=0 && strstr(variab1,dijethlt_name[jk]) && strlen(variab1)-strlen(dijethlt_name[jk])<5)) {
           const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltPrescaleProvider_.prescaleValuesInDetail(iEvent,iSetup,variab1));
            preL1 = prescalesInDetail.first[0].second;  //acesses the L1 prescale
            preHLT = prescalesInDetail.second;     // acesses the HLT prescale
//           if(preL1<=0){preL1=1;}      //skip the l1 prescale if it give negative or zero
             prescale = preL1 * preHLT;
             compres[jk] = prescale;
     //cout << " By index " << triggerPrescales->getPrescaleForIndex(ij) << "  L1 : "<< preL1 << "  HLT :" << preHLT << "  L1*HLT "  << compres[jk] <<endl;
     // cout << "prescale check " <<   compres[jk] << endl; 
     }
    }
  }//calculation of prescale
  
  //Second condition:  select events with two reconstructed jets above a certain (low) threshold
  if ((!ak4PFJets.isValid()) ||  ak4PFJets->size() <2) return;
  if ((*ak4PFJets)[0].pt()<30.0 || (*ak4PFJets)[1].pt()<30.0)  return ;
  if (fabs((*ak4PFJets)[0].eta())>2.5 || fabs((*ak4PFJets)[1].eta())>2.5)  return ;
  if (!isInEtaRange[0] && (aveleadingpt<30)) return;
  
  
  
  
  
  //const edm::TriggerNames & trigName = iEvent.triggerNames(*trigRes);
  /*
    bool singleMu =false;
    if( !trigRes.failedToGet() ) {
    int N_Triggers = trigRes->size();
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
    
    if (trigRes.product()->accept(i_Trig)) {
    TString TrigPath =trigName.triggerName(i_Trig);
    //    cout << "Passed trigger Path : " << TrigPath <<endl;    
    //  if(TrigPath.Index("HLT_IsoMu24_v") >=0) {cout << " HLT_IsoMu24_v " << " Trigger Found "<< endl ;}
    for (unsigned ij=0; ij<nMuHLTmx; ij++){ 
    if ((strstr(TrigPath,muhlt_name[ij])) && strlen(TrigPath)-strlen(muhlt_name[ij])<5) {
    cout <<"Passed SingleMu trigger Name : "<< TrigPath <<" || Trigger Bits :"<< trigRes.product()->accept(i_Trig) << endl;
    singleMu = true;
    nSingleMutig++;
    
    
    }
    } //nMuHLTmx
    
    }
    }
    } //!trigRes.failedToGet()
  */
  
  //if(!singleMu) return;
  //if(!singleMu) {cout << "great "<< endl;}
  
  /*
  //Fill event selected event
  for (int iet=0; iet<njetetamn; iet++) {
  for (int jk=0; jk<nDiJetHLTmx; jk++) {
  //      if(aveleadingpt>jethlt_thr[jk-1]) {
  hlt_dijettrg_all_evt[jk][iet]->Fill(aveleadingpt, compres[jk]);
  // }
  }
  } //Fill the all valid event 
  */
const edm::TriggerNames & trigName = iEvent.triggerNames(*trigRes);
bool HLTtrg = false;
 if(!trigRes.failedToGet()) {
       int N_Triggers = trigRes->size();
       for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigRes.product()->accept(i_Trig)) {
         TString TrigPath =trigName.triggerName(i_Trig);
//     cout << "TrigPath  fire name : "<< TrigPath<< endl;
         for (unsigned ij=0; ij<nDiJetHLTmx; ij++){
//       if ((strstr(variab1,dijethlt_name[ij])) && strlen(variab1)-strlen(dijethlt_name[ij])<5) {
       if ((strstr(TrigPath,dijethlt_name[ij])) && strlen(TrigPath)-strlen(dijethlt_name[ij])<5) {
          trgpas[ij] = true ;
          HLTtrg = true;
//      cout <<"Dijet Trigger Name : " << TrigPath << "  Trigger Bits : " << trigRes.product()->accept(i_Trig) <<" trgpass[]="<<  trgpas[ij] <<endl;
               }
            } //nDiJjetHLTmx
         }
      }
   } //!trigRes.failedToGet()

if(!HLTtrg) return;
 
  
  
  //Third Condition : match the two jets with the trigger objects in eta-phi space
  bool delta_object=false;
  double object_pt[2];
  object_pt[0]=-1;
  object_pt[1]=-1;
  bool obj1 = false;
  bool obj2 = false;

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { //Trigger object standalone
    // if(obj.pt() < 10.0 && obj.eta() > 2.4) continue;              //look for HLT objects with Pt > 10 GeV only and within eta  2.4 
    obj.unpackFilterLabels(iEvent,*trigRes); 
    obj.unpackPathNames(names);
    
    for (unsigned ih = 0; ih < obj.filterLabels().size(); ++ih){
      std::string objfillabel = obj.filterLabels()[ih];
      variab2 = objfillabel.c_str();
      for (int jk=0; jk<nDiJetHLTmx; jk++) { 
	  if (strstr(variab2,dijethlt_label[jk]) && strlen(variab2)-strlen(dijethlt_label[jk])<5){ 
        if(trgpas[jk]){
//	  cout << "Trigger object name, ID pt, eta, phi: " << variab2 <<","<< obj.filterIds()[ih]<<", " << obj.pt()<<", "<<obj.eta()<<", "<<obj.phi() << endl;
	  
	  //-------------------------------------------------------- 
	  double dr1 = dR((*ak4PFJets)[0].eta(), (*ak4PFJets)[0].phi(), obj.eta(), obj.phi());
	  double dr2 = dR((*ak4PFJets)[1].eta(), (*ak4PFJets)[1].phi(), obj.eta(), obj.phi());
	  if (dr1 <0.3){object_pt[0]=obj.pt();
                    obj1 = true;}
	  if (dr2 <0.3){object_pt[1]=obj.pt();
                      obj2 = true; }
//	  cout <<"Reco1 PT :"<< (*ak4PFJets)[0].pt()<<" Reco2 PT :"<< (*ak4PFJets)[1].pt()<<" ; Object PT " << obj.pt() << " ; " << (*ak4PFJets)[0].eta() <<" ; " << obj.eta()<< " ; "               
//	       << (*ak4PFJets)[0].phi()<< " ; " << obj.phi() << endl;
//	  cout << "dR = " << dr1 << "  :  " << dr2<< endl;
	  if( dr1<.3 || dr2 <.3 ) {delta_object = true; }
	} //trgpass    
       }//label matched
      }//jk<nDiJetHLTmx
    }//ih < obj.filterLabels().size();
  }//Trigger object standalone
  
  if(!delta_object) return;  
  if(!obj1 || !obj2) return;
  if(object_pt[0]< 0 || object_pt[1] < 0) return;

  // cout <<"Delta Condtion: " << delta_object <<  obj1 << obj2 << endl;
  
  
  double objectPt_ave = (object_pt[0]+object_pt[1])/2;
  //mimic the DijetPFJet trigger decision (average of the two trigger objects pT above a certain threshold) 
  if(objectPt_ave < 30 ) return;
  
 // cout <<"Reco JetHT : " << aveleadingpt << "  HLT jetAve : " << objectPt_ave << " Trig fired : " << ntrig << " Jets No :" << ak4PFJets->size() <<endl; 
  ntrig++;
  
  
  //Fill event selected event
  for (int iet=0; iet<njetetamn; iet++) {
    for (int jk=0; jk<nDiJetHLTmx; jk++) {
      //     if(aveleadingpt>jethlt_thr[jk]) {
      hlt_dijettrg_all_evt[jk][iet]->Fill(aveleadingpt, compres[jk]);
//     cout << "Filled PS : " << compres[jk] << endl;
        //    }
    }
  } //Fill the all reco jet 
  
  
  for (int iet=0; iet<njetetamn; iet++) {
    for (int jk=0; jk<nDiJetHLTmx; jk++) {
      if(objectPt_ave>=jethlt_thr[jk]) {
	hlt_dijettrg_fired[jk][iet]->Fill(aveleadingpt, compres[jk]);
      }
    }
  } //Fill the HLT trigger fired*/
  //
  
  
  //- plot the efficiency (events where the diet emulation method fires divided by the total number of events) as a function of the average pT at the reconstruction level.
  
  
  /*
  if (trigRes.isValid()) {
    unsigned int size = trigRes->size();
    for(unsigned ij = 0; ij != size; ij++) {
      if (trigRes->accept(ij)) {
	std::string name = triggerNames.triggerName(ij);
	variab2 = name.c_str();
	for (unsigned ij=0; ij<nDiJetHLTmx; ij++){
	  if ((strstr(variab2,dijethlt_name[ij])) && strlen(variab2)-strlen(dijethlt_name[ij])<5) {
	    //   if (!trigRes->accept(ij)) continue; 
	    cout <<"Dijet trigger Check : " <<trigRes->accept(ij) <<"  Name : "<< variab2 <<endl;
	    
	    for (int iet=0; iet<njetetamn; iet++) {
	      //    hlt_dijettrg_fired[ij][iet]->Fill(aveleadingpt, compres[ij]);
	      //     cout << " New Trigger Found for event number : " << nevt << endl;       
	      //       hlt_dijettrg_all_evt[ij-1][iet]->Fill(aveleadingpt, compres[jk]);
	    }
	  }
	  
	} //nDiJjetHLTmx
      }//size
    } 
  } // trigRes.isValid()
  */
/*
  const edm::TriggerNames & trigName = iEvent.triggerNames(*trigRes);
  if( !trigRes.failedToGet() ) {
    int N_Triggers = trigRes->size();
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigRes.product()->accept(i_Trig)) {
	TString TrigPath =trigName.triggerName(i_Trig);
	for (unsigned ij=0; ij<nDiJetHLTmx; ij++){
	  if ((strstr(TrigPath,dijethlt_name[ij])) && strlen(TrigPath)-strlen(dijethlt_name[ij])<5) {
	    for (int iet=0; iet<njetetamn; iet++) {
	      //       hlt_dijettrg_fired[ij][iet]->Fill(aveleadingpt, compres[ij]);
	      cout <<"Dijet Trigger Name : " << TrigPath << "  Trigger Bits : " << trigRes.product()->accept(i_Trig) <<endl;
	      
	    }
	  }
	} //nDiJjetHLTmx
	//	 if(TrigPath.Index("HLT_Mu3_PFJet200DeepCSV_1p59_v") >=0) passHLT_Mu3_PFJet200DeepCSV_1p59=true; 
	//	 if(TrigPath.Index("HLT_Mu3_L1SingleJet180_v") >=0)passHLT_Mu3_L1SingleJet180=true;
	//	 if(TrigPath.Index("HLT_PFJet200DeepCSV_1p59_v") >=0)passHLT_PFJet200DeepCSV_1p59=true;
	//Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
      }
    }
  } //!trigRes.failedToGet()
  */
  
  
  
} // end of event


// ------------ method called once each job just before starting event loop  ------------
void
Triggereffi::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
Triggereffi::endJob()
{
  cout << "Total Event =" << nevt << endl;
  cout << "Total  triggerd event  =" << ntrig << endl;
  // cout << "Single Mu trigger  : "<<nSingleMutig<< endl; 
}

void
Triggereffi::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
// Initialize hltConfig

  bool changed(true);
  if (hltPrescaleProvider_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
    hltConfig.dump("Triggers");
    hltConfig.dump("PrescaleTable");

    for (unsigned int ij=0; ij<nDiJetHLTmx; ij++) {
      l1pres[ij] = hltpres[ij]=-7;
    }

       } else {
         }
}



//Delta R function
double dR(double eta1, double phi1, double eta2, double phi2) {
  double DR=0.0;
  DR=sqrt(pow((eta1- eta2),2) +pow(PhiInRange(phi1 - phi2),2));
  //  cout << "DRRRRR " << eta1 << phi1 << eta2 << phi2 <<endl;
  return DR;
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


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Triggereffi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  
  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Triggereffi);
