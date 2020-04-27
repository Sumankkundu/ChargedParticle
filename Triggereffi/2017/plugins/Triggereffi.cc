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
// Original Author:  Tanmay Sarkar
//         Created:  Wed, 11 Jan 2017 06:26:05 GMT
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

//#include "Test/QCDEventShape/plugins/EventShape_vector.h"


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

using namespace edm;
using namespace reco;
using namespace std;
using namespace CLHEP;
using namespace trigger;
using namespace math;

const int nHLTmx=480;                      //Total Number of Trigger
static const  int nMuHLTmx=16;     // muon trigger number
//static const  int nMuHLTmx=10;     // muon trigger number
static const int nDiJetHLTmx=8;   //DiJet Trigger no.
static const int nJetHLTmx=8;   //DiJet Trigger no.

const int njetetamn=1;                     // one eta space is choosen 

const char* jethlt_name[nJetHLTmx]={"HLT_PFJet80_v",
				    "HLT_PFJet140_v",
				    "HLT_PFJet200_v",
				    "HLT_PFJet260_v",
				    "HLT_PFJet320_v",
				    "HLT_PFJet400_v",
				    "HLT_PFJet450_v",
				    "HLT_PFJet500_v"};

const char* dijethlt_name[nDiJetHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};                           // Di-Jet Trigger Name

const char* muhlt_name[nMuHLTmx]={"HLT_IsoMu17_eta2p1_",
"HLT_IsoMu20_eta2p1_",
"HLT_IsoMu24_eta2p1_",
"HLT_IsoMu20_",
"HLT_IsoMu27_",
"HLT_IsoTkMu20_eta2p1_",
"HLT_IsoTkMu20_",
"HLT_IsoTkMu24_eta2p1_",
"HLT_IsoTkMu27_",
"HLT_Mu20_",
"HLT_Mu27_",
"HLT_Mu45_eta2p1_",
"HLT_Mu50_",
"HLT_TkMu20_",
"HLT_TkMu24_eta2p1_",
"HLT_TkMu27_"};// Muon trigger name
/*"HLT_IsoMu17_eta2p1_", //63
				  "HLT_IsoMu20_", //72
				  "HLT_IsoMu24_eta2p1_", //77
				  "HLT_IsoMu27_", //78
				  "HLT_IsoTkMu20_", //79
				  "HLT_IsoTkMu20_eta2p1_", //80
				  "HLT_IsoTkMu24_eta2p1_", //81
				  "HLT_IsoTkMu27_", //82
				  "HLT_IsoMu18_",   //449
				  "HLT_IsoMu22_"
*/
//};   //453     // Muon trigger name



double l1Pt[nDiJetHLTmx] = {0,35,60,90,120,170,170,170};
double jethlt_thr[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
double leadingPtThreshold[nDiJetHLTmx]={60,80,140,200,260,320,400,500};
const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
double etarange[njetetamn] ={2.4};
bool trgpas[nDiJetHLTmx];//={0,0,0,0,0,0,0,0};
// class declaration
//

int l1pres[nDiJetHLTmx], hltpres[nDiJetHLTmx], compres[nDiJetHLTmx];    // For jet
int l1mupres[nMuHLTmx], hltmupres[nMuHLTmx];                           // For Muons

//double deltaR(double eta1, double phi1, double eta2, double phi2);
//static double PhiInRange(const double& phi);

int hlt_list[nHLTmx];
int hltmu_list[nHLTmx]; //[nMuHLTmx];
int hltjet_list[nHLTmx]; //[nDiJetHLTmx];
int hltdijet_list[nHLTmx]; //[nDiJetHLTmx];
int hltdijet_mimic[nDiJetHLTmx]; //[nDiJetHLTmx];

std::map<std::string, bool> fired;
int nevt;


// DeltaR function
template <class T, class U> static double deltaR(const T& t, const U& u);
double deltaR(double eta1, double phi1, double eta2, double phi2);
int sbitx(unsigned ival, int ibit);



struct triggervar{
  HepLorentzVector trg4v;
  bool            both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
};
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
//

class Triggereffi : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit Triggereffi(const edm::ParameterSet&);
  ~Triggereffi();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //    TH1F* hlt_dijet[nDiJetHLTmx];
  //    TH1F* hlt_dijetref[nDiJetHLTmx];
  
  //Dijet trigger efficiency
  TH1F* hlt_dijettag[nDiJetHLTmx][njetetamn];
  TH1F* hlt_dijetprob[nDiJetHLTmx][njetetamn];
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
  generator1_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("evtinfo"))),
  MuonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),  //Edited
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
   //usesResource("TFileService");
  // theRootFileName = iConfig.getUntrackedParameter<string>("RootFileName");
  edm::Service<TFileService> fs;
  theHLTTag = iConfig.getUntrackedParameter<string>("HLTTag", "HLT");
  char name[200];
  char title[200];
  
  
  /*   for (int ij=0; ij<nDiJetHLTmx; ij++) {
  // for (int jk=0; jk<1; jk++) {
  sprintf(name, "hlt_dijetref_%i", ij);
  sprintf(title, "dijetref probed P_T : (%s) |i#eta|<", jethlt_name[ij]);
  hlt_dijetref[ij] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
  hlt_dijetref[ij]->Sumw2();
  sprintf(name, "hlt_dijetprob_%i", ij);
  sprintf(title, "dijet probed P_T : (%s) |i#eta|<", jethlt_name[ij]);
  hlt_dijet[ij] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
  hlt_dijet[ij]->Sumw2();
  //}
  }*/
  
  
  //define histogram
  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
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
  }
  
  
  //Define Histogram for efficiency
  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
    for (int unsigned jk=0; jk<njetetamn; jk++) {
      
      sprintf(name, "hlt_dijet_effi_%i_%i", ij, jk);
      sprintf(title, "dijet trigger fired: (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      
      hlt_dijettrg_fired[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettrg_fired[ij][jk]->Sumw2();
    }
  }
  
  for (int unsigned ij=0; ij<nDiJetHLTmx; ij++) {
    for (int unsigned jk=0; jk<njetetamn; jk++) {
      
      sprintf(name, "hlt_dijet_all_evt_%i_%i", ij, jk);
      sprintf(title, "dijet trigger All event: (%s) |i#eta|<%g", dijethlt_name[ij], etarange[jk]);
      
      hlt_dijettrg_all_evt[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettrg_all_evt[ij][jk]->Sumw2();
    }
  }
  
  
  
  nevt=0;  //event counter
  
  
}//End constractor


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
  
  //  nevt++;
  //  if (nevt%5000==1) cout <<"TriggerEfficiency::analyze run number = "<<nevt<<endl;
  if (nevt%500==1) cout <<"TriggerEfficiency::analyze run number = "<<nevt<<endl;
  cout << endl ;
  cout << endl ;
  cout <<"TriggerEfficiency::analyze run number = "<< nevt <<endl;
 nevt++ ;
  const char* variab1;
  const char* variab2;
  double aveleadingpt =0;
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerBits_, trigRes);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  
  const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  bool isInEtaRange[njetetamn]={0};
  
  edm::Handle<pat::JetCollection> ak4PFJets;
  iEvent.getByToken(jetSrcToken_, ak4PFJets);
  
  
  
  //Frist Condition : Check for  the Trigger Name For Mu and Jet only
  if (!trigRes.isValid()) return;
  
  bool mutrig = false;
  bool dijettrig = false;
  
 if (trigRes.isValid()) {
    unsigned int size = trigRes->size();
//    cout << "trigger result size = " << size << endl;
//    edm::TriggerNames triggerNames = iEvent.triggerNames(*trigRes);
    
    for(unsigned ij = 0; ij != size; ij++) {
           int ihlt =  trigRes->accept(ij);
        if (ihlt <1) continue;

        std::string name = names.triggerName(ij);
//       if(nevt==0)  {cout <<"inclu:hltobject "<<" "<<ij<<" "<<name<<endl;  }
          variab1 = name.c_str();
         cout << " Trigger fried : " << variab1  << endl;
	 for (unsigned kl=0; kl<nMuHLTmx; kl++) {
	  if ((strstr(variab1,muhlt_name[kl])) && strlen(variab1)-strlen(muhlt_name[kl])<5) {
	    
               cout << "muon trigger found :" << trigRes->accept(ij) << variab1 <<endl;
                 mutrig=true; break;
	  
            }
	}
	
	for (unsigned kl=0; kl<nDiJetHLTmx; kl++) {
	  if ((strstr(variab1,dijethlt_name[kl])) && strlen(variab1)-strlen(dijethlt_name[kl])<5) {
               	cout << " Dijet trigger found :" << variab1 << endl;
                   dijettrig = true; break;
	  }
	}
    }    // loop over triggerResult size
       
  } //if (trigRes.isValid())
  
  if (!mutrig) return;  
  if (!dijettrig) return;
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(MuonToken_,muons);
  
  //  edm::Handle<edm::View<reco::MuonCollection> > muons;
  //  edm::Handle<edm::View<reco::Muon> > muons;
  //  iEvent.getByLabel("muons",muons);
  
  //  edm::Handle<reco::PFJetCollection> PFJets;
  //  iEvent.getByLabel("ak4PFJets", PFJets); 
  //  iEvent.getByLabel("ak4PFJetsL1FastL2L3Residual", PFJets);
  //  if ((!muons.isValid()) || (!PFJets.isValid())) return;
  
  unsigned int musz=muons->size();
  if(musz <= 0) return;
  cout << "Muons Trigger Ok with muon  Numbers : " <<musz<< endl;
  
  
  //Calculate average Pt 
  if (ak4PFJets.isValid() &&  ak4PFJets->size()>1) {
    for (int iet=0; iet<njetetamn; iet++) {
      isInEtaRange[iet] = true;
    }
    
    for (int ij=0; ij<2; ij++) {
      for (int iet=0; iet<njetetamn; iet++) {
	if (abs((*ak4PFJets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
      }
      double NHF = (*ak4PFJets)[ij].neutralHadronEnergyFraction();
      double NEMF = (*ak4PFJets)[ij].neutralEmEnergyFraction();
      double CHF = (*ak4PFJets)[ij].chargedHadronEnergyFraction();
      double CEMF = (*ak4PFJets)[ij].chargedEmEnergyFraction();
      int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
      int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
      int CHM = (*ak4PFJets)[ij].chargedMultiplicity();
      bool looseJetID =false;
      if(abs((*ak4PFJets)[ij].eta())<=3.0){
	if( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs((*ak4PFJets)[ij].eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) ||
						      
						      abs((*ak4PFJets)[ij].eta())>2.4) ) looseJetID =true;
      } else {
	if( (NEMF<0.90 && NumNeutralParticles>10) ) looseJetID =true;
      }
      
      if (abs((*ak4PFJets)[ij].eta())>3.0) {looseJetID = false;}
      if ((*ak4PFJets)[ij].pt()<30.0) {looseJetID = false;}
      
      if (looseJetID) {
	aveleadingpt +=(*ak4PFJets)[ij].pt();
      } else {
	aveleadingpt -=100000;
      }
    }
    aveleadingpt /=2.0;
    
  } // calculation of average pt
  
  if( aveleadingpt < 0) return;
//  nevt++;
  
  
  //Preslace Factor Calculation
  // trgpas[nDiJetHLTmx]={false};
  for (int jk=-1; jk<nDiJetHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str();
      if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) ||
	  (jk>=0 && strstr(variab1,dijethlt_name[jk]) && strlen(variab1)-strlen(dijethlt_name[jk])<5)) {
	if (jk>=0) {
	  compres[jk] = triggerPrescales->getPrescaleForIndex(ij);
	  //       if (trigRes->accept(ij)) {trgpas[jk] = true;}
	  
	  
	}
      }
    }
  }//calculation of prescale
  
  
  
  
  
  //Second condition:  select events with two reconstructed jets above a certain (low) threshold
  if ((!ak4PFJets.isValid()) ||  ak4PFJets->size() <2) return;
  if ((*ak4PFJets)[0].pt()<30.0 || (*ak4PFJets)[1].pt()<30.0)  return ;
  if (fabs((*ak4PFJets)[0].eta())>2.4 || fabs((*ak4PFJets)[1].eta())>2.4)  return ;
  if (!isInEtaRange[0] && (aveleadingpt<30)) return;
  
  cout << "1st Jet PT=" << (*ak4PFJets)[0].pt() <<" : 2nd Jet PT= " << (*ak4PFJets)[1].pt()<<endl;
  cout << "1st Jet eta=" << fabs((*ak4PFJets)[0].eta()) <<" : 2nd Jet eta= " << fabs((*ak4PFJets)[1].eta())<<endl;
  
 

//Fill event selected event
 for (int iet=0; iet<njetetamn; iet++) {
    for (int jk=0; jk<nDiJetHLTmx; jk++) {
      //   if(aveleadingpt>jethlt_thr[jk]) {
      hlt_dijettrg_all_evt[jk][iet]->Fill(aveleadingpt, compres[jk]);
      // }
    }
  } //Fill the all valid event */
  
  
  
  
 //Third Condition : match the two jets with the trigger objects in eta-phi space
 bool delta_object=false;
 
 for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
   obj.unpackPathNames(names);
//   std::cout << "\t   Collection: " << obj.collection() << std::endl;
   //---------------------------------------------
        std::vector<std::string> pathNamesAll = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
     //   std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";  


      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
            bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
       //     std::cout << "   " << pathNamesAll[h];
            if (isBoth) std::cout << "(L,3)";
            if (isL3 && !isBoth) std::cout << "(*,3)";
            if (isLF && !isBoth) std::cout << "(L,*)";
            if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        }


      for (unsigned ih = 0, n = pathNamesAll.size(); ih < n; ++ih) {
              variab2 = pathNamesAll[ih].c_str();
               cout << "obejct pathNamesAll "  << " Name "<< variab2 << endl;
     
          for (int jk=0; jk<nDiJetHLTmx; jk++) {
                if (strstr(variab2,dijethlt_name[jk]) && strlen(variab2)-strlen(dijethlt_name[jk])<5){
	 
         	cout<< "MATCHED Trigger Object ****************"<< "jjjk= "<< jk << " Name "<< variab2 << endl;
	     // for (unsigned h = 0; h < obj.filterIds().size(); ++h) {cout << "ids " << obj.filterIds()[h]<<endl;}
          	 //  for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {cout << "Label " << obj.filterLabels()[h]<< endl;}
	 
	 //-------------------------------------------------------- 
	 delta_object = true;
	 double dr1 = deltaR((*ak4PFJets)[0].eta(), (*ak4PFJets)[0].phi(), obj.eta(), obj.phi());
	 double dr2 = deltaR((*ak4PFJets)[1].eta(), (*ak4PFJets)[1].phi(), obj.eta(), obj.phi());
         
         double p1 = (*ak4PFJets)[0].phi();
         double p2 = obj.phi();
         float  n1 = (*ak4PFJets)[0].eta();
         float  n2 = obj.eta();
         auto   dp = std::abs(p1 - p2);
         if (dp > float(M_PI)){dp -= float(2 * M_PI);}
         double dr1test = sqrt(pow((n1-n2), 2) +pow(dp,2));
	 
        cout <<"Jets eta =" <<(*ak4PFJets)[0].eta() << " : "<<(*ak4PFJets)[1].eta() << endl;
        cout <<"Jets phi =" <<(*ak4PFJets)[0].phi() << " : "<<(*ak4PFJets)[1].phi() << endl;
        cout <<"Object eta =" << obj.eta() << ": Object Phi =" << obj.phi() << endl;
	   cout << "deltaR = " << dr1 << "  :  " << dr2<<" -------DELTA R test=" << dr1test <<endl;
	 if( dr1>=.5 || dr2 >=.5 ) {delta_object = false;
	   //   cout << delta_object << endl;
	   
	 }
       }//if (strstr(variab2,dijethlt_name[jk]) && strlen(variab2)-strlen(dijethlt_name[jk])<5){
     } //nDiJetHLTmx      
   }
 }//object loop
  //  if(!delta_object) return;  
 cout <<"Delta Condtion: " << delta_object << endl;
 
 
 
 //Fourth condition : mimic the DijetPFJet trigger decision (average of the two trigger objects pT above a certain threshold)
 if (trigRes.isValid()) {
   //  unsigned int size = trigRes->size();
//   edm::TriggerNames triggerNames = iEvent.triggerNames(*trigRes);
 
  for (unsigned jk=0; jk<nDiJetHLTmx; jk++) {
     hltdijet_mimic[jk]=-1;
   }
   
   for (unsigned ij=0; ij<nDiJetHLTmx; ij++){
     for (unsigned jk=0; jk<nHLTmx; jk++) {
      std::string name = names.triggerName(ij);
   //    std::string name = triggerNames.triggerName(jk);
       variab2 = name.c_str();
       if ((strstr(variab2,dijethlt_name[ij])) && strlen(variab2)-strlen(dijethlt_name[ij])<5) {
	 hltdijet_mimic[ij]=ij; break;}
     }
   }
   
 }
 
 
 
 
 //................................
 // trgpas[nDiJetHLTmx]={0,0,0,0,0,0,0,0};
 /*  bool trg_prev=false;
 // trgpas[nDiJetHLTmx]={false};
 for (int jk=-1; jk<nDiJetHLTmx; jk++) {
 for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
 std::string name = names.triggerName(ij);
 variab1 = name.c_str();
 if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) ||
 (jk>=0 && strstr(variab1,dijethlt_name[jk]) && strlen(variab1)-strlen(dijethlt_name[jk])<5)) {
 if (jk>=0) {
 compres[jk] = triggerPrescales->getPrescaleForIndex(ij);
 if (trigRes->accept(ij)) {trgpas[jk] = true;}
 
 //   cout<< "OK1"<< variab1 <<endl;
 trg_prev = trigRes->accept(ij);           
 
 
 if(trg_prev){
 for (int iet=0; iet<njetetamn; iet++) {
 if (isInEtaRange[iet] && (aveleadingpt>30)) {
 
 //		if(aveleadingpt>=jethlt_thr[jk]) {hlt_dijettrg_fired[jk][iet]->Fill(aveleadingpt, compres[jk]);
 //		  cout << " New Trigger Found for event number : " << nevt << endl;
 }
 } 
 }
 
 trg_prev = trigRes->accept(ij);
 break;
 } else {
 trg_prev = trigRes->accept(ij);
 break;
 }
 }
 }
 } */
 //for (int iet=0; iet<njetetamn; iet++) {
 //   for (int jk=0; jk<nDiJetHLTmx; jk++) {
 //     if(aveleadingpt>jethlt_thr[jk]) {hlt_dijettrg_all_evt[jk][iet]->Fill(aveleadingpt, compres[jk]); 
 //     }
 //   }
 // } //Fill the all valid event */
 
 
 
 
 
 //Fill the fired event
 for (int iet=0; iet<njetetamn; iet++) {
   for (int jk=0; jk<nDiJetHLTmx; jk++) {
     if(aveleadingpt>=jethlt_thr[jk] && hltdijet_mimic[jk]>=0) {hlt_dijettrg_fired[jk][iet]->Fill(aveleadingpt, compres[jk]);
 //      cout << " New Trigger Found for event number : " << nevt << endl;
      }
   }
 }
 
 
 
 
 
 cout<< "end ievent loop     ...........................   "<<endl;
 cout <<endl;
 
}//end of event loop


// ------------ method called once each job just before starting event loop  ------------



void 
Triggereffi::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Triggereffi::endJob() 
{
  cout << "Total Evnet =" << nevt << endl;
  //Devide by the Total event
  /*for(int ij=0; ij<nDiJetHLTmx; ij++){
    for (int jk=0; jk<njetetamn; jk++) {
    int tmpnbn = hlt_dijettrg_effi[ij][jk]->GetNbinsX();
    for (int ix=0; ix<tmpnbn; ix++) {
    hlt_dijettrg_effi[ij][jk]->SetBinContent(ix+1, hlt_dijettrg_effi[ij][jk]->GetBinContent(ix+1)/nevt);
    }
    }
    }*/
  
}
/*
  void
  Triggereffi::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
  {
  cout << "Write test 4 = ok " << endl;
  bool changed(true);
  if (hltPrescaleProvider_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
  hltConfig.dump("Triggers");
  hltConfig.dump("PrescaleTable");
  
  } else {
  }
  
  }
*/



/*
double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1- eta2, 2) +pow(PhiInRange(phi1 - phi2),2));
  
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

*/



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Triggereffi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Triggereffi);
