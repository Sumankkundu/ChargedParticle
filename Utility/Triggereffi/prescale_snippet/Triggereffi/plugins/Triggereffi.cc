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
double jethlt_thr[nDiJetHLTmx]={40,60,80,140,200,260,320,400,500};
double leadingPtThreshold[nDiJetHLTmx]={40,60,80,140,200,260,320,400,500};
const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
double etarange[njetetamn] ={2.4};
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
//  char name[200];
 // char title[200];
  
  
  
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
  
  
  
  if (nevt%500==1) cout <<"TriggerEfficiency::analyze run number = "<<nevt<<endl;
  const char* variab1;
 // const char* variab2;
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
 // bool mutrig=false;
 // bool dijtrig=false;
  edm::TriggerNames triggerNames = iEvent.triggerNames(*trigRes);   //Get The trigger name for TriggerResult Input
  
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
      int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
     // int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
      int CHM = (*ak4PFJets)[ij].chargedMultiplicity();
      bool TightJetID =false;
      if(abs((*ak4PFJets)[ij].eta())<=2.7){
                     if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && (abs((*ak4PFJets)[ij].eta())<=2.4 && CHF>0 && CHM>0 )) TightJetID =true;
      } else {
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
  

  //Prescale Factor Calculation
  int preL1(-1), preHLT(-1); //prescale(-1);

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
            cout << " By index " << triggerPrescales->getPrescaleForIndex(ij) << "    L1 : "<< preL1 << "   HLT :" << preHLT <<endl;
     }
    }
  }//calculation of prescale
   if (!isInEtaRange[0] && (aveleadingpt<30)) return;
 
  
  
  
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
