
//CM energy 13Tev
//use the fastjet for jet finding

//Author : S.k.kundu
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"                   // This is the minimal interface needed to access FastJet.
#include "fastjet/ClusterSequenceArea.hh"         //clusterSequence with the area support
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"                                 // for histrograming
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"                     // ROOT, for saving file.
#include "TTree.h"                     //for Tree file 
#include "TROOT.h"
#include <vector>

#include "EventShape_vector.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"


using namespace std;
using namespace fastjet;
using namespace Pythia8;

using namespace std;
using namespace CLHEP;

///////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
 

  int const ntype = 3;         // Jet & Charage particles , all Particles
  const int nHLTmx = 8; //hist_averjet12_pt Range
  const int njetetamn = 1;  //eta value used 2.4
  
  float lowtrig = 83;
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  char histname[100], name[100], title[100], Axisname[100] , RecoBinName[100], GenBinName[100];
  
   double leadingPtThreshold[nHLTmx+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger
   double htbins[nHLTmx+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0};
   double recohtbins[nHLTmx+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0};
   
//   TDirectory *TUnfoldBinng2D =new TDirectoryFile("analyzeBasicPat2D","TUnfoldBinning 1D  for 2D Historgams");


   Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass



    int Mnbinsx0[nvar]={0,0,0,16,0,0,0,0,
                  0,16,0,0,0,0,0,7,
                  0,0,12,0,0,0,0,0,
                  13,0,0,0,0,0,0,0};

    int  Mnbinsx1[nvar]={0,0,0,9,0,0,0,0,
                  0,9,0,0,0,0,0,4,
                  0,0,8,0,0,0,0,0,
                  7,0,0,0,0,0,0,0};


   //Gen Level
   const int Mnmxbins=16;
   double Mbinrngs0[nvar][Mnmxbins+1] = {{},{},{},
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

   double Mbinrngs1[nvar][Mnmxbins+1]={{},{},{},
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
   
   TH1* h_genvar_1D[ntype][nusedvar][nHLTmx]; //Gen
   
for(int ity=0; ity <ntype; ity++){
     for(int ivar=0; ivar < nusedvar ; ivar++){
       for(int ipt = 0 ; ipt < nHLTmx ; ipt++){        
   
       sprintf(histname, "gen_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]);    sprintf(title, "gen type %i pt%i eta0 %i", ity, ipt, var[ivar]);
       h_genvar_1D[ity][ivar][ipt] = new TH1D(histname, title, (ity==0) ? Mnbinsx0[var[ivar]] : Mnbinsx1[var[ivar]], (ity==0) ? Mbinrngs0[var[ivar]] : Mbinrngs1[var[ivar]] );
       }
     }
   }

  // Create Pythia instance and set it up to generate hard QCD processes
  Pythia pythia;                          // Generator
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("HardQCD:all = on");  // 
  pythia.readString("PhaseSpace:pTHatMin = 15");
  pythia.readString("PhaseSpace:pTHatMax = 7000.0");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.0");
  pythia.readString("Tune:pp = 14");
  pythia.readString("Tune:ee = 7");
  pythia.readString("Beams:eCM = 13000.");

//----------------------------CP5 Tune
//in argumant 
//-------------------------------------
  pythia.readFile(argv[1]);
  pythia.init();
  
  TFile* outFile = new TFile(argv[2], "RECREATE");
  int nEvent = pythia.mode("Main:numberOfEvents");
  //int nEvent = 100;
  //int nEvent = argi[1];
   TDirectory *Hist_dir =new TDirectoryFile("analyzeBasicPat","1D Historgams");

////////////////////////////////////////////////////////////////////////
  
//Define and initialaize an Tree file 
  TTree *Event_tree = new TTree("Event_tree","Event_Information");
  Double_t PX, PY, PZ, P_T, E_Jet, eta, phi, weight; //Define the Variable for branch
  //Double_t weight;  
  vector<Double_t> Px, Py, Pt, Pz, Eta, Phi, Ejet;
  Int_t  Njets , sjtevn=0 ,djtevn =0,mjtevn=0;
  
  Event_tree->Branch("Px", "vector<Double_t>", &Px);
  Event_tree->Branch("Py", "vector<Double_t>", &Py);
  Event_tree->Branch("Pz", "vector<Double_t>", &Pz);
  Event_tree->Branch("Pt", "vector<Double_t>", &Pt);
  Event_tree->Branch("Eta", "vector<Double_t>", &Eta);
  Event_tree->Branch("Phi", "vector<Double_t>", &Phi);
  Event_tree->Branch("Ejet", "vector<Double_t>", &Ejet);
  Event_tree->Branch("weight", &weight, "weight/D");
  Event_tree->Branch("Njets", &Njets, "Njets/I");


  TH1D *hist_inclusive_pt = new TH1D("hist_inclusive_pt","PT_of_all_jets",400,20,2020.0);  hist_inclusive_pt->Sumw2();
  TH1D *hist_inclusive_eta = new TH1D("hist_inclusive_eta","Eta_of all_Jets",100,-2.5, 2.5);   hist_inclusive_eta->Sumw2();
  TH1D *hist_inclusive_phi = new TH1D("hist_inclusive_phi","phi_of_all_Jets",100,-3.14,3.14);  hist_inclusive_phi->Sumw2();
  TH1D *hist_averjet12_pt = new TH1D("genjetallave_pt_0","HT2", 400,20.0 ,2020.0);  hist_averjet12_pt->Sumw2();
  TH1D *hist_jet1_pt = new TH1D("genjet1_pt_0","PT_of_1st_Jet",400,20.0,2020);  hist_jet1_pt->Sumw2();
  TH1D *hist_jet1_eta = new TH1D("genjet1_eta","Eta of 1st Jet",100,-2.5, 2.5);  hist_jet1_eta->Sumw2();
  TH1D *hist_jet1_phi = new TH1D("genchg1_phi","#phi_{gencharge_jet1}",100,-M_PI, M_PI);  hist_jet1_phi->Sumw2();
  
  TH1D *hist_jet2_pt = new TH1D("genjet2_pt_0","PT_of_second_leading_jet",400,20.0,2020.0);  hist_jet2_pt->Sumw2();
  TH1D *hist_jet2_eta =  new TH1D("genjet2_eta","second_jet_eta",100,-2.5, 2.5);  hist_jet2_eta->Sumw2();
  TH1D *hist_jet2_phi = new TH1D("genchg2_phi","#phi_{gencharge_jet2}",100,-M_PI, M_PI);  hist_jet2_phi->Sumw2();
  
  TH1D *hist_deltapt12 = new TH1D("genjetdpt_0","Delta_PT_of_ist_and Second_jets",100,20.0,500.0);  hist_deltapt12->Sumw2();
  TH1D *hist_detalphi12 = new TH1D("genjetdphi_0","Delta_phi_of_ist_and_2nd_jets",100,-3.14,3.14);  hist_detalphi12->Sumw2();
  TH1D *hist_comppt1pt2 = new TH1D("genjetpt2bypt1_0","Component_of_2nd_jet_on_leading_jet",60,0.0,1.0);  hist_comppt1pt2->Sumw2();
  
  TH1D *Numjet = new TH1D("njets_0","njets_0",30,0.0,30);  Numjet->Sumw2();
  TH1D *hist_ncharged = new TH1D("ncharges_0","Number_of_charged particles",600,0.0,600);  hist_ncharged->Sumw2();
  
  double R = 0.4;                  //Jet Radius 
  
  JetDefinition jet_def(antikt_algorithm, R);              //Choose the jet difinition in :jet_def
  vector<PseudoJet> fjInputs;   //Define Fastjet input :fjInput
  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {                          /*start event loop*/
  if (!pythia.next()) continue;
   
   fjInputs.resize(0);   //Reset Fastjet input 

   for (int i = 0 ; i < pythia.event.size() ; ++i){

     if (!pythia.event[i].isFinal())        continue;                                                                      //allow final state particle only
     fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(),pythia.event[i].pz(), pythia.event[i].e());
     if (pythia.event[i].isCharged()){ particle.set_user_index(-1);                                                        // Put index -1 for charged particles 
       }else if (pythia.event[i].isNeutral()){ particle.set_user_index(pythia.event[i].isNeutral());};                     //Put index +1 for neutral particles
    
    //fjInputs.push_back(PseudoJet(pythia.event[i].px(), pythia.event[i].py(),pythia.event[i].pz(), pythia.event[i].e()));            
    fjInputs.push_back(particle);            
      
   } 
   
   //check the event  
   if (fjInputs.size() == 0){
     cout << "Error:event with no final state particles"<<endl;
     continue; 
   }
    
    //Run Fastjet algorithm//**********************************************************************************
   ClusterSequence clust_seq(fjInputs, jet_def);  // run the jet clustering with the above jet 
   
   //sorted the jet with the minimum pT
   double ptmin =30.0;
   Selector jet_selector = SelectorAbsEtaMax(4.7);
   vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(clust_seq.inclusive_jets(ptmin)));


   Njets =0;
   int ihtbin=-1;
   int numjet = inclusive_jets.size();                    //number of jets in any events 
   int nchar = 0;                    //number of charged particles
   weight=0.0;


////////////////////////////////////////////
     if(inclusive_jets.size()<2) continue ;
     if((inclusive_jets[0].perp()) < 30 || fabs(inclusive_jets[0].eta()) >2.4) continue;
     if((inclusive_jets[1].perp()) < 30 || fabs(inclusive_jets[1].eta()) >2.4) continue;
     if(((inclusive_jets[0].perp() + inclusive_jets[1].perp())/2) < lowtrig  ) continue;

        double  ht2 = (inclusive_jets[0].perp() + inclusive_jets[1].perp())/2;
        weight = pythia.info.weight();

        for(Int_t iipt =0 ; iipt <nHLTmx+1 ; iipt++){                                                     // Loop to check the HT range
          if(ht2 > htbins[iipt] && ht2 < htbins[iipt+1]){ihtbin = iipt;
            };
	}
       //       cout << " Ht2: "<< ht2 << "   save at "<< ihtbin << "th phasesape"  <<endl;  


        hist_averjet12_pt->Fill(ht2,weight);
        hist_jet1_pt->Fill(inclusive_jets[0].perp(),weight);
        hist_jet2_pt->Fill(inclusive_jets[1].perp(),weight);

        hist_jet1_eta->Fill(inclusive_jets[0].eta(),weight);
        hist_jet2_eta->Fill(inclusive_jets[1].eta(),weight);
        Numjet->Fill(inclusive_jets.size(),weight);
        hist_comppt1pt2->Fill(((inclusive_jets[1].perp() * sin(inclusive_jets[0].phi_std() - inclusive_jets[1].phi_std()))/inclusive_jets[0].perp(), weight));
        hist_deltapt12->Fill((inclusive_jets[0].perp() - inclusive_jets[1].perp()),weight);
        hist_detalphi12->Fill(inclusive_jets[0].phi_std() - inclusive_jets[1].phi_std(),weight);

    

    vector<HepLorentzVector> gen_object[ntype];

    
//    cout << endl <<"Event " << iEvent  << " nJets " << inclusive_jets.size() ;
    for (int ijet =0; ijet < inclusive_jets.size(); ijet++){ //Jet Loop

     hist_inclusive_pt->Fill(inclusive_jets[ijet].perp(),weight);
     hist_inclusive_eta->Fill(inclusive_jets[ijet].eta(),weight);
     hist_inclusive_phi->Fill(inclusive_jets[ijet].phi_std(),weight);

      HepLorentzVector tmp4v(inclusive_jets[ijet].px(), inclusive_jets[ijet].py(), inclusive_jets[ijet].pz(), inclusive_jets[ijet].e());                 // create a temporary fourvector

      gen_object[0].push_back(tmp4v);

      vector<PseudoJet> constituents = inclusive_jets[ijet].constituents();
  //    cout << " Jet : "<<  ijet <<" nPartiles " << constituents.size();
      if(constituents.size()<1)continue;

      for (unsigned j = 0; j <constituents.size(); j++) { // All particle loop
              
             HepLorentzVector tmp4v_all(constituents[j].px(), constituents[j].py() , constituents[j].pz(), constituents[j].e());  //create  temporary 4-vector for all particles
             gen_object[2].push_back(tmp4v_all);

             if(constituents[j].user_index() ==-1){  // check for Charged Particles

             nchar++; 
             HepLorentzVector tmp4v_char(constituents[j].px(), constituents[j].py() , constituents[j].pz(), constituents[j].e());  // create 4-vector fpr Charged particles
             gen_object[1].push_back(tmp4v_char);

                      }
         	}
          }   //-------------------------------------Jet Loop
     
 //calculate Eventshape
  if(ihtbin>=0){
      for(int ity =0; ity< ntype ; ity++){
  //    genvar.clear();
       vector<double> genvar;

      EventShape_vector  genevtshape(gen_object[ity], 2.4, 0, 2, 1);
      genvar =  genevtshape.getEventShapes();

      for( int ifill = 0; ifill < nusedvar ; ifill ++){
          h_genvar_1D[ity][ifill][ihtbin]->Fill(genvar[var[ifill]],weight);
                }
            } 
         }

 hist_ncharged->Fill(nchar,weight);
  } //--------------------------------end of the event loop

 /*****************************************************************************/
 pythia.stat();
 
// Event_tree->Write();
Hist_dir->cd();
 hist_inclusive_pt->Write();
 hist_inclusive_eta->Write();
 hist_inclusive_phi->Write();
 hist_jet1_pt->Write();
 hist_jet2_pt->Write();
 hist_jet1_eta->Write();
 hist_jet2_eta->Write();
 hist_jet1_phi->Write();
 hist_jet2_phi->Write();
 
 hist_averjet12_pt->Write();
 
 hist_deltapt12->Write();
 hist_detalphi12->Write();
 hist_comppt1pt2->Write();
 
 Numjet->Write();
 hist_ncharged->Write();

for(int ity=0; ity <ntype; ity++){
   for(int ivar=0; ivar < nusedvar ; ivar++){
     for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
       h_genvar_1D[ity][ivar][ipt]->Write();
     }
   }
}



 delete outFile;
 return 0;

} //end main program
