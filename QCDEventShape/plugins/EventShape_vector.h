#ifndef EventShape_vector_h
#define EventShape_vector_h

/*  \class EventShape
*
*  Class that, given the 4-momenta of the objects in the event, 
*  allows to calculate event shapes (see below)
*
*  Authors: Matthias Weber                 Date: July/October 2007
*
*/

#include <vector>
#include <iostream>
#include <cmath>
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using namespace std;
using namespace CLHEP;
const unsigned int nevtvar=33;
class EventShape_vector {

public:

  //constructor taking as argument vectors of 
  //HepLorentzVector of input Objects 
  //in the event to calculate the event shapes 
  // eta_central : Maximum absolute rapoidity
  //by rapidity the choice between using the rapidity 
  //or the pseudorapidity Eta is given
  // 0: pseudorapidity eta is used (recommended)
  // 1: rapidity is used - please note that the central region
  // is defined in terms of the real rapidity in this case
  // so you should use it only in the case of 
  // massless input objects
  //the choices are then defined by SetMethod(rapidity);

  // nmn : Returns the number of jets withing |eta_central|
  //       -ve number if all jets are in same hemesphere
  // Alog : +1 for KT alogorith
  //      : -1 for anti-Kt algorithm
  // Rad : Cone radius

  EventShape_vector(vector<HepLorentzVector> four_vector, double eta_central, int rapidity, int nmn, double Rad=0.4, int algo=1);


  //Destructor
~EventShape_vector(){};

 std::vector<double> getEventShapes(); 
//returns the values of the event shapes

//============================================================================
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//============================================================================
//beware if you have less than three (central) input objects the 
//directly global (central) three-jet resolution thresholds will be set 
//to -1.0 per default
//in the case of no central input object the central variables will 
//be set to -1.0 per default
//-1.0 is outside of values of all event shape variables so if you have 
//events which give -1.0 as value for an event shape,they had 
//less input momenta than needed for the calculations
//============================================================================

//===============================================================
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//transverse thrust values are the tau values, tau = 1-thrust
//in the following order 
//00. directly global transverse thrust
//here the value of tau = 1 - thrust is returned
//01. directly global thrust minor
//02. directly global three-jet resolution threshold
//03. central transverse thrust
//04. central transverse thrust with exponentially suppressed forward term 
//her the value tau + exp is returned here
//05. central transverse thrust with recoil term
//here the value of tau + rec is returned
//06. central thrust minor
//07. central thrust minor with exponentially suppressed forward term
//08. central thrust minor with recoil term
//09. central total jet mass
//10. central total jet mass with exponentially suppressed forward term
//11. central total jet mass with recoil term
//12. central heavy jet mass
//13. central heavy jet mass with exponentially suppressed forward term
//14. central heavy jet mass with recoil term
//15. central three-jet resolution threshold
//16. central three-jet resolution threshold with exponentially suppressed forward term
//17. central three-jet resolution threshold with recoil term
//18. central total jet broadening 
//19. central total jet broadening with exponentially suppressed forward term
//20. central total jet broadening with recoil term
//21. central wide jet broadening 
//22. central wide jet broadening with exponentially suppressed forward term
//23. central wide jet broadening with recoil term
//24. central total transverse jet mass
//25. central total transverse jet mass with exponentially suppressed forward term
//26. central total transverse jet mass with recoil term
//27. central heavy transverse jet mass
//28. central heavy transverse jet mass with exponentially suppressed forward term
//29. central heavy transverse jet mass with recoil term
//30. Transverse Spheticity
//31. Transvesse C-paramter
//32. Number of jets in central region
//=============================================================================
 void getSphericity(double* );
 //return transverse sphericity and C-parameter

 std::vector<double> getThrustAxis(); //returns the global thrust axis Nx, Ny, Nz=0
 // std::vector<double> getThrustAxis_C(); //returns the central thrust axis Nx, Ny, Nz=0 

 void SetThrustDir(std::vector<double>, double val=0.9); //Store input direction of thrust vector for PFCandidates

 //eta_c: choice of the central region
 //recommended: the two hardest jets should be within the central region

 // void setAlgo(int alg) { m_algo = alg; }
 void setRadius(double rad) { m_radSq= rad*rad; }

 void SetEtac(double eta_central){
   eta_c = eta_central;
 }


 //wether to take the rapidity y (rap==1)  or the pseudorapidity eta (rap==0)
 void SetMethod(int rapidity){
   rap = rapidity;
 }

 private:

 int calculate();
 std::vector<HepLorentzVector> Object_4v;
 std::vector<double> testval;

 std::vector<double> Object_P;
 std::vector<double> Object_Pt;
 std::vector<double> Object_E;
 std::vector<double> Object_Phi;
 std::vector<double> Object_Eta;
 std::vector<double> EventShapes;
 std::vector<double> ThrustAxis;
 // std::vector<double> ThrustAxis_c;
 std::vector<double> iniDir;
 double mncosthe;

 double eta_c;
 int rap;
 int y3_recom;
 unsigned int nmnjet;
 // int m_algo;
 double m_radSq;
 //returns the difference in phi between two vectors
 double DeltaPhi(double, double);

 double max (double, double);

 double min (double, double);

 //the lorentz scalar product
 double lorentz_sp (std::vector<double>, std::vector<double>);

 //calculates the three-jet resolutions
 double three_jet_res(std::vector<HepLorentzVector>, int, int);

 //calculates the thrust axis and the tau values
 std::vector<double> Thrust_calculate(std::vector<HepLorentzVector>);

};

#endif


