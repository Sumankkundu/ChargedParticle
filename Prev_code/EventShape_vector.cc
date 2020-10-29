#include "home/suman/QCD/EventShape_vector.h"

using namespace std;

using std::vector;
using std::cout;
using std::endl;

/*  \class EventShape
*
*  Class that, given the 4-momenta of the objects in the event, 
*  allows to calculate event shapes 
*
*  Authors: Matthias Weber                 Date: July/October 2007
*
*/

EventShape_vector::EventShape_vector(vector<HepLorentzVector> four_vector, double eta_central, int rapidity, int nmn, double rad, int alg) : Object_4v(four_vector), eta_c(eta_central),rap(rapidity), nmnjet(nmn), m_radSq(rad*rad) // , m_algo(alg)
{
  iniDir.push_back(-100.0);
  iniDir.push_back(-100.0);
  mncosthe=-1.;
}


vector<double> EventShape_vector::getEventShapes(){
  this->calculate();
  return EventShapes;
}

vector<double> EventShape_vector::getThrustAxis(){
  //  this->calculate();
  return ThrustAxis;
}

void EventShape_vector::SetThrustDir(std::vector<double> var, double val){
  iniDir = var;
  mncosthe=val;
}

int EventShape_vector::calculate(){
  unsigned int length = (unsigned int) Object_4v.size();
  if (!Object_P.empty()){
    Object_P.clear();
    Object_Pt.clear();
    Object_Eta.clear();
    Object_Phi.clear();
    EventShapes.clear();
    ThrustAxis.clear();
  }
  for(unsigned int ij = 0; ij < Object_4v.size(); ij++){
    Object_P.push_back(0.);
    Object_Pt.push_back(0.);
    Object_Eta.push_back(0.);
    Object_Phi.push_back(0.);
  }

  for(unsigned int ij = 0; ij < nevtvar; ij++){
    EventShapes.push_back(-50.);
  }

  EventShapes.push_back(double(Object_4v.size()));

  for(unsigned int ij = 0; ij < 3; ij++){
    ThrustAxis.push_back(0.);
  }

  for(unsigned int jk =0; jk<length; jk++){
    Object_P[jk] = Object_4v[jk].rho();
    Object_Pt[jk] = Object_4v[jk].perp();

    if(Object_P[jk]>Object_4v[jk].e()  + 1E-4){
      cout << "ERROR!!! Object " << jk <<" has P = " << Object_P[jk] << " which is bigger than E = " << Object_4v[jk].e()  <<" of total length "<< length<<endl;
      return 0;
    }
    //to prevent a division by zero
    if(rap==0){
      Object_Eta[jk] = Object_4v[jk].eta();
    } else if (rap==1) {
      Object_Eta[jk] = Object_4v[jk].rapidity();
    }
    
    if((rap!=0)&&(rap!=1)){
      cout<<"ERROR!!!, The choice to use the rapidity y or the pseudorapidity eta is not set correctly! Change that please!"<<endl;
      return 0;}
    Object_Phi[jk] = Object_4v[jk].phi();
  }
  


  //GMA make logarithm of its
  //||||||||||||||||||||||||||||||||||||||||||||||
  //here the central event shape variables begin||
  //||||||||||||||||||||||||||||||||||||||||||||||

  std::vector<HepLorentzVector> Object_v4_in;
  std::vector<double> Object_Px_in;
  std::vector<double> Object_Py_in;
  std::vector<double> Object_Pz_in;
  std::vector<double> Object_Pt_in;
  std::vector<double> Object_E_in;
  std::vector<double> Object_Et_in;
  std::vector<double> Object_Eta_in;
  std::vector<double> Object_Px_out;
  std::vector<double> Object_Py_out;
  std::vector<double> Object_Pz_out;
  std::vector<double> Object_E_out;
  std::vector<double> Object_Pt_out;
  std::vector<double> Object_Eta_out;

   if (!Object_Px_in.empty()){
     Object_v4_in.clear();
     Object_Px_in.clear();
     Object_Py_in.clear();
     Object_Pz_in.clear();
     Object_Pt_in.clear();
     Object_E_in.clear();
     Object_Et_in.clear();
     Object_Eta_in.clear();
     Object_Px_out.clear();
     Object_Py_out.clear();
     Object_Pz_out.clear();
     Object_Pt_out.clear();
     Object_E_out.clear();
     Object_Eta_out.clear();
   }
  unsigned int nin = 0;
  for(unsigned int ij=0;ij<length;ij++){
    if(fabs(Object_Eta[ij])<eta_c){
      Object_v4_in.push_back(Object_4v[ij]);
      Object_Px_in.push_back(Object_4v[ij].px());
      Object_Py_in.push_back(Object_4v[ij].py());
      Object_Pz_in.push_back(Object_4v[ij].pz());
      Object_E_in.push_back(Object_4v[ij].e());
      Object_Pt_in.push_back(Object_4v[ij].perp());
      Object_Et_in.push_back(sqrt((pow(Object_4v[ij].e(),2)*pow(Object_4v[ij].perp(),2))/(pow(Object_4v[ij].perp(),2)+pow(Object_4v[ij].pz(),2))));
      Object_Eta_in.push_back(Object_Eta[ij]);
      nin+=1;
    } else {
      Object_Px_out.push_back(Object_4v[ij].px());
      Object_Py_out.push_back(Object_4v[ij].py());
      Object_Pz_out.push_back(Object_4v[ij].pz());
      Object_E_out.push_back(Object_4v[ij].e());
      Object_Pt_out.push_back(Object_4v[ij].perp());
      Object_Eta_out.push_back(Object_Eta[ij]);
    }
  }

  if(Object_Px_in.size()!=nin){cout<<"ERROR!!! wrong dimension of in momenta"<<endl;}

  unsigned int nout = length - nin; 
  //  if(nin==0){
  if(nin<nmnjet){
    //     cout<<"WARNING!! Number momenta in the central region is "<<nin<<"  Central event shapes won't be calculated and set = -1.0 per default!"<<endl;
    for(unsigned int i=0; i<nevtvar;i++){
      EventShapes[i]=-50.0;
    }
  }

  EventShapes[nevtvar-1] = nin;
  
  //  if (nin>0){
  if (nin>=nmnjet){
    double P_sum_c = 0; //GMA
    double Pt_sum_c = 0;
    double Eta_cw=0;
    double RecTerm =0;
    double px_sum_in = 0;
    double py_sum_in = 0;
    for(unsigned int ij=0;ij<nin;ij++){
      Pt_sum_c+=Object_Pt_in[ij];
      P_sum_c +=sqrt(pow(Object_Pt_in[ij],2.) + pow(Object_Pz_in[ij], 2.0)); //GMA
      Eta_cw+=Object_Pt_in[ij]*Object_Eta_in[ij];
      px_sum_in+=Object_Px_in[ij];
      py_sum_in+=Object_Py_in[ij];
    }
    Eta_cw=Eta_cw/Pt_sum_c;
    RecTerm=sqrt(pow(px_sum_in,2)+pow(py_sum_in,2))/Pt_sum_c;
    
    double ExpTerm =0;
    for(unsigned int ij=0;ij<nout;ij++){
      ExpTerm+=Object_Pt_out[ij]*exp(-fabs(Object_Eta_out[ij]-Eta_cw));
     }
    ExpTerm = ExpTerm/Pt_sum_c;
    
    
    //the central global transverse thrust centrthr is calculated
    double centrthr=0;

    std::vector<double> Thrust_central = Thrust_calculate(Object_v4_in);

    for(unsigned int ij=0; ij<3; ij++){
      ThrustAxis[ij]=Thrust_central[ij];
      //      cout <<"ij "<< ThrustAxis[ij]<<" "<<Thrust_central[ij]<<endl;
    }
    
    //the variable which gets resummed is not thrust 
    //but tau = 1 - thrust - see calculation
    centrthr=Thrust_central[3];
    EventShapes[3] = centrthr;

    //the central transverse thrust with exponentially
    //suppressed forward term
    //centrthr is again centau = 1 - centrthr
    EventShapes[4]=centrthr+ExpTerm;
    
    //the central transverse thrust with 
    //recoil term
    //centrthr is again centau = 1 - centrthr
    EventShapes[5]=centrthr+RecTerm;
    
    //the central thrust minor centhrmin
    double centhrmin =0;
    //rotate the coordinate system around the beam axis that
    //the thrust axis is the new y'-Axis - the projections are
    //simply the new y-values then
    double alpha_c=atan2(ThrustAxis[1],ThrustAxis[0]);
    for(unsigned int ij=0; ij<nin; ij++){
      centhrmin+= fabs(-sin(alpha_c)*Object_Px_in[ij]+cos(alpha_c)*Object_Py_in[ij]);
    }
    centhrmin=centhrmin/Pt_sum_c;
    EventShapes[6] = centhrmin;

    //the central thrust minor with exponentially
    //suppressed forward term centhrminwexp
    
    EventShapes[7]=centhrmin+ExpTerm;
    
    //the central thrust minor with 
    //recoil term centhrminwr
    EventShapes[8]=centhrmin+RecTerm;
  
    //central jet masses
    //define two jet masses in region U and D
    double cenjm_up=0;
    double cenjm_down=0;
    double dot_product =0;
    
    std::vector<double> up_sum;
    std::vector<double> down_sum;
    for(unsigned int ij=0; ij<4; ij++){
      up_sum.push_back(0.);
      down_sum.push_back(0.);
    }
    for(unsigned int ij=0; ij<nin; ij++){
      dot_product = Object_Px_in[ij]*ThrustAxis[0]+Object_Py_in[ij]*ThrustAxis[1];
      if(dot_product>=0){
	up_sum[0]+=Object_Px_in[ij];
	up_sum[1]+=Object_Py_in[ij];
	up_sum[2]+=Object_Pz_in[ij];
	up_sum[3]+=Object_E_in[ij];
      }else{
	down_sum[0]+=Object_Px_in[ij];
	down_sum[1]+=Object_Py_in[ij];
	down_sum[2]+=Object_Pz_in[ij];
	down_sum[3]+=Object_E_in[ij];
      }
    }
    cenjm_up = lorentz_sp(up_sum,up_sum)/pow(P_sum_c,2.); //GMA pow(Pt_sum_c,2);
    cenjm_down = lorentz_sp(down_sum,down_sum)/pow(P_sum_c,2.); //GMA pow(Pt_sum_c,2);

    
    //central total jet mass centotjm
    double centotjm=0;
    centotjm = cenjm_up+cenjm_down;
    
    EventShapes[9]=centotjm;
    
    //central total jet mass with exponentially
    //suppressed forward term centotjmwexp
    double centotjmwexp=0;
    centotjmwexp = centotjm+ExpTerm;
    
    EventShapes[10]=centotjmwexp;
    
    //central total jet mass with recoil term
    //centotjmwr
    double centotjmwr=0;
    centotjmwr = centotjm+RecTerm;
    
    EventShapes[11]=centotjmwr;
    
    //central heavy jet mass cenheavjm
    double cenheavjm=0;
    cenheavjm = max(cenjm_up,cenjm_down);
    
    EventShapes[12]=cenheavjm;
    
    //central heavy jet mass with exponentially
    //suppressed forward term cenheavjmwexp
    double cenheavjmwexp=0;
    cenheavjmwexp = cenheavjm+ExpTerm;
    
    EventShapes[13]=cenheavjmwexp;
    
    //central heavy jet mass with recoil term
    //cenheavjmwr
 
    double cenheavjmwr=0;
    cenheavjmwr = cenheavjm+RecTerm;
    
    EventShapes[14]=cenheavjmwr;
    
    double centrjm_up=0;
    double centrjm_down=0;	
    std::vector<double> upsum;
    std::vector<double> downsum;
    for(unsigned int ij=0; ij<3;ij++){
      upsum.push_back(0.);
      downsum.push_back(0.);
    }
    for(unsigned int ij=0; ij<nin; ij++){
      dot_product = Object_Px_in[ij]*ThrustAxis[0]+Object_Py_in[ij]*ThrustAxis[1];
      if(dot_product>=0){
	upsum[0]+=Object_Px_in[ij];
	upsum[1]+=Object_Py_in[ij];
	upsum[2]+=Object_Et_in[ij];
      }else{
	downsum[0]+=Object_Px_in[ij];
	downsum[1]+=Object_Py_in[ij];
	downsum[2]+=Object_Et_in[ij];
      }
    }
    centrjm_up= lorentz_sp(upsum,upsum)/pow(Pt_sum_c,2);
    centrjm_down= lorentz_sp(downsum,downsum)/pow(Pt_sum_c,2);	
    double centottrjm = centrjm_up + centrjm_down;
    double cenheavtrjm = max(centrjm_up,centrjm_down);

    EventShapes[24] = centottrjm;
    EventShapes[25] = centottrjm + ExpTerm;
    EventShapes[26] = centottrjm + RecTerm;

    EventShapes[27] = cenheavtrjm;
    EventShapes[28] = cenheavtrjm + ExpTerm;
    EventShapes[29] = cenheavtrjm + RecTerm;
    

    //central three-jet resolution threshold
    double ceny3=0;

    //central three-jet resolution threshold with
    //exponentially suppressed forward term
    double ceny3wexp=0;

    //central three-jet resolution threshold
    //with recoil term
    double ceny3wr=0;

    //central three-jet resolution threshold for antikt
    double aceny3=0,  aceny3wexp=0,  aceny3wr = 0;

    if(nin<3){
      //cout<<"WARNING: less than three central input momenta!! The central three-jet resolution threshold variables will be set to -1.0 per default!"<<endl; 
      ceny3 =-1.0;
      ceny3wexp =-1.0;
      ceny3wr =-1.0;

      aceny3= aceny3wexp = aceny3wr = -1;

    } else {
      ceny3 = three_jet_res(Object_v4_in, rap, 1);
      //this time we need to add the squares of the recoil- and exp-terms
      ceny3wexp = ceny3+pow(ExpTerm,2);
      ceny3wr = ceny3+pow(RecTerm,2);

      aceny3 = three_jet_res(Object_v4_in, rap, -1);
      //this time we need to add the squares of the recoil- and exp-terms
      aceny3wexp = aceny3+pow(ExpTerm,2);
      aceny3wr = aceny3+pow(RecTerm,2);

    }

   //cout<<"ceny3==="<<ceny3<<"  ceny3wexp==="<<ceny3wexp<<"  ceny3wr=="<<ceny3wr<<endl;
    //cout<<"aceny3==="<<aceny3<<"  aceny3wexp==="<<aceny3wexp<<"  aceny3wr=="<<aceny3wr<<endl;

    EventShapes[15] = ceny3;
    EventShapes[16] = ceny3wexp;
    EventShapes[17] = ceny3wr;

    EventShapes[0] = aceny3;
    EventShapes[1] = aceny3wexp;
    EventShapes[2] = aceny3wr;

    //the central jet broadenings in the up and down
    //region
    double cenbroad_up=0;
    double cenbroad_down=0;

    double eta_up=0;
    unsigned int num_up=0;
    double eta_down =0;
    unsigned int num_down =0;
    double phi_temp =0;
    double phi_up_aver =0;
    double phi_down_aver =0;
    double Pt_sum_up =0;
    double Pt_sum_down =0;
    double dot_product_b =0;
    vector<double> phi_up;
    vector<double> phi_down;
    double py_rot =0;
    double px_rot =0;

    for(unsigned int ij=0; ij<4; ij++){
      up_sum.push_back(0.);
      down_sum.push_back(0.);
    }
    for(unsigned int ij=0; ij<nin; ij++){
      dot_product_b =sqrt(Object_Px_in[ij]*ThrustAxis[0]+Object_Py_in[ij]*ThrustAxis[1]);
      if(dot_product_b>=0){
	Pt_sum_up+=Object_Pt_in[ij];
	//rotate the coordinate system so that
	//the central thrust axis is e_x
        px_rot = cos(alpha_c)*Object_Px_in[ij]+sin(alpha_c)*Object_Py_in[ij];
        py_rot = - sin(alpha_c)*Object_Px_in[ij]+cos(alpha_c)*Object_Py_in[ij];
        //calculate the eta and phi in the rotated system
        eta_up+=Object_Pt_in[ij]*Object_Eta_in[ij];
        phi_temp =atan2(py_rot,px_rot);
	
	if(phi_temp>M_PI/2){
	  phi_temp = phi_temp - M_PI/2;}
	if (phi_temp<-M_PI/2){
	  phi_temp = phi_temp + M_PI/2;
	  }
	phi_up.push_back(phi_temp);
	phi_up_aver+=Object_Pt_in[ij]*phi_temp;
	num_up+=1;
      }else{
	eta_down+=Object_Pt_in[ij]*Object_Eta_in[ij];
	Pt_sum_down+=Object_Pt_in[ij];
	px_rot = cos(alpha_c)*Object_Px_in[ij]+sin(alpha_c)*Object_Py_in[ij];
	py_rot = - sin(alpha_c)*Object_Px_in[ij]+cos(alpha_c)*Object_Py_in[ij];
	phi_temp =atan2(py_rot,px_rot);
	if(phi_temp>M_PI/2){
          //if phi is bigger than pi/2 in the new system calculate 
          //the difference to the thrust axis 
	  phi_temp = M_PI -phi_temp;}
	if (phi_temp<-M_PI/2){
          //if phi is smaller than 
	  phi_temp = -M_PI-phi_temp;}
	phi_down.push_back(phi_temp);
        //calculate the pt-weighted phi
	phi_down_aver+=Object_Pt_in[ij]*phi_temp;	
	num_down+=1;
      }
    }

    if (num_up!=0){
      eta_up=eta_up/Pt_sum_up;
      phi_up_aver=phi_up_aver/Pt_sum_up;}
    if(num_down!=0){
      eta_down = eta_down/Pt_sum_down;
      phi_down_aver=phi_down_aver/Pt_sum_down;}

    unsigned int index_up=0;
    unsigned int index_down=0;

    
    for(unsigned int ij=0; ij<nin; ij++){
      dot_product_b =Object_Px_in[ij]*ThrustAxis[0]+Object_Py_in[ij]*ThrustAxis[1];
      if(dot_product_b>=0){
        //calculate the broadenings of the regions with the rotated system
        //and the pt-weighted average of phi in the rotated system
	cenbroad_up+=Object_Pt_in[ij]*sqrt(pow(Object_Eta_in[ij]-eta_up,2)+pow(DeltaPhi(phi_up[index_up],phi_up_aver),2));

	
	//	if (nin>=2) cout <<"up   "<< i<<" "<<cenbroad_up<<" "<<Object_Pt_in[ij]<<" "<<Object_Eta_in[ij]<<" "<<eta_up<<" "<<phi_up[index_up]<<" "<<phi_up_aver<<" "<<2*Pt_sum_c<<" "<<(cenbroad_up+cenbroad_down)/(2*Pt_sum_c)<<endl;

	index_up+=1;
      }else{
	cenbroad_down+=Object_Pt_in[ij]*sqrt(pow(Object_Eta_in[ij]-eta_down,2)+pow(DeltaPhi(phi_down[index_down],phi_down_aver),2));
	
	//	if (nin>=2) cout <<"down "<<i<<" "<<cenbroad_down<<" "<<Object_Pt_in[ij]<<" "<<Object_Eta_in[ij]<<" "<<eta_down<<" "<<phi_down[index_down]<<" "<<phi_down_aver<<" "<<2*Pt_sum_c<<" "<<(cenbroad_up+cenbroad_down)/(2*Pt_sum_c)<<endl;

	index_down+=1;
      }
    }
    if (index_up == 0 || index_down ==0) EventShapes[nevtvar-1] *=-1.;
    
    cenbroad_up=cenbroad_up/(2*Pt_sum_c);
    cenbroad_down=cenbroad_down/(2*Pt_sum_c);

    //    if (cenbroad_up+cenbroad_down>.0001 && nin==2) cout <<"nin================== "<< nin<<" "<< cenbroad_up<<" "<<cenbroad_down<<" "<<EventShapes[30]<<endl;

    
    //central total jet broadening
    double centotbroad=0;
    centotbroad = cenbroad_up+cenbroad_down;

    EventShapes[18]=centotbroad;

    //central total jet broadening with exponentially
    //suppressed forward term
    double centotbroadwexp=0;
    centotbroadwexp =centotbroad + ExpTerm;

    EventShapes[19]=centotbroadwexp;

    //central total jet broadening with
    //recoil term
    double centotbroadwr=0;
    centotbroadwr = centotbroad + RecTerm;

    EventShapes[20]=centotbroadwr;

    //central wide jet broadening
    double cenwidbroad=0;
    cenwidbroad = max(cenbroad_up,cenbroad_down);

    EventShapes[21]=cenwidbroad;

    //central wide jet broadening with exponentially
    //suppressed forward term
    double cenwidbroadwexp=0;
    cenwidbroadwexp =cenwidbroad + ExpTerm;

    EventShapes[22]=cenwidbroadwexp;

    //central wide jet broadening with
    //recoil term
    double cenwidbroadwr=0;
    cenwidbroadwr = cenwidbroad + RecTerm;

    EventShapes[23]=cenwidbroadwr;

    //GMA added transvers sphericity and C-paramters on 2nd August 2010
    
    double a1=0, a2=0, b1=0, b2=0, c1=0, c2=0, d1=0, d2=0;
    //    if (Object_Px_in.size() == Object_Py_in.size() && Object_Px_in.size() == Object_Pz_in.size()) {
    for (int unsigned ij=0; ij<Object_Px_in.size(); ij++) {
      a1 += Object_Px_in[ij]*Object_Px_in[ij];
      b1 += Object_Px_in[ij]*Object_Py_in[ij];
      d1 += Object_Py_in[ij]*Object_Py_in[ij];
      
      a2 += Object_Px_in[ij]*Object_Px_in[ij]/Object_Pt_in[ij];
      b2 += Object_Px_in[ij]*Object_Py_in[ij]/Object_Pt_in[ij];
      d2 += Object_Py_in[ij]*Object_Py_in[ij]/Object_Pt_in[ij];
    }
    c1 = b1; c2 = b2;
    
    double b2m4ac = sqrt((a1 -d1)*(a1-d1) + 4*b1*c1);
    double val1 = ((a1 + d1) + b2m4ac)/2.;
    double val2 = ((a1 + d1) - b2m4ac)/2.;
    
    
    EventShapes[30] = 2*min(val1,val2)/(max(1.e-12,val1 + val2));
    
    //    cout <<"val[0]x "<< Object_Px_in.size()<<" "<<EventShapes[30]<<" "<<val1/(val1 + val2)<<" "<< val2/(val1 + val2)<<endl; 
    
    b2m4ac = sqrt((a2 -d2)*(a2-d2) + 4*b2*c2);
    val1 = ((a2 + d2) + b2m4ac)/2.;
    val2 = ((a2 + d2) - b2m4ac)/2.;
    
    double sum = max(1.e-12,val1+val2);
    val1 /=sum;
    val2 /=sum;
    EventShapes[31] = 4*val1*val2;
    //    cout <<"val[1]x "<< EventShapes[31]<<" "<<val1<<" "<< val2<<" "<<sum<<endl; 
    //  }

    //    for (int ij=3; ij<32; ij++) {
    for (int ij=0; ij<30; ij++) {

      if (EventShapes[ij] < 1.e-20) EventShapes[ij] = 1.e-20;
      EventShapes[ij] = log(EventShapes[ij]);
      //	cout <<"ij "<< ij <<" "<<EventShapes[ij]<<endl;
    }
  }
  
  //  if (EventShapes[30]<-40) {
  //    cout <<"==========================================="<<endl;
  //    for (int ij=0; ij<Object_Px.size(); ij++) {
  //      cout << Object_Px[ij]<<" "<<Object_Py[ij]<<" "<<Object_Pz[ij]<<" "<<Object_P[ij]<<endl;
  //    }
  //  }
  
  return 1;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||
//here functions begin||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||

double EventShape_vector::DeltaPhi(double phi1, double phi2)
{
  double delta_phi = fabs(phi2 - phi1);
  if (delta_phi > M_PI){ 
    delta_phi = 2*M_PI-delta_phi;
  } 
  return delta_phi;
} 

double EventShape_vector::max (double a, double b)
{ if(a>=b){return a;}
 else{return b;}
}

double EventShape_vector::min (double a, double b)
{ if(a<=b){return a;}
 else{return b;}
}

//returns the scalar product between two 4 momenta
double EventShape_vector::lorentz_sp(std::vector<double> a, std::vector<double> b){
  unsigned int dim = (unsigned int) a.size();
  if(a.size()!=b.size()){
    cout<<"ERROR!!! Dimension of input vectors are different! Change that please!"<<endl;
    return 0;}
  else{double l_dot_product=a[dim-1]*b[dim-1];
  for(unsigned int ij=0; ij<dim-1; ij++){
    l_dot_product-=a[ij]*b[ij];
  }
  return l_dot_product;
  }
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//function which calculates the three-jet resolution thresholds|||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
double EventShape_vector::three_jet_res(vector<HepLorentzVector> In_Object_v4, int rap, int m_algo){

  unsigned int y3_length = (unsigned int)In_Object_v4.size();
  
  std::vector<double> In_Object_P;
  std::vector<double> In_Object_Pt;
  std::vector<double> In_Object_Eta;
  std::vector<double> In_Object_Phi;

 if (!In_Object_P.empty()){
    In_Object_P.clear();
    In_Object_Pt.clear();
    In_Object_Eta.clear();
    In_Object_Phi.clear();
  }

  for(unsigned int ij = 0; ij < y3_length; ij++){
    In_Object_P.push_back(0.);
    In_Object_Pt.push_back(0.);
    In_Object_Eta.push_back(0.);
    In_Object_Phi.push_back(0.);
  }
  double theta_y3_1st=0;
  for(unsigned int jk =0; jk<y3_length; jk++){
    In_Object_P[jk] = In_Object_v4[jk].rho();
    In_Object_Pt[jk] =  In_Object_v4[jk].perp();
    //calculates the pseudorapidity
    //to prevent a division by zero
    if(rap==0){
      if (fabs(In_Object_v4[jk].pz()) > 1E-5) {
	theta_y3_1st = atan(In_Object_Pt[jk]/(In_Object_v4[jk].pz()));
      }
      else {theta_y3_1st = M_PI/2;}
      if(theta_y3_1st<0.){theta_y3_1st = theta_y3_1st + M_PI;}
      In_Object_Eta[jk] = - log(tan(0.5*theta_y3_1st));
    }
    //calculates the real rapidity
    if(rap==1){
      In_Object_Eta[jk]=0.5*log((In_Object_v4[jk].e()+In_Object_v4[jk].pz())/(In_Object_v4[jk].e()-In_Object_v4[jk].pz()));
    }
    In_Object_Phi[jk] = In_Object_v4[jk].phi();
  }

  //the three-jet resolution 
  //threshold y3
  double y3 =0;

  //vector which will be filled with the 
  //minimum of the distances
  double max_dmin_temp=0;

  double max_dmin =0;

  //distance input object jk, beam
  double distance_jB=0;
  double distance_jB_min=0;
  //distance of input object jk to l
  double distance_jk=0;
  double distance_jk_min =0;
  //as we search the minimum of the distances
  //give them values which are for sure higher
  //than those we evaluate first in the for-loups


  unsigned int index_jB=0;

  unsigned int index_j_jk = 0;
  unsigned int index_k_jk =0;

  //to decide later if the minmum is a jB or jk
  int decide_jB=-1;

  std::vector<double> input_Pt;
  std::vector<double> input_Px;
  std::vector<double> input_Py;
  std::vector<double> input_Pz;
  std::vector<double> input_P;
  std::vector<double> input_E;
  std::vector<double> input_Phi;
  std::vector<double> input_Eta;

   if (!input_Pt.empty()){
     input_Pt.clear();
     input_Px.clear();
     input_Px.clear();
     input_Pz.clear();
     input_P.clear();
     input_E.clear();
     input_Phi.clear();
     input_Eta.clear();
   }

  for(unsigned int j=0; j<y3_length;j++){
    input_Pt.push_back(In_Object_Pt[j]);
    input_Px.push_back(In_Object_v4[j].px());
    input_Py.push_back(In_Object_v4[j].py());
    input_Pz.push_back(In_Object_v4[j].pz());
    input_P.push_back(In_Object_P[j]);
    input_E.push_back(In_Object_v4[j].e());
    input_Phi.push_back(In_Object_Phi[j]);
    input_Eta.push_back(In_Object_Eta[j]);
  }
  if(y3_length<3){
    //cout<<"WARNING: less than three input momenta!! The directly global three-jet resolution threshold will be set to -1.0 per default!"<<endl; 
    y3 =-1.0;}
  else{
    unsigned int rest = y3_length;
    for(unsigned int i = 0; i<y3_length; i++){
      //make the minima at the initialization step
      //of each louping bigger than the first values
      distance_jB_min = pow(input_Pt[0],2*m_algo) + 10000;
      //DELTA PHIs wanted not the pure difference
      distance_jk_min=min(pow(input_Pt[1],2*m_algo),pow(input_Pt[0],2*m_algo))*(pow(input_Eta[1]-input_Eta[0],2)+pow(DeltaPhi(input_Phi[1],input_Phi[0]),2)) + 10;
      //do the procedure only until we have only 2 objects left anymore
      if(rest>2){
	for(unsigned int j=0; j<rest;j++){
	  //calculate the distance between object j and the beam
	  distance_jB = pow(input_Pt[j],2*m_algo);
	  if(distance_jB<distance_jB_min){
	    distance_jB_min = distance_jB;
	    index_jB = j;}	
	  if(j>0){
	    for(unsigned int k=0; k<j;k++){
	      //calculate the distance in delta eta and delta phi
	      //between object i and object j
	      distance_jk = min(pow(input_Pt[j],2*m_algo),pow(input_Pt[k],2*m_algo))*(pow(input_Eta[j]-input_Eta[k],2)+pow(DeltaPhi(input_Phi[j],input_Phi[k]),2))/(m_radSq);
	      if(distance_jk<distance_jk_min){
		distance_jk_min = distance_jk;
		index_j_jk = j;	      
		index_k_jk =k;
	      }
	    }
	    //cout<<"distance jk min "<<distance_jk_min<<endl;
	  }
	}
	//decide if the minimum is from a jB 
	//or jk combination
	if(distance_jk_min<distance_jB_min){
	  max_dmin_temp = max(distance_jk_min,max_dmin_temp);
	  decide_jB = 0;
	}else{
	  max_dmin_temp = max(distance_jB_min,max_dmin_temp);
	  decide_jB=1;}	
	//if we have only three jets left calculate
	//the maxima of the dmin's
	//if the minimum is a jB eliminate the input object
	if(decide_jB==1){
	  //if index_jB is the last one nothing is to do
	  if(index_jB!=rest-1){
	    for(unsigned int ij=index_jB; ij<rest-1;ij++){
	      input_Pt[ij]=input_Pt[i+1];
	      input_Phi[ij]=input_Phi[i+1];
	      input_Eta[ij]=input_Eta[i+1];
	      input_Px[ij]=input_Px[i+1];
	      input_Py[ij]=input_Py[i+1];
	      input_Pz[ij]=input_Pz[i+1];
	      input_E[ij]=input_E[i+1];
	    }
	  }
	}
	//if the minimum is a jk combine both input objects
	if(decide_jB==0){
	  input_Px[index_k_jk] = input_Px[index_k_jk]+input_Px[index_j_jk];
	  input_Py[index_k_jk] = input_Py[index_k_jk]+input_Py[index_j_jk];
	  input_Pz[index_k_jk] = input_Pz[index_k_jk]+input_Pz[index_j_jk];
	  input_E[index_k_jk] = input_E[index_k_jk]+input_E[index_j_jk];
	  input_P[index_k_jk]=sqrt(pow(input_Px[index_k_jk],2)+pow(input_Py[index_k_jk],2)+pow(input_Pz[index_k_jk],2));
	  //calculate the pt, eta and phi of the new combined momenta k_jk
	  input_Pt[index_k_jk] = sqrt(pow(input_Px[index_k_jk],2)+pow(input_Py[index_k_jk],2));
	  //in the case of pseudorapidity
	  if(rap==0){
	    double theta_new =0;
	    if (fabs(input_Pz[index_k_jk]) > 1E-5){
	      theta_new = atan(input_Pt[index_k_jk]/(input_Pz[index_k_jk]));
	    }
	    else {theta_new = M_PI/2;}
	    if(theta_new<0.){theta_new = theta_new + M_PI;}
	    input_Eta[index_k_jk] = - log(tan(0.5*theta_new));
	  }
	  //in the real rapidity y is wanted 
	  if(rap==1){
	    input_Eta[index_k_jk]=0.5*log((input_E[index_k_jk]+input_Pz[index_k_jk])/(input_E[index_k_jk]-input_Pz[index_k_jk]));
	  }
	  input_Phi[index_k_jk] = atan2( input_Py[index_k_jk],input_Px[index_k_jk]);
	  if(index_j_jk!=rest-1){
	    for(unsigned int ij=index_j_jk; ij<rest-1; ij++){
	      input_Pt[ij]=input_Pt[ij+1];
	      input_Phi[ij]=input_Phi[ij+1];
	      input_Eta[ij]=input_Eta[ij+1];
	      input_Px[ij]=input_Px[ij+1];
	      input_Py[ij]=input_Py[ij+1];
	      input_Pz[ij]=input_Pz[ij+1];
	      input_E[ij]=input_E[ij+1];
	    }
	  }
	}
      }
      if(rest==3){
	max_dmin=max_dmin_temp;
      }
      //cout<<"rest "<<rest<<endl;
      rest=rest-1;
    }
  }
  
  double Et2=0;
  Et2= input_Pt[0]+input_Pt[1];

/*
 if (m_algo==1) { 
    y3=max_dmin/pow(Et2,2);
  cout<< "Kt- ET= "<<Et2<<" dmin= "<<max_dmin<<endl;
  } else if (max_dmin>0) {
    y3=pow(Et2,-2)/max_dmin;
 cout<< "AKt- ET=  "<<m_algo<<"  " <<Et2<<" dmin= "<<max_dmin<<endl;
  }
  return y3;
 // */
  y3=max_dmin/pow(Et2,2);

  return y3;
}

//|||||||||||||||||||||||||||||||||||||||||
//function which calculates the thrusts||||
//|||||||||||||||||||||||||||||||||||||||||

std::vector<double> EventShape_vector::Thrust_calculate (std::vector<HepLorentzVector> Input_4v){
  std::vector<double> Input_Px;
  std::vector<double> Input_Py;
  double thrustmax_calc =0;
  double temp_calc =0;
  unsigned int length_thrust_calc =0;
  std::vector<double> ThrustValues;
  std::vector<double> Thrust_Axis_calc;
  std::vector<double> p_thrust_max_calc;
  std::vector<double> p_dec_1_calc;
  std::vector<double> p_dec_2_calc;
  std::vector<double> p_pt_beam_calc;

  if (!ThrustValues.empty()){
    ThrustValues.clear();
    Thrust_Axis_calc.clear();
    p_thrust_max_calc.clear();
    p_dec_1_calc.clear();
    p_dec_2_calc.clear();
    p_pt_beam_calc.clear();
  }
  
  for(unsigned int ij = 0; ij < 3; ij++){
    p_pt_beam_calc.push_back(0.);
    p_dec_1_calc.push_back(0.);
    p_dec_2_calc.push_back(0.);
    p_thrust_max_calc.push_back(0.);
    Thrust_Axis_calc.push_back(0.);
  }
  
  for(unsigned int ij =0;ij<4;ij++){
    ThrustValues.push_back(0.);
  }

  //  cout <<"lenght "<< Input_4v.size()<<endl;
  length_thrust_calc = Input_4v.size();
  for (unsigned jk=0; jk<length_thrust_calc; jk++) {
    Input_Px.push_back(Input_4v[jk].px());
    Input_Py.push_back(Input_4v[jk].py());
  }

  float Pt_sum_calc =0;
  
  Hep2Vector inivec(-100, -100);
  if (mncosthe>0) {
    inivec.setX(iniDir[0]);
    inivec.setY(iniDir[1]);
  }
    //   Hep2Vector inivec(iniDir[0], iniDir[1]);

  for(unsigned int jk=0;jk<length_thrust_calc;jk++){
    bool isColl=false;  //Don't consider this particle, if it is far away from thrust direction
    if (mncosthe>0) {
      Hep2Vector parvec(Input_Px[jk], Input_Py[jk]);
      double costhe =cos(inivec.angle(parvec));
      //      cout<<"costheta "<<parvec<<" "<<inivec<<" "<< costhe<<" "<<mncosthe<<endl;
      if (abs(costhe)< mncosthe) {
 	isColl =true;
      }
    }

    if (isColl) continue;
    
    Pt_sum_calc+=sqrt(pow(Input_Px[jk],2)+pow(Input_Py[jk],2)); 
    for(unsigned int ij = 0; ij < 3; ij++){
      p_thrust_max_calc[ij]=0;
    }
    //get a vector perpendicular to the beam axis and 
    //perpendicular to the momentum of particle k
    //per default beam axis b = (0,0,1)   
    p_pt_beam_calc[0] = Input_Py[jk]*1; 
    p_pt_beam_calc[1] = - Input_Px[jk]*1;
    p_pt_beam_calc[2] = 0.; // GMA p_pt_beam_calc[3] = 0.;
    for(unsigned int ij=0; ij<length_thrust_calc; ij++){
      if(ij!=jk){
	if((Input_Px[ij]*p_pt_beam_calc[0]+Input_Py[ij]*p_pt_beam_calc[1])>=0){
	  p_thrust_max_calc[0]= p_thrust_max_calc[0]+Input_Px[ij];
	  p_thrust_max_calc[1]= p_thrust_max_calc[1]+Input_Py[ij];
	}
	else{
	  p_thrust_max_calc[0]= p_thrust_max_calc[0]-Input_Px[ij];
	  p_thrust_max_calc[1]= p_thrust_max_calc[1]-Input_Py[ij];
	}
      }
    }
    p_dec_1_calc[0]=p_thrust_max_calc[0]+Input_Px[jk];
    p_dec_1_calc[1]=p_thrust_max_calc[1]+Input_Py[jk];
    p_dec_1_calc[2]=0;
    p_dec_2_calc[0]=p_thrust_max_calc[0]-Input_Px[jk];
    p_dec_2_calc[1]=p_thrust_max_calc[1]-Input_Py[jk];
    p_dec_2_calc[2]=0;
    temp_calc = pow(p_dec_1_calc[0],2)+pow(p_dec_1_calc[1],2);

    if(temp_calc>thrustmax_calc){
      thrustmax_calc =temp_calc;
      for(unsigned int ij=0; ij<3; ij++){
	Thrust_Axis_calc[ij]=p_dec_1_calc[ij]/sqrt(thrustmax_calc);
      }
    }
    temp_calc = pow(p_dec_2_calc[0],2)+pow(p_dec_2_calc[1],2);
    if(temp_calc>thrustmax_calc){
      thrustmax_calc =temp_calc;
      for(unsigned int ij=0; ij<3; ij++){
	Thrust_Axis_calc[ij]=p_dec_2_calc[ij]/sqrt(thrustmax_calc);
      }
    }
  }

  for (unsigned int ij=0; ij<3; ij++){
    ThrustValues[ij]=Thrust_Axis_calc[ij];
  }

  double thrust_calc=0;
  thrust_calc = sqrt(thrustmax_calc)/Pt_sum_calc;
  
  //the variable which gets resummed is not the thrust
  //but tau=1-thrust
  ThrustValues[3]=1.-thrust_calc; 
  
  if (ThrustValues[3] < 1.e-20) ThrustValues[3] = 1.e-20;

  return ThrustValues;
}

void EventShape_vector::getSphericity(double* val) {
  
  int nin=0;
  double a1=0, a2=0, b1=0, b2=0, c1=0, c2=0, d1=0, d2=0;

  std::vector<double> Object_Px;
  std::vector<double> Object_Py;
  std::vector<double> Object_Pz;
  //  cout <<"lenght "<<  Object_4v.size()<<endl;
  for (unsigned jk=0; jk<Object_4v.size(); jk++) {
    if (abs(Object_4v[jk].eta())>eta_c) continue;
    Object_Px.push_back( Object_4v[jk].px());
    Object_Py.push_back( Object_4v[jk].py());
    Object_Pz.push_back( Object_4v[jk].pz());
  }
 

  if (Object_Px.size() == Object_Py.size()) {
    for (unsigned int ij=0; ij<Object_Px.size(); ij++) {
      double pt = sqrt(Object_Px[ij]*Object_Px[ij] + Object_Py[ij]*Object_Py[ij]);
      double mag= sqrt(Object_Pz[ij]*Object_Pz[ij] + pt*pt);
      double eta = -log(tan(acos(Object_Pz[ij]/mag)/2.));
      if (abs(eta) < eta_c) {
	nin++;
	a1 += Object_Px[ij]*Object_Px[ij];
	b1 += Object_Px[ij]*Object_Py[ij];
	d1 += Object_Py[ij]*Object_Py[ij];
	
	a2 += Object_Px[ij]*Object_Px[ij]/pt;
	b2 += Object_Px[ij]*Object_Py[ij]/pt;
	d2 += Object_Py[ij]*Object_Py[ij]/pt;     
      }
    }
    c1 = b1; c2 = b2;
    
    double b2m4ac = sqrt((a1 -d1)*(a1-d1) + 4*b1*c1);
    double val1 = ((a1 + d1) + b2m4ac)/2.;
    double val2 = ((a1 + d1) - b2m4ac)/2.;
    
    
    val[0] = log(max(1.e-20, 2*min(val1,val2)/(max(1.e-12,val1 + val2))));
    

    b2m4ac = sqrt((a2 -d2)*(a2-d2) + 4*b2*c2);
    val1 = ((a2 + d2) + b2m4ac)/2.;
    val2 = ((a2 + d2) - b2m4ac)/2.;
    
    double sum = max(1.e-12,val1+val2);
    val1 /=sum;
    val2 /=sum;
    val[1] = log (max(1.e-20, 4*val1*val2));

  }
}

  /*
double sphericity[2];
getSpehticity(spheticity);



*/
