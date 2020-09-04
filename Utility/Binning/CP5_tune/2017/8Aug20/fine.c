void fine(){

//For Jets
static const int nHLTmx = 8;
static const int nvar = 32;
//static const int nmxbins = 18;

static const int nmxbins=26;
int nbinsx0[nvar]={0,0,0,26,0,0,
                   0,0,0,20,0,0,
                   0,0,0,12,0,0,
                   20,0,0,0,0,0,
                   22,0,0,0,0,0,0,0};


//For charge Particles
int nbinsx1[nvar]={0,0,0,20,0,0,
                   0,0,0,16,0,0,
                   0,0,0,8,0,0,
                   22,0,0,0,0,0,
                   14,0,0,0,0,0,0,0};
//

static const int rnmxbins=26;
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


//For Gen Level
//for Jets
double binrngs0[nvar][nmxbins+1] ={{},{},{},
                                  {-6.38, -6.08, -5.79, -5.52, -5.26, -5.01, -4.77, -4.53, -4.3, -4.07, -3.85, -3.64, -3.43, -3.23, -3.04, -2.85, -2.67, -2.5, -2.34, -2.19, -2.05, -1.92, -1.8, -1.68, -1.57, -1.47, -1.38},//, -1.3, -1.23, -1.17, -1.11 },//30
                                  {},{},{},{},{},
                                  {-4.69, -4.45, -4.22, -4, -3.78, -3.56, -3.35, -3.14, -2.93, -2.72, -2.52, -2.32, -2.13, -1.94, -1.76, -1.59, -1.42, -1.26, -1.11, -0.97, -0.85},// -0.74, -0.64, -0.55,-0.47, -0.4, -0.33, -0.27, -0.22 },   //28 bins
                                  {},{},{},{},{},
                                  {/*-6.17,*/ -5.68, -5.24, -4.83, -4.45, -4.1, -3.77, -3.45, -3.15, -2.86, -2.58, -2.31, -2.04, -1.78},//, -1.52}, //14 bin
                                  {},{},
                                  {/*-4.69, -4.44, -4.18, -3.91, -3.65,*/ -3.39, -3.14, -2.9, -2.67, -2.45, -2.25, -2.06, -1.88, -1.72, -1.57, -1.43, -1.3, -1.18, -1.07, -0.97, -0.88, -0.79, -0.71, -0.63, -0.56, -0.49},//, -0.43, -0.37, -0.32, -0.27, -0.22, -0.17, -0.13}, //32 Bins
                                  {},{},{},{},{},
                                  {-5.58, -5.29, -5.01, -4.75, -4.5, -4.26, -4.03, -3.8, -3.58, -3.37, -3.16, -2.96, -2.77, -2.59, -2.41, -2.24, -2.08, -1.93, -1.8, -1.68, -1.57, -1.47, -1.37},//, -1.28, -1.2},//24 bins
                                  {},{},{},{},{},{},{}};


//-----------------------------------------------------------------------------------
//For charge particle
double binrngs1[nvar][nmxbins+1] ={{},{},{},
                                  {-5.27, -5.02, -4.77, -4.53, -4.29, -4.05, -3.81, -3.57, -3.33, -3.09, -2.85, -2.62, -2.39, -2.17, -1.97, -1.78, -1.61, -1.46, -1.33, -1.22, -1.13}, //20 bins
                                  {},{},{},{},{},
                                  {-4.57, -4.21, -3.86, -3.52, -3.19, -2.87, -2.56, -2.26, -1.98, -1.71, -1.46, -1.23, -1.02, -0.83, -0.66, -0.51, -0.38},// -0.27, -0.18}, // 18 bins
                                  {},{},{},{},{},
                                  {-6.71, -6, -5.29, -4.62, -3.99, -3.39, -2.8, -2.18, -1.49},// 8 bins
                                  {},{},
                                  {-4.55, -4.25, -3.99, -3.76, -3.55, -3.35, -3.16, -2.97, -2.78, -2.59, -2.4, -2.21, -2.02, -1.83, -1.64, -1.45, -1.27, -1.09, -0.92, -0.76, -0.62, -0.49, -0.38},//22
                                  {},{},{},{},{},
                                  {-5.39, -4.96, -4.57, -4.21, -3.87, -3.54, -3.22, -2.91, -2.61, -2.33, -2.06, -1.81, -1.58, -1.37, -1.18},//14
                                  {},{},{},{},{},{},{}};

double binrngs0_reco[nvar][nmxbins*2];

cout << " Jets" << endl;

for(int ijvar=0; ijvar < nvar  ; ijvar++){
	if(ijvar==3 ||ijvar==9||ijvar==15||ijvar==18||ijvar==24){
    cout << "var name " << ijvar<< endl;
//		for(int ijpt=0; ijpt < nHLTmx  ; ijpt++){
     
	cout << "{";
     
	for(int ijbin=0; ijbin <=  nbinsx0[ijvar]  ; ijbin++){


	cout << binrngs0[ijvar][ijbin] ;
	if(ijbin+1<= nbinsx0[ijvar]){cout << "," ;}
        if(ijbin+1<= nbinsx0[ijvar]) {cout << binrngs0[ijvar][ijbin]+(binrngs0[ijvar][ijbin+1]-binrngs0[ijvar][ijbin])/2 << ",";}
//	binrngs0_reco[ijvar][ijpt][ijbin*2]=binrngs0[ijvar][ijpt][ijbin];
//	binrngs0_reco[ijvar][ijpt][ijbin+1]=(binrngs0[nvar][nHLTmx][ijbin])+((binrngs0[ijvar][ijpt][ijbin+1])-(binrngs0[nvar][nHLTmx][ijbin]))/2;}
      }
    cout << "},"<< endl;
//   }
 }
}


cout << " Charge Particle" << endl ;
for(int ijvar=0; ijvar < nvar  ; ijvar++){
        if(ijvar==3 ||ijvar==9||ijvar==15||ijvar==18||ijvar==24){
    cout << "var name " << ijvar<< endl;
     //           for(int ijpt=0; ijpt < nHLTmx  ; ijpt++){

        cout << "{";

        for(int ijbin=0; ijbin <=  nbinsx1[ijvar]  ; ijbin++){


        cout << binrngs1[ijvar][ijbin];
        if(ijbin+1<= nbinsx1[ijvar]){cout << "," ;}

	if(ijbin+1<= nbinsx1[ijvar]) {cout << binrngs1[ijvar][ijbin]+(binrngs1[ijvar][ijbin+1]-binrngs1[ijvar][ijbin])/2 << ",";}

//      binrngs0_reco[ijvar][ijpt][ijbin*2]=binrngs0[ijvar][ijpt][ijbin];
//      binrngs0_reco[ijvar][ijpt][ijbin+1]=(binrngs0[nvar][nHLTmx][ijbin])+((binrngs0[ijvar][ijpt][ijbin+1])-(binrngs0[nvar][nHLTmx][ijbin]))/2;}
      }
    cout << "},"<< endl;
  /// }
 }
}

for(int ijvar=0; ijvar < nvar  ; ijvar++){
        if(ijvar==3 ||ijvar==9||ijvar==15||ijvar==18||ijvar==24){
    cout << "var name " << ijvar<< endl;
      cout << "{" ;      
//    for(int ijpt=0; ijpt < nHLTmx  ; ijpt++){
                 cout << nbinsx0[ijvar]*2 << "," ;

//		}
	 cout << "}"<< endl;
	}
}

for(int ijvar=0; ijvar < nvar  ; ijvar++){
        if(ijvar==3 ||ijvar==9||ijvar==15||ijvar==18||ijvar==24){
    cout << "var name " << ijvar<< endl;
      cout << "{" ;
  //  for(int ijpt=0; ijpt < nHLTmx  ; ijpt++){
                 cout << nbinsx1[ijvar]*2 << "," ;

  //              }
         cout << "}"<< endl;
        }
}



}//end code
