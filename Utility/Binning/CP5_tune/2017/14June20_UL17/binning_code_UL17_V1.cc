//For Jets

int nbinsx0[nvar][nHLTmx]={{},{},{},
                   {13,13,13,15,15,17,16,18},//3
                   {},{},{},{},{},
                   {10,10,11,12,12,13,13,12},//9
                   {},{},{},{},{},
                   {4,4,5,5,6,6,6,7},//15
                   {},{},
                   {8,7,7,8,9,11,11,12},//18
                   {},{},{},{},{},
                   {8,8,10,11,11,12,13,13},//24
                   {},{},{},{},{},{},{}};

//For charge Particles
int nbinsx1[nvar][nHLTmx]={{},{},{},
                   {11,12,12,12,13,12,11,9},//3
                   {},{},{},{},{},
                   {8,8,8,8,9,9,8,7},//9
                   {},{},{},{},{},
                   {4,4,4,4,4,4,4,4},//15
                   {},{},
                   {9,9,9,9,9,9,8,7},//18
                   {},{},{},{},{},
                   {6,7,8,8,8,8,8,7},//24
                   {},{},{},{},{},{},{}};
///////////////////Modify /////////////////////////////
//for Jets
double binrngs0[nvar][nHLTmx][nmxbins+1] ={
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //Thrust 3
        {{-6.71, -5.87, -5.18, -4.59, -4.07, -3.6, -3.17, -2.77, -2.4, -2.05, -1.73, -1.43, -1.16, -0.91},
{-6.71, -5.76, -5.09, -4.53, -4.04, -3.58, -3.15, -2.74, -2.35, -1.98, -1.65, -1.36, -1.12, -0.92},
{-6.71, -5.78, -5.14, -4.61, -4.14, -3.7, -3.28, -2.87, -2.48, -2.11, -1.77, -1.47, -1.22, -1.02},
{-6.71, -5.94, -5.33, -4.8, -4.32, -3.87, -3.44, -3.03, -2.64, -2.28, -1.94, -1.64, -1.38, -1.16, -0.97, -0.82},
{-6.71, -5.99, -5.4, -4.88, -4.4, -3.95, -3.52, -3.11, -2.72, -2.36, -2.03, -1.74, -1.48, -1.27, -1.1, -0.96},
{-6.71, -6.04, -5.46, -4.94, -4.46, -4, -3.57, -3.16, -2.77, -2.41, -2.09, -1.8, -1.55, -1.35, -1.18, -1.05, -0.95, -0.87},
{-6.71, -6.04, -5.45, -4.92, -4.43, -3.97, -3.54, -3.14, -2.77, -2.43, -2.12, -1.84, -1.6, -1.39, -1.21, -1.06, -0.94},
{-6.71, -6.06, -5.48, -4.95, -4.46, -4, -3.57, -3.18, -2.82, -2.49, -2.2, -1.94, -1.71, -1.51, -1.34, -1.2, -1.09, -1, -0.92}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoT 9
        {{-4.99, -4.44, -3.89, -3.35, -2.82, -2.3, -1.81, -1.35, -0.93, -0.56, -0.24},
{-4.99, -4.43, -3.88, -3.35, -2.83, -2.33, -1.85, -1.4, -0.99, -0.62, -0.29 },
{-4.99, -4.44, -3.93, -3.44, -2.96, -2.49, -2.04, -1.61, -1.22, -0.87, -0.56, -0.31},
{-4.99, -4.47, -3.98, -3.51, -3.05, -2.61, -2.18, -1.77, -1.39, -1.05, -0.75, -0.5, -0.3},
{-4.99, -4.49, -4.02, -3.57, -3.13, -2.7, -2.29, -1.89, -1.52, -1.18, -0.89, -0.64, -0.44},
{-4.99, -4.49, -4.02, -3.57, -3.14, -2.72, -2.32, -1.94, -1.59, -1.27, -0.99, -0.74, -0.54, -0.38},
{-4.99, -4.49, -4.03, -3.6, -3.18, -2.78, -2.39, -2.02, -1.67, -1.35, -1.06, -0.81, -0.61, -0.45 },
{-5.08, -4.6, -4.15, -3.72, -3.3, -2.9, -2.51, -2.14, -1.8, -1.49, -1.21, -0.97, -0.77 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
     //Y23
        {{-4.35, -3.46, -2.52, -1.66, -0.98 },
{-4.35, -3.39, -2.43, -1.51, -0.8 },
{-5.17, -4.19, -3.17, -2.25, -1.47, -0.83 },
{-5.17, -4.16, -3.22, -2.4, -1.69, -1.07 },
{-5.58, -4.67, -3.78, -2.97, -2.26, -1.62, -1.01 },
{-5.99, -5.02, -4.11, -3.31, -2.61, -1.98, -1.41 },
{-6.08, -5.02, -4.17, -3.44, -2.78, -2.19, -1.66 },
{-6.89, -5.82, -4.93, -4.17, -3.5, -2.89, -2.33, -1.79 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //BT
        {{-2.54, -2.25, -1.96, -1.67, -1.38, -1.09, -0.81, -0.55, -0.32 },
{-2.54, -2.2, -1.86, -1.53, -1.22, -0.92, -0.65, -0.4 },
{-2.54, -2.19, -1.83, -1.49, -1.18, -0.91, -0.67, -0.46 },
{-2.54, -2.18, -1.81, -1.47, -1.19, -0.95, -0.74, -0.55, -0.37 },
{-2.94, -2.56, -2.18, -1.82, -1.5, -1.22, -0.98, -0.78, -0.61, -0.46 },
{-3.35, -2.93, -2.51, -2.12, -1.77, -1.46, -1.19, -0.96, -0.77, -0.61, -0.47, -0.35 },
{-3.35, -2.9, -2.47, -2.08, -1.74, -1.45, -1.2, -0.99, -0.82, -0.67, -0.54, -0.43 },
{-3.76, -3.27, -2.8, -2.37, -2, -1.69, -1.43, -1.21, -1.02, -0.85, -0.69, -0.55, -0.41 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
       //RhoTT 24
        {{-5.08, -4.31, -3.7, -3.15, -2.62, -2.11, -1.63, -1.2, -0.85 },
{-5.08, -4.28, -3.68, -3.14, -2.63, -2.14, -1.68, -1.27, -0.93 },
{-5.08, -4.46, -3.92, -3.41, -2.91, -2.43, -1.98, -1.59, -1.27, -1.04, -0.89 },
{-5.08, -4.47, -3.94, -3.45, -2.98, -2.53, -2.11, -1.73, -1.41, -1.16, -0.97, -0.85 },
{-5.08, -4.52, -4.01, -3.54, -3.09, -2.66, -2.26, -1.9, -1.58, -1.31, -1.09, -0.92 },
{-5.08, -4.55, -4.06, -3.59, -3.14, -2.72, -2.33, -1.98, -1.67, -1.41, -1.2, -1.04, -0.92 },
{-5.49, -4.9, -4.38, -3.9, -3.45, -3.03, -2.64, -2.28, -1.95, -1.67, -1.43, -1.23, -1.07, -0.95 },
{-5.49, -4.96, -4.45, -3.97, -3.52, -3.11, -2.73, -2.39, -2.08, -1.8, -1.54, -1.3, -1.09, -0.89 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}}};


//-----------------------------------------------------------------------------------
//For charge particle
double binrngs1[nvar][nHLTmx][nmxbins+1] ={
//double binrngs1[nvar][nHLTmx][37] ={
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //Thrust 3
        {{-5.29, -4.57, -4.1, -3.7, -3.33, -2.96, -2.59, -2.21, -1.84, -1.5, -1.22, -1 },
{-5.29, -4.65, -4.2, -3.8, -3.42, -3.03, -2.63, -2.23, -1.84, -1.49, -1.21, -1.02, -0.91 },
{-5.49, -4.9, -4.44, -4.02, -3.61, -3.19, -2.76, -2.33, -1.91, -1.54, -1.25, -1.05, -0.93 },
{-5.54, -5, -4.53, -4.09, -3.65, -3.2, -2.75, -2.31, -1.9, -1.54, -1.26, -1.06, -0.93 },
{-5.65, -5.11, -4.62, -4.15, -3.68, -3.21, -2.74, -2.28, -1.86, -1.51, -1.24, -1.06, -0.94, -0.87 },
{-5.65, -5.11, -4.59, -4.08, -3.56, -3.05, -2.55, -2.09, -1.69, -1.37, -1.13, -0.97, -0.87 },
{-5.66, -5.08, -4.52, -3.96, -3.4, -2.85, -2.33, -1.86, -1.48, -1.2, -1.02, -0.9 },
{-5.54, -4.89, -4.22, -3.54, -2.87, -2.26, -1.74, -1.33, -1.05, -0.87 }},
	{{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoT 9
        {{-5.8, -4.59, -3.69, -2.94, -2.29, -1.71, -1.19, -0.73, -0.33 },
{-5.8, -4.65, -3.81, -3.1, -2.46, -1.87, -1.33, -0.86, -0.46 },
{-5.8, -4.74, -3.98, -3.32, -2.7, -2.1, -1.54, -1.05, -0.65 },
{-6.62, -5.4, -4.55, -3.83, -3.16, -2.53, -1.93, -1.38, -0.9, -0.51 },
{-6.62, -5.43, -4.59, -3.87, -3.21, -2.58, -1.98, -1.43, -0.96, -0.58 },
{-6.62, -5.47, -4.61, -3.86, -3.16, -2.49, -1.86, -1.3, -0.83, -0.47 },
{-6.62, -5.45, -4.57, -3.79, -3.05, -2.34, -1.68, -1.11, -0.66 },
{-5.8, -4.9, -3.88, -2.95, -2.2, -1.58, -1.02, -0.46 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
//        {{},{},{},{},{},{},{},{}}, //extra
        //Y23 15
        {{-6.89, -5.41, -3.93, -2.62, -1.23 },
{-6.89, -5.42, -3.94, -2.67, -1.28 },
{-6.89, -5.47, -4.09, -2.88, -1.61 },
{-6.89, -5.47, -4.12, -2.94, -1.72 },
{-6.89, -5.58, -4.28, -3.09, -1.9 },
{-6.89, -5.55, -4.24, -3.02, -1.78 },
{-6.08, -4.91, -3.57, -2.32, -0.85 }}, 
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //BT
        {{-3.35, -2.91, -2.56, -2.24, -1.93, -1.61, -1.29, -0.97, -0.66, -0.39 },
{-3.76, -3.23, -2.84, -2.5, -2.17, -1.83, -1.48, -1.13, -0.79, -0.49 },
{-3.76, -3.28, -2.9, -2.55, -2.2, -1.83, -1.45, -1.07, -0.73, -0.46 },
{-3.76, -3.31, -2.92, -2.55, -2.17, -1.79, -1.4, -1.03, -0.69, -0.41 },
{-4.17, -3.64, -3.19, -2.77, -2.36, -1.94, -1.53, -1.14, -0.79, -0.51 },
{-4.02, -3.51, -3.05, -2.6, -2.15, -1.7, -1.27, -0.89, -0.59, -0.38 },
{-3.79, -3.26, -2.74, -2.23, -1.74, -1.28, -0.88, -0.56, -0.32 },
{-3.57, -2.92, -2.31, -1.8, -1.36, -0.96, -0.55, -0.08 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoTT 24
        {{-5.49, -4.21, -3.36, -2.67, -2.07, -1.52, -1.01 },
{-5.49, -4.29, -3.52, -2.86, -2.25, -1.68, -1.19, -0.81 },
{-5.89, -4.67, -3.89, -3.22, -2.6, -2.03, -1.51, -1.08, -0.76 },
{-5.89, -4.81, -4.06, -3.4, -2.78, -2.2, -1.67, -1.23, -0.9 },
{-5.89, -4.88, -4.13, -3.46, -2.83, -2.24, -1.71, -1.27, -0.95 },
{-5.89, -4.95, -4.18, -3.48, -2.84, -2.25, -1.72, -1.27, -0.91 },
{-5.89, -4.96, -4.13, -3.38, -2.69, -2.07, -1.53, -1.07, -0.7 },
{-5.89, -4.93, -4.13, -3.43, -2.78, -2.15, -1.52, -0.86 }},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}}};


 int nbinsx[nvar]={120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120};



double endxt[nvar]={8.0, 8.0, 6.0, 8.0, 6.0,  5.0,
                    4.0, 4.0,  2.0, 7.0, 7.0, 4.0,
                    7.0, 7.0, 4.0, 8.0, 8.0, 6.0,
                    5.0, 5.0, 4.0, 5.0, 5.0, 4.0,
                    10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                    0.0, 0.0};
double startx[nvar]={2., 0.5, 1., 1, 0.0, 2.,
                     0., 0., 1., 0.2, -1., -1.,
                     0., 0., 1., 1.0, 2., 0.5,
                     0.2, -0.5, -1., 0., -1., -1.,
                     1.0, 0., -0.5, 0., 0., -0.5,
                    -1., -1.};
double endx[nvar]={10.0, 10.0, 10.0, 8.0, 12.0, 12.0,
                   19.0, 19.0, 10.0, 6.2, 10.0, 10.0,
                   13.0, 13.0,  10.0, 8.0, 10.0, 6.0,
                   6.0, 5.0, 5.0, 8.0, 8.0, 8.0,
                   7.0, 8.0, 8.0, 12.0, 12.0, 8.0,
                   0.0, 0.0};

//----------------------for reco rebin
static const int  rnmxbins =74;
int rnbinsx0[nvar][nHLTmx]={{},{},{},
                   {26,26,26,30,30,34,32,36},//3
                   {},{},{},{},{},
                   {20,20,22,24,24,26,26,24},//9
                   {},{},{},{},{},
                   {8,8,10,10,12,12,12,14},//15
                   {},{},
                   {16,14,14,16,18,22,22,24},//18
                   {},{},{},{},{},
                   {16,16,20,22,22,24,26,26},//24
                   {},{},{},{},{},{},{}};

//For charge Particles
int rnbinsx1[nvar][nHLTmx]={{},{},{},
                   {22,24,24,24,26,24,22,18},//3
                   {},{},{},{},{},
                   {16,16,16,16,18,18,16,14},//9
                   {},{},{},{},{},
                   {8,8,8,8,8,8,8,8},//15
                   {},{},
                   {18,18,18,18,18,18,16,14},//18
                   {},{},{},{},{},
                   {12,14,16,16,16,16,16,14},//24
                   {},{},{},{},{},{},{}};
///////////////////Modify /////////////////////////////
//for Jets
double rbinrngs0[nvar][nHLTmx][rnmxbins+1] ={
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //Thrust 3
        {{-6.71,-6.29,-5.87,-5.525,-5.18,-4.885,-4.59,-4.33,-4.07,-3.835,-3.6,-3.385,-3.17,-2.97,-2.77,-2.585,-2.4,-2.225,-2.05,-1.89,-1.73,-1.58,-1.43,-1.295,-1.16,-1.035,-0.91},
{-6.71,-6.235,-5.76,-5.425,-5.09,-4.81,-4.53,-4.285,-4.04,-3.81,-3.58,-3.365,-3.15,-2.945,-2.74,-2.545,-2.35,-2.165,-1.98,-1.815,-1.65,-1.505,-1.36,-1.24,-1.12,-1.02,-0.92},
{-6.71,-6.245,-5.78,-5.46,-5.14,-4.875,-4.61,-4.375,-4.14,-3.92,-3.7,-3.49,-3.28,-3.075,-2.87,-2.675,-2.48,-2.295,-2.11,-1.94,-1.77,-1.62,-1.47,-1.345,-1.22,-1.12,-1.02},
{-6.71,-6.325,-5.94,-5.635,-5.33,-5.065,-4.8,-4.56,-4.32,-4.095,-3.87,-3.655,-3.44,-3.235,-3.03,-2.835,-2.64,-2.46,-2.28,-2.11,-1.94,-1.79,-1.64,-1.51,-1.38,-1.27,-1.16,-1.065,-0.97,-0.895,-0.82},
{-6.71,-6.35,-5.99,-5.695,-5.4,-5.14,-4.88,-4.64,-4.4,-4.175,-3.95,-3.735,-3.52,-3.315,-3.11,-2.915,-2.72,-2.54,-2.36,-2.195,-2.03,-1.885,-1.74,-1.61,-1.48,-1.375,-1.27,-1.185,-1.1,-1.03,-0.96},
{-6.71,-6.375,-6.04,-5.75,-5.46,-5.2,-4.94,-4.7,-4.46,-4.23,-4,-3.785,-3.57,-3.365,-3.16,-2.965,-2.77,-2.59,-2.41,-2.25,-2.09,-1.945,-1.8,-1.675,-1.55,-1.45,-1.35,-1.265,-1.18,-1.115,-1.05,-1,-0.95,-0.91,-0.87},
{-6.71,-6.375,-6.04,-5.745,-5.45,-5.185,-4.92,-4.675,-4.43,-4.2,-3.97,-3.755,-3.54,-3.34,-3.14,-2.955,-2.77,-2.6,-2.43,-2.275,-2.12,-1.98,-1.84,-1.72,-1.6,-1.495,-1.39,-1.3,-1.21,-1.135,-1.06,-1,-0.94},
{-6.71,-6.385,-6.06,-5.77,-5.48,-5.215,-4.95,-4.705,-4.46,-4.23,-4,-3.785,-3.57,-3.375,-3.18,-3,-2.82,-2.655,-2.49,-2.345,-2.2,-2.07,-1.94,-1.825,-1.71,-1.61,-1.51,-1.425,-1.34,-1.27,-1.2,-1.145,-1.09,-1.045,-1,-0.96,-0.92}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoT 9
        {{-4.99,-4.715,-4.44,-4.165,-3.89,-3.62,-3.35,-3.085,-2.82,-2.56,-2.3,-2.055,-1.81,-1.58,-1.35,-1.14,-0.93,-0.745,-0.56,-0.4,-0.24},
{-4.99,-4.71,-4.43,-4.155,-3.88,-3.615,-3.35,-3.09,-2.83,-2.58,-2.33,-2.09,-1.85,-1.625,-1.4,-1.195,-0.99,-0.805,-0.62,-0.455,-0.29},
{-4.99,-4.715,-4.44,-4.185,-3.93,-3.685,-3.44,-3.2,-2.96,-2.725,-2.49,-2.265,-2.04,-1.825,-1.61,-1.415,-1.22,-1.045,-0.87,-0.715,-0.56,-0.435,-0.31},
{-4.99,-4.73,-4.47,-4.225,-3.98,-3.745,-3.51,-3.28,-3.05,-2.83,-2.61,-2.395,-2.18,-1.975,-1.77,-1.58,-1.39,-1.22,-1.05,-0.9,-0.75,-0.625,-0.5,-0.4,-0.3},
{-4.99,-4.74,-4.49,-4.255,-4.02,-3.795,-3.57,-3.35,-3.13,-2.915,-2.7,-2.495,-2.29,-2.09,-1.89,-1.705,-1.52,-1.35,-1.18,-1.035,-0.89,-0.765,-0.64,-0.54,-0.44},
{-4.99,-4.74,-4.49,-4.255,-4.02,-3.795,-3.57,-3.355,-3.14,-2.93,-2.72,-2.52,-2.32,-2.13,-1.94,-1.765,-1.59,-1.43,-1.27,-1.13,-0.99,-0.865,-0.74,-0.64,-0.54,-0.46,-0.38},
{-4.99,-4.74,-4.49,-4.26,-4.03,-3.815,-3.6,-3.39,-3.18,-2.98,-2.78,-2.585,-2.39,-2.205,-2.02,-1.845,-1.67,-1.51,-1.35,-1.205,-1.06,-0.935,-0.81,-0.71,-0.61,-0.53,-0.45},
{-5.08,-4.84,-4.6,-4.375,-4.15,-3.935,-3.72,-3.51,-3.3,-3.1,-2.9,-2.705,-2.51,-2.325,-2.14,-1.97,-1.8,-1.645,-1.49,-1.35,-1.21,-1.09,-0.97,-0.87,-0.77}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
     //Y23
        {{-4.35,-3.905,-3.46,-2.99,-2.52,-2.09,-1.66,-1.32,-0.98},
{-4.35,-3.87,-3.39,-2.91,-2.43,-1.97,-1.51,-1.155,-0.8},
{-5.17,-4.68,-4.19,-3.68,-3.17,-2.71,-2.25,-1.86,-1.47,-1.15,-0.83},
{-5.17,-4.665,-4.16,-3.69,-3.22,-2.81,-2.4,-2.045,-1.69,-1.38,-1.07},
{-5.58,-5.125,-4.67,-4.225,-3.78,-3.375,-2.97,-2.615,-2.26,-1.94,-1.62,-1.315,-1.01},
{-5.99,-5.505,-5.02,-4.565,-4.11,-3.71,-3.31,-2.96,-2.61,-2.295,-1.98,-1.695,-1.41},
{-6.08,-5.55,-5.02,-4.595,-4.17,-3.805,-3.44,-3.11,-2.78,-2.485,-2.19,-1.925,-1.66},
{-6.89,-6.355,-5.82,-5.375,-4.93,-4.55,-4.17,-3.835,-3.5,-3.195,-2.89,-2.61,-2.33,-2.06,-1.79}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //BT
        {{-2.54,-2.395,-2.25,-2.105,-1.96,-1.815,-1.67,-1.525,-1.38,-1.235,-1.09,-0.95,-0.81,-0.68,-0.55,-0.435,-0.32},
{-2.54,-2.37,-2.2,-2.03,-1.86,-1.695,-1.53,-1.375,-1.22,-1.07,-0.92,-0.785,-0.65,-0.525,-0.4},
{-2.54,-2.365,-2.19,-2.01,-1.83,-1.66,-1.49,-1.335,-1.18,-1.045,-0.91,-0.79,-0.67,-0.565,-0.46},
{-2.54,-2.36,-2.18,-1.995,-1.81,-1.64,-1.47,-1.33,-1.19,-1.07,-0.95,-0.845,-0.74,-0.645,-0.55,-0.46,-0.37},
{-2.94,-2.75,-2.56,-2.37,-2.18,-2,-1.82,-1.66,-1.5,-1.36,-1.22,-1.1,-0.98,-0.88,-0.78,-0.695,-0.61,-0.535,-0.46},
{-3.35,-3.14,-2.93,-2.72,-2.51,-2.315,-2.12,-1.945,-1.77,-1.615,-1.46,-1.325,-1.19,-1.075,-0.96,-0.865,-0.77,-0.69,-0.61,-0.54,-0.47,-0.41,-0.35},
{-3.35,-3.125,-2.9,-2.685,-2.47,-2.275,-2.08,-1.91,-1.74,-1.595,-1.45,-1.325,-1.2,-1.095,-0.99,-0.905,-0.82,-0.745,-0.67,-0.605,-0.54,-0.485,-0.43},
{-3.76,-3.515,-3.27,-3.035,-2.8,-2.585,-2.37,-2.185,-2,-1.845,-1.69,-1.56,-1.43,-1.32,-1.21,-1.115,-1.02,-0.935,-0.85,-0.77,-0.69,-0.62,-0.55,-0.48,-0.41}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
       //RhoTT 24
        {{-5.08,-4.695,-4.31,-4.005,-3.7,-3.425,-3.15,-2.885,-2.62,-2.365,-2.11,-1.87,-1.63,-1.415,-1.2,-1.025,-0.85},
{-5.08,-4.68,-4.28,-3.98,-3.68,-3.41,-3.14,-2.885,-2.63,-2.385,-2.14,-1.91,-1.68,-1.475,-1.27,-1.1,-0.93},
{-5.08,-4.77,-4.46,-4.19,-3.92,-3.665,-3.41,-3.16,-2.91,-2.67,-2.43,-2.205,-1.98,-1.785,-1.59,-1.43,-1.27,-1.155,-1.04,-0.965,-0.89},
{-5.08,-4.775,-4.47,-4.205,-3.94,-3.695,-3.45,-3.215,-2.98,-2.755,-2.53,-2.32,-2.11,-1.92,-1.73,-1.57,-1.41,-1.285,-1.16,-1.065,-0.97,-0.91,-0.85},
{-5.08,-4.8,-4.52,-4.265,-4.01,-3.775,-3.54,-3.315,-3.09,-2.875,-2.66,-2.46,-2.26,-2.08,-1.9,-1.74,-1.58,-1.445,-1.31,-1.2,-1.09,-1.005,-0.92},
{-5.08,-4.815,-4.55,-4.305,-4.06,-3.825,-3.59,-3.365,-3.14,-2.93,-2.72,-2.525,-2.33,-2.155,-1.98,-1.825,-1.67,-1.54,-1.41,-1.305,-1.2,-1.12,-1.04,-0.98,-0.92},
{-5.49,-5.195,-4.9,-4.64,-4.38,-4.14,-3.9,-3.675,-3.45,-3.24,-3.03,-2.835,-2.64,-2.46,-2.28,-2.115,-1.95,-1.81,-1.67,-1.55,-1.43,-1.33,-1.23,-1.15,-1.07,-1.01,-0.95},
{-5.49,-5.225,-4.96,-4.705,-4.45,-4.21,-3.97,-3.745,-3.52,-3.315,-3.11,-2.92,-2.73,-2.56,-2.39,-2.235,-2.08,-1.94,-1.8,-1.67,-1.54,-1.42,-1.3,-1.195,-1.09,-0.99,-0.89}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}}};


//-----------------------------------------------------------------------------------
//For charge particle
double rbinrngs1[nvar][nHLTmx][rnmxbins+1] ={
//double binrngs1[nvar][nHLTmx][37] ={
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //Thrust 3 
       {{-5.29,-4.93,-4.57,-4.335,-4.1,-3.9,-3.7,-3.515,-3.33,-3.145,-2.96,-2.775,-2.59,-2.4,-2.21,-2.025,-1.84,-1.67,-1.5,-1.36,-1.22,-1.11,-1},
{-5.29,-4.97,-4.65,-4.425,-4.2,-4,-3.8,-3.61,-3.42,-3.225,-3.03,-2.83,-2.63,-2.43,-2.23,-2.035,-1.84,-1.665,-1.49,-1.35,-1.21,-1.115,-1.02,-0.965,-0.91},
{-5.49,-5.195,-4.9,-4.67,-4.44,-4.23,-4.02,-3.815,-3.61,-3.4,-3.19,-2.975,-2.76,-2.545,-2.33,-2.12,-1.91,-1.725,-1.54,-1.395,-1.25,-1.15,-1.05,-0.99,-0.93},
{-5.54,-5.27,-5,-4.765,-4.53,-4.31,-4.09,-3.87,-3.65,-3.425,-3.2,-2.975,-2.75,-2.53,-2.31,-2.105,-1.9,-1.72,-1.54,-1.4,-1.26,-1.16,-1.06,-0.995,-0.93},
{-5.65,-5.38,-5.11,-4.865,-4.62,-4.385,-4.15,-3.915,-3.68,-3.445,-3.21,-2.975,-2.74,-2.51,-2.28,-2.07,-1.86,-1.685,-1.51,-1.375,-1.24,-1.15,-1.06,-1,-0.94,-0.905,-0.87},
{-5.65,-5.38,-5.11,-4.85,-4.59,-4.335,-4.08,-3.82,-3.56,-3.305,-3.05,-2.8,-2.55,-2.32,-2.09,-1.89,-1.69,-1.53,-1.37,-1.25,-1.13,-1.05,-0.97,-0.92,-0.87},
{-5.66,-5.37,-5.08,-4.8,-4.52,-4.24,-3.96,-3.68,-3.4,-3.125,-2.85,-2.59,-2.33,-2.095,-1.86,-1.67,-1.48,-1.34,-1.2,-1.11,-1.02,-0.96,-0.9},
{-5.54,-5.215,-4.89,-4.555,-4.22,-3.88,-3.54,-3.205,-2.87,-2.565,-2.26,-2,-1.74,-1.535,-1.33,-1.19,-1.05,-0.96,-0.87}},
	{{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoT 9
  {{-5.8,-5.195,-4.59,-4.14,-3.69,-3.315,-2.94,-2.615,-2.29,-2,-1.71,-1.45,-1.19,-0.96,-0.73,-0.53,-0.33},
{-5.8,-5.225,-4.65,-4.23,-3.81,-3.455,-3.1,-2.78,-2.46,-2.165,-1.87,-1.6,-1.33,-1.095,-0.86,-0.66,-0.46},
{-5.8,-5.27,-4.74,-4.36,-3.98,-3.65,-3.32,-3.01,-2.7,-2.4,-2.1,-1.82,-1.54,-1.295,-1.05,-0.85,-0.65},
{-6.62,-6.01,-5.4,-4.975,-4.55,-4.19,-3.83,-3.495,-3.16,-2.845,-2.53,-2.23,-1.93,-1.655,-1.38,-1.14,-0.9},
{-6.62,-6.025,-5.43,-5.01,-4.59,-4.23,-3.87,-3.54,-3.21,-2.895,-2.58,-2.28,-1.98,-1.705,-1.43,-1.195,-0.96,-0.77,-0.58},
{-6.62,-6.045,-5.47,-5.04,-4.61,-4.235,-3.86,-3.51,-3.16,-2.825,-2.49,-2.175,-1.86,-1.58,-1.3,-1.065,-0.83,-0.65,-0.47},
{-6.62,-6.035,-5.45,-5.01,-4.57,-4.18,-3.79,-3.42,-3.05,-2.695,-2.34,-2.01,-1.68,-1.395,-1.11,-0.885,-0.66},
{-5.8,-5.35,-4.9,-4.39,-3.88,-3.415,-2.95,-2.575,-2.2,-1.89,-1.58,-1.3,-1.02,-0.74,-0.46}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
//        {{},{},{},{},{},{},{},{}}, //extra
        //Y23 15
        {{-6.89,-6.15,-5.41,-4.67,-3.93,-3.275,-2.62,-1.925,-1.23},
{-6.89,-6.155,-5.42,-4.68,-3.94,-3.305,-2.67,-1.975,-1.28},
{-6.89,-6.18,-5.47,-4.78,-4.09,-3.485,-2.88,-2.245,-1.61},
{-6.89,-6.18,-5.47,-4.795,-4.12,-3.53,-2.94,-2.33,-1.72},
{-6.89,-6.235,-5.58,-4.93,-4.28,-3.685,-3.09,-2.495,-1.9},
{-6.89,-6.22,-5.55,-4.895,-4.24,-3.63,-3.02,-2.4,-1.78},
{-6.08,-5.495,-4.91,-4.24,-3.57,-2.945,-2.32,-1.585,-0.85},
{-6.08,-5.495,-4.91,-4.24,-3.57,-2.945,-2.32,-1.585,-0.85}}, 
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //BT
    {{-3.35,-3.13,-2.91,-2.735,-2.56,-2.4,-2.24,-2.085,-1.93,-1.77,-1.61,-1.45,-1.29,-1.13,-0.97,-0.815,-0.66,-0.525,-0.39},
{-3.76,-3.495,-3.23,-3.035,-2.84,-2.67,-2.5,-2.335,-2.17,-2,-1.83,-1.655,-1.48,-1.305,-1.13,-0.96,-0.79,-0.64,-0.49},
{-3.76,-3.52,-3.28,-3.09,-2.9,-2.725,-2.55,-2.375,-2.2,-2.015,-1.83,-1.64,-1.45,-1.26,-1.07,-0.9,-0.73,-0.595,-0.46},
{-3.76,-3.535,-3.31,-3.115,-2.92,-2.735,-2.55,-2.36,-2.17,-1.98,-1.79,-1.595,-1.4,-1.215,-1.03,-0.86,-0.69,-0.55,-0.41},
{-4.17,-3.905,-3.64,-3.415,-3.19,-2.98,-2.77,-2.565,-2.36,-2.15,-1.94,-1.735,-1.53,-1.335,-1.14,-0.965,-0.79,-0.65,-0.51},
{-4.02,-3.765,-3.51,-3.28,-3.05,-2.825,-2.6,-2.375,-2.15,-1.925,-1.7,-1.485,-1.27,-1.08,-0.89,-0.74,-0.59,-0.485,-0.38},
{-3.79,-3.525,-3.26,-3,-2.74,-2.485,-2.23,-1.985,-1.74,-1.51,-1.28,-1.08,-0.88,-0.72,-0.56,-0.44,-0.32},
{-3.57,-3.245,-2.92,-2.615,-2.31,-2.055,-1.8,-1.58,-1.36,-1.16,-0.96,-0.755,-0.55,-0.315,-0.08}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        //RhoTT 24
        {{-5.49,-4.85,-4.21,-3.785,-3.36,-3.015,-2.67,-2.37,-2.07,-1.795,-1.52,-1.265,-1.01},
{-5.49,-4.89,-4.29,-3.905,-3.52,-3.19,-2.86,-2.555,-2.25,-1.965,-1.68,-1.435,-1.19,-1,-0.81},
{-5.89,-5.28,-4.67,-4.28,-3.89,-3.555,-3.22,-2.91,-2.6,-2.315,-2.03,-1.77,-1.51,-1.295,-1.08,-0.92,-0.76},
{-5.89,-5.35,-4.81,-4.435,-4.06,-3.73,-3.4,-3.09,-2.78,-2.49,-2.2,-1.935,-1.67,-1.45,-1.23,-1.065,-0.9},
{-5.89,-5.385,-4.88,-4.505,-4.13,-3.795,-3.46,-3.145,-2.83,-2.535,-2.24,-1.975,-1.71,-1.49,-1.27,-1.11,-0.95},
{-5.89,-5.42,-4.95,-4.565,-4.18,-3.83,-3.48,-3.16,-2.84,-2.545,-2.25,-1.985,-1.72,-1.495,-1.27,-1.09,-0.91},
{-5.89,-5.425,-4.96,-4.545,-4.13,-3.755,-3.38,-3.035,-2.69,-2.38,-2.07,-1.8,-1.53,-1.3,-1.07,-0.885,-0.7},
{-5.89,-5.41,-4.93,-4.53,-4.13,-3.78,-3.43,-3.105,-2.78,-2.465,-2.15,-1.835,-1.52,-1.19,-0.86}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}},
        {{},{},{},{},{},{},{},{}}};



