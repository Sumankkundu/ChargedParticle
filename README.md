# Chargedparticle
Study of Event Shape Variables with the charged particles inside jets 

Analysis code based on MINIAOD Format
# Data Used: 
```
/JetHT/Run2017F-17Nov2017-v1/MINIAOD
/JetHT/Run2017E-17Nov2017-v1/MINIAOD
/JetHT/Run2017D-17Nov2017-v1/MINIAOD
/JetHT/Run2017C-17Nov2017-v1/MINIAOD
/JetHT/Run2017B-17Nov2017-v1/MINIAOD
```
# MC Used :  
Flat PYTHIA sample
```
file dataset=/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
```
```
git clone https://github.com/Sumankkundu/Chargedparticle.git
```
# How to run :
```
cmsrel  CMSSW_10_2_18
cd  CMSSW_10_2_18/src
cmsenv
mkdir Test
mkedanlzr QCDEventShape
cd QCDEventShape
cd plugin
```
**In plugin put those files: 
               QCDEventShape.cc(Replace existing one) , EventShape_vector.cc(Copy)& EventShape_vector.h(copy) **
 ```              
 cd test 
 ```
 **In test 
   put those .py files from the  git folder $Chargedparticle/QCDEventShape/test/ **
 ``` 
 cd QCDEventShape
 scramv1 b    **this will complie 
 cd test 
 cmsRun Run_QCD_test_miaod_v2_76x_mc_cfg.py
```


# TO Do List 

- [x] Reco distribution : Done
- [ ] Basic Kinematics  : not done
- [ ] Trigger           : not done 
- [ ] Unfolding         : not done
