# ESVs with Charged particle
Study of Event Shape Variables with the charged particles inside jets 

## Trigger use 
HLT_DiPFJetAve40_v, HLT_DiPFJetAve60_v ,HLT_DiPFJetAve80_v, HLT_DiPFJetAve140_v, HLT_DiPFJetAve200_v, HLT_DiPFJetAve260_v, HLT_DiPFJetAve320_v, HLT_DiPFJetAve400_v, HLT_DiPFJetAve500_v


Analysis code based on MINIAOD Format
## 2017 
**Luminosity : 41.53 /fb**
### Data Used:
```
~~/JetHT/Run2017F-17Nov2017-v1/MINIAOD
/JetHT/Run2017E-17Nov2017-v1/MINIAOD
/JetHT/Run2017D-17Nov2017-v1/MINIAOD
/JetHT/Run2017C-17Nov2017-v1/MINIAOD
/JetHT/Run2017B-17Nov2017-v1/MINIAOD ~~
New Data
/JetHT/Run2017F-31Mar2018-v1/MINIAOD
/JetHT/Run2017E-31Mar2018-v1/MINIAOD
/JetHT/Run2017D-31Mar2018-v1/MINIAOD
/JetHT/Run2017C-31Mar2018-v1/MINIAOD
/JetHT/Run2017B-31Mar2018-v1/MINIAOD
We have plan to use UL DATA 
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'

```
### MC Used  
Flat PYTHIA sample

Global tag: 94X_mc2017_realistic_v14
```
/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
/QCD_Pt-15to7000_TuneCP5_Flat2017_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
Flat Herwig++

```
### JetID    
https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2017

## 2018
### Data Used
```
JetHT/Run2018A-17Sep2018-v1/MINIAOD
/JetHT/Run2018B-17Sep2018-v1/MINIAOD
/JetHT/Run2018C-17Sep2018-v1/MINIAOD
/JetHT/Run2018D-PromptReco-v2/MINIAOD
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
```
### MC Used
```
To be Update
```
### JetID
https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
``
```
git clone https://github.com/Sumankkundu/Chargedparticle.git
```

### How to run :
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


### TO Do List 

- [x] Reco distribution : Done
- [ ] Basic Kinematics  : not done
- [ ] Trigger           : not done 
- [ ] Unfolding         : not done
