# ESVs with Charged particles inside jets
Study of Event Shape Variables with the charged particles inside jets 

## Trigger used
HLT_DiPFJetAve40_v, HLT_DiPFJetAve60_v ,HLT_DiPFJetAve80_v, HLT_DiPFJetAve140_v, HLT_DiPFJetAve200_v, HLT_DiPFJetAve260_v, HLT_DiPFJetAve320_v, HLT_DiPFJetAve400_v, HLT_DiPFJetAve500_v


Analysis code based on MINIAOD Format
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
## 2017 
**Luminosity : 41.53 /fb**
### Data Used:
# ReReco
```
/JetHT/Run2017F-31Mar2018-v1/MINIAOD
/JetHT/Run2017E-31Mar2018-v1/MINIAOD
/JetHT/Run2017D-31Mar2018-v1/MINIAOD
/JetHT/Run2017C-31Mar2018-v1/MINIAOD
/JetHT/Run2017B-31Mar2018-v1/MINIAOD

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
````
# UL17
```
 /JetHT/Run2017B-09Aug2019_UL2017-v1/MINIAOD
/JetHT/Run2017C-09Aug2019_UL2017-v1/MINIAOD
/JetHT/Run2017D-09Aug2019_UL2017-v1/MINIAOD
/JetHT/Run2017E-09Aug2019_UL2017-v1/MINIAOD
/JetHT/Run2017F-09Aug2019_UL2017-v1/MINIAOD

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

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
Tight Jet 
https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2017

## 2018
### Data Used(To be Updated)
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
**In QCDEventShape folder, put the code in plugin,test from the folders,
https://github.com/Sumankkundu/ChargedParticle/tree/master/QCDEventShape/2017/MC

 Plugin : QCDEventShape.cc(Replace existing one) , EventShape_vector.cc(Copy)& EventShape_vector.h(copy) 
 test : All python files
 **
 
 ```
 cd QCDEventShape
 scramv1 b    
 cd test 
 cmsRun Run_QCD_test_miaod_v2_76x_mc_cfg.py
```


### TO Do List 

- [x] Reco distribution : Done
- [x] Basic Kinematics  : Done(17.04.2020)
- [x] Trigger           : Done (17.04.2020)
- [x] 1DUnfolding      :Partially done(10.07.19) 
- [x] 2DUnfolding      : not done

