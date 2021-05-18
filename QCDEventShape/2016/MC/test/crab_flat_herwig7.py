#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

config.General.requestName ='ESV_CH3_herwig7_Flat'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_MiniAOD_102x_mc_cfg.py'
#config.JobType.psetName = options.cfg
#config.JobType.maxMemoryMB = 9000 # Default is 2500 : Max I have used is 13000
#config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000
#config.JobType.numCores = 4


config.JobType.inputFiles= [
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunB_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunC_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunD_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunA_V5_DATA_Uncertainty_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunB_V5_DATA_Uncertainty_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunC_V5_DATA_Uncertainty_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2018/AK4PFCHS_Summer19UL/Summer19UL18_RunD_V5_DATA_Uncertainty_AK4PFchs.txt"
]





config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring15DR74-HFscaleFlat10to30Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'EventBased'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 10  # for Automatic must be 180-2700 range
config.Data.unitsPerJob = 1  #For Filebased or Lumibased
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/%s/' % (sukundu)
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_UL2018_Herwig7_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
