#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

#config.General.requestName ='ESVQCD_UL_Ptbinned_3200toinf_tuneCP5_bin_16July20'
config.General.requestName ='ESV_CP5_PY82017UL_Flat_2D'
#config.General.requestName ='ESV_Pt_C5_PY82017_249_14PU94X_20APr_CMSSW94X'
#config.General.requestName ='ESV_Pt_TuneCUETP8M1_PY82017_14PU_GT94X_LJEC_23Apr_CMSSW94X'
#config.General.requestName ='ESV_Pt_TuneCUETP8M1_PY82017_OLPU_GT94X_LJEC_23Apr_CMSSW94X'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_106x_mc_cfg.py'

#config.JobType.psetName = options.cfg
#config.JobType.maxMemoryMB = 9000 # Default is 2500 : Max I have used is 13000
#config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000
#config.JobType.numCores = 4

config.JobType.inputFiles= [
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_RunB_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_RunC_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_RunD_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_RunE_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/AK4PFCHS_Summer19UL/Summer19UL17_RunF_V5_DATA_UncertaintySources_AK4PFchs.txt"
]


config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6_ext2-v2/MINIAODSIM'
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
config.Data.outputDatasetTag = 'MC_UL17_Py8_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
