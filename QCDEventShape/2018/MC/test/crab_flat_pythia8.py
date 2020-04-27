from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName ='ESVQCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8_Autumn2018_2020Mar29'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_94x_mc_cfg.py'

config.JobType.inputFiles= ["/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_V7_DATA_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_V7_DATA_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_V7_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_V7_MC_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_RunD_V19_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_RunC_V19_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_RunB_V19_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2018/AK4PFCHS/Autumn18_RunA_V19_DATA_UncertaintySources_AK4PFchs.txt"
]

#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIAutumn18MiniAOD-FlatPU0to70_102X_upgrade2018_realistic_v15_ext4-v1/MINIAODSIM'
config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIIAutumn18MiniAOD-FlatPU0to70RAW_HEM_102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM'
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
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_analysis_Pythia8_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
