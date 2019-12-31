from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName ='QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_76x_mc_cfg.py'

config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_analysis_Pythia8_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
