from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'tutorial_May2015_MC_analysis_herwihplusplus'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_mc_cfg.py'

config.Data.inputDataset ='/QCD_Pt-35toInf_fwdJet_bwdJet_TuneCUETHS1_13TeV-herwigpp/RunIIFall15MiniAODv2-PU25nsData2015v1_castor_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'MC_analysis_Pythia8_flat'

config.Site.storageSite ='T2_IN_TIFR'
