from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName ='QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_2017_FileBased'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_76x_mc_cfg.py'

config.JobType.inputFiles= ["/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_V2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_V2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_17Nov2017F_V6_DATA_UncertaintySources_AK4PFchs.txt","/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt","/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_17Nov2017_V6_MC_UncertaintySources_AK4PFchs.txt" ]


#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-NoPU_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset ='/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring15DR74-HFscaleFlat10to30Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'EventBased'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 10  # for Automatic must be 180-2700 range
config.Data.unitsPerJob = 15  #For Filebased or Lumibased
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_analysis_Pythia8_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
