from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName ='ESV_PY8_Flat_CP5_PU_Old_Code_14PU_20April_CMSSW94X'
#config.General.requestName ='ESV_PY8_Flat_CUETP8M1_PU_Old_Code_14PU_19April_CMSSW94X'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_miaod_v2_94x_mc_cfg.py'

config.JobType.inputFiles= ["/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_V2_MC/Fall17_V2_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_V2_MC/Fall17_V2_MC_SF_AK4PFchs.txt", 
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_UncertaintySources_AK4PFchs.txt"
]

#config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall17MiniAODv2-NoPU_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat2017_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
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
#config.Data.outLFNDirBase = '/store/user/sukundu' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_analysis_Pythia8_flat'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
