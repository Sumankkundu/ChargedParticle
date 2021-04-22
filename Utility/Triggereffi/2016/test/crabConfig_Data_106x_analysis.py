from CRABClient.UserUtilities import config
config = config()


#config.General.requestName = 'Trigger_2016UL_B_v1'
#config.General.requestName = 'Trigger_2016UL_B_v2'
#config.General.requestName = 'Trigger_2016UL_C'
#config.General.requestName = 'Trigger_2016UL_D'
#config.General.requestName = 'Trigger_2016UL_E'
#config.General.requestName = 'Trigger_2016UL_F'
#config.General.requestName = 'Trigger_2016UL_G'
config.General.requestName = 'Trigger_2016UL_H'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'
config.JobType.inputFiles= [
#"/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt"
]

#config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-21Feb2020_UL2016-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016G-21Feb2020_UL2016-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2016H-21Feb2020_UL2016-v1/MINIAOD'


config.Data.inputDBS = 'global'

#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.unitsPerJob = 200

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'


#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
