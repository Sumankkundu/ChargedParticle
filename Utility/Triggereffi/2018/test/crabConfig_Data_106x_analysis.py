from CRABClient.UserUtilities import config
config = config()


#config.General.requestName = 'Trigger_2018UL_A'
#config.General.requestName = 'Trigger_2018UL_B'
#config.General.requestName = 'Trigger_2018UL_C'
config.General.requestName = 'Trigger_2018UL_D'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'
config.JobType.inputFiles= [
#"/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt"
]
#config.Data.inputDataset = '/JetHT/Run2018A-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018C-12Nov2019_UL2018_rsb-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2018D-12Nov2019_UL2018_rsb-v1/MINIAOD'



#config.Data.inputDataset = ' /JetHT/Run2018D-UL2018_MiniAODv2-v1/MINIAOD'

config.Data.inputDBS = 'global'

#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.unitsPerJob = 200

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

#config.Data.runRange = '246908-260627' # '193093-194075'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
