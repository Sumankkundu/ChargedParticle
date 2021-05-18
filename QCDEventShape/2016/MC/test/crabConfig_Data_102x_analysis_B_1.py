#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016B_v1'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016B_v2'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016C'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016D'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016E'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016F'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016G'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2016H'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_102x_data_cfg.py'

config.JobType.inputFiles= [
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Summer19UL16APV_RunBCDEF_V5_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Summer19UL16APV_RunBCD_V7_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Summer19UL16APV_RunEF_V7_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Spring16_25nsV10_MC_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2016/AK4PFCHS_Summer19UL/Summer16_25nsV5_MC_Uncertainty_AK4PFchs.txt",
]



config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset ='/JetHT/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset ='/JetHT/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD'

#Post VFP
#config.Data.inputDataset ='/JetHT/Run2016F-21Feb2020_UL2016-v1/MINIAOD'
#config.Data.inputDataset ='/JetHT/Run2016G-21Feb2020_UL2016-v1/MINIAOD'
#config.Data.inputDataset ='/JetHT/Run2016H-21Feb2020_UL2016-v1/MINIAOD'



config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.lumiMask ='HLT_13TeV_UL2016_Collisions16_GoldenJSON.txt'


#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
