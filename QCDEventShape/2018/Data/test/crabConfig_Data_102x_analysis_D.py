#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2018A'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2018B'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2018C'
config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2018D'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_102x_data_cfg.py'

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

#config.Data.inputDataset = '/JetHT/Run2018A-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = ''
#config.Data.inputDataset = ''

#config.Data.inputDataset = '/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset ='/JetHT/Run2018D-UL2018_MiniAODv2-v1/MINIAOD'


config.Data.inputDataset = '/JetHT/Run2018D-12Nov2019_UL2018_rsb-v1/MINIAOD'


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.lumiMask ='HLT_13TeV_UL2018_Collisions18_GoldenJSON.txt'


#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
