from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017F-UL_14Oct20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017E-UL_14Oct20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017D-UL_14Oct20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017C-UL_14Oct20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017B-UL_14Oct20'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_102x_data_cfg.py'

config.JobType.inputFiles= ["/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/Fall17_V2_MC/Fall17_V2_MC_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/Fall17_V2_MC/Fall17_V2_MC_SF_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_RunB_V4_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_RunC_V4_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_RunD_V4_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_RunE_V4_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_RunF_V4_DATA_UncertaintySources_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_JRV2_DATA_PtResolution_AK4PFchs.txt",
"/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/DataMC17UL_V4/Summer19UL17_JRV2_DATA_SF_AK4PFchs.txt",
]

#config.JobType.inputFiles= ["/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt"]


#config.Data.inputDataset = '/JetHT/Run2017E-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017B-31Mar2018-v1/MINIAOD'


#config.Data.inputDataset = '/JetHT/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2017F-09Aug2019_UL2017-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v4/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW_10_2_18/Uncertainty2017/cert/HLT_diPFJetAve_2017_ReReco_v1.txt'
config.Data.lumiMask = '/afs/cern.ch/work/s/sukundu/private/ESV_charge_CMSSW/Uncertainty2017/cert/HLT_13TeV_UL2017_Collisions17_GoldenJSON.txt'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.lumiMask ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.runRange = '246908-260627' # '193093-194075'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
