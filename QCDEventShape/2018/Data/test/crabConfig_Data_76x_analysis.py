from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017F-31Mar2018-v1_6Mar20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017E-31Mar2018-v1_6Mar20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017D-31Mar2018-v1_6Mar20'
#config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017C-31Mar2018-v1_6Mar20'
config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017B-31Mar2018-v1_6Mar20'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_76x_data_cfg.py'

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

#config.JobType.inputFiles= ["/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt"]


#config.Data.inputDataset = '/JetHT/Run2017F-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017E-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017B-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-17Nov2017-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2017F-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017E-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-31Mar2018-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2017B-31Mar2018-v1/MINIAOD'


#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v4/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v3/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2015C_25ns-05Oct2015-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/MINIAOD';

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.lumiMask ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.runRange = '246908-260627' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
