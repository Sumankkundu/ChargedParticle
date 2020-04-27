from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '2015D_Event_Shape_Data_analysis_25ns_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_data-cfg.py'

#config.Data.inputDataset = '/JetHT/Run2015C_25ns-05Oct2015-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v4/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2015D-05Oct2015-v1/MINIAOD';
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
config.Data.runRange = '254231-257599' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.Site.storageSite = 'T2_IN_TIFR'
