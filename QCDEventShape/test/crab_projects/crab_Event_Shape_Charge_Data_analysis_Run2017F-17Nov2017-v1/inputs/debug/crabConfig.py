from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = True
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'Event_Shape_Charge_Data_analysis_Run2017F-17Nov2017-v1'
config.section_('JobType')
config.JobType.psetName = 'Run_QCD_test_76x_data_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_V2_MC_PtResolution_AK4PFchs.txt', '/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_V2_MC_SF_AK4PFchs.txt', '/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall17_17Nov2017F_V6_DATA_UncertaintySources_AK4PFchs.txt', '/afs/cern.ch/work/s/sukundu/private/ESV_Charge_Jet/CMSSW_9_4_6/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt']
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/JetHT/Run2017F-17Nov2017-v1/MINIAOD'
config.Data.publication = False
config.Data.unitsPerJob = 15
config.Data.splitting = 'LumiBased'
config.Data.inputDBS = 'global'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.outLFNDirBase = '/store/user/sukundu/'
config.section_('Site')
config.Site.storageSite = 'T2_IN_TIFR'
config.section_('User')
config.section_('Debug')
