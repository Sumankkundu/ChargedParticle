import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B6F8B6-90F1-E311-B72C-0025905A6092.root'
#'/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_12Apr2018_94X_mc2017_realistic_v14-v1/120000/FE40878F-EDEA-E811-85E5-0025905A611C.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/F3A1AC28-A59D-A843-A675-33AE077F75F0.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/C6F0958E-BFB8-7C43-B7BD-02CAC1B592C1.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/A24FE704-61D8-BC4F-81EB-C32148AA7D90.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/3F3D177E-B73F-BF46-AAE9-6589D4476573.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/329292B7-ACB8-F144-ABA5-A7348ECB5B0E.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/1D7395E7-1A4C-9247-B0E1-5192C2B5FFBF.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/18A9FF14-EB08-1C4E-97A9-F9CC6CA7CEDE.root',
#'/store/relval/CMSSW_10_6_14_Pyt8240BugFix/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/0DFAB61C-8207-A74C-A348-5D11770C16EB.root',
#'/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/7E86CBA6-9842-E811-B26F-0CC47AA53D6E.root',
#'/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/7E65A585-9C42-E811-BA22-0025901ABB72.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/FADB5AF4-89AA-A54A-8203-DF285FCB5F5B.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/EAAD438A-5B7C-6C46-8D94-A3CBFF3FAC54.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/E5E1D2B5-FE94-6E4C-A00F-E7D49FF3CBCB.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/E3DC4BED-D9FE-AA42-B2EA-FE2016BCEEFB.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/E265A236-EC64-104F-8A38-614F7942AF12.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/C0D7C8F9-14EC-5840-8D5B-EA33167390DD.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/B7EFD98B-6DCF-B14F-A491-9C4F1A25E62E.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/AA98982C-15EB-7B43-BF1E-858C0F1E5A47.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/89A3A3C6-090A-4847-A324-4D4CAE0E8C02.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/72064DDF-4849-D94A-9B30-6BBCC549C572.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/6195B1D0-5C13-D840-A6A6-084477B0140D.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/44ACA960-FAAD-8A4B-8630-FE9A55BA8F89.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/4153783F-CAD2-134D-8EC7-18707038EB4A.root',
'/store/relval/CMSSW_10_6_14/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_106X_mc2017_realistic_v7_HS-v1/10000/0BDE025E-81CB-D141-8A00-560212DCDC99.root',
 )
#eventsToSkip = cms.untracked.VEventRange('1:1950-1:2000'),
#eventsToSkip = cms.untracked.EventRange('1:1950-1:2000'),
#eventRanges = cms.untracked.VEventRange('1:1000-1:2000'),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mc2017_realistic_v7')
#process.GlobalTag = GlobalTag(process.GlobalTag,'74X_mcRun2_asymptotic_realisticBS_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'MCRUN2_74_V9')
from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("HLTrigger.HLTcore.hltPrescaleRecorder_cfi")


process.MessageLogger = cms.Service("MessageLogger",  
 cout = cms.untracked.PSet(  
  default = cms.untracked.PSet( ## kill all messages in the log  
 
   limit = cms.untracked.int32(0)  
  ),  
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted  
 
   limit = cms.untracked.int32(-1)  
  )  
 ),  
 categories = cms.untracked.vstring('FwkJob'), 
 destinations = cms.untracked.vstring('cout')  
)  





#ak5 PF & Gen Jets
#from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
#from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
#from RecoMET.METProducers.PFMET_cfi import pfMet


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')


# Select candidates that would pass CHS requirements
#process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#makes chs ak5 jets   (instead of ak4 that are default in miniAOD )
#process.ak5PFJetsCHS = ak5PFJets.clone(src = 'chs')

process.TFileService=cms.Service("TFileService",
    #fileName=cms.string("Test_MC2017_PY8Bugfix.root")
    fileName=cms.string("Test_MC2017_PY8_old.root")
)
print "test1"
# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
# Fix POWHEG if buggy (this PDF set will also appear on output,
# so only two more PDF sets can be added in PdfSetNames if not "")
#FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
#GenTag = cms.untracked.InputTag("genParticles"),
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring(
#               "CT10nlo.LHgrid"
                 "PDF4LHC15_100.LHgrid"
                , "MSTW2008nlo68cl.LHgrid"
                , "NNPDF30_100.LHgrid"
            )
      )

process.analyzeBasicPat = cms.EDAnalyzer("QCDEventShape",
#       photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
#       electronSrc = cms.untracked.InputTag("cleanPatElectrons"),
#       muonSrc = cms.untracked.InputTag("cleanPatMuons"),
#       tauSrc = cms.untracked.InputTag("cleanPatTaus"),
#        r = cms.EventRange('1:1000-1:2000'),
        jetSrc = cms.InputTag("slimmedJets"),
        metSrc = cms.InputTag("slimmedMETs"),
        genSrc = cms.untracked.InputTag("packedGenParticles"),
        pfSrc = cms.InputTag("packedPFCandidates"),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
        objects = cms.InputTag("selectedPatTrigger"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        bsSrc = cms.InputTag("offlineBeamSpot"),
        genjetSrc = cms.InputTag("slimmedGenJets"),
        pileupSrc =cms.InputTag("slimmedAddPileupInfo"),
        ak5pfJetSrc = cms.InputTag("ak5PFJets"),
        ak5genJetSrc = cms.InputTag("ak5GenJets"),
        evtinfo =cms.InputTag("generator"),
        rho = cms.InputTag('fixedGridRhoAll'),
        LHEEventProductInputTag   = cms.InputTag('externalLHEProducer'),
        LHERunInfoProductInputTag = cms.InputTag('externalLHEProducer'),
        PDFCTEQWeightsInputTag   = cms.InputTag('pdfWeights:CT14'),
        PDFMMTHWeightsInputTag   = cms.InputTag('pdfWeights:MMHT2014lo68cl'),
        PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights:NNPDF30'),
        #ak5PFJetCHSSrc    = cms.InputTag("ak5PFJetsCHS")
	RootFileName = cms.untracked.string('pythia8_test_13tev.root'),
	GenJET =  cms.untracked.bool(True),
	HistFill = cms.untracked.bool(True),
	MonteCarlo =  cms.untracked.bool(True),
	ParticleLabel =  cms.untracked.bool(False),
	Reconstruct =cms.untracked.bool(True),
#  EtaRange =  cms.untracked.double(5.0),
#  PtThreshold = cms.untracked.double(12.0),
  	EtaRange =  cms.untracked.double(3.0),
  	PtThreshold = cms.untracked.double(55.0), #effective is 21
  	LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81       
#        scaleFactorsFile = cms.FileInPath('xxCondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('xxCondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        scaleFactorsFile = cms.FileInPath('Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),


 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



