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
 #'/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/164247BB-6577-E411-A063-00259073E50A.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt50nsRecodebug_MCRUN2_74_V9A-v1/10000/06C3281B-CE01-E511-B1B1-002590593920.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25nsRecodebug_MCRUN2_74_V9-v1/50000/4887CFD0-7303-E511-9CAA-001D09FDD7B9.root'
# 'file:/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_4_1/src/Test/QCDEventShape/test/0E4CEBFE-ECFB-E411-9F0C-842B2B29273C.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/2E68E42D-EDFB-E411-8027-001E67397CC9.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/6E887F0F-EDFB-E411-875B-BCAEC54B303A.root',
#'/store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/EA06BB15-EDFB-E411-9672-003048F010F4.root',
#'/store/mc/RunIISpring15DR74/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/50000/1E726D13-0A0D-E511-BB09-782BCB6E134C.root',
#'/store/mc/RunIISpring15DR74/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/1274D4AF-12FC-E411-9A65-0025B31E330A.root'
#'/store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/521/00000/0C2A2945-6E29-E511-A6F9-02163E013576.root'
#'/store/mc/RunIISpring15DR74/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/80000/027FF217-8729-E511-B60B-782BCB408BD2.root',
#'/store/mc/RunIISpring15DR74/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/80000/126ACEC6-0728-E511-B60B-0025905B8572.root'
'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/AsymptFlat0to50bx50Reco_MCRUN2_74_V9A-v3/10000/0E70AE5E-7B08-E511-83A7-02163E00EACB.root'
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#ak5 PF & Gen Jets
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from RecoMET.METProducers.PFMET_cfi import pfMet


process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
process.ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')


# Select candidates that would pass CHS requirements
#process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#makes chs ak5 jets   (instead of ak4 that are default in miniAOD )
#process.ak5PFJetsCHS = ak5PFJets.clone(src = 'chs')




process.TFileService=cms.Service("TFileService",
    fileName=cms.string("Test_MCQCD.root")
)

process.analyzeBasicPat = cms.EDAnalyzer("QCDEventShape",
#       photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
#       electronSrc = cms.untracked.InputTag("cleanPatElectrons"),
#       muonSrc = cms.untracked.InputTag("cleanPatMuons"),
#       tauSrc = cms.untracked.InputTag("cleanPatTaus"),
        jetSrc = cms.InputTag("slimmedJets"),
        metSrc = cms.InputTag("slimmedMETs"),
        genSrc = cms.untracked.InputTag("packedGenParticles"),
        pfSrc = cms.InputTag("packedPFCandidates"),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
        objects = cms.InputTag("selectedPatTrigger"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        genjetSrc = cms.InputTag("slimmedGenJets"),
        ak5pfJetSrc = cms.InputTag("ak5PFJets"),
        ak5genJetSrc = cms.InputTag("ak5GenJets"),
        evtinfo =cms.InputTag("generator"),
        #ak5PFJetCHSSrc    = cms.InputTag("ak5PFJetsCHS")
	RootFileName = cms.untracked.string('pythia8_test_13tev.root'),
	GenJET =  cms.untracked.bool(True),
	HistFill = cms.untracked.bool(True),
	MonteCarlo =  cms.untracked.bool(True),
	ParticleLabel =  cms.untracked.bool(False),
	Reconstruct =cms.untracked.bool(False),
#  EtaRange =  cms.untracked.double(5.0),
#  PtThreshold = cms.untracked.double(12.0),
  	EtaRange =  cms.untracked.double(3.0),
  	PtThreshold = cms.untracked.double(55.0), #effective is 21
  	LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81       

 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
#process.p = cms.Path(process.analyzeBasicPat)
process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



