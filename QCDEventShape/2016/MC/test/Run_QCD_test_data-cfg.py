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
#'/store/data/Run2015C/JetHT/MINIAOD/PromptReco-v1/000/253/808/00000/3E6025B5-7340-E511-A8B7-02163E01440E.root'
#'/store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/161/00000/CACA9422-9C26-E511-AD2E-02163E0117FF.root'
#'/store/data/Run2015C_25ns/JetHT/MINIAOD/05Oct2015-v1/60000/1E4C4F31-7374-E511-875C-0025905A6138.root'
'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/750/00000/28938773-BD72-E511-A479-02163E01432A.root'
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/0CE8F23E-3B6C-E511-B68A-02163E013744.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/36DC8060-3B6C-E511-BC73-02163E0143DD.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/50A3A073-3B6C-E511-A997-02163E0144CD.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/709AEE6D-3B6C-E511-A634-02163E0145B8.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/7884066E-3B6C-E511-8733-02163E014218.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/B02ABB3F-3B6C-E511-9CDD-02163E0128DC.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/B8D31E3F-3B6C-E511-936E-02163E014552.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/174/00000/6E079318-CC6C-E511-B435-02163E01410C.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/175/00000/64BD481A-FA6C-E511-AFBC-02163E0135AC.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/175/00000/B89DC103-FA6C-E511-B5A0-02163E011D7E.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/024E87D4-776D-E511-BE59-02163E01375C.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/082E00CB-776D-E511-9549-02163E01468E.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/1497C06A-986D-E511-9463-02163E01422B.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/14D0A3EB-506D-E511-8F0B-02163E013641.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/2CC51A3E-706D-E511-B139-02163E0143D6.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/3AC7D0E8-776D-E511-A88E-02163E01463A.root'
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_P_V56::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_R_44_V11::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'74X_dataRun2_Prompt_v1')
process.GlobalTag = GlobalTag(process.GlobalTag,'74X_dataRun2_Prompt_v4')
from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("HLTrigger.HLTcore.hltPrescaleRecorder_cfi")

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
    fileName=cms.string("Test_QCD_data.root")
)
print "test1"
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
        bsSrc = cms.InputTag("offlineBeamSpot"),
        genjetSrc = cms.InputTag("slimmedGenJets"),
        pileupSrc =cms.InputTag("slimmedAddPileupInfo"),
        ak5pfJetSrc = cms.InputTag("ak5PFJets"),
        ak5genJetSrc = cms.InputTag("ak5GenJets"),
        evtinfo =cms.InputTag("generator"),
        #ak5PFJetCHSSrc    = cms.InputTag("ak5PFJetsCHS")
	RootFileName = cms.untracked.string('pythia8_test_13tev.root'),
	GenJET =  cms.untracked.bool(False),
	HistFill = cms.untracked.bool(True),
	MonteCarlo =  cms.untracked.bool(False),
	ParticleLabel =  cms.untracked.bool(False),
	Reconstruct =cms.untracked.bool(True),
#  EtaRange =  cms.untracked.double(5.0),
#  PtThreshold = cms.untracked.double(12.0),
  	EtaRange =  cms.untracked.double(3.0),
  	PtThreshold = cms.untracked.double(55.0), #effective is 21
  	LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81       

 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



