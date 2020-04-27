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
#'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/AsymptFlat0to50bx50Reco_MCRUN2_74_V9A-v3/10000/0E70AE5E-7B08-E511-83A7-02163E00EACB.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25nsRecodebug_MCRUN2_74_V9-v1/50000/02787365-4F03-E511-BAD2-842B2B182385.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/60000/161AEAEB-7E0C-E511-9181-0025B3E05C60.root'
#'/store/mc/RunIISpring15DR74/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/221B185A-3502-E511-AC72-0025905C4262.root'
#'/store/mc/RunIISpring15MiniAODv2/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/F4E78EB7-9272-E511-823C-0026B9278678.root'
#'/store/mc/RunIISpring15MiniAODv2/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/042F0EE3-9B71-E511-98C7-B083FECFC6ED.root'
#'/store/mc/RunIIFall15MiniAODv2/QCD_Pt-35toInf_fwdJet_bwdJet_TuneCUETHS1_13TeV-herwigpp/MINIAODSIM/PU25nsData2015v1_castor_76X_mcRun2_asymptotic_v12-v1/60000/4E65EB29-E64F-E611-BF1B-D8D385AE8848.root'
'root://cms-xrd-global.cern.ch///store/mc/RunIIFall15MiniAODv1/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/20E80D49-5AA5-E511-84A7-0025905A60D0.root'
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
process.GlobalTag = GlobalTag(process.GlobalTag,'76X_mcRun2_asymptotic_RunIIFall15DR76_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'74X_mcRun2_asymptotic_realisticBS_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'MCRUN2_74_V9')
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
    fileName=cms.string("Test_MC_QCD.root")
)
#print "test1"
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
        rho = cms.InputTag('fixedGridRhoAll'),
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
#        scaleFactorsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
#print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



