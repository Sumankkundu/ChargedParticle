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
'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FEA2ED14-5CDF-E711-ACA6-02163E012AF0.root',
'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FE211553-36DF-E711-BAB7-02163E019BD0.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FE155D02-00DF-E711-BA34-02163E011A55.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/FE08F446-63DF-E711-A338-A4BF0112BCF8.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/750/00000/28938773-BD72-E511-A479-02163E01432A.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/0075E33B-3B6C-E511-BCC8-02163E01455C.root'
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/0CE8F23E-3B6C-E511-B68A-02163E013744.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/36DC8060-3B6C-E511-BC73-02163E0143DD.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/50A3A073-3B6C-E511-A997-02163E0144CD.root'
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/301A497D-70B0-E511-9630-002590D0AFA8.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7210C351-67B0-E511-A34C-7845C4FC37AF.root'
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/745E2A4F-67B0-E511-9DA3-0090FAA57620.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7E46D250-67B0-E511-BB96-0025905C3E66.root'
 )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v2')
process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v6')
from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

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
    fileName=cms.string("Test_Data_QCD_char_2017.root")
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
        rho = cms.InputTag('fixedGridRhoAll'),
        LHEEventProductInputTag   = cms.InputTag('externalLHEProducer'),
        LHERunInfoProductInputTag = cms.InputTag('externalLHEProducer'),
        PDFCTEQWeightsInputTag   = cms.InputTag('pdfWeights:CT14'),
        PDFMMTHWeightsInputTag   = cms.InputTag('pdfWeights:MMHT2014lo68cl'),
        PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights:NNPDF30'),
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
#        scaleFactorsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'), 
#        scaleFactorsFile = cms.FileInPath('Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
#        scaleFactorsFile = cms.FileInPath('Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),

      
 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)



