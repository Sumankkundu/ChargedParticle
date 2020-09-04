import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

#DataSet  ReReco  A,B,C : Promt : D   
#Datatyp='A'
#Datatyp='B'
#Datatyp='C'
Datatyp='D'


process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B6F8B6-90F1-E311-B72C-0025905A6092.root'
#'/store/data/Run2017D/JetHT/MINIAOD/31Mar2018-v1/90000/E4402B21-F739-E811-BA23-0CC47A4C8E1E.root',
#'/store/data/Run2017B/JetHT/MINIAOD/17Nov2017-v1/60000/DA587060-EED8-E711-80CC-0CC47A4D7678.root',
#'/store/data/Run2018A/JetHT/MINIAOD/17Sep2018-v1/60000/FF52A3F2-0FA6-6F42-AAA6-9EED28B21F3D.root',
#'/store/data/Run2018A/JetHT/MINIAOD/17Sep2018-v1/60000/FEF101AF-E88F-4D42-8E05-FBFEF0F3F88E.root',
#'/store/data/Run2018A/JetHT/MINIAOD/17Sep2018-v1/60000/FBC99C48-0AF1-5E4B-ADBE-A2CE1FC065E9.root',
#'/store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/60000/96A26888-EC64-5048-A1E6-0DFB816583ED.root',
#'/store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/60000/973859C5-8DB5-AA4F-A54C-5F93ACC00FB4.root',
#'/store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/00000/0022AC1D-1DBC-C949-8910-EE2885633502.root',
#'/store/data/Run2018C/JetHT/MINIAOD/17Sep2018-v1/80000/DDC38B74-3A1C-BF4B-9B01-11A3A6A4078A.root',
#'/store/data/Run2018D/JetHT/MINIAOD/PromptReco-v2/000/325/175/00000/D8D6F54F-8CD2-FA4F-8DEF-92253BE78358.root',
'/store/data/Run2018D/JetHT/MINIAOD/PromptReco-v2/000/325/172/00000/FAA909FC-24D3-0F40-B817-9E98320AE6F5.root',
#'/store/data/Run2018D/JetHT/MINIAOD/PromptReco-v2/000/325/172/00000/C1856928-49D0-2143-B3BB-9394444DEAE9.root',
'/store/data/Run2018D/JetHT/MINIAOD/PromptReco-v2/000/325/172/00000/AEA34B67-5386-C44E-8ACD-27960F3E756D.root',
#'/store/data/Run2018C/JetHT/MINIAOD/17Sep2018-v1/80000/D764371C-8A54-E846-BC3C-AB4DFB30205E.root',
#'/store/data/Run2018C/JetHT/MINIAOD/17Sep2018-v1/80000/D3CB5A51-1FAA-BC48-BAB6-4DED7307DB96.root',
#'/store/data/Run2018C/JetHT/MINIAOD/17Sep2018-v1/80000/AC59982E-1178-1F46-A9DE-FAC659B7F66D.root',

 )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_P_V56::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_R_44_V11::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v11')
if Datatyp=='D':
        process.GlobalTag = GlobalTag(process.GlobalTag,'102X_dataRun2_Prompt_v15')
else:
	process.GlobalTag = GlobalTag(process.GlobalTag,'102X_dataRun2_v12')

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
    fileName=cms.string("Test_QCD_data.root")
)
print "test1"
process.analyzeBasicPat = cms.EDAnalyzer("Triggereffi",
#       photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
#       electronSrc = cms.untracked.InputTag("cleanPatElectrons"),
#        muonSrc = cms.untracked.InputTag("cleanPatMuons"),
        #muonSrc = cms.InputTag("cleanPatMuons"),
        muonSrc = cms.InputTag("slimmedMuons"),
#       tauSrc = cms.untracked.InputTag("cleanPatTaus"),
        jetSrc = cms.InputTag("slimmedJets"),
        metSrc = cms.InputTag("slimmedMETs"),
        genSrc = cms.untracked.InputTag("packedGenParticles"),
        pfSrc = cms.InputTag("packedPFCandidates"),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
        objects = cms.InputTag("slimmedPatTrigger"),
#      objects = cms.InputTag("selectedPatTrigger"),
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



