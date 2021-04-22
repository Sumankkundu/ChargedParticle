import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2017D/JetHT/MINIAOD/09Aug2019_UL2017-v1/00000/03A89500-ED1F-3B4F-AA84-33D55D96D9B9.root',
'/store/data/Run2017D/JetHT/MINIAOD/09Aug2019_UL2017-v1/00000/4EAE67FC-658E-1647-B91F-A978348DA634.root',
'/store/data/Run2017D/JetHT/MINIAOD/09Aug2019_UL2017-v1/00000/56454910-080C-714B-9B0B-D76B3AD61185.root',
'/store/data/Run2017D/JetHT/MINIAOD/09Aug2019_UL2017-v1/00000/9A2F7287-D3C7-EB42-B8CC-D4F5F536BFC4.root',
'/store/data/Run2017E/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/016BE62D-1105-AC4F-8A58-59BD14326D8B.root',
'/store/data/Run2017E/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/01BB9E36-70E0-D64A-8164-87AEA03925B2.root',
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7E46D250-67B0-E511-BB96-0025905C3E66.root'
#'/store/data/Run2017D/JetHT/MINIAOD/31Mar2018-v1/90000/E4402B21-F739-E811-BA23-0CC47A4C8E1E.root',
#'/store/data/Run2017B/JetHT/MINIAOD/17Nov2017-v1/60000/DA587060-EED8-E711-80CC-0CC47A4D7678.root',
 )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
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
#process.GlobalTag = GlobalTag(process.GlobalTag,'102X_dataRun2_v12')
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v33')
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



