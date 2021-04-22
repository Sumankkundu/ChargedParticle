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
#'/store/data/Run2017F/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FD6F6891-A1F9-534E-9B4C-6D49F8620D0B.root',
'/store/data/Run2016B/JetHT/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/10000/1DEECE7D-4F11-7943-BD9F-B41CBF5CA860.root',
'/store/data/Run2016B/JetHT/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/10000/27E0797F-88D8-D544-9887-B79DC448A686.root',
'/store/data/Run2016B/JetHT/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/10000/2EC520CC-F255-BA44-AC5A-7DEFC8DB9052.root',
'/store/data/Run2016B/JetHT/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/10000/30A082C5-4362-4C46-B42E-1032707D113C.root',
'/store/data/Run2016B/JetHT/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/10000/31FC85C7-80C7-354A-B6BF-556018F02160.root',

#'/store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/006D4A04-6DF4-7941-B54F-12CEF66AEE20.root',
#'/store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/00DC6D7F-C300-504F-A3BC-497C7D8FDA98.root',
#'/store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/0101796A-F847-3340-8AEA-0617AA884CA5.root',
#'/store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/010EF86B-7A36-EF4D-B96E-EF1D390FA68A.root',
#'/store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/0141D902-DDEF-C240-8F77-A770EBB7BE86.root',
)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v11')
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v32')
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
    fileName=cms.string("Test_Data_UL16.root")
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



