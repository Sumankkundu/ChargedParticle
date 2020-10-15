import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B6F8B6-90F1-E311-B72C-0025905A6092.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FEDF66BB-75F5-6746-A6E0-CE47C31F7F6C.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FE41CF02-0D7E-914E-B610-BE747691D499.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FE3755A8-7960-BB48-8B2C-31478A80E7C5.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/280000/00F73378-FCD0-6847-B007-C90878390A8A.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/20000/086E8B6F-28D0-AE47-A07B-AEF644006C7D.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/20000/3D81A6E0-C31E-254B-B49D-C264BA44C71A.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/20000/188FCC60-4351-0F4A-9FE0-9D4AE90546D5.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FFDCB325-AA00-ED43-A45A-EE4704788533.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FFD92D26-BC76-D54E-916D-3326443A581E.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/0879FBF1-5326-F94C-B12A-8E88EEBF639C.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/0859D23C-3B96-F94D-85C5-2A277F6704FE.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/084EEF62-871E-1E49-A9FA-EEC55E72FE35.root',
#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/FB6B1B36-AF31-BB4F-8B24-37063D7D7CD9.root',
'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/FA0DBDDB-4DC6-5A44-934B-B3100EB43C5E.root',
'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/F199462C-AD3F-CD41-A707-168EBF5BA2B9.root',
'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/E98B458E-FF8B-1F41-993E-05D6AC948556.root',
#'/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/7E86CBA6-9842-E811-B26F-0CC47AA53D6E.root',
#'/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/7E65A585-9C42-E811-BA22-0025901ABB72.root',
 )
#eventsToSkip = cms.untracked.VEventRange('1:1950-1:2000'),
#eventsToSkip = cms.untracked.EventRange('1:1950-1:2000'),
#eventRanges = cms.untracked.VEventRange('1:1000-1:2000'),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_mc2017_realistic_v17')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_mc2017_realistic_v14')
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mc2017_realistic_v7')
#process.GlobalTag = GlobalTag(process.GlobalTag,'76X_mcRun2_asymptotic_RunIIFall15DR76_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'76X_mcRun2_asymptotic_v12')
#process.GlobalTag = GlobalTag(process.GlobalTag,'74X_mcRun2_asymptotic_realisticBS_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'MCRUN2_74_V9')
from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("HLTrigger.HLTcore.hltPrescaleRecorder_cfi")


process.options = cms.untracked.PSet(

)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)



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
    fileName=cms.string("Test_MC_QCD.root")
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



