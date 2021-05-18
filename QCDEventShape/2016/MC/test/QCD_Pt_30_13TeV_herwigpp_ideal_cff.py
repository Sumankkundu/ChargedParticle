# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/QCD_Pt_30_13TeV_herwigpp_cff.py -s GEN --conditions auto:mc --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 10 --no_exec --python_filename QCD_Pt_30_13TeV_herwigpp_ideal_cff.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/Generator/python/QCD_Pt_30_13TeV_herwigpp_cff.py nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('QCD_Pt_30_13TeV_herwigpp_cff_py_GEN.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("ThePEGGeneratorFilter",
    configFiles = cms.vstring(),
    crossSection = cms.untracked.double(-1.0),
    dataLocation = cms.string('${HERWIGPATH}'),
    eventHandlers = cms.string('/Herwig/EventHandlers'),
    filterEfficiency = cms.untracked.double(1.0),
    generatorModule = cms.string('/Herwig/Generators/LHCGenerator'),
    hwpp_basicSetup = cms.vstring('create ThePEG::RandomEngineGlue /Herwig/RandomGlue', 
        'set /Herwig/Generators/LHCGenerator:RandomNumberGenerator /Herwig/RandomGlue', 
        'set /Herwig/Generators/LHCGenerator:NumberOfEvents 10000000', 
        'set /Herwig/Generators/LHCGenerator:DebugLevel 1', 
        'set /Herwig/Generators/LHCGenerator:UseStdout 0', 
        'set /Herwig/Generators/LHCGenerator:PrintEvent 0', 
        'set /Herwig/Generators/LHCGenerator:MaxErrors 10000'),
    hwpp_cm_13TeV = cms.vstring('set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 13000.0', 
        'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV'),
    hwpp_cmsDefaults = cms.vstring('+hwpp_basicSetup', 
        '+hwpp_setParticlesStableForDetector'),
    hwpp_pdf_CTEQ6L1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Common', 
        '+hwpp_ue_EE5C'),
    hwpp_pdf_CTEQ6L1_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Common', 
        '+hwpp_ue_CUETHS1'),
    hwpp_pdf_CTEQ6L1_Common = cms.vstring('create ThePEG::LHAPDF /Herwig/Partons/cmsPDFSet ThePEGLHAPDF.so', 
        'set /Herwig/Partons/cmsPDFSet:PDFName cteq6ll.LHpdf', 
        'set /Herwig/Partons/cmsPDFSet:RemnantHandler /Herwig/Partons/HadronRemnants', 
        'set /Herwig/Particles/p+:PDF /Herwig/Partons/cmsPDFSet', 
        'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/cmsPDFSet'),
    hwpp_pdf_CTEQ6L1_Hard = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_Common', 
        '+hwpp_ue_EE5C'),
    hwpp_pdf_CTEQ6L1_Hard_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_Common', 
        '+hwpp_ue_CUETHS1'),
    hwpp_pdf_CTEQ6L1_Hard_Common = cms.vstring('create ThePEG::LHAPDF /Herwig/Partons/cmsHardPDFSet ThePEGLHAPDF.so', 
        'set /Herwig/Partons/cmsHardPDFSet:PDFName cteq6ll.LHpdf', 
        'set /Herwig/Partons/cmsHardPDFSet:RemnantHandler /Herwig/Partons/HadronRemnants'),
    hwpp_pdf_CTEQ6LL = cms.vstring('+hwpp_pdf_CTEQ6L1'),
    hwpp_pdf_CTEQ6LL_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_CUETHS1'),
    hwpp_pdf_CTEQ6LL_Hard = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard'),
    hwpp_pdf_CTEQ6LL_Hard_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_CUETHS1'),
    hwpp_setParticlesStableForDetector = cms.vstring('set /Herwig/Particles/mu-:Stable Stable', 
        'set /Herwig/Particles/mu+:Stable Stable', 
        'set /Herwig/Particles/Sigma-:Stable Stable', 
        'set /Herwig/Particles/Sigmabar+:Stable Stable', 
        'set /Herwig/Particles/Lambda0:Stable Stable', 
        'set /Herwig/Particles/Lambdabar0:Stable Stable', 
        'set /Herwig/Particles/Sigma+:Stable Stable', 
        'set /Herwig/Particles/Sigmabar-:Stable Stable', 
        'set /Herwig/Particles/Xi-:Stable Stable', 
        'set /Herwig/Particles/Xibar+:Stable Stable', 
        'set /Herwig/Particles/Xi0:Stable Stable', 
        'set /Herwig/Particles/Xibar0:Stable Stable', 
        'set /Herwig/Particles/Omega-:Stable Stable', 
        'set /Herwig/Particles/Omegabar+:Stable Stable', 
        'set /Herwig/Particles/pi+:Stable Stable', 
        'set /Herwig/Particles/pi-:Stable Stable', 
        'set /Herwig/Particles/K+:Stable Stable', 
        'set /Herwig/Particles/K-:Stable Stable', 
        'set /Herwig/Particles/K_S0:Stable Stable', 
        'set /Herwig/Particles/K_L0:Stable Stable'),
    hwpp_ue_EE5C = cms.vstring('+hwpp_ue_EE5CEnergyExtrapol', 
        'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes', 
        'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49', 
        'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80', 
        'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30', 
        'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes', 
        'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes', 
        'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2'),
    hwpp_ue_EE5CEnergyExtrapol = cms.vstring('set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power', 
        'set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV', 
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33', 
        'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV'),
    parameterSets = cms.vstring('hwpp_cmsDefaults', 
        'hwpp_cm_13TeV', 
        'hwpp_ue_EE5C', 
        'hwpp_pdf_CTEQ6L1', 
        'processParameters'),
    processParameters = cms.vstring('insert /Herwig/MatrixElements/SimpleQCD:MatrixElements[0] /Herwig/MatrixElements/MEQCD2to2', 
        'set /Herwig/Cuts/JetKtCut:MinKT 30*GeV', 
        'set /Herwig/Cuts/QCDCuts:MHatMin 0.0*GeV', 
        'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0'),
    repository = cms.string('HerwigDefaults.rpo'),
    run = cms.string('LHC')
)
process.qcdevt = cms.EDAnalyzer("QCDEventShape",
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
        Trigger = cms.untracked.bool(False),
        GenJET =  cms.untracked.bool(True),
        HistFill = cms.untracked.bool(True),
        MonteCarlo =  cms.untracked.bool(True),
        ParticleLabel =  cms.untracked.bool(False),
        Reconstruct =cms.untracked.bool(False),
        EtaRange =  cms.untracked.double(3.0),
        PtThreshold = cms.untracked.double(55.0), #effective is 21
        LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81       


 )

#process.TFileService=cms.Service("TFileService",
#    fileName=cms.string("Test_MC_QCD.root")
#)

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
process.p = cms.Path(process.qcdevt)
# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


