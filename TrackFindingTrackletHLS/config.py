import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process( "DTCIR" )

process.load( 'Configuration.Geometry.GeometryExtended2026D49Reco_cff' )
process.load( 'Configuration.Geometry.GeometryExtended2026D49_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_cff' )
process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
# process.load( 'Configuration.StandardSequences.L1TrackTrigger_cff' )
process.load( 'L1Trigger.TrackTrigger.TrackTrigger_cff' )
process.load( 'SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff' )
process.load( 'L1Trigger.TrackerDTC.ProducerED_cff' )
process.load( 'L1Trigger.TrackerDTC.ProducerES_cff' )
process.L1TrackTrigger=cms.Sequence(process.TrackTriggerClustersStubs*process.TrackTriggerAssociatorClustersStubs)
process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(True)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag( process.GlobalTag, 'auto:phase2_realistic', '' )


process.load( "FWCore.MessageLogger.MessageLogger_cfi" )

options = VarParsing.VarParsing( 'analysis' )

#--- Specify input MC
options.register( 'inputMC', 
  '/store/relval/CMSSW_11_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/A65BAFE2-1EB8-7448-899C-B96FAB2781A7.root',
  VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed"
)

#--- Specify number of events to process.
options.register( 'Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Number of Events to analyze" )

#--- Specify whether to output a GEN-SIM-DIGI-RAW dataset containing the DTCStub Collections.
options.register( 'outputDataset' ,0 ,VarParsing.VarParsing.multiplicity.singleton ,VarParsing.VarParsing.varType.int, "Create GEN-SIM-DIGI-RAW dataset containing DTCStub Collections" )

options.parseArguments()

#--- input and output

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring( options.inputMC ),
  secondaryFileNames = cms.untracked.vstring(),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

process.TFileService = cms.Service( "TFileService", fileName = cms.string( "Hist.root" ) )

process.Timing = cms.Service( "Timing", summaryOnly = cms.untracked.bool( True ) )

#--- Load code that produces DTCStubs
# load Track Trigger Configuration
process.load( 'L1Trigger.TrackerDTC.ProducerES_cff' )
# load code that produces DTCStubs
process.load( 'L1Trigger.TrackerDTC.ProducerED_cff' )
# load code that analyzes DTCStubs
process.load( 'L1Trigger.TrackerDTC.Analyzer_cff' )
process.dtc = cms.Path( process.TrackerDTCProducer )#* process.TrackerDTCAnalyzer )

#--- Load code that runs IR
process.load( 'L1Trigger.TrackFindingTrackletHLS.IR_cff' )

process.ir = cms.Path( process.IRProducer * process.IRAnalyzer )

# Write output dataset
process.load( 'Configuration.EventContent.EventContent_cff' )

process.writeDataset = cms.OutputModule(
"PoolOutputModule",
splitLevel = cms.untracked.int32( 0 ),
eventAutoFlushCompressedSize = cms.untracked.int32( 5242880 ),
outputCommands = process.RAWSIMEventContent.outputCommands,
fileName = cms.untracked.string( 'output_dataset.root' ), ## ADAPT IT ##
dataset  = cms.untracked.PSet(
filterName = cms.untracked.string(''),
dataTier   = cms.untracked.string( 'GEN-SIM' )
)
)
process.writeDataset.outputCommands.append( 'keep  *_*_*_*' )
process.writeDataset.outputCommands.append( 'drop  *_TrackerDTCProducer_StubAccepted_*' )
process.pd = cms.EndPath( process.writeDataset )

process.schedule = cms.Schedule( process.dtc, process.ir )