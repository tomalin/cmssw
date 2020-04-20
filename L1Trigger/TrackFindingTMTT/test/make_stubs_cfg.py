import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("STUBS")

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

options = VarParsing.VarParsing ('analysis')

def getTxtFile(txtFileName): 
  return os.environ['CMSSW_BASE']+'/src/L1Trigger/TrackFindingTMTT/data/'+txtFileName

#--- Specify input MC
inputMCtxt = getTxtFile('MCsamples/1110/RelVal/TTbar/localRAL/PU200.txt')
options.register('inputMC', inputMCtxt, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")

#--- Specify number of events to process.
options.register('Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"Number of Events to analyze")

options.parseArguments()

#--- input and output
list = FileUtils.loadListFromFile(options.inputMC)
readFiles = cms.untracked.vstring(*list)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            )

process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('output_PU200.root'), ## ADAPT IT ##
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)
process.RAWSIMoutput.outputCommands.append('keep  *_*_ClusterAccepted_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_StubAccepted_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)


process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")
process.p = cms.Path(process.TrackTriggerClustersStubs * process.TrackTriggerAssociatorClustersStubs)

