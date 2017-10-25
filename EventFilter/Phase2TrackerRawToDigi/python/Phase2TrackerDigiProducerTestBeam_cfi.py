import FWCore.ParameterSet.Config as cms

Phase2TrackerDigiProducerTestBeam = cms.EDProducer(
    'Phase2TrackerDigiProducerTestBeam',
    ProductLabel = cms.InputTag("rawDataCollector")
)
