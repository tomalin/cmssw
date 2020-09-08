import FWCore.ParameterSet.Config as cms

IRProducer_params = cms.PSet (
	InputTagTTDTC	       = cms.InputTag  ( "TrackerDTCProducer", "StubAccepted" ),   # 
	ProductLabel	       = cms.string    ( "IRMemories" ),                                    #
	linkWord_is2sBit = cms.uint32( 16 ),
	linkWord_nLayersOffset = cms.uint32( 17 ),
	linkWord_layerBarrelEndcapOffset = cms.uint32( 0 ),
	linkWord_layerOrDiskNumberOffset = cms.uint32( 1 ),
	maxNStubsPerDTC = cms.uint32( 256 )

)