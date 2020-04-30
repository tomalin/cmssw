import FWCore.ParameterSet.Config as cms

IRProducer_params = cms.PSet (
	linkWord_is2sBit = cms.uint32( 17 ),
	linkWord_hasFirstBarrelLayerBit = cms.uint32( 18 ),
	linkWord_layerBarrelEndcapOffset = cms.uint32( 3 ),
	maxNStubsPerDTC = cms.uint32( 256 )

)