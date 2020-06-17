import FWCore.ParameterSet.Config as cms

#---------------------------------------------------------------------------------------------------------
# This describes the IR Stub processing
#---------------------------------------------------------------------------------------------------------

from L1Trigger.TrackFindingTrackletHLS.IR_Defaults_cfi import IRProducer_params

# from L1Trigger.TrackerDTC.Producer_Defaults_cfi import TrackerDTCProducer_params
# from L1Trigger.TrackerDTC.Format_Hybrid_cfi import TrackerDTCFormat_params
# from L1Trigger.TrackerDTC.Analyzer_Defaults_cfi import TrackerDTCAnalyzer_params

IRProducer = cms.EDProducer('IRProducer',
			IRProducer_params,
			# TrackerDTCProducer_params,
			# TrackerDTCFormat_params,
			# TrackerDTCAnalyzer_params
	)

IRAnalyzer = cms.EDAnalyzer('IRAnalyzer'
	)