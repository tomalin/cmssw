// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "CondFormats/SiPhase2TrackerObjects/interface/TrackerDetToDTCELinkCablingMap.h"
#include "CondFormats/DataRecord/interface/TrackerDetToDTCELinkCablingMapRcd.h"
#include "DataFormats/L1TrackTrigger/interface/TTDTC.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/TrackerDTC/interface/Settings.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithm_official.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithmRecord.h"

// HLS arbitrary precision types
#include "ap_fixed.h"

// Track finding HLS modules
// #include "Constants.hh"
#include "InputRouterTop.h"

#include <array>
#include <iostream>

//
// class declaration
//

class IRProducer : public edm::stream::EDProducer<> {
   public:
      explicit IRProducer(const edm::ParameterSet&);
      ~IRProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static constexpr int linkWordNBits_ = 18;
      static constexpr int stubWordNBits_ = 39;

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT <TTDTC> tokenDTC_;
      // helper DTC class that stores configuration of DTC
      trackerDTC::Settings settings_;
      // These sources are only (95% confidence) required by the DTC Settings class
      // An ESSource containing the DTC Settings info is being developed
      // So will no longer need these at some point
      edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> getTokenTrackerGeometry_;
      edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> getTokenTrackerTopology_;
      edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> getTokenCablingMap_;
      edm::ESGetToken<TTStubAlgorithm<Ref_Phase2TrackerDigi_>, TTStubAlgorithmRecord> getTokenTTStubAlgorithm_;

      const TrackerGeometry* trackerGeometry_;
      const TrackerTopology* trackerTopology_;
      const TrackerDetToDTCELinkCablingMap* ttCablingMap_;
      edm::ESHandle<TTStubAlgorithm<Ref_Phase2TrackerDigi_>> handleTTStubAlgorithm_;

      const unsigned int linkWord_is2sBit_;
      const unsigned int linkWord_hasFirstBarrelLayerBit_;
      const unsigned int linkWord_layerBarrelEndcapOffset_;
      const unsigned int maxNStubsPerDTC_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
IRProducer::IRProducer(const edm::ParameterSet& iConfig) :
  tokenDTC_( consumes< TTDTC >( edm::InputTag( "TrackerDTCProducer", "StubAccepted" ) ) ),
  settings_(iConfig),
  linkWord_is2sBit_( iConfig.getParameter<unsigned int>( "linkWord_is2sBit" ) ),
  linkWord_hasFirstBarrelLayerBit_( iConfig.getParameter<unsigned int>( "linkWord_hasFirstBarrelLayerBit" ) ),
  linkWord_layerBarrelEndcapOffset_( iConfig.getParameter<unsigned int>( "linkWord_layerBarrelEndcapOffset" ) ),
  maxNStubsPerDTC_( iConfig.getParameter<unsigned int>( "maxNStubsPerDTC" ) )
{
  getTokenTTStubAlgorithm_ = esConsumes<TTStubAlgorithm<Ref_Phase2TrackerDigi_>, TTStubAlgorithmRecord, edm::Transition::BeginRun>( settings_.inputTagTTStubAlgorithm() );
  getTokenTrackerGeometry_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>( settings_.inputTagTrackerGeometry() );
  getTokenTrackerTopology_ = esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>( settings_.inputTagTrackerTopology() );
  getTokenCablingMap_ = esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd, edm::Transition::BeginRun>( settings_.inputTagCablingMap() );
}


IRProducer::~IRProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
IRProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const std::vector<std::vector<int>> layerIdEncodings = settings_.hybrid()->layerIdEncodings();
  edm::Handle< TTDTC > handleDTC;
  iEvent.getByToken< TTDTC >( tokenDTC_, handleDTC );
  for ( const int& region : handleDTC->tfpRegions() ) {
    for ( const int& channel : handleDTC->tfpChannels() ) {

      const int thisDtcId{ handleDTC->dtcId( region, channel ) };
      const std::vector<int>* thisDtcLayerEncodings{ &layerIdEncodings.at( channel % settings_.numDTCsPerRegion() ) };

      std::bitset<linkWordNBits_> linkWord;

      bool readsFromFirstBarrelLayer = thisDtcLayerEncodings->at(0) == settings_.offsetLayerId();
      linkWord |= ( readsFromFirstBarrelLayer << ( linkWord_hasFirstBarrelLayerBit_ - 1 ) );

      bool is2S = !handleDTC->psModlue( thisDtcId );
      linkWord |= ( is2S << ( linkWord_is2sBit_ - 1 ) );

      for ( size_t encodedLayer = 0; encodedLayer < thisDtcLayerEncodings->size(); ++encodedLayer ) {
        int decodedLayer = thisDtcLayerEncodings->at( encodedLayer );
        bool isBarrelLayer = decodedLayer < settings_.offsetLayerDisks();
        if ( !isBarrelLayer ) {
          decodedLayer -= settings_.offsetLayerDisks();
        }
        linkWord |= ( ( isBarrelLayer << linkWord_layerBarrelEndcapOffset_ ) ) << ( settings_.numLayers() * encodedLayer ) ;
        linkWord |= ( decodedLayer ) << ( settings_.numLayers() * encodedLayer ) ;
      }

      ap_uint< linkWordNBits_ > linkWord_ap = linkWord.to_ulong();

      // std::array< ap_uint< stubWordNBits_ >, size_t(maxNStubsPerDTC_) > stubsForIR{};
      std::vector< ap_uint< stubWordNBits_ > > stubsForIR( maxNStubsPerDTC_ );
      std::cout << "================================" << std::endl;
      std::cout << "DTC word : " << linkWord << std::endl;
      const TTDTC::Stream& streamFromDTC{ handleDTC->stream( region, channel ) };

      if ( streamFromDTC.size() > stubsForIR.size() ) {
        // Make into exception
        std::cout << "You might crash soon : " << streamFromDTC.size() << " " << stubsForIR.size() << std::endl;
      }

      for ( size_t stubIndex = 0; stubIndex < streamFromDTC.size(); ++stubIndex ) {
      // for ( const TTDTC::Frame& stub : handleDTC->stream( region, channel ) ) {
        const TTDTC::Frame& stub{ streamFromDTC[ stubIndex ] };
        // Shorten stub word
        // IR expects 39 bit word, DTC produces 39 bit word and then pads with 0s up to 64 bits
        std::bitset< stubWordNBits_ > stubWordForIR{ ( stub.second<< ( TTBV::S - stubWordNBits_) ).to_string() };
        // Reverse order of word for IR
        // Current implementation of IR assumes different endianess to DTC (DTC follows twiki recommendations)
        auto stubWordForIR_string = stubWordForIR.to_string();
        std::reverse(stubWordForIR_string.begin(), stubWordForIR_string.end());
        std::bitset< stubWordNBits_ > stubWordForIR_reversed{ stubWordForIR_string };
        // Convert to ap_uint and store stub
        stubsForIR[ stubIndex ] = ap_uint< stubWordNBits_ >( stubWordForIR_reversed.to_ullong() );


      }

      for ( const auto& stub : stubsForIR ) {
        std::cout << std::bitset<stubWordNBits_>( stub ) << std::endl;
      }
      
      // Have stubs in an array, and word describing where the stubs from this DTC are coming from
      // Now just call the HLS IR top level function
      // stubsForIR.data() should give pointer to first element of array of vector
    }
  }
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
IRProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
IRProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
IRProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  trackerGeometry_ = &iSetup.getData( getTokenTrackerGeometry_ );
  trackerTopology_ = &iSetup.getData( getTokenTrackerTopology_ );
  ttCablingMap_ = &iSetup.getData( getTokenCablingMap_ );

  settings_.setTrackerGeometry( &iSetup.getData( getTokenTrackerGeometry_ ) );
  settings_.setTrackerTopology( &iSetup.getData( getTokenTrackerTopology_ ) );
  settings_.setCablingMap( &iSetup.getData( getTokenCablingMap_ ) );
  settings_.setTTStubAlgorithm( iSetup.getHandle( getTokenTTStubAlgorithm_ ) );
  settings_.beginRun();
}
 
// ------------ method called when ending the processing of a run  ------------
/*
void
IRProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
IRProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
IRProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IRProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IRProducer);
