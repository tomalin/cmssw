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
// #include "L1Trigger/TrackerDTC/interface/Settings.h"
#include "L1Trigger/TrackerDTC/interface/Setup.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithm_official.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithmRecord.h"

// Output product
#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"
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
      // helper class to store DTC configuration
      trackerDTC::Setup setup_;
      // Setup token
      edm::ESGetToken<trackerDTC::Setup, trackerDTC::SetupRcd> esGetToken_;
      // Input token
      edm::EDGetTokenT<TTDTC> tokenDTC_;
      // Product token
      edm::EDPutTokenT<TTIRMemory> putTokenTTIRMemory_;


      // helper DTC class that stores configuration of DTC
      // trackerDTC::Settings settings_;
      // These sources are only (95% confidence) required by the DTC Settings class
      // An ESSource containing the DTC Settings info is being developed
      // So will no longer need these at some point
      // edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> getTokenTrackerGeometry_;
      // edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> getTokenTrackerTopology_;
      // edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> getTokenCablingMap_;
      // edm::ESGetToken<TTStubAlgorithm<Ref_Phase2TrackerDigi_>, TTStubAlgorithmRecord> getTokenTTStubAlgorithm_;

      // const TrackerGeometry* trackerGeometry_;
      // const TrackerTopology* trackerTopology_;
      // const TrackerDetToDTCELinkCablingMap* ttCablingMap_;
      // edm::ESHandle<TTStubAlgorithm<Ref_Phase2TrackerDigi_>> handleTTStubAlgorithm_;

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
  tokenDTC_( consumes< TTDTC >( edm::InputTag( iConfig.getParameter<edm::InputTag>( "InputTagTTDTC" ) ) ) ),
  putTokenTTIRMemory_( produces< TTIRMemory >( iConfig.getParameter<std::string>( "ProductLabel" ) ) ),
  // settings_(iConfig),
  linkWord_is2sBit_( iConfig.getParameter<unsigned int>( "linkWord_is2sBit" ) ),
  linkWord_hasFirstBarrelLayerBit_( iConfig.getParameter<unsigned int>( "linkWord_hasFirstBarrelLayerBit" ) ),
  linkWord_layerBarrelEndcapOffset_( iConfig.getParameter<unsigned int>( "linkWord_layerBarrelEndcapOffset" ) ),
  maxNStubsPerDTC_( iConfig.getParameter<unsigned int>( "maxNStubsPerDTC" ) )
{
  // book ES product
  esGetToken_ = esConsumes<trackerDTC::Setup, trackerDTC::SetupRcd, edm::Transition::BeginRun>();

  // getTokenTTStubAlgorithm_ = esConsumes<TTStubAlgorithm<Ref_Phase2TrackerDigi_>, TTStubAlgorithmRecord, edm::Transition::BeginRun>( settings_.inputTagTTStubAlgorithm() );
  // getTokenTrackerGeometry_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>( settings_.inputTagTrackerGeometry() );
  // getTokenTrackerTopology_ = esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>( settings_.inputTagTrackerTopology() );
  // getTokenCablingMap_ = esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd, edm::Transition::BeginRun>( settings_.inputTagCablingMap() );
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

  TTIRMemory product( setup_ );

  const std::vector<std::set<int>> layerIdEncodings = setup_.encodingsLayerId();
  edm::Handle< TTDTC > handleDTC;
  iEvent.getByToken< TTDTC >( tokenDTC_, handleDTC );
  for ( const int& region : handleDTC->tfpRegions() ) {
    for ( const int& channel : handleDTC->tfpChannels() ) {

      if ( region != 0 || channel != 0 ) continue;

      const int thisDtcId{ setup_.dtcId( region, channel ) };
      const std::set<int>* thisDtcLayerEncodings{ &layerIdEncodings.at( channel % setup_.numDTCsPerRegion() ) };

      std::bitset<linkWordNBits_> linkWord;

      // bool readsFromFirstBarrelLayer = thisDtcLayerEncodings->at(0) == setup_.offsetLayerId();
      bool readsFromFirstBarrelLayer = ( thisDtcLayerEncodings->find( setup_.offsetLayerId() ) != thisDtcLayerEncodings->end() );
      linkWord |= ( readsFromFirstBarrelLayer << ( linkWord_hasFirstBarrelLayerBit_ - 1 ) );

      bool is2S = !setup_.psModule( thisDtcId );
      linkWord |= ( is2S << ( linkWord_is2sBit_ - 1 ) );

      size_t encodedLayer = 0;
      for ( auto decodedLayer : *thisDtcLayerEncodings ) {
        bool isBarrelLayer = decodedLayer < setup_.offsetLayerDisks();
        if ( !isBarrelLayer ) {
          decodedLayer -= setup_.offsetLayerDisks();
        }
        linkWord |= ( ( isBarrelLayer << linkWord_layerBarrelEndcapOffset_ ) ) << ( setup_.numLayers() * encodedLayer ) ;
        linkWord |= ( decodedLayer ) << ( setup_.numLayers() * encodedLayer ) ;

        ++encodedLayer;
      }

      ap_uint< linkWordNBits_ > linkWord_ap = linkWord.to_ulong();

      // std::array< ap_uint< stubWordNBits_ >, size_t(maxNStubsPerDTC_) > stubsForIR{};
      std::vector< ap_uint< stubWordNBits_ > > stubsForIR( maxNStubsPerDTC_ );
      // std::cout << "================================" << std::endl;
      // std::cout << "Region, channel : " << region << " " << channel << std::endl;
      // std::cout << "DTC word : " << linkWord << std::endl;
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
        // // MODIFIED DTC stub to reverse word
        // // As should only reverse order of variables, not just reverse all bits
        // // Reverse order of word for IR
        // // Current implementation of IR assumes different endianess to DTC (DTC follows twiki recommendations)
        // auto stubWordForIR_string = stubWordForIR.to_string();
        // std::reverse(stubWordForIR_string.begin(), stubWordForIR_string.end());
        // std::bitset< stubWordNBits_ > stubWordForIR_reversed{ stubWordForIR_string };

        // Convert to ap_uint and store stub
        // stubsForIR[ stubIndex ] = ap_uint< stubWordNBits_ >( stubWordForIR_reversed.to_ullong() );
        stubsForIR[ stubIndex ] = ap_uint< stubWordNBits_ >( stubWordForIR.to_ullong() );


      }

      // for ( const auto& stub : stubsForIR ) {
      //   std::cout << std::bitset<stubWordNBits_>( stub ) << std::endl;
      // }

      // Have stubs in an array, and word describing where the stubs from this DTC are coming from
      // Now just call the HLS IR top level function
      // stubsForIR.data() should give pointer to first element of array of vector

      StubsBarrelPS hBarrelPS;
      StubsDiskPS hDiskPS;

      StubsBarrel2S hBarrel2S;
      StubsDisk2S hDisk2S;

      // Simulate IR here for now
      // Remove non-valid stubs
      // Remove valid bit and layerID from stubs words
      // Calculate coarse phi region
      // Format into "IputStub" format (output form IR)
      for ( const auto& stub : stubsForIR ) {

        std::bitset< stubWordNBits_ > stub_bv( stub );
        // std::cout << stub_bv << std::endl;

        // Remove stubs that aren't valid, and remove valid bit
        bool isStubValid = stub_bv.test( stubWordNBits_ - 1 );
        // std::cout << "Is stub valid : " << isStubValid << std::endl;
        if ( !isStubValid ) continue;

        // Get layer bits
        auto stub_layerID = ( stub_bv>>( 36 ) &= std::bitset<stubWordNBits_>(0x3) );
        // std::bitset< stubWordNBits_ > stubWordForIR{ ( stub.second<< ( TTBV::S - stubWordNBits_) ).to_string() };
        // std::cout << "Layer bits : " << stub_bv.test(37) << " " << stub_bv.test(36) << " " << stub_layerID << std::endl;
        // bool isBarrelLayer = thisDtcLayerEncodings->at( stub_layerID.to_ulong() ) < setup_.offsetLayerDisks();
        std::set<int>::iterator it = thisDtcLayerEncodings->begin();
        std::advance(it, stub_layerID.to_ulong());
        bool isBarrelLayer = ( *it < setup_.offsetLayerDisks() );

        // Get phi bits
        int nBitsPhi{ 14 };
        int lsbPhi{ 20 };
        int nBitsCoarsePhi{ 2 };
        if ( readsFromFirstBarrelLayer && ( stub_layerID == 0 ) ) {
          nBitsCoarsePhi = 3;
        }
        if ( is2S ) {
          if ( isBarrelLayer ) {
            nBitsPhi = 17;
            lsbPhi = 16;
          }
          else {
            nBitsPhi = 14;
            lsbPhi = 15;           
          }
        }
        else if ( !isBarrelLayer ) {
            nBitsPhi = 14;
            lsbPhi = 20;          
        }

        // auto stub_phi = ( stub_bv>>( lsbPhi - 1 ) &= std::bitset<stubWordNBits_>( std::pow( 2, nBitsPhi ) - 1) );
        // std::bitset< stubWordNBits_ > stubWordForIR{ ( stub.second<< ( TTBV::S - stubWordNBits_) ).to_string() };
        // std::cout << "Phi bits : " << stub_phi << std::endl;
        // std::cout << "Info : " << readsFromFirstBarrelLayer << " " << isBarrelLayer << " " << is2S << std::endl;
        auto phiRegion = ( stub_bv>>( lsbPhi - 1 + nBitsPhi - nBitsCoarsePhi ) &= std::bitset<stubWordNBits_>( std::pow( 2, nBitsCoarsePhi ) - 1) );
        // std::cout << "Coarse phi bits : " << phiRegion << std::endl;

        std::bitset< TTIRMemory::stubWordIRNBits_ > stubWordFromIR{ ( stub_bv &= std::bitset<stubWordNBits_>().set() ).to_string() };
        // std::cout << "Stub word from IR : " << stubWordFromIR << std::endl;

        if ( isBarrelLayer ) {
          if ( is2S ) {
            if ( stub_layerID.to_ulong() == 0 ) {
              hBarrel2S.m1[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrel2S.m1[phiRegion.to_ulong()].getEntries(0) );
            }
            else if ( stub_layerID.to_ulong() == 1 ) {
              hBarrel2S.m2[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrel2S.m2[phiRegion.to_ulong()].getEntries(0) );             
            }
            else {
              hBarrel2S.m3[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrel2S.m3[phiRegion.to_ulong()].getEntries(0) );
            }
          }
          else {
            if ( stub_layerID.to_ulong() == 0 ) {
              hBarrelPS.m1[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrelPS.m1[phiRegion.to_ulong()].getEntries(0) );
            }
            else if ( stub_layerID.to_ulong() == 1 ) {
              hBarrelPS.m2[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrelPS.m2[phiRegion.to_ulong()].getEntries(0) );             
            }
            else {
              hBarrelPS.m3[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hBarrelPS.m3[phiRegion.to_ulong()].getEntries(0) );
            }
          }
        }
        else { 
          if ( is2S ) {
            // hDisk2S;
            if ( stub_layerID.to_ulong() == 0 ) {
              hDisk2S.m1[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDisk2S.m1[phiRegion.to_ulong()].getEntries(0) );
            }
            else if ( stub_layerID.to_ulong() == 1 ) {
              hDisk2S.m2[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDisk2S.m2[phiRegion.to_ulong()].getEntries(0) );             
            }
            else if ( stub_layerID.to_ulong() == 2 ) {
              hDisk2S.m3[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDisk2S.m3[phiRegion.to_ulong()].getEntries(0) );
            }     
            else if ( stub_layerID.to_ulong() == 3 ) {
              hDisk2S.m4[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDisk2S.m4[phiRegion.to_ulong()].getEntries(0) );
            }     
            else {
              hDisk2S.m5[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDisk2S.m5[phiRegion.to_ulong()].getEntries(0) );
            }     

          }
          else {
            // hDiskPS;
            if ( stub_layerID.to_ulong() == 0 ) {
              hDiskPS.m1[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDiskPS.m1[phiRegion.to_ulong()].getEntries(0) );
            }
            else if ( stub_layerID.to_ulong() == 1 ) {
              hDiskPS.m2[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDiskPS.m2[phiRegion.to_ulong()].getEntries(0) );             
            }
            else if ( stub_layerID.to_ulong() == 2 ) {
              hDiskPS.m3[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDiskPS.m3[phiRegion.to_ulong()].getEntries(0) );             
            }
            else if ( stub_layerID.to_ulong() == 3 ) {
              hDiskPS.m4[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDiskPS.m4[phiRegion.to_ulong()].getEntries(0) );             
            }
            else {
              hDiskPS.m5[phiRegion.to_ulong()].write_mem( 0, ap_uint< TTIRMemory::stubWordIRNBits_ >( stubWordFromIR.to_ullong() ), hDiskPS.m5[phiRegion.to_ulong()].getEntries(0) );             
            }

          }

        }


      }

      // if ( !is2S && stubsForIR.size() > 0 ) {
      //   std::cout << "Writing this stub : " << std::bitset<stubWordNBits_>(stubsForIR[ 0 ]) << std::endl;
      //   hBarrelPS.m1[0].write_mem(0, InputStub<BARRELPS>( stubsForIR[ 0 ] ), 0 );
      // }

      product.setIRMemory( thisDtcId, hBarrelPS );
      // for( size_t cRegion=0; cRegion < 8 ; cRegion++) {
      //   std::cout << "Number of entries : " << cRegion << " : " << hBarrelPS.m1[cRegion].getEntries(0) << std::endl;
      //   for ( size_t iStub = 0; iStub < hBarrelPS.m1[cRegion].getEntries(0); ++iStub ) {
      //     std::cout << "Word : " << std::bitset< TTIRMemory::stubWordIRNBits_> ( hBarrelPS.m1[cRegion].read_mem(0,iStub).raw() ) << std::endl;
      //   }
      // }

    }
  }

  iEvent.emplace( putTokenTTIRMemory_, std::move( product ) );
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
  setup_ = iSetup.getData(esGetToken_);

  // trackerGeometry_ = &iSetup.getData( getTokenTrackerGeometry_ );
  // trackerTopology_ = &iSetup.getData( getTokenTrackerTopology_ );
  // ttCablingMap_ = &iSetup.getData( getTokenCablingMap_ );

  // settings_.setTrackerGeometry( &iSetup.getData( getTokenTrackerGeometry_ ) );
  // settings_.setTrackerTopology( &iSetup.getData( getTokenTrackerTopology_ ) );
  // settings_.setCablingMap( &iSetup.getData( getTokenCablingMap_ ) );
  // settings_.setTTStubAlgorithm( iSetup.getHandle( getTokenTTStubAlgorithm_ ) );
  // settings_.beginRun();
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
