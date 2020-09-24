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

// IRT
#include "L1Trigger/TrackFindingTrackletHLS/interface/TxtFileWriter.h"
#include <utility>
#include <map>
#include <mutex>


// LUTs from HLS
// Need to store these/access via and ESSource?
// Or provided from wiring map?

static const std::vector< ap_uint<20> > kLinkAssignmentTable = {
 0x500b9,
 0x3000b,
 0x3000d,
 0x5006d,
 0x50082,
 0x500a4,
 0x60843,
 0x8a623,
 0x40045, // 0x20005,
 0x60a62,
 0x40047,
 0x40087,
 0x500b9,
 0x3000b,
 0x3000d,
 0x5006d,
 0x50082,
 0x500a4,
 0x60843,
 0x8a623,
 0x40045, // 0x20005,
 0x60a62,
 0x40047,
 0x40087
};

// LUT with phi corrections to the nominal radius. Only used by layers.
// Values are determined by the radius and the bend of the stub.
static constexpr int kPhiCorrtable_L1[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL1.txt"
;
static constexpr int kPhiCorrtable_L2[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL2.txt"
;
static constexpr int kPhiCorrtable_L3[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL3.txt"
;
static constexpr int kPhiCorrtable_L4[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL4.txt"
;
static constexpr int kPhiCorrtable_L5[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL5.txt"
;
static constexpr int kPhiCorrtable_L6[] = 
#include "emData/MemPrints/Tables/VMPhiCorrL6.txt"
;


//
// class declaration
//

class IRProducer : public edm::stream::EDProducer<> {
   public:
      explicit IRProducer(const edm::ParameterSet&);
      ~IRProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static constexpr int linkWordNBits_ = 20;
      static constexpr int stubWordNBits_ = 38;
      static constexpr int irStubWordNBits_ = 36;
      static constexpr int maxIRMemories = 20; 

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      ap_uint< linkWordNBits_ > getLinkWord( const unsigned int region, const unsigned int channel );

      TTIRMemory::TTStubRefs getTTStubRefsForIRStubs( const TTIRMemory::IRMemories& memories, const TTDTC::Stream& dtcStubs, const std::vector<int> encodedLayers );


      // ----------member data ---------------------------
      // helper class to store DTC configuration
      trackerDTC::Setup setup_;
      // Setup token
      edm::ESGetToken<trackerDTC::Setup, trackerDTC::SetupRcd> esGetToken_;
      // Input token
      edm::EDGetTokenT<TTDTC> tokenDTC_;
      // Product token
      edm::EDPutTokenT<TTIRMemory> putTokenTTIRMemory_;


      const unsigned int linkWord_is2sBit_;
      const unsigned int linkWord_nLayersOffset_;
      const unsigned int linkWord_layerBarrelEndcapOffset_;
      const unsigned int linkWord_layerOrDiskNumberOffset_;
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
  linkWord_nLayersOffset_( iConfig.getParameter<unsigned int>( "linkWord_nLayersOffset" ) ),
  linkWord_layerBarrelEndcapOffset_( iConfig.getParameter<unsigned int>( "linkWord_layerBarrelEndcapOffset" ) ),
  linkWord_layerOrDiskNumberOffset_( iConfig.getParameter<unsigned int>( "linkWord_layerOrDiskNumberOffset" ) ),
  maxNStubsPerDTC_( iConfig.getParameter<unsigned int>( "maxNStubsPerDTC" ) )
{
  // book ES product
  esGetToken_ = esConsumes<trackerDTC::Setup, trackerDTC::SetupRcd, edm::Transition::BeginRun>();
}


IRProducer::~IRProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

ap_uint< IRProducer::linkWordNBits_ > IRProducer::getLinkWord( const unsigned int region, const unsigned int channel ) {

  std::bitset<linkWordNBits_> linkWord{0};

  const int thisDtcId{ setup_.dtcId( region, channel ) };

  const std::vector<int>* thisDtcLayerEncodings{ &setup_.encodingLayerId( channel ) };

  linkWord |= ( thisDtcLayerEncodings->size() << linkWord_nLayersOffset_ );

  bool is2S = !setup_.psModule( thisDtcId );
  linkWord |= ( is2S << ( linkWord_is2sBit_ ) );

  size_t encodedLayer = 0;
  for ( auto decodedLayer : *thisDtcLayerEncodings ) {
    bool isBarrelLayer = decodedLayer < setup_.offsetLayerDisks();
    if ( !isBarrelLayer ) {
      decodedLayer -= setup_.offsetLayerDisks();
    }
    linkWord |= ( ( isBarrelLayer << linkWord_layerBarrelEndcapOffset_ ) ) << ( setup_.hybridNumLayers() * encodedLayer );
    linkWord |= ( decodedLayer << linkWord_layerOrDiskNumberOffset_ ) << ( setup_.hybridNumLayers() * encodedLayer );
    ++encodedLayer;
  }

  ap_uint< linkWordNBits_ > linkWord_ap = linkWord.to_ulong();

  return linkWord_ap;
}

TTIRMemory::TTStubRefs IRProducer::getTTStubRefsForIRStubs( const TTIRMemory::IRMemories& memories, const TTDTC::Stream& dtcStubs,  const std::vector<int> encodedLayers ) {

  TTIRMemory::TTStubRefs allMatchedTTStubs;

  std::vector<unsigned int> layersOfMemories;
  unsigned int decodedLayer{0};
  for ( const auto& layer : encodedLayers ) {
    if ( layer == setup_.offsetLayerId() ) {
      layersOfMemories.insert( layersOfMemories.end(), 8, decodedLayer );
    }
    else {
      layersOfMemories.insert( layersOfMemories.end(), 4, decodedLayer );
    }
    ++decodedLayer;
  }

  // Add exception, check that layersOfMemories has same size as memories
  if ( layersOfMemories.size() > memories.size() ) {
    std::cout << "You might crash soon : " << layersOfMemories.size() << " " << memories.size() << std::endl;
  }

  unsigned int memoryIndex{0};
  for ( const auto& memory : memories ) {

    std::vector< TTStubRef > matchedTTStubs;
    matchedTTStubs.reserve( memory.getEntries(0) );

    const unsigned int encodedLayer{ layersOfMemories[memoryIndex] };
    ++memoryIndex;
    for ( unsigned int iMem = 0; iMem < memory.getEntries(0); ++iMem ) {
      const std::bitset< irStubWordNBits_ > stubWord { memory.read_mem(0, iMem).raw() };


      auto dtcFrame = std::find_if( dtcStubs.begin(), dtcStubs.end(),
                                  [&](const TTDTC::Frame f) {
                                    const std::bitset< TTBV::S > mask{3};
                                    const std::bitset< 2 > dtcEncodedLayer{ ( ( f.second >> 36 ) & mask ).to_ulong()  };
                                    const bool sameEncodedLayer = ( encodedLayer == dtcEncodedLayer.to_ulong() );

                                    const std::bitset< irStubWordNBits_ > dtcWord { ( f.second<< ( TTBV::S - irStubWordNBits_) ).to_string() }; 
                                    return ( sameEncodedLayer && dtcWord == stubWord );
                                  }
                                  );
      if ( dtcFrame != dtcStubs.end() ) {
        matchedTTStubs.emplace_back( dtcFrame->first );
      }
      else {
        matchedTTStubs.emplace_back( TTStubRef() );
      }
    }
    allMatchedTTStubs.emplace_back( matchedTTStubs );
  }
  return allMatchedTTStubs;
}

// ------------ method called to produce the data  ------------
void
IRProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  TTIRMemory product( setup_ );

  edm::Handle< TTDTC > handleDTC;
  iEvent.getByToken< TTDTC >( tokenDTC_, handleDTC );

  // Process stubs in each region and channel within that region
  for ( const int& region : handleDTC->tfpRegions() ) {
    for ( const int& channel : handleDTC->tfpChannels() ) {

      // Get the DTC link ID
      // This is 'reconstructed' from working out the link word
      ap_uint< linkWordNBits_ > linkWord_ap = getLinkWord( region, channel );
      // And then finding the location of the link word in the link assignment table
      auto linkWordInMap = std::find( kLinkAssignmentTable.begin(), kLinkAssignmentTable.end(), linkWord_ap );
      int linkID{-1};
      if ( linkWordInMap != kLinkAssignmentTable.end() ) {
        linkID = linkWordInMap - kLinkAssignmentTable.begin();
      }

      // Container for input stubs to the IR
      std::vector< ap_uint< stubWordNBits_ > > stubsForIR( maxNStubsPerDTC_, 0 );

      // Get the stubs from the DTC
      const TTDTC::Stream& streamFromDTC{ handleDTC->stream( region, channel ) };

      // Make this check into an exception
      if ( streamFromDTC.size() > stubsForIR.size() ) {
        std::cout << "You might crash soon : " << streamFromDTC.size() << " " << stubsForIR.size() << std::endl;
      }

      // Prepare the DTC stubs for the IR
      for ( size_t stubIndex = 0; stubIndex < streamFromDTC.size(); ++stubIndex ) {

        const TTDTC::Frame& stub{ streamFromDTC[ stubIndex ] };

        // Need to modify the stub format to match what the IR expects
        // Should work to remove the need for this and synchronise the stub data formats
        // Set not valid stubs to 0x0000 instead of using valid bit
        if ( stub.first.isNull() ) {
          stubsForIR[ stubIndex ] = ap_uint< stubWordNBits_ >( 0 );
        }
        else {
          // Shorten stub word
          // IR expects 38 bit word, DTC produces 38 bit word and then pads with 0s up to 64 bits
          std::bitset< stubWordNBits_ > stubWordForIR{ ( stub.second<< ( TTBV::S - stubWordNBits_) ).to_string() };         
    
          // Convert to ap_uint and store stub
          stubsForIR[ stubIndex ] = ap_uint< stubWordNBits_ >( stubWordForIR.to_ullong() );
        }

      }

      // Have stubs in an array, and word/index describing where the stubs from this DTC are coming from
      // Now call the HLS IR top level functions
      InputStubMemory<TRACKER> hMemories[maxIRMemories];
      InputRouterTop( linkID,
                      kLinkAssignmentTable.data(),
                      kPhiCorrtable_L1 ,
                      kPhiCorrtable_L2,
                      kPhiCorrtable_L3,
                      kPhiCorrtable_L4,
                      kPhiCorrtable_L5,
                      kPhiCorrtable_L6,
                      stubsForIR.data(),
                      hMemories);

      TTIRMemory::IRMemories outputMemories( std::begin(hMemories), std::end(hMemories) );
      TTIRMemory::TTStubRefs outputTTStubs;

      // Store IR memories with index based on tfp region and channel, rather than ID of input DTC stream
      const unsigned int outputIndex( region * handleDTC->tfpChannels().size() + channel );
      product.setIRMemory( outputIndex, outputMemories );
      product.setTTStubs( outputIndex, getTTStubRefsForIRStubs( outputMemories, streamFromDTC, setup_.encodingLayerId( channel ) ) );
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

  // Write txt files (only once if there are multiple threads).
  const TxtFileWriter txtFileWriter(setup_,
				    //				    linkWordNBits_,
				    linkWord_is2sBit_,
				    linkWord_nLayersOffset_,
				    linkWord_layerBarrelEndcapOffset_,
				    linkWord_layerOrDiskNumberOffset_); 
  static std::once_flag doOnceA, doOnceB;
  std::call_once(doOnceA, &TxtFileWriter::writeDtcPhiRange, txtFileWriter, "dtcPhiRange.txt");
  std::call_once(doOnceB, &TxtFileWriter::writeDtcLinkWords, txtFileWriter, "dtcLinkWords.txt");
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
