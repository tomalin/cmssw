#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/GenericHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"

#include "L1Trigger/TrackerDTC/interface/Setup.h"

#include <TProfile.h>
#include <TH2F.h>


using namespace std;
using namespace edm;

class IRAnalyzer : public one::EDAnalyzer<one::WatchRuns, one::SharedResources> {
public:
  IRAnalyzer(const ParameterSet& iConfig);
  void beginJob() override {}
  void beginRun(const Run& iEvent, const EventSetup& iSetup) override;
  void analyze(const Event& iEvent, const EventSetup& iSetup) override;
  void endRun(const Run& iEvent, const EventSetup& iSetup) override {}
  void endJob() override;

private:
  // ed input tokens
  EDGetTokenT<TTIRMemory> getTokenTTIRMemory_;

  // helper class to store DTC configuration
  trackerDTC::Setup setup_;
  // Setup token
  edm::ESGetToken<trackerDTC::Setup, trackerDTC::SetupRcd> esGetToken_;
  // Input token
  edm::EDGetTokenT<TTDTC> tokenDTC_;


  // Histograms
  TProfile* profStubsPerIR_;
  TProfile* profTTStubsPerIR_;

  TH1F* hisStubsPerCoarseIRRegion_;

  TH2F* hisRZStubs_;
  TH2F* hisRPhiStubs_;

  unsigned int nValidDTCStubs_;
  unsigned int nValidIRStubs_;
};

IRAnalyzer::IRAnalyzer(const ParameterSet& iConfig) :
  tokenDTC_( consumes< TTDTC >( edm::InputTag( iConfig.getParameter<edm::InputTag>( "InputTagTTDTC" ) ) ) ),
  nValidDTCStubs_{0},
  nValidIRStubs_{0}
{
  usesResource("TFileService");
  getTokenTTIRMemory_ = consumes<TTIRMemory>( InputTag( "IRProducer", "IRMemories" ) );

  // book ES product
  esGetToken_ = esConsumes<trackerDTC::Setup, trackerDTC::SetupRcd, edm::Transition::BeginRun>();

}

void IRAnalyzer::beginRun(const Run& iEvent, const EventSetup& iSetup) {
  setup_ = iSetup.getData(esGetToken_);

  // book histograms
  Service<TFileService> fs;
  TFileDirectory dir;
  dir = fs->mkdir("Profiles");
  const unsigned int nIRs = setup_.numRegions() * setup_.numOverlappingRegions() * setup_.numDTCsPerRegion();
  profStubsPerIR_ = dir.make<TProfile>("NStubs", ";", nIRs, 0.5, nIRs+0.5);
  profTTStubsPerIR_ = dir.make<TProfile>("NTTStubs", ";", nIRs, 0.5, nIRs+0.5);
  hisStubsPerCoarseIRRegion_ = dir.make<TH1F>("NStubsPerCoarseIRRegion",";",150,0,150);
  hisRZStubs_ = dir.make<TH2F>("RZ Stubs", ";;", 400, -300, 300, 400, 0., 120);
  hisRPhiStubs_ = dir.make<TH2F>("RPhi Stubs", ";;", 400, -120, 120, 400, -120, 120);

}

void IRAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // Read IR Memories
  Handle<TTIRMemory> handleTTIRMemory;
  iEvent.getByToken<TTIRMemory>(getTokenTTIRMemory_, handleTTIRMemory);
  // Get DTC stubs
  edm::Handle< TTDTC > handleDTC;
  iEvent.getByToken< TTDTC >( tokenDTC_, handleDTC );

  for ( const int& region : handleDTC->tfpRegions() ) {
    for ( const int& channel : handleDTC->tfpChannels() ) {

      unsigned int nValidDTCStubsThisIR{0};
      unsigned int nValidIRStubsThisIR{0};

      const unsigned int irIndex(region * handleDTC->tfpChannels().size() + channel);

      // Get the IR memories containing the output stubs (from the IR)
      TTIRMemory::IRMemories memories{ handleTTIRMemory->IRMemory( irIndex ) };

      // Get the stubs from the DTC (input to IR)
      // Count number of valid stubs 
      const TTDTC::Stream& streamFromDTC{ handleDTC->stream( region, channel ) };
      for ( const auto& frame: streamFromDTC ) {
        if ( frame.first.isNonnull() ) {
          ++nValidDTCStubsThisIR;
        }
      }
      
      // Count number of stubs in IR memories
      unsigned int nStubs{ 0 };
      for ( const auto& memory : memories ) {
        const unsigned int nStubsMem( memory.getEntries(0) );
        nStubs += nStubsMem;
        hisStubsPerCoarseIRRegion_->Fill( nStubsMem );
      }

      profStubsPerIR_->Fill(irIndex, nStubs);

      // Count number of valid ttStubs
      TTIRMemory::TTStubRefs ttStubVectors{ handleTTIRMemory->TTStubs( irIndex ) };
      unsigned int nTTStubs{ 0 };
      for ( const auto& ttStubs : ttStubVectors ) {
        for ( const auto& ttStub: ttStubs ) {
          if ( ttStub.isNonnull() ) {
            nTTStubs += 1;

            ++nValidIRStubsThisIR;

            const GlobalPoint& ttPos = setup_.stubPos(ttStub);
            hisRZStubs_->Fill(ttPos.z(), ttPos.perp());
            hisRPhiStubs_->Fill(ttPos.x(), ttPos.y());
          }
        }
      }
      
      profTTStubsPerIR_->Fill(irIndex, nTTStubs);

      nValidDTCStubs_ += nValidDTCStubsThisIR;
      nValidIRStubs_ += nValidIRStubsThisIR;
    }
  }
}

void IRAnalyzer::endJob() {
  cout << "=============================================================" << endl;
  cout << "                         IR  SUMMARY                         " << endl;

  cout << "Number of valid DTC stubs : " << nValidDTCStubs_ << endl;
  cout << "Number of valid IR stubs : " << nValidIRStubs_ << endl;
  cout << "Fraction of valid DTC stubs in IR memories : " << (double) nValidIRStubs_ / (double) nValidDTCStubs_ << endl;
  cout << "=============================================================" << endl;

}

DEFINE_FWK_MODULE(IRAnalyzer);