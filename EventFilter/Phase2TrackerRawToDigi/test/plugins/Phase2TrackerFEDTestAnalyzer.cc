//---------------------------------------------------
// Prints info about RAW data file.
//---------------------------------------------------

// system includes
#include <utility>
#include <vector>

// user includes
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDRawChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDZSChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#define LOGPRINT edm::LogPrint("Phase2TrackerFEDTestAnalyzer")

/**
   @class Phase2TrackerFEDTestAnalyzer 
   @brief Analyzes contents of FED_test_ collection
*/

class Phase2TrackerFEDTestAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  typedef std::pair<uint16_t, uint16_t> Fed;
  typedef std::vector<Fed> Feds;
  typedef std::vector<uint16_t> Channels;
  typedef std::map<uint16_t, Channels> ChannelsMap;

  Phase2TrackerFEDTestAnalyzer(const edm::ParameterSet&);
  ~Phase2TrackerFEDTestAnalyzer() override {}

  void beginJob() override {}
  void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override {}
  void endJob() override {}

private:
  const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
  const edm::EDGetTokenT<FEDRawDataCollection> token_;
  const Phase2TrackerCabling* cabling_ = nullptr;
};

using namespace std;
using namespace Phase2Tracker;

// -----------------------------------------------------------------------------
//
Phase2TrackerFEDTestAnalyzer::Phase2TrackerFEDTestAnalyzer(const edm::ParameterSet& pset)
    : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
      token_(consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"))) {
  LogDebug("Phase2TrackerFEDTestAnalyzer") << "[Phase2TrackerFEDTestAnalyzer::" << __func__ << "]"
                                           << "Constructing object...";
}
// -----------------------------------------------------------------------------
//
void Phase2TrackerFEDTestAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& es) {
  // fetch cabling from event setup
  cabling_ = &es.getData(ph2CablingESToken_);
}
// -----------------------------------------------------------------------------
//
void Phase2TrackerFEDTestAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  // Retrieve FEDRawData collection
  edm::Handle<FEDRawDataCollection> buffers;
  event.getByToken(token_, buffers);

  vector<int> feds = cabling_->listFeds();

  // Loop over DTCs
  for (int fedIndex : feds) {
    const FEDRawData& fed = buffers->FEDData(fedIndex);
    if (fed.size() == 0)
      continue;
    // Check which DTC inputs are connected to a module.
    vector<bool> connectedInputs = cabling_->connectedInputs(fedIndex);
    // construct buffer
    Phase2Tracker::Phase2TrackerFEDBuffer* buffer = nullptr;
    buffer = new Phase2Tracker::Phase2TrackerFEDBuffer(fed.data(), fed.size(), connectedInputs);

    LOGPRINT << " -------------------------------------------- ";
    LOGPRINT << " buffer debug ------------------------------- ";
    LOGPRINT << " -------------------------------------------- ";
    LOGPRINT << " buffer size : " << buffer->bufferSize();
    LOGPRINT << " fed id      : " << fedIndex;
    LOGPRINT << " -------------------------------------------- ";
    LOGPRINT << " tracker header debug ------------------------";
    LOGPRINT << " -------------------------------------------- ";

    Phase2TrackerFEDHeader tr_header = buffer->trackerHeader();
    LOGPRINT << " Version  : " << hex << setw(2) << (int)tr_header.getDataFormatVersion();
    LOGPRINT << " Mode     : " << hex << setw(2) << (int)tr_header.getDebugMode();
    LOGPRINT << " Type     : " << hex << setw(2) << (int)tr_header.getEventType();
    LOGPRINT << " Readout  : " << hex << setw(2) << (int)tr_header.getReadoutMode();
    LOGPRINT << " Status   : " << hex << setw(16) << (int)tr_header.getGlibStatusCode();
    LOGPRINT << " FE stat  : ";
    for (int i = 15; i >= 0; i--) {
      if ((tr_header.frontendStatus())[i]) {
        LOGPRINT << "1";
      } else {
        LOGPRINT << "0";
      }
    }
    LOGPRINT << endl;
    LOGPRINT << " Nr CBC   : " << hex << setw(16) << (int)tr_header.getNumberOfCBC() << endl;
    LOGPRINT << " FE/Chip status : ";
    vector<Phase2TrackerFEDFEDebug> all_fe_debug = tr_header.CBCStatus();
    vector<Phase2TrackerFEDFEDebug>::iterator FE_it;
    for (FE_it = all_fe_debug.begin(); FE_it < all_fe_debug.end(); FE_it++) {
      if (FE_it->IsOn()) {
        LOGPRINT << " FE L1ID: " << endl;
        LOGPRINT << "    " << hex << setw(4) << FE_it->getFEL1ID()[0] << dec << endl;
        LOGPRINT << "    " << hex << setw(4) << FE_it->getFEL1ID()[1] << dec << endl;
        for (int i = 0; i < 16; i++) {
          LOGPRINT << " Chip Error" << hex << setw(1) << FE_it->getChipError(i) << dec << endl;
          LOGPRINT << " Chip L1ID " << hex << setw(4) << FE_it->getChipL1ID(i) << dec << endl;
          LOGPRINT << " Chip PA   " << hex << setw(4) << FE_it->getChipPipelineAddress(i) << dec << endl;
        }
      }
    }
    LOGPRINT << endl;
    LOGPRINT << " -------------------------------------------- " << endl;
    LOGPRINT << " Payload  ----------------------------------- " << endl;
    LOGPRINT << " -------------------------------------------- " << endl;

    // loop channels
    int ichan = 0;
    for (int ife = 0; ife < 16; ife++) {
      for (int icbc = 0; icbc < 16; icbc++) {
        const Phase2TrackerFEDChannel& channel = buffer->channel(ichan);
        if (channel.length() > 0) {
          LOGPRINT << dec << " reading channel : " << icbc << " on FE " << ife;
          LOGPRINT << dec << " with length  : " << (int)channel.length();
          Phase2TrackerFEDRawChannelUnpacker unpacker = Phase2TrackerFEDRawChannelUnpacker(channel);
          while (unpacker.hasData()) {
            LOGPRINT << (unpacker.stripOn() ? "1" : "_");
            unpacker++;
          }
          LOGPRINT << "\n";
        }
        ichan++;
      }
    }  // end loop on channels
  }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Phase2TrackerFEDTestAnalyzer);
