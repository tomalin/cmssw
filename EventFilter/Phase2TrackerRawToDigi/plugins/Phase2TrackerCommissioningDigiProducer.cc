#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerCommissioningDigi.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace Phase2Tracker {
  typedef edm::DetSet<Phase2TrackerCommissioningDigi> condata_map;
class Phase2TrackerCommissioningDigiProducer : public edm::one::EDProducer<edm::one::WatchRuns> {
  public:
    /// constructor
    Phase2TrackerCommissioningDigiProducer(const edm::ParameterSet& pset);
    /// default constructor
    ~Phase2TrackerCommissioningDigiProducer() override = default;
    void beginRun(const edm::Run& run, const edm::EventSetup& es)  override;
    void endRun(const edm::Run& run, const edm::EventSetup& es)  override {}
    void produce(edm::Event& ev, const edm::EventSetup& es) override;

  private:
    const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
    edm::EDGetTokenT<FEDRawDataCollection> token_;
    const Phase2TrackerCabling* cabling_ = nullptr;
  };

Phase2TrackerCommissioningDigiProducer::Phase2TrackerCommissioningDigiProducer(const edm::ParameterSet& pset) :
  ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()) {
    produces<condata_map>("ConditionData");
    token_ = consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"));
  }

  void Phase2TrackerCommissioningDigiProducer::beginRun(const edm::Run& run, const edm::EventSetup& es) {
    // fetch cabling from event setup
    cabling_ = &es.getData(ph2CablingESToken_);
}

  void Phase2TrackerCommissioningDigiProducer::produce(
                                                       edm::Event& event,
                                                       const edm::EventSetup& es)  {
    // Retrieve FEDRawData collection
    edm::Handle<FEDRawDataCollection> buffers;
    event.getByToken(token_, buffers);

    // Analyze strip tracker FED buffers in data
    std::vector<int> feds = cabling_->listFeds();
    
    // Loop over DTCs
    for (int fedIndex : feds) {

      // Check which DTC inputs are connected to a module.
      std::vector<bool> connectedInputs = cabling_->connectedInputs(fedIndex);      
    
      // reading
      const FEDRawData& fed = buffers->FEDData(fedIndex);
      if (fed.size() == 0)
        continue;
      // construct buffer
      Phase2Tracker::Phase2TrackerFEDBuffer buffer(fed.data(), fed.size(), connectedInputs);
      std::map<uint32_t, uint32_t> cond_data = buffer.conditionData();

      // print cond data for debug
      LogTrace("Phase2TrackerCommissioningDigiProducer") << "--- Condition data debug ---" << std::endl;
      std::map<uint32_t, uint32_t>::const_iterator it;
      //for (it = cond_data.begin(); it != cond_data.end(); it++) {
      for (auto& [key, value] : cond_data) {
        LogTrace("Phase2TrackerCommissioningDigiProducer")
            << "key: " << std::hex << std::setw(8) << std::setfill('0') << key << " value: " << std::hex << std::setw(8)
            << std::setfill('0') << value << " (hex) " << std::dec << value << " (dec) " << std::endl;
      }
      LogTrace("Phase2TrackerCommissioningDigiProducer") << "----------------------------" << std::endl;
      // store it into digis
      condata_map* cond_data_digi = new condata_map(fedIndex);
      for (it = cond_data.begin(); it != cond_data.end(); it++) {
        cond_data_digi->push_back(Phase2TrackerCommissioningDigi(it->first, it->second));
      }
      std::unique_ptr<condata_map> cdd(cond_data_digi);
      event.put(std::move(cdd), "ConditionData");
    }
  }
}  // namespace Phase2Tracker
#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerCommissioningDigiProducer Phase2TrackerCommissioningDigiProducer;
DEFINE_FWK_MODULE(Phase2TrackerCommissioningDigiProducer);
