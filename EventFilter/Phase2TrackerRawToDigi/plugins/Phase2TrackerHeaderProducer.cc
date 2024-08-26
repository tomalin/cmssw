// -*- C++ -*-
//
// Package:    EventFilter/Phase2TrackerRawToDigi/Phase2TrackerHeaderProducer
// Class:      Phase2TrackerHeaderProducer
//
/**\class Phase2TrackerHeaderProducer Phase2TrackerHeaderProducer.cc EventFilter/Phase2TrackerRawToDigi/plugins/Phase2TrackerHeaderProducer.cc

 Description: Producer for the phase 2 tracker header digi

*/
//
// Original Author:  Jerome De Favereau De Jeneret
//         Created:  Mon, 01 Sep 2014 08:42:31 GMT
//

#include <memory>
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerHeaderDigi.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace std;

namespace Phase2Tracker {

  typedef std::vector<Phase2TrackerHeaderDigi> header_map;

  class Phase2TrackerHeaderProducer : public edm::stream::EDProducer<> {
  public:
    explicit Phase2TrackerHeaderProducer(const edm::ParameterSet&);
    ~Phase2TrackerHeaderProducer() override = default;
    void beginRun(const edm::Run& run, const edm::EventSetup& es)  override;
    void produce(edm::Event& event, const edm::EventSetup& es) override;

  private:
    const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
    edm::EDGetTokenT<FEDRawDataCollection> token_;
    const Phase2TrackerCabling* cabling_ = nullptr;
  };

  Phase2Tracker::Phase2TrackerHeaderProducer::Phase2TrackerHeaderProducer(const edm::ParameterSet& pset) :
    ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>())
  {
    produces<header_map>("TrackerHeader");
    token_ = consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"));
  }


  void Phase2TrackerHeaderProducer::beginRun(const edm::Run& run, const edm::EventSetup& es) {
    // fetch cabling from event setup
    cabling_ = &es.getData(ph2CablingESToken_);
  }

  void Phase2Tracker::Phase2TrackerHeaderProducer::produce(edm::Event& event, const edm::EventSetup& es) {
    // Retrieve FEDRawData collection
    edm::Handle<FEDRawDataCollection> buffers;
    event.getByToken(token_, buffers);

    // fill collection
    std::unique_ptr<header_map> hdigis(new header_map);

    // Analyze strip tracker FED buffers in data
    std::vector<int> feds = cabling_->listFeds();
    
    // Loop over DTCs
    for (int fedIndex : feds) {
    
      const FEDRawData& fed = buffers->FEDData(fedIndex);
      if (fed.size() == 0)
        continue;
      // Check which DTC inputs are connected to a module.
      std::vector<bool> connectedInputs = cabling_->connectedInputs(fedIndex);
      // construct buffer
      Phase2Tracker::Phase2TrackerFEDBuffer buffer(fed.data(), fed.size(), connectedInputs);
      Phase2TrackerHeaderDigi head_digi = Phase2TrackerHeaderDigi(buffer.trackerHeader());
      // store digis
      hdigis->push_back(head_digi);
    }
    event.put(std::move(hdigis), "TrackerHeader");
  }
}  // namespace Phase2Tracker

#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerHeaderProducer Phase2TrackerHeaderProducer;
DEFINE_FWK_MODULE(Phase2TrackerHeaderProducer);
