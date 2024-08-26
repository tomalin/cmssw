#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerStub.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDStubChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ext/algorithm>

using namespace std;

namespace Phase2Tracker {

  class Phase2TrackerStubProducer : public edm::stream::EDProducer<> {
  public:
    Phase2TrackerStubProducer(const edm::ParameterSet& pset);
    ~Phase2TrackerStubProducer() override = default;
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
    edm::EDGetTokenT<FEDRawDataCollection> token_;
    const Phase2TrackerCabling* cabling_ = nullptr;
    const TrackerGeometry* tGeom_ = nullptr;
    const TrackerTopology* tTopo_ = nullptr;
    std::map<int, std::pair<int, int>> stackMap_;
    DetIdCollection detids_;
  };

  Phase2TrackerStubProducer::Phase2TrackerStubProducer(const edm::ParameterSet& pset)
      : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
        geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
        topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
        token_(consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"))) {
    // define product
    produces<edmNew::DetSetVector<Phase2TrackerStub>>("Stubs");
  }

  void Phase2TrackerStubProducer::beginRun(const edm::Run& run, const edm::EventSetup& es) {
    // fetch cabling from event setup
    cabling_ = &es.getData(ph2CablingESToken_);
    tGeom_ = &es.getData(geomToken_);
    tTopo_ = &es.getData(topoToken_);

    for (auto iu = tGeom_->detUnits().begin(); iu != tGeom_->detUnits().end(); ++iu) {
      unsigned int detId_raw = (*iu)->geographicalId().rawId();
      DetId detId = DetId(detId_raw);
      if (detId.det() == DetId::Detector::Tracker) {
        // build map of upper and lower for each module
        if (tTopo_->isLower(detId) != 0) {
          stackMap_[tTopo_->stack(detId)].first = detId;
        }
        if (tTopo_->isUpper(detId) != 0) {
          stackMap_[tTopo_->stack(detId)].second = detId;
        }
      }
    }  // end loop on detunits
  }

  void Phase2TrackerStubProducer::produce(edm::Event& event, const edm::EventSetup& es) {
    std::unique_ptr<edmNew::DetSetVector<Phase2TrackerStub>> stubs(new edmNew::DetSetVector<Phase2TrackerStub>());

    // Retrieve FEDRawData collection
    edm::Handle<FEDRawDataCollection> buffers;
    event.getByToken(token_, buffers);

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
      // Skip FED if buffer is not a valid tracker FEDBuffer
      if (buffer.isValid() == 0) {
        LogTrace("Phase2TrackerStubProducer") << "[Phase2Tracker::Phase2TrackerStubProducer::" << __func__ << "]: \n";
        LogTrace("Phase2TrackerStubProducer") << "Skipping invalid buffer for FED nr " << fedIndex << endl;
        continue;
      }
      // loop on channels
      // TODO : check loops consistency
      for (int ife = 0; ife < MAX_FE_PER_FED; ife++) {
        // container for this module's digis
        std::vector<Phase2TrackerStub> festubs;
        const Phase2TrackerFEDChannel& channel = buffer.stubchannel(ife);
        if (channel.length() > 0) {
#ifdef EDM_ML_DEBUG
          LogTrace("Phase2TrackerStubProducer") << "[Phase2Tracker::Phase2TrackerStubProducer::" << __func__ << "]: \n";
          LogTrace("Phase2TrackerStubProducer") << "Reading stubs for FED: " << fedIndex << ", FE: " << ife << endl;
#endif

          // get fedid from cabling
          const Phase2TrackerModule mod = cabling_->findFedCh(std::make_pair(fedIndex, ife));
          uint32_t detid = mod.getDetid();
          // create appropriate unpacker
          if (channel.dettype() == DET_Son2S) {
            Phase2TrackerFEDStubSon2SChannelUnpacker unpacker = Phase2TrackerFEDStubSon2SChannelUnpacker(channel);
            while (unpacker.hasData()) {
#ifdef EDM_ML_DEBUG
              LogTrace("Phase2TrackerStubProducer")
                  << "Stub X: " << unpacker.stubX() << " Bend: " << int(unpacker.Bend()) << endl;
#endif
              festubs.push_back(Phase2TrackerStub(unpacker.stubX(), unpacker.Bend()));
              unpacker++;
            }
          } else if (channel.dettype() == DET_PonPS) {
            Phase2TrackerFEDStubPonPSChannelUnpacker unpacker = Phase2TrackerFEDStubPonPSChannelUnpacker(channel);
            while (unpacker.hasData()) {
#ifdef EDM_ML_DEBUG
              LogTrace("Phase2TrackerStubProducer")
                  << "Stub X: " << unpacker.stubX() << ", Stub Y: " << unpacker.stubY()
                  << " Bend: " << int(unpacker.Bend()) << endl;
#endif
              festubs.push_back(Phase2TrackerStub(unpacker.stubX(), unpacker.Bend(), unpacker.stubY()));
              unpacker++;
            }
          }
          if (detid > 0) {
            std::vector<Phase2TrackerStub>::iterator it;
            edmNew::DetSetVector<Phase2TrackerStub>::FastFiller spct(*stubs, detid + 1);
            for (it = festubs.begin(); it != festubs.end(); it++) {
              spct.push_back(*it);
            }
          }  // end if detid > 0 (?)
        }    // end reading stubs channel
      }      // end loop on FE
    }        // en loop on FED
    // store stubs
    event.put(std::move(stubs), "Stubs");
  }
}  // namespace Phase2Tracker

#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerStubProducer Phase2TrackerStubProducer;
DEFINE_FWK_MODULE(Phase2TrackerStubProducer);
