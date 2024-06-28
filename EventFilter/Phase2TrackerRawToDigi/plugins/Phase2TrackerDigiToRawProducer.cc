#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
//#include "DataFormats/FEDRawData/src/fed_header.h"
//#include "DataFormats/FEDRawData/src/fed_trailer.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerDigiToRaw.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include <sstream>
#include <iomanip>
#include <ext/algorithm>

using namespace std;
namespace Phase2Tracker {

  class Phase2TrackerDigiToRawProducer : public edm::stream::EDProducer<> {
  public:
    /// constructor
    Phase2TrackerDigiToRawProducer(const edm::ParameterSet& pset);
    /// default constructor
    ~Phase2TrackerDigiToRawProducer() override = default;
    void beginRun(edm::Run const&, edm::EventSetup const&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
    const edm::EDGetTokenT<edmNew::DetSetVector<Phase2TrackerCluster1D>> token_;
    const Phase2TrackerCabling* cabling_ = nullptr;
    const TrackerTopology* tTopo_ = nullptr;
    const TrackerGeometry* tGeom_ = nullptr;
    std::map<int, std::pair<int, int>> stackMap_;
  };

  Phase2TrackerDigiToRawProducer::Phase2TrackerDigiToRawProducer(const edm::ParameterSet& pset)
      : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
        geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
        topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
        token_(
            consumes<edmNew::DetSetVector<Phase2TrackerCluster1D>>(pset.getParameter<edm::InputTag>("ProductLabel"))) {
    produces<FEDRawDataCollection>();
  }

  void Phase2TrackerDigiToRawProducer::beginRun(edm::Run const& run, edm::EventSetup const& es) {
    cabling_ = &es.getData(ph2CablingESToken_);
    tGeom_ = &es.getData(geomToken_);
    tTopo_ = &es.getData(topoToken_);

    // // Build list of detids for dummy cabling: this should go in another producer"
    // // FIXME: remove this and replace by an external package
    // int fedid = 0;
    // int channel = 0;
    // for (auto iu = tGeom_->detUnits().begin(); iu != tGeom_->detUnits().end(); ++iu) {
    //   unsigned int detId_raw = (*iu)->geographicalId().rawId();
    //   DetId detId = DetId(detId_raw);
    //   if (detId.det() == DetId::Detector::Tracker) {
    //     if ( tTopo_->isLower(detId) != 0 ) {
    //       std::cout << tTopo_->stack(detId) << " " << fedid << " " << channel << std::endl;
    //       channel += 1;
    //       if (channel == 72) {
    //         channel = 0;
    //         fedid += 1;
    //       }
    //     }
    //   }
    // }

    // build map of upper and lower for each module
    for (auto iu = tGeom_->detUnits().begin(); iu != tGeom_->detUnits().end(); ++iu) {
      unsigned int detId_raw = (*iu)->geographicalId().rawId();
      DetId detId = DetId(detId_raw);
      if (detId.det() == DetId::Detector::Tracker) {
        if (tTopo_->isLower(detId) != 0) {
          stackMap_[tTopo_->stack(detId)].first = detId;
        }
        if (tTopo_->isUpper(detId) != 0) {
          stackMap_[tTopo_->stack(detId)].second = detId;
        }
      }
    }  // end loop on detunits
  }

  void Phase2TrackerDigiToRawProducer::produce(edm::Event& event, const edm::EventSetup& es) {
    std::unique_ptr<FEDRawDataCollection> buffers(new FEDRawDataCollection);
    edm::Handle<edmNew::DetSetVector<Phase2TrackerCluster1D>> digis_handle;
    event.getByToken(token_, digis_handle);
    Phase2TrackerDigiToRaw raw_producer(cabling_, tGeom_, tTopo_, stackMap_, digis_handle, 1);
    raw_producer.buildFEDBuffers(buffers);
    event.put(std::move(buffers));
  }
}  // namespace Phase2Tracker

#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerDigiToRawProducer Phase2TrackerDigiToRawProducer;
DEFINE_FWK_MODULE(Phase2TrackerDigiToRawProducer);
