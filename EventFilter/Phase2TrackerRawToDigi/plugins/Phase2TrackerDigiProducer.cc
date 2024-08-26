#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDRawChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDZSChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace Phase2Tracker {

  class Phase2TrackerDigiProducer : public edm::stream::EDProducer<> {
  public:
    /// constructor
    Phase2TrackerDigiProducer(const edm::ParameterSet& pset);
    /// default constructor
    ~Phase2TrackerDigiProducer() override = default;
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    std::map<int, std::pair<int, int>> stackMap_;
    const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
    const edm::EDGetTokenT<FEDRawDataCollection> token_;
    const edm::EDPutTokenT<edm::DetSetVector<Phase2TrackerDigi>> putToken_;
    const edm::EDPutTokenT<edmNew::DetSetVector<Phase2TrackerCluster1D>> putTokenSparsified_;
    const TrackerGeometry* tkGeom_ = nullptr;
    const TrackerTopology* tTopo_ = nullptr;
    const Phase2TrackerCabling* cabling_ = nullptr;
    class Registry {
    public:
      /// constructor
      Registry(uint32_t aDetid, uint16_t firstStrip, size_t indexInVector, uint16_t numberOfDigis)
          : detid(aDetid), first(firstStrip), index(indexInVector), length(numberOfDigis) {}
      /// < operator to sort registries
      bool operator<(const Registry& other) const {
        return (detid != other.detid ? detid < other.detid : first < other.first);
      }
      /// public data members
      uint32_t detid;
      uint16_t first;
      size_t index;
      uint16_t length;
    };
  };
}  // namespace Phase2Tracker

#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerDigiProducer Phase2TrackerDigiProducer;
DEFINE_FWK_MODULE(Phase2TrackerDigiProducer);

using namespace std;

namespace Phase2Tracker {

  Phase2TrackerDigiProducer::Phase2TrackerDigiProducer(const edm::ParameterSet& pset)
      : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
        geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
        topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
        token_(consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"))),
        putToken_(produces<edm::DetSetVector<Phase2TrackerDigi>>("Unsparsified")),
        putTokenSparsified_(produces<edmNew::DetSetVector<Phase2TrackerCluster1D>>("Sparsified")) {}

  void Phase2TrackerDigiProducer::beginRun(const edm::Run& run, const edm::EventSetup& es) {
    // fetch cabling from event setup
    cabling_ = &es.getData(ph2CablingESToken_);

    // FIXME: build map of stacks to compensate for missing trackertopology methods
    tkGeom_ = &es.getData(geomToken_);
    tTopo_ = &es.getData(topoToken_);

    for (auto iu = tkGeom_->detUnits().begin(); iu != tkGeom_->detUnits().end(); ++iu) {
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

  void Phase2TrackerDigiProducer::produce(edm::Event& event, const edm::EventSetup& es) {
    std::vector<Registry> proc_work_registry_;
    std::vector<Phase2TrackerDigi> proc_work_digis_;

    auto clusters = std::make_unique<edmNew::DetSetVector<Phase2TrackerCluster1D>>();

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
        LogTrace("Phase2TrackerDigiProducer") << "[Phase2Tracker::Phase2TrackerDigiProducer::" << __func__ << "]: \n";
        LogTrace("Phase2TrackerDigiProducer") << "Skipping invalid buffer for FED nr " << fedIndex << endl;
        continue;
      }

#ifdef EDM_ML_DEBUG
      std::ostringstream ss;
      ss << " -------------------------------------------- " << endl;
      ss << " buffer debug ------------------------------- " << endl;
      ss << " -------------------------------------------- " << endl;
      ss << " buffer size : " << buffer.bufferSize() << endl;
      ss << " fed id      : " << fedIndex << endl;
      ss << " -------------------------------------------- " << endl;
      ss << " tracker header debug ------------------------" << endl;
      ss << " -------------------------------------------- " << endl;
      LogTrace("Phase2TrackerDigiProducer") << ss.str();
      ss.clear();
      ss.str("");
#endif

      Phase2TrackerFEDHeader tr_header = buffer.trackerHeader();

#ifdef EDM_ML_DEBUG
      ss << " Version  : " << hex << setw(2) << (int)tr_header.getDataFormatVersion() << endl;
      ss << " Mode     : " << hex << setw(2) << tr_header.getDebugMode() << endl;
      ss << " Type     : " << hex << setw(2) << (int)tr_header.getEventType() << endl;
      ss << " Readout  : " << hex << setw(2) << tr_header.getReadoutMode() << endl;
      ss << " Condition Data : " << (tr_header.getConditionData() ? "Present" : "Absent") << "\n";
      ss << " Data Type      : " << (tr_header.getDataType() ? "Real" : "Fake") << "\n";
      ss << " Status   : " << hex << setw(16) << (int)tr_header.getGlibStatusCode() << endl;
      ss << " FE stat  : ";
      for (int i = MAX_FE_PER_FED - 1; i >= 0; i--) {
        if ((tr_header.frontendStatus())[i]) {
          ss << "1";
        } else {
          ss << "0";
        }
      }
      ss << endl;
      ss << " Nr CBC   : " << hex << setw(16) << (int)tr_header.getNumberOfCBC() << endl;
      ss << " FE/Chip status : ";
      std::vector<Phase2TrackerFEDFEDebug> all_fe_debug = tr_header.CBCStatus();
      std::vector<Phase2TrackerFEDFEDebug>::iterator FE_it;
      for (FE_it = all_fe_debug.begin(); FE_it < all_fe_debug.end(); FE_it++) {
        if (FE_it->IsOn()) {
          ss << " FE L1ID: " << endl;
          ss << "    " << hex << setw(4) << FE_it->getFEL1ID()[0] << dec << endl;
          ss << "    " << hex << setw(4) << FE_it->getFEL1ID()[1] << dec << endl;
          for (int i = 0; i < 16; i++) {
            ss << " Chip Error" << hex << setw(1) << FE_it->getChipError(i) << dec << endl;
            ss << " Chip L1ID " << hex << setw(4) << FE_it->getChipL1ID(i) << dec << endl;
            ss << " Chip PA   " << hex << setw(4) << FE_it->getChipPipelineAddress(i) << dec << endl;
          }
        }
      }
      LogTrace("Phase2TrackerDigiProducer") << ss.str();
      ss.clear();
      ss.str("");
      ss << " -------------------------------------------- " << endl;
      ss << " Payload  ----------------------------------- " << endl;
      ss << " -------------------------------------------- " << endl;
#endif
      // check readout mode
      //std::cout<<"READOUT MODE "<<(tr_header.getReadoutMode() == READOUT_MODE_PROC_RAW) <<" "<<(tr_header.getReadoutMode() == READOUT_MODE_ZERO_SUPPRESSED) <<std::endl;
      if (tr_header.getReadoutMode() == READOUT_MODE_PROC_RAW) {
        // N.B. RAW mode only exists for 2S modules.
        // loop channels (one per CBC chip)
        int ichan = -1;
        for (int ife = 0; ife < MAX_FE_PER_FED; ife++) {
          for (int icbc = 0; icbc < MAX_CBC_PER_FE; icbc++) {
            ichan++;
            const Phase2TrackerFEDChannel& channel = buffer.channel(ichan);
            if (channel.length() > 0) {  // DTC input channel active
              // get detid from cabling
              const Phase2TrackerModule& mod = cabling_->findFedCh(std::make_pair(fedIndex, ife));
              if (not mod.connected())
                throw cms::Exception("Phase2TrackerDigiProducer") << "Active DTC input is disconnected in cabling map";
              uint32_t detid = mod.getDetid();

#ifdef EDM_ML_DEBUG
              ss << dec << " id from cabling : " << detid << endl;
              ss << dec << " reading channel : " << icbc << " on FE " << ife;
              ss << dec << " with length  : " << (int)channel.length() << endl;
#endif
              // container for this channel's digis
              std::vector<Phase2TrackerDigi> stripsTop;
              std::vector<Phase2TrackerDigi> stripsBottom;

              // unpacking data
              // TODO : set Y position in digi as a function of the side of the CBC

              // IAN: Think this sets y-position (=iside) and avoids x-overflow.
              // (CBC <= 8 is y-end 0 & >=8 is y-end 1 of module).
              int icbc8 = icbc % (MAX_CBC_PER_FE / 2);
              int iside = icbc / (MAX_CBC_PER_FE / 2);

              Phase2TrackerFEDRawChannelUnpacker unpacker = Phase2TrackerFEDRawChannelUnpacker(channel);
              while (unpacker.hasData()) {
                if (unpacker.stripOn()) {
                  if (unpacker.stripIndex() % 2) {
                    stripsTop.push_back(
                        Phase2TrackerDigi((int)(STRIPS_PER_CBC * icbc8 + unpacker.stripIndex()) / 2, iside));
#ifdef EDM_ML_DEBUG
                    ss << "t";
#endif
                  } else {
                    stripsBottom.push_back(
                        Phase2TrackerDigi((int)(STRIPS_PER_CBC * icbc8 + unpacker.stripIndex()) / 2, iside));
#ifdef EDM_ML_DEBUG
                    ss << "b";
#endif
                  }
                } else {
#ifdef EDM_ML_DEBUG
                  ss << "_";
#endif
                }
                unpacker++;
              }
#ifdef EDM_ML_DEBUG
              ss << endl;
              LogTrace("Phase2TrackerDigiProducer") << ss.str();
              ss.clear();
              ss.str("");
#endif

              // store beginning and end of this digis for this detid and add this registry to the list
              // and store data
              // FIXME : detid scheme should be taken from topology / geometry
              Registry regItemTop(
                  stackMap_[detid].second, STRIPS_PER_CBC * icbc / 2, proc_work_digis_.size(), stripsTop.size());
              proc_work_registry_.push_back(regItemTop);
              proc_work_digis_.insert(proc_work_digis_.end(), stripsTop.begin(), stripsTop.end());
              Registry regItemBottom(
                  stackMap_[detid].first, STRIPS_PER_CBC * icbc / 2, proc_work_digis_.size(), stripsBottom.size());
              proc_work_registry_.push_back(regItemBottom);
              proc_work_digis_.insert(proc_work_digis_.end(), stripsBottom.begin(), stripsBottom.end());
            }
          }
        }  // end loop on channels

      } else if (tr_header.getReadoutMode() == READOUT_MODE_ZERO_SUPPRESSED) {
        // loop channels (4 per DTC input channel)
        int ichan = -1;
        for (int ife = 0; ife < MAX_FE_PER_FED; ife++) {  // Loop inputs of this DTC

          // container for this module's digis
          std::vector<Phase2TrackerCluster1D> clustersTop;
          std::vector<Phase2TrackerCluster1D> clustersBottom;

          // loop over concentrators (4 virtual concentrators in case of PS)
          // N.B. Digis in each DTC input channel assumed ordered in raw data
          // by lower-left, lower-right, upper-left, upper-right sensor.
          // During Raw2Digi conversion, Phase2TrackerFEDBuffer profits
          // from this to create 4 channels corresponding to
          // p-left, p-right, s-left, s-right for each DTC input.
          // (p-left & p-right have zero clusters for 2S module).

          uint32_t detid = 0;
          for (int iconc = 0; iconc < 4; iconc++) {
            ichan++;
            const Phase2TrackerFEDChannel& channel = buffer.channel(ichan);
            std::cout << "CHANNEL LENGTH " << fedIndex << " " << ife << " " << ichan << " " << channel.length() << " "
                      << tr_header.frontendStatus()[ife] << std::endl;
            if (channel.length() > 0) {  // DTC input channel active

              // get detid from cabling
              const Phase2TrackerModule& mod = cabling_->findFedCh(std::make_pair(fedIndex, ife));
              if (not mod.connected())
                throw cms::Exception("Phase2TrackerDigiProducer") << "Active DTC input is disconnected in cabling map";
              detid = mod.getDetid();

              std::cout << "FED channel " << fedIndex << " " << ife << " " << iconc << " ichan " << ichan << " detid "
                        << (uint64_t)detid << " length " << channel.length() << " "
                        << " status " << tr_header.frontendStatus()[ife] << " type " << int(channel.dettype())
                        << std::endl;

#ifdef EDM_ML_DEBUG
              ss << dec << " id from cabling : " << detid << endl;
              ss << dec << " reading channel : " << iconc << " on FE " << ife;
              ss << dec << " with length  : " << (int)channel.length() << endl;
#endif
              // create appropriate unpacker
              if (channel.dettype() == DET_Son2S) {
                Phase2TrackerFEDZSSon2SChannelUnpacker unpacker = Phase2TrackerFEDZSSon2SChannelUnpacker(channel);
                while (unpacker.hasData()) {
                  unpacker.Merge();
#ifdef EDM_ML_DEBUG
                  ss << std::dec << " Son2S " << (int)unpacker.clusterX() << " " << (int)unpacker.clusterSize() << " "
                     << (int)unpacker.chipId() << endl;
#endif
                  if (unpacker.rawX() % 2) {
                    clustersTop.push_back(
                        Phase2TrackerCluster1D(unpacker.clusterX(), unpacker.clusterY(), unpacker.clusterSize()));
                  } else {
                    clustersBottom.push_back(
                        Phase2TrackerCluster1D(unpacker.clusterX(), unpacker.clusterY(), unpacker.clusterSize()));
                  }
                  unpacker++;
                }
              } else if (channel.dettype() == DET_SonPS) {
                Phase2TrackerFEDZSSonPSChannelUnpacker unpacker = Phase2TrackerFEDZSSonPSChannelUnpacker(channel);
                while (unpacker.hasData()) {
                  unpacker.Merge();
#ifdef EDM_ML_DEBUG
                  ss << std::dec << " SonPS " << (int)unpacker.clusterX() << " " << (int)unpacker.clusterSize() << " "
                     << (int)unpacker.chipId() << endl;
#endif
                  clustersTop.push_back(Phase2TrackerCluster1D(
                      unpacker.clusterX(), unpacker.clusterY(), unpacker.clusterSize(), unpacker.threshold()));
                  unpacker++;
                }
              } else if (channel.dettype() == DET_PonPS) {
                Phase2TrackerFEDZSPonPSChannelUnpacker unpacker = Phase2TrackerFEDZSPonPSChannelUnpacker(channel);
                while (unpacker.hasData()) {
                  unpacker.Merge();
#ifdef EDM_ML_DEBUG
                  ss << std::dec << " PonPS " << (int)unpacker.clusterX() << " " << (int)unpacker.clusterSize() << " "
                     << (int)unpacker.clusterY() << " " << (int)unpacker.chipId() << endl;
#endif
                  clustersBottom.push_back(
                      Phase2TrackerCluster1D(unpacker.clusterX(), unpacker.clusterY(), unpacker.clusterSize()));
                  unpacker++;
                }
              }
#ifdef EDM_ML_DEBUG
              ss << endl;
              LogTrace("Phase2TrackerDigiProducer") << ss.str();
              ss.clear();
              ss.str("");
#endif
            }  // end reading left/right/top/bottom sensor.
          }    // end loop on channels

          // Store clusters of this FED.
          std::vector<Phase2TrackerCluster1D>::iterator it;
          {
            // outer detid is defined as inner detid + 1 or module detid + 2
            edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller spct(*clusters, stackMap_[detid].second);
            for (it = clustersTop.begin(); it != clustersTop.end(); it++) {
              spct.push_back(*it);
            }
          }
          {
            edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller spcb(*clusters, stackMap_[detid].first);
            for (it = clustersBottom.begin(); it != clustersBottom.end(); it++) {
              spcb.push_back(*it);
            }
          }

        }  // end loop on FE
      } else {
        throw cms::Exception("Phase2TrackerDigiProducer") << "Unknown readout mode";
      }
    }

    // Sort and store Unsparsified digis. N.B. Here, digis are from all DTCs.
    std::sort(proc_work_registry_.begin(), proc_work_registry_.end());
    std::vector<edm::DetSet<Phase2TrackerDigi>> sorted_and_merged;
    std::vector<Registry>::iterator it = proc_work_registry_.begin(), it2 = it + 1, end = proc_work_registry_.end();
    while (it < end) {
      sorted_and_merged.push_back(edm::DetSet<Phase2TrackerDigi>(it->detid));
      std::vector<Phase2TrackerDigi>& digis = sorted_and_merged.back().data;
      // first count how many digis we have
      size_t len = it->length;
      for (it2 = it + 1; (it2 != end) && (it2->detid == it->detid); ++it2) {
        len += it2->length;
      }
      // reserve memory
      digis.reserve(len);
      // push them in
      for (it2 = it + 0; (it2 != end) && (it2->detid == it->detid); ++it2) {
        digis.insert(digis.end(), &proc_work_digis_[it2->index], &proc_work_digis_[it2->index + it2->length]);
      }
      it = it2;
    }
    //edm::DetSetVector<Phase2TrackerDigi> proc_raw_dsv( sorted_and_merged, true );
    //pr->swap( proc_raw_dsv );
    //event.put(std::unique_ptr<edm::DetSetVector<Phase2TrackerDigi>>(pr), "Unsparsified");

    //Is this correct?..Jerome did the above
    event.emplace(putToken_, sorted_and_merged, true);

    // store Sparsified Digis
    event.put(putTokenSparsified_, std::move(clusters));
  }
}  // namespace Phase2Tracker
