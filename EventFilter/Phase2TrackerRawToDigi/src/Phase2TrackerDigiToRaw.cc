#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerDigiToRaw.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/DetSetVector.h"

namespace Phase2Tracker {
  const int MAX_NP = 31;  // max P clusters per concentrator i.e. per side
  const int MAX_NS = 31;  // same for S clusters

  std::pair<int, int> SortExpandAndLimitClusters(std::vector<StackedClus>& scluss, int max_ns, int max_np) {
    std::vector<StackedClus> processed;
    // number of clusters allowed : P-left, P-right, S-left, S-right
    int roomleft[4] = {max_np, max_np, max_ns, max_ns};
    // fill left and right vectors, expand big clusters
    for (auto sclus = scluss.begin(); sclus < scluss.end(); sclus++) {
      std::vector<StackedClus> parts = sclus->splitDigi();
      for (auto id = parts.begin(); id < parts.end(); id++) {
        if (roomleft[id->getSideType()] > 0) {
          processed.push_back(*id);
          roomleft[id->getSideType()] -= 1;
        }
      }
    }
    // Sort vector
    std::sort(processed.begin(), processed.end());
    // replace input vector with truncated one
    scluss.swap(processed);
    // return number of S and P clusters respectively
    return std::make_pair(2 * max_ns - roomleft[2] - roomleft[3], 2 * max_np - roomleft[0] - roomleft[1]);
  }

  Phase2TrackerDigiToRaw::Phase2TrackerDigiToRaw(const Phase2TrackerCabling* cabling,
                                                 const TrackerGeometry* tGeom,
                                                 const TrackerTopology* tTopo,
                                                 std::map<int, std::pair<int, int>> stackMap,
                                                 edm::Handle<DetSetVecClus> cluss_handle,
                                                 int mode)
      : cabling_(cabling),
        tTopo_(tTopo),
        tGeom_(tGeom),
        stackMap_(stackMap),
        clusshandle_(cluss_handle),
        mode_(mode),
        FedDaqHeader_(0, 0, 0, DAQ_EVENT_TYPE_SIMULATED),  // TODO : add L1ID
        FedDaqTrailer_(0, 0) {
    FedHeader_.setDataFormatVersion(2);
    FedHeader_.setDebugMode(SUMMARY);
    FedHeader_.setEventType((uint8_t)0x04);
  }

  void Phase2TrackerDigiToRaw::buildFEDBuffers(std::unique_ptr<FEDRawDataCollection>& rcollection) {
    // store all clusters for a given fedid
    std::vector<DetSetVecClus::const_iterator> cluss_t;
    // note which DTC input channels have one cluster in this event
    std::vector<bool> festatus(72, false);

    std::vector<int> feds = cabling_->listFeds();
    for (int fedIndex : feds) {                         // Loop over DTCs
      for (int ife = 0; ife < MAX_FE_PER_FED; ife++) {  // Loop inputs of this DTC

        // get detid from cabling
        const Phase2TrackerModule& mod = cabling_->findFedCh(std::make_pair(fedIndex, ife));
        if (not mod.connected())
          continue;  // Skip unconnected DTC inputs.
        int detidStack = mod.getDetid();
        // detid of first plane is detid of module + 1
        // FIXME (when it is implemented) : we should get this from the topology / geometry
        unsigned int detid = stackMap_[detidStack].first;
        // end of fixme

        DetSetVecClus::const_iterator cluss;
        cluss = clusshandle_->find(detid);
        if (cluss != clusshandle_->end()) {
          cluss_t.push_back(cluss);
          festatus[ife] = true;
        }
        cluss = clusshandle_->find(tTopo_->partnerDetId(detid));
        if (cluss != clusshandle_->end()) {
          cluss_t.push_back(cluss);
          festatus[ife] = true;
        }
        std::cout << "FED CHANNEL " << fedIndex << " " << ife << " " << festatus[ife] << std::endl;
      }
      // save buffer
      FedHeader_.setFrontendStatus(festatus);

      // Write to buffer the clusters of one DTC.
      // The vector cluss_t contains one entry for each silicon sensor that contains at
      // least one cluster and is connected to this DTC.
      std::vector<uint64_t> fedbuffer = makeBuffer(cluss_t);

      std::cout << "BUFFERMADE " << cluss_t.size() << " ";
      for (bool aa : festatus)
        std::cout << aa;
      std::cout << std::endl;

      FEDRawData& frd = rcollection->FEDData(fedIndex);
      int size = fedbuffer.size() * 8;
      frd.resize(size);
      memcpy(frd.data(), &fedbuffer[0], size);
      festatus.assign(72, false);
      cluss_t.clear();
    }
  }

  // Build buffer for a single DTC, consisting of FEDDAQHeader & FEDHeader,
  // and for each DTC input channel that receives at least one cluster,
  // a FEDHeaderSparsified and all that channels clusters.

  std::vector<uint64_t> Phase2TrackerDigiToRaw::makeBuffer(std::vector<DetSetVecClus::const_iterator> cluss) {
    uint64_t bitindex = 0;
    int moduletype = -1;
    // IAN note: buffer size continually increased. Very inefficient. Should calculate size of buffer. Added call to reserve() to help.
    std::vector<uint64_t> fedbuffer;
    constexpr unsigned int sizeGuess = 1000;
    fedbuffer.reserve(sizeGuess);
    // add daq header
    fedbuffer.push_back(*(uint64_t*)FedDaqHeader_.data());
    bitindex += 64;
    // add fed header
    uint8_t* feh = FedHeader_.data();
    fedbuffer.push_back(*(uint64_t*)feh);
    fedbuffer.push_back(*(uint64_t*)(feh + 8));
    bitindex += 128;
    // loop on sensors in tracker modules read by this DTC.
    std::vector<DetSetVecClus::const_iterator>::const_iterator sensor;
    for (sensor = cluss.begin(); sensor != cluss.end(); sensor++) {
      // get id of sensor
      unsigned int detid = (*sensor)->detId();
      TrackerGeometry::ModuleType det_type = tGeom_->getDetectorType(detid);
      if (det_type == TrackerGeometry::ModuleType::Ph2PSP or det_type == TrackerGeometry::ModuleType::Ph2PSS) {
        moduletype = 1;
      } else if (det_type == TrackerGeometry::ModuleType::Ph2SS) {
        moduletype = 0;
      } else {
        // FIXME: raise exception : we should never be here
      }
      // container for clusters in stacked module, to be sorted afterwards
      std::vector<StackedClus> digs_stack;
      DetSetClus::const_iterator it;
      // pair sensors if there are clusters for both
      if (tTopo_->isLower((*sensor)->detId()) == 1) {
        // clusters for inner plane (P in case of PS)
        for (it = (*sensor)->begin(); it != (*sensor)->end(); it++) {
          digs_stack.push_back(StackedClus(it, LAYER_INNER, moduletype));
        }
        // If next digi is the corresponding outer plane : join them
        if ((sensor + 1) != cluss.end() and (int)((*(sensor + 1))->detId()) == (int)tTopo_->partnerDetId(detid)) {
          sensor++;
          for (it = (*sensor)->begin(); it != (*sensor)->end(); it++) {
            digs_stack.push_back(StackedClus(it, LAYER_OUTER, moduletype));
          }
        }
      } else {
        // clusters from outer plane (S in case of PS)
        for (it = (*sensor)->begin(); it != (*sensor)->end(); it++) {
          digs_stack.push_back(StackedClus(it, LAYER_OUTER, moduletype));
        }
      }
      // here we:
      // - sort all clusters -- order according to StackedClus::operator<(),
      //      so lower-left, lower-right, upper-left, upper-right sensor.
      // - divide big clusters into 8-strips parts
      // - count clusters on each side/layer (concentrator)
      // - remove truncated clusters

      // IAN: This writes one header per active DTC input channel,
      // rather than one per CIC chip, as in latest DTC DAQ format.

      std::sort(digs_stack.begin(), digs_stack.end());
      std::pair<int, int> nums = SortExpandAndLimitClusters(digs_stack, MAX_NS, MAX_NP);
      // - write appropriate header
      writeFeHeaderSparsified(fedbuffer, bitindex, moduletype, nums.second, nums.first);
      // - write the clusters
      std::vector<StackedClus>::iterator its;
      //std::cout<<"=== CHECK ORDER === "<<std::endl;
      for (its = digs_stack.begin(); its != digs_stack.end(); its++) {
        //  std::cout<<"    CLUS "<<" "<<its->getModuleType()<<" "<<its->getSide()<<" "<<its->getLayer()<<" "<<its->getChipId()<<std::endl;
        writeCluster(fedbuffer, bitindex, *its);
      }

    }  // end stack (FE) loop
    // add daq trailer
    fedbuffer.push_back(*(uint64_t*)FedDaqTrailer_.data());
    return fedbuffer;
  }

  void Phase2TrackerDigiToRaw::writeFeHeaderSparsified(
      std::vector<uint64_t>& buffer, uint64_t& bitpointer, int modtype, int np, int ns) {
    uint8_t length = 0;
    uint16_t header = ((uint16_t)ns & 0x3F);
    // module type switch
    if (modtype == 1) {
      header |= ((uint16_t)np & 0x3F) << 6;
      header |= ((uint16_t)modtype & 0x01) << 12;
      length = 13;
    } else {
      header |= ((uint16_t)modtype & 0x01) << 6;
      length = 7;
    }
    write_n_at_m(buffer, length, bitpointer, header);
    bitpointer += length;
  }

  // layer = 0 for inner, 1 for outer (-1 if irrelevant)
  void Phase2TrackerDigiToRaw::writeCluster(std::vector<uint64_t>& buffer,
                                            uint64_t& bitpointer,
                                            const StackedClus& digi) {
    if (digi.getModuleType() == 0) {
      // 2S module
      writeSCluster(buffer, bitpointer, digi, false);
    } else {
      // PS module
      if (digi.getLayer() == LAYER_INNER) {
        writePCluster(buffer, bitpointer, digi);
      } else {
        writeSCluster(buffer, bitpointer, digi, true);
      }
    }
  }

  void Phase2TrackerDigiToRaw::writeSCluster(std::vector<uint64_t>& buffer,
                                             uint64_t& bitpointer,
                                             const StackedClus& digi,
                                             bool threshold) {
    int csize = 15;
    uint16_t scluster = (digi.getChipId() & 0x0F) << 11;
    scluster |= (digi.getRawX() & 0xFF) << 3;
    scluster |= ((digi.getSizeX() - 1) & 0x07);
    if (threshold) {
      csize += 1;
      scluster <<= 1;
      scluster |= (digi.getThreshold() & 0x01);
    }
    write_n_at_m(buffer, csize, bitpointer, scluster);
    bitpointer += csize;
// debug
#ifdef EDM_ML_DEBUG
    std::ostringstream ss;
    ss << "S chip: " << digi.getChipId() << " digiX: " << digi.getDigiX() << " raw size: " << digi.getSizeX()
       << " digiY: " << digi.getDigiY() << " Layer: " << digi.getLayer();
    LogTrace("Phase2TrackerDigiProducer") << ss.str();
    ss.clear();
    ss.str("");
#endif
  }

  void Phase2TrackerDigiToRaw::writePCluster(std::vector<uint64_t>& buffer,
                                             uint64_t& bitpointer,
                                             const StackedClus& digi) {
    uint32_t pcluster = (digi.getChipId() & 0x0F) << 14;
    pcluster |= (digi.getRawX() & 0x7F) << 7;
    pcluster |= (digi.getRawY() & 0x0F) << 3;
    pcluster |= ((digi.getSizeX() - 1) & 0x07);
    write_n_at_m(buffer, 18, bitpointer, pcluster);
    bitpointer += 18;
// debug
#ifdef EDM_ML_DEBUG
    std::ostringstream ss;
    ss << "P chip: " << digi.getChipId() << " digiX: " << digi.getDigiX() << " raw size: " << digi.getSizeX()
       << " digiY: " << digi.getDigiY();
    LogTrace("Phase2TrackerDigiProducer") << ss.str();
    ss.clear();
    ss.str("");
#endif
  }
}  // namespace Phase2Tracker
