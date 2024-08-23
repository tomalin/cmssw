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

  std::pair<int, int> SortExpandAndLimitClusters(std::vector<StackedDigi>& digis, int max_ns, int max_np) {
    std::vector<StackedDigi> processed;
    // number of clusters allowed : P-left, P-right, S-left, S-right
    int roomleft[4] = {max_np, max_np, max_ns, max_ns};
    // fill left and right vectors, expand big clusters
    for (auto dig = digis.begin(); dig < digis.end(); dig++) {
      std::vector<StackedDigi> parts = dig->splitDigi();
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
    digis.swap(processed);
    // return number of S and P clusters respectively
    return std::make_pair(2 * max_ns - roomleft[2] - roomleft[3], 2 * max_np - roomleft[0] - roomleft[1]);
  }

  Phase2TrackerDigiToRaw::Phase2TrackerDigiToRaw(const Phase2TrackerCabling* cabling,
                                                 const TrackerGeometry* tGeom,
                                                 const TrackerTopology* tTopo,
                                                 std::map<int, std::pair<int, int>> stackMap,
                                                 edm::Handle<edmNew::DetSetVector<Phase2TrackerCluster1D>> digis_handle,
                                                 int mode)
      : cabling_(cabling),
        tTopo_(tTopo),
        tGeom_(tGeom),
        stackMap_(stackMap),
        digishandle_(digis_handle),
        mode_(mode),
        FedDaqHeader_(0, 0, 0, DAQ_EVENT_TYPE_SIMULATED),  // TODO : add L1ID
        FedDaqTrailer_(0, 0) {
    FedHeader_.setDataFormatVersion(2);
    FedHeader_.setDebugMode(SUMMARY);
    FedHeader_.setEventType((uint8_t)0x04);
  }

  void Phase2TrackerDigiToRaw::buildFEDBuffers(std::unique_ptr<FEDRawDataCollection>& rcollection) {
    // store digis for a given fedid
    std::vector<edmNew::DetSet<Phase2TrackerCluster1D>> digis_t;
    // store active connections for a given fedid
    std::vector<bool> festatus(72, false);
    // iterate on all possible channels
    Phase2TrackerCabling::cabling conns = cabling_->orderedConnections(0);
    Phase2TrackerCabling::cabling::const_iterator iconn = conns.begin(), end = conns.end(), icon2;
    while (iconn != end) {
      unsigned int fedid = (*iconn)->getCh().first;
      // Loop over subset of modules belonging to same DTC.
      for (icon2 = iconn; icon2 != end && (*icon2)->getCh().first == fedid; icon2++) {
        int detidStack =  (*icon2)->getDetid() ;        
        // FIXME (when proper cabling exists) : because we use test cabling, we have some detids set to 0 : we should ignore them
        // Note from Ian. DummyCablingTxt_cfi.py sets detid = 0 for any unconnected DTC channels. If we skip them here, why does it bother?
        if (detidStack == 0)
          continue;
        if (stackMap_.find(detidStack) == stackMap_.end()) edm::LogError("Phase2TrackerDigiProducer") << "DetID "<< detidStack << " in cabling map is unknown to Tracker geometry";
        // detid of first plane is detid of module + 1
        // FIXME (when it is implemented) : we should get this from the topology / geometry
        unsigned int detid = stackMap_[detidStack].first;
        // end of fixme
        unsigned int fedCh = (*icon2)->getCh().second;
        edmNew::DetSetVector<Phase2TrackerCluster1D>::const_iterator digis;
        digis = digishandle_->find(detid);
        if (digis != digishandle_->end()) {
          digis_t.push_back(*digis);
          festatus[fedCh] = true;
        }
        digis = digishandle_->find(tTopo_->partnerDetId(detid));
        if (digis != digishandle_->end()) {
          digis_t.push_back(*digis);
          festatus[fedCh] = true;
        }
      }
      // save buffer
      FedHeader_.setFrontendStatus(festatus);
      // write digis of one DTC to buffer
      std::vector<uint64_t> fedbuffer = makeBuffer(digis_t);
      FEDRawData& frd = rcollection->FEDData(fedid);
      int size = fedbuffer.size() * 8;
      frd.resize(size);
      memcpy(frd.data(), &fedbuffer[0], size);
      festatus.assign(72, false);
      digis_t.clear();
      // advance connections pointer
      iconn = icon2;
    }
  }

  // Build buffer for a single DTC,
  // consisting of FEDDAQHeader, FEDHeader & for each DTC input channel
  // the feHeaderSparsified & all the clusters.

  std::vector<uint64_t> Phase2TrackerDigiToRaw::makeBuffer(std::vector<edmNew::DetSet<Phase2TrackerCluster1D>> digis) {
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
    // looping on stacked modules in this DTC.
    std::vector<edmNew::DetSet<Phase2TrackerCluster1D>>::const_iterator stack;
    for (stack = digis.begin(); stack != digis.end(); stack++) {
      // get id of stack
      unsigned int detid = stack->detId();
      TrackerGeometry::ModuleType det_type = tGeom_->getDetectorType(detid);
      if (det_type == TrackerGeometry::ModuleType::Ph2PSP or det_type == TrackerGeometry::ModuleType::Ph2PSS) {
        moduletype = 1;
      } else if (det_type == TrackerGeometry::ModuleType::Ph2SS) {
        moduletype = 0;
      } else {
        // FIXME: raise exception : we should never be here
      }
      // container for digis in stacked module, to be sorted afterwards
      std::vector<StackedDigi> digs_stack;
      edmNew::DetSet<Phase2TrackerCluster1D>::const_iterator it;
      // pair modules if there are digis for both
      if (tTopo_->isLower(stack->detId()) == 1) {
        // digis for inner plane (P in case of PS)
        if ((stack + 1) != digis.end() and (int)(stack + 1)->detId() == (int)tTopo_->partnerDetId(detid)) {
          // next digi is the corresponding outer plane : join them
          for (it = stack->begin(); it != stack->end(); it++) {
            digs_stack.push_back(StackedDigi(it, LAYER_INNER, moduletype));
          }
          stack++;
          for (it = stack->begin(); it != stack->end(); it++) {
            digs_stack.push_back(StackedDigi(it, LAYER_OUTER, moduletype));
          }
        } else {
          // next digi is from another module, only use this one
          for (it = stack->begin(); it != stack->end(); it++) {
            digs_stack.push_back(StackedDigi(it, LAYER_INNER, moduletype));
          }
        }
      } else {
        // digis from outer plane (S in case of PS)
        for (it = stack->begin(); it != stack->end(); it++) {
          digs_stack.push_back(StackedDigi(it, LAYER_OUTER, moduletype));
        }
      }
      // here we:
      // - sort all digis -- order according to StackedDigi::operator<(),
      //      so lower-left, lower-right, upper-left, upper-right sensor.
      // - divide big clusters into 8-strips parts
      // - count digis on each side/layer (concentrator)
      // - remove truncated digis

      // IAN: This seems to be writing one header per DTC input channel,
      // rather than one per CIC chip, as in latest DTC DAQ format.
      
      std::sort(digs_stack.begin(), digs_stack.end());
      std::pair<int, int> nums = SortExpandAndLimitClusters(digs_stack, MAX_NS, MAX_NP);
      // - write appropriate header
      writeFeHeaderSparsified(fedbuffer, bitindex, moduletype, nums.second, nums.first);
      // - write the digis
      std::vector<StackedDigi>::iterator its;
      //std::cout<<"=== CHECK ORDER === "<<std::endl;
      //for (its = digs_stack.begin(); its != digs_stack.end(); its++) {
      //  std::cout<<"    CLUS "<<" "<<its->getModuleType()<<" "<<its->getSide()<<" "<<its->getLayer()<<" "<<its->getChipId()<<std::endl;
      //  writeCluster(fedbuffer, bitindex, *its);
      //}

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
  void Phase2TrackerDigiToRaw::writeCluster(std::vector<uint64_t>& buffer, uint64_t& bitpointer, StackedDigi digi) {
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
                                             StackedDigi digi,
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

  void Phase2TrackerDigiToRaw::writePCluster(std::vector<uint64_t>& buffer, uint64_t& bitpointer, StackedDigi digi) {
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
