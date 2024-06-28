#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDZSChannelUnpacker.h"

namespace Phase2Tracker {

  Phase2TrackerFEDZSChannelUnpacker::Phase2TrackerFEDZSChannelUnpacker(const Phase2TrackerFEDChannel& channel)
      : data_(channel.data()) {
    currentOffset_ = channel.offset() * 8 + channel.bitoffset();
  }

  void Phase2TrackerFEDZSChannelUnpacker::Merge() {
    clusterx_ = unMergedX();
    clustery_ = unMergedY();
    clustersize_ = unMergedSize();
    while (gluedToNextCluster()) {
      advance();
      clustersize_ += unMergedSize();
    }
  }

  Phase2TrackerFEDZSChannelUnpacker& Phase2TrackerFEDZSChannelUnpacker::operator++() {
    currentOffset_ += clusterdatasize_;
    clustersLeft_--;
    return (*this);
  }

  Phase2TrackerFEDZSChannelUnpacker& Phase2TrackerFEDZSChannelUnpacker::advance() {
    currentOffset_ += clusterdatasize_;
    clustersLeft_--;
    return (*this);
  }

  Phase2TrackerFEDZSChannelUnpacker& Phase2TrackerFEDZSChannelUnpacker::operator++(int) {
    ++(*this);
    return *this;
  }

  int Phase2TrackerFEDZSSon2SChannelUnpacker::unMergedX() const {
    uint8_t id = chipId();
    if (id > 7)
      id -= 8;
    return (STRIPS_PER_CBC * id + rawX()) / 2;
  }

  int Phase2TrackerFEDZSSon2SChannelUnpacker::unMergedY() const { return (chipId() > 7) ? 1 : 0; }

  bool Phase2TrackerFEDZSSon2SChannelUnpacker::gluedToNextCluster() const {
    uint8_t size = rawSize();
    if (clustersLeft_ <= 1 or next().unMergedY() != unMergedY() or next().Plane() != Plane())
      return false;
    if (next().unMergedX() == unMergedX() + size) {
      if (size == 8 or (unMergedX() + size) % (STRIPS_PER_CBC / 2) == 0)
        return true;
      // else : two glued clusters split for not apparent reason
    }
    return false;
  }

  Phase2TrackerFEDZSSon2SChannelUnpacker Phase2TrackerFEDZSSon2SChannelUnpacker::next() const {
    Phase2TrackerFEDZSSon2SChannelUnpacker next(*this);
    next.advance();
    return next;
  }

  int Phase2TrackerFEDZSSonPSChannelUnpacker::unMergedX() const {
    uint8_t id = chipId();
    if (id > 7)
      id -= 8;
    return PS_ROWS * id + rawX();
  }

  int Phase2TrackerFEDZSSonPSChannelUnpacker::unMergedY() const { return (chipId() > 7) ? 1 : 0; }

  Phase2TrackerFEDZSSonPSChannelUnpacker Phase2TrackerFEDZSSonPSChannelUnpacker::next() const {
    Phase2TrackerFEDZSSonPSChannelUnpacker next(*this);
    next.advance();
    return next;
  }

  bool Phase2TrackerFEDZSSonPSChannelUnpacker::gluedToNextCluster() const {
    uint8_t size = rawSize();
    if (clustersLeft_ <= 1 or next().unMergedY() != unMergedY()) {
      return false;
    }
    if (next().unMergedX() == unMergedX() + size) {
      if (size == 8 or (unMergedX() + size) % PS_ROWS == 0) {
        return true;
      }
    }
    return false;
  }

  int Phase2TrackerFEDZSPonPSChannelUnpacker::unMergedX() const {
    uint8_t id = chipId();
    if (id > 7)
      id -= 8;
    return PS_ROWS * id + rawX();
  }

  int Phase2TrackerFEDZSPonPSChannelUnpacker::unMergedY() const {
    return (chipId() > 7) ? (rawY() + PS_COLS / 2) : rawY();
  }

  Phase2TrackerFEDZSPonPSChannelUnpacker Phase2TrackerFEDZSPonPSChannelUnpacker::next() const {
    Phase2TrackerFEDZSPonPSChannelUnpacker next(*this);
    next.advance();
    return next;
  }

  bool Phase2TrackerFEDZSPonPSChannelUnpacker::gluedToNextCluster() const {
    uint8_t size = rawSize();
    if (clustersLeft_ <= 1 or next().unMergedY() != unMergedY()) {
      return false;
    }
    if (next().unMergedX() == unMergedX() + size) {
      if (size == 8 or (unMergedX() + size) % PS_ROWS == 0)
        return true;
    }
    return false;
  }

}  // namespace Phase2Tracker
