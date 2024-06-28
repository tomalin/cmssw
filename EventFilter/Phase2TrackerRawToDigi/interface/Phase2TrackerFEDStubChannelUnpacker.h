#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDStubChannelUnpacker_H  // {
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDStubChannelUnpacker_H

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include <stdint.h>

namespace Phase2Tracker {

  class Phase2TrackerFEDStubChannelUnpacker {
  public:
    Phase2TrackerFEDStubChannelUnpacker(const Phase2TrackerFEDChannel& channel);
    bool hasData() const { return (clustersLeft_ > 0); }
    Phase2TrackerFEDStubChannelUnpacker& operator++();
    Phase2TrackerFEDStubChannelUnpacker& operator++(int);
    virtual int stubX() const = 0;
    virtual int stubY() const = 0;

  protected:
    // raw methods (bitwise operations)
    virtual uint8_t rawX() const = 0;
    virtual uint8_t chipId() const = 0;
    const uint8_t* data_;
    uint16_t currentOffset_;  // caution : this is in bits, not bytes
    int clustersLeft_;
    uint8_t clusterdatasize_;
  };

  inline Phase2TrackerFEDStubChannelUnpacker::Phase2TrackerFEDStubChannelUnpacker(const Phase2TrackerFEDChannel& channel)
      : data_(channel.data()) {
    currentOffset_ = channel.offset() * 8 + channel.bitoffset();
  }

  inline Phase2TrackerFEDStubChannelUnpacker& Phase2TrackerFEDStubChannelUnpacker::operator++() {
    currentOffset_ += clusterdatasize_;
    clustersLeft_--;
    return (*this);
  }

  inline Phase2TrackerFEDStubChannelUnpacker& Phase2TrackerFEDStubChannelUnpacker::operator++(int) {
    ++(*this);
    return *this;
  }

  class Phase2TrackerFEDStubSon2SChannelUnpacker : public Phase2TrackerFEDStubChannelUnpacker {
  public:
    Phase2TrackerFEDStubSon2SChannelUnpacker(const Phase2TrackerFEDChannel& channel)
        : Phase2TrackerFEDStubChannelUnpacker(channel) {
      clusterdatasize_ = STUBS_SIZE_2S;
      clustersLeft_ = channel.length() * 8 / clusterdatasize_;
    }
    inline uint8_t chipId() const { return (uint8_t)read_n_at_m_l2r(data_, 4, currentOffset_); }
    inline uint8_t rawX() const { return (uint8_t)read_n_at_m_l2r(data_, 8, currentOffset_ + 4); }
    inline uint8_t Bend() const { return (uint8_t)read_n_at_m_l2r(data_, 4, currentOffset_ + 12); }
    int stubX() const;
    int stubY() const;
  };

  int Phase2TrackerFEDStubSon2SChannelUnpacker::stubX() const {
    uint8_t id = chipId();
    if (id > 7)
      id -= 8;
    return (STRIPS_PER_CBC * id + rawX() - 2);
  }

  int Phase2TrackerFEDStubSon2SChannelUnpacker::stubY() const { return (chipId() > 7) ? 1 : 0; }

  class Phase2TrackerFEDStubPonPSChannelUnpacker : public Phase2TrackerFEDStubChannelUnpacker {
  public:
    Phase2TrackerFEDStubPonPSChannelUnpacker(const Phase2TrackerFEDChannel& channel)
        : Phase2TrackerFEDStubChannelUnpacker(channel) {
      clusterdatasize_ = STUBS_SIZE_PS;
      clustersLeft_ = channel.length() * 8 / clusterdatasize_;
    }
    inline uint8_t chipId() const { return (uint8_t)read_n_at_m_l2r(data_, 4, currentOffset_); }
    inline uint8_t rawX() const { return (uint8_t)read_n_at_m_l2r(data_, 8, currentOffset_ + 4); }
    inline uint8_t Bend() const { return (uint8_t)read_n_at_m_l2r(data_, 4, currentOffset_ + 12); }
    inline uint8_t rawY() const { return (uint8_t)read_n_at_m_l2r(data_, 4, currentOffset_ + 16); }
    int stubX() const;
    int stubY() const;
  };

  int Phase2TrackerFEDStubPonPSChannelUnpacker::stubX() const {
    uint8_t id = chipId();
    if (id > 7)
      id -= 8;
    return PS_ROWS * id + rawX();
  }

  int Phase2TrackerFEDStubPonPSChannelUnpacker::stubY() const {
    return (chipId() > 7) ? (rawY() + PS_COLS / 2) : rawY();
  }

}  // namespace Phase2Tracker

#endif  // } end def EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDStubChannelUnpacker_H
