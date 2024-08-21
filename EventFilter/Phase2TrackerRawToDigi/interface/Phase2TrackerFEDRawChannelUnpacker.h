#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDRawChannelUnpacker_H  // {
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDRawChannelUnpacker_H

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include <cstdint>

namespace Phase2Tracker {

  // unpacker for RAW CBC data
  // each bit of the channel is related to one strip
  class Phase2TrackerFEDRawChannelUnpacker {
  public:
    Phase2TrackerFEDRawChannelUnpacker(const Phase2TrackerFEDChannel& channel);
    uint8_t stripIndex() const { return currentStrip_; }
    bool stripOn() const {
      return bool(static_cast<uint8_t>(read_n_at_m_L2R(data_, 1, currentOffset_ * 8 + currentStrip_)));
    }
    bool hasData() const { return valuesLeft_; }
    Phase2TrackerFEDRawChannelUnpacker& operator++();
    Phase2TrackerFEDRawChannelUnpacker& operator++(int);

  private:
    const uint8_t* data_;
    uint16_t currentOffset_;
    uint8_t currentStrip_;
    uint16_t valuesLeft_;
  };  // end of Phase2TrackerFEDRawChannelUnpacker

}  // namespace Phase2Tracker

#endif  // } end def EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDRawChannelUnpacker_H
