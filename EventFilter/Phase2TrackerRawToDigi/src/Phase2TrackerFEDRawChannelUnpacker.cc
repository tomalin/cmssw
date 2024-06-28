#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDRawChannelUnpacker.h"

namespace Phase2Tracker {

  Phase2TrackerFEDRawChannelUnpacker::Phase2TrackerFEDRawChannelUnpacker(const Phase2TrackerFEDChannel& channel)
      : data_(channel.data()),
        currentOffset_(channel.offset()),
        currentStrip_(0),
        valuesLeft_((channel.length()) * 8 - STRIPS_PADDING) {}

  Phase2TrackerFEDRawChannelUnpacker& Phase2TrackerFEDRawChannelUnpacker::operator++() {
    currentStrip_++;
    valuesLeft_--;
    return (*this);
  }

  Phase2TrackerFEDRawChannelUnpacker& Phase2TrackerFEDRawChannelUnpacker::operator++(int) {
    ++(*this);
    return *this;
  }

}  // namespace Phase2Tracker
