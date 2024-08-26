#ifndef EventFilter_Phase2TrackerRawToDigi_FEDChannel_H
#define EventFilter_Phase2TrackerRawToDigi_FEDChannel_H

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include <cstdint>

// In zero-suppressed mode, each DTC input channel organised by raw2digi
// conversion into 4 sub-channels, for p-left, p-right, s-left, s-right sensors.
// This class represents a single one of these sub-channels.
// (Though in RAW data mode, this class represents a single CBC chip).

namespace Phase2Tracker {

  class Phase2TrackerFEDBuffer;

  class Phase2TrackerFEDChannel {
  public:
    // Give location of FED buffer, offset to this channel & size of this channel (in bytes).
    Phase2TrackerFEDChannel(const uint8_t* const data,
                            const size_t offset,
                            const uint16_t length,
                            const uint8_t bitoffset = 0,
                            const DET_TYPE dettype = UNUSED)
        : data_(data), offset_(offset), length_(length), bitoffset_(bitoffset), dettype_(dettype) {}

    //gets length from first 2 bytes (assuming normal FED channel)
    Phase2TrackerFEDChannel(const uint8_t* const data, const size_t offset);
    // Location of FED buffer (bytes)
    const uint8_t* data() const { return data_; }
    // Offset of this channel inside buffer (bytes)
    size_t offset() const { return offset_; }
    //  Length of data for this channel (bytes). Is zero if channel not connected.
    uint16_t length() const { return length_; }
    uint16_t bitoffset() const { return bitoffset_; }
    DET_TYPE dettype() const { return dettype_; }

  private:
    friend class Phase2TrackerFEDBuffer;
    //third byte of channel data for normal FED channels
    uint8_t packetCode() const;
    const uint8_t* data_;
    size_t offset_;
    uint16_t length_;
    uint16_t bitoffset_;
    DET_TYPE dettype_;
  };  // end Phase2TrackerFEDChannel class

}  // namespace Phase2Tracker

#endif  // } end def EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDChannel_H
