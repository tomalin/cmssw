#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Phase2Tracker {

  // implementation of Phase2TrackerFEDBuffer
Phase2TrackerFEDBuffer::Phase2TrackerFEDBuffer(const uint8_t* fedBuffer, const size_t fedBufferSize, const std::vector<bool>& connectedInputs)
  : buffer_(fedBuffer), bufferSize_(fedBufferSize), valid_(1) {
    LogTrace("Phase2TrackerFEDBuffer") << "[Phase2Tracker::Phase2TrackerFEDBuffer::" << __func__ << "] "
                                       << "\n";
    LogTrace("Phase2TrackerFEDBuffer") << "content of buffer with size: " << int(fedBufferSize) << std::endl;
    for (size_t i = 0; i < fedBufferSize; i += 8) {
      uint64_t word = read64(i, buffer_);
      LogTrace("Phase2TrackerFEDBuffer") << " word " << std::setfill(' ') << std::setw(2) << i / 8 << " | " << std::hex
                                         << std::setw(16) << std::setfill('0') << word << std::dec << std::endl;
    }
    LogTrace("Phase2TrackerFEDBuffer") << std::endl;
    // reserve all channels to avoid vector reservation updates (should be 16x16 in our case)
    channels_.reserve(MAX_FE_PER_FED * MAX_CBC_PER_FE);
    // first 64 bits word is for DAQ header
    daqHeader_ = FEDDAQHeader(buffer_);
    // last 64 bit word is daq trailer
    daqTrailer_ = FEDDAQTrailer(buffer_ + bufferSize_ - 8);
    // tracker header follows daq header
    trackerHeader_ = Phase2TrackerFEDHeader(buffer_ + 8);
    valid_ = trackerHeader_.isValid();
    // get pointer to payload
    payloadPointer_ = getPointerToPayload();
    // fill list of Phase2TrackerFEDChannels and get pointers to trigger and comissioning data
    findChannels( connectedInputs );
  }

  Phase2TrackerFEDBuffer::~Phase2TrackerFEDBuffer() {}

  void Phase2TrackerFEDBuffer::findChannels(const std::vector<bool>& connectedInputs) {
    // Each FED can be connectd to up to 72 frontends (read from header).
    // Each fronted can be connected to up to 16 CBC.
    // In unsparsified (raw) mode, a header of 16bits tells which CBC are activated on this FE. One channel corresponds to one CBC chip.
    // In sparsified (ZS) mode, the header tells how many P and S clusters to expect for this FE. Each cluster contains the CBC ID at the beginning. One channel corresponds to one end of one sensor in a tracker module.

    // offset of beginning of current channel
    size_t offsetBeginningOfChannel = 0;

    // Note which DTC input channels received at least one cluster,
    // and so have entries in the RAW data buffer.
    const std::vector<bool> status = trackerHeader_.frontendStatus();

    if (readoutMode() == READOUT_MODE_PROC_RAW) {
      for (int ife = 0; ife < MAX_FE_PER_FED; ife++) { // Loop inputs of one DTC
        // Fill channels & advance pointer to end of channel if ...
        // IAN: Need to decide if want (1) or (2) ...
        // (1) DTC input is connected to a tracker module
        // if (connectedInputs[ife]) {
        // (2) DTC input receives at least one digi in this event.
        if (status[ife]) {
        // read CBC status bits and advance pointer to after them
          uint16_t cbc_status =
              static_cast<uint32_t>(read_n_at_m_L2R(payloadPointer_, 16, offsetBeginningOfChannel * 8));
          offsetBeginningOfChannel += MAX_CBC_PER_FE / 8;  // counts bytes
          for (int i = 0; i < MAX_CBC_PER_FE; i++) {
            // if CBC is ON, fill channel and advance pointer. else, push back empty channel
            if ((cbc_status >> i) & 0x1) {
              // Warning: STRIPS_PADDING+STRIPS_PER_CBC should always be an entire number of bytes
              channels_.push_back(Phase2TrackerFEDChannel(
                  payloadPointer_, offsetBeginningOfChannel, (STRIPS_PADDING + STRIPS_PER_CBC) / 8));
              offsetBeginningOfChannel += (STRIPS_PADDING + STRIPS_PER_CBC) / 8; // bytes
            } else {
              // Dead CBC. Fill with null channel.
              channels_.push_back(Phase2TrackerFEDChannel(nullptr, 0, 0));
            }
          }
        } else {
          // Disconected DTC channel. Fill MAX_CBC_PER_FE null channels
          channels_.insert(channels_.end(), size_t(MAX_CBC_PER_FE), Phase2TrackerFEDChannel(nullptr, 0, 0));
        }
      }
      
    } else if (readoutMode() == READOUT_MODE_ZERO_SUPPRESSED) {
      // save current bit index
      int bitOffset = 0;

      // loop on input channels of this DTC
      for (int ife = 0; ife < MAX_FE_PER_FED; ife++) { // Loop inputs of one DTC

        // Fill channels & advance pointer to end of channel if ...
        // IAN: Need to decide if want (1) or (2) ...
        // (1) DTC input is connected to a tracker module
        // if (connectedInputs[ife]) {
        // (2) DTC input receives at least one digi in this event.
        if (status[ife]) {
          // Read FE header :
          // FE header is : P/S (1bit) + #P clusters (5bits if P, 0 otherwise) + #S clusters (5bits)
          // read type of module (2S/PS)
          uint8_t mod_type = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 1, bitOffset));
          // note the index of the next channel to fill
          int ichan = channels_.size();
          // add the proper number of channels for this FE (4 channels per module, 1 for each side of S and P)
          channels_.insert(channels_.end(), size_t(4), Phase2TrackerFEDChannel(payloadPointer_, 0, 0));
          // read header
          uint8_t num_p, num_s;
          int s_cluster_size_bits = Son2S_CLUSTER_SIZE_BITS;
          if (mod_type == 0) {
            num_p = 0;
            num_s = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 6, bitOffset + 1));
            bitOffset += 7;
          } else {
            // s clusters on PS modules have an extra MIPS threshold bit
            s_cluster_size_bits = SonPS_CLUSTER_SIZE_BITS;
            num_p = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 6, bitOffset + 1));
            num_s = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 6, bitOffset + 7));
            bitOffset += 13;
          }
          // Loop over clusters in this stacked module, for p & s type sensors,
          int iCBC;
          int iOffset = bitOffset;
          // 0 = P-left-end, 1 = P-right-end, 2 = S-left-end, 3 = S-right-end
          int chansize_0 = 0, chansize_1 = 0, chansize_2 = 0, chansize_3 = 0;
          for (int i = 0; i < num_p; i++) {
            iCBC = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 4, bitOffset + 14));
            if (iCBC < 8) { // which end of stacked module
              chansize_0 += P_CLUSTER_SIZE_BITS;
            } else {
              chansize_1 += P_CLUSTER_SIZE_BITS;
            }
            bitOffset += P_CLUSTER_SIZE_BITS;
          }
          for (int i = 0; i < num_s; i++) {
            iCBC = static_cast<uint8_t>(read_n_at_m_L2R(payloadPointer_, 4, bitOffset + 12));
            if (iCBC < 8) {
              chansize_2 += s_cluster_size_bits;
            } else {
              chansize_3 += s_cluster_size_bits;
            }
            bitOffset += s_cluster_size_bits;
          }
          
          // N.B. Raw2Digi conversion assumes digis on any given DTC input
          // channel sorted (by StackedClus::operator<()) in order
          // lower-left, lower-right, upper-left, upper-right sensor.
          // where lower is p-type in PS modules.
          // This order assumed here.
          
          // p sensor (left & right end) -- size 0 if 2S module.
          channels_[ichan + 0] =
              Phase2TrackerFEDChannel(payloadPointer_, iOffset / 8, (chansize_0 + 8 - 1) / 8, iOffset % 8, DET_PonPS);
          iOffset += chansize_0;
          channels_[ichan + 1] =
              Phase2TrackerFEDChannel(payloadPointer_, iOffset / 8, (chansize_1 + 8 - 1) / 8, iOffset % 8, DET_PonPS);
          iOffset += chansize_1;
          // s sensor (left & right end)
          DET_TYPE det_type = (mod_type == 0) ? DET_Son2S : DET_SonPS;
          channels_[ichan + 2] =
              Phase2TrackerFEDChannel(payloadPointer_, iOffset / 8, (chansize_2 + 8 - 1) / 8, iOffset % 8, det_type);
          iOffset += chansize_2;
          channels_[ichan + 3] =
              Phase2TrackerFEDChannel(payloadPointer_, iOffset / 8, (chansize_3 + 8 - 1) / 8, iOffset % 8, det_type);
          iOffset += chansize_3;
        } else {
          // Inactive DTC input channel. Add 4 null channels for it.
          channels_.insert(channels_.end(), size_t(4), Phase2TrackerFEDChannel(nullptr, 0, 0));
        }
      }
      // compute byte offset for payload
      offsetBeginningOfChannel = (bitOffset + 8 - 1) / 8;
    }
    
    // round the offset to the next 64 bits word
    // IAN: No option to switch off reading stub data???
    //         Why doesn't this crash if it's not present?
    std::cout<<"READING STUB DATA???"<<std::endl;
    
    int words64 = (offsetBeginningOfChannel + 8 - 1) / 8;  // size in 64 bit
    int payloadSize = words64 * 8;                         // size in bytes
    triggerPointer_ = (uint8_t*)(payloadPointer_ + payloadSize);

    // register stub channels
    // save current bit index
    int bitOffset = 0;
    for (int ife = 0; ife < MAX_FE_PER_FED; ife++) { // Loop inputs of one DTC
      // if the current fronted is on, fill channels and advance pointer to end of channel
      if (connectedInputs[ife]) {
        // read FE stub header (6 bits)
        uint8_t nstubs = static_cast<uint8_t>(read_n_at_m_L2R(triggerPointer_, 5, bitOffset));
        uint8_t modtype = static_cast<uint8_t>(read_n_at_m_L2R(triggerPointer_, 1, bitOffset + 5));
        bitOffset += 6;
        // set module type and data size
        DET_TYPE det_type = (modtype == 0) ? DET_Son2S : DET_PonPS;
        int stub_size = (modtype == 0) ? STUBS_SIZE_2S : STUBS_SIZE_PS;
        // add channel
        stub_channels_.push_back(Phase2TrackerFEDChannel(
            triggerPointer_, bitOffset / 8, (stub_size * nstubs + 8 - 1) / 8, bitOffset % 8, det_type));
        bitOffset += stub_size * nstubs;
      } else {
        // else add a null channel, don't advance the channel pointer
        stub_channels_.push_back(Phase2TrackerFEDChannel(triggerPointer_, 0, 0));
      }
    }
    // trigger size rounded up to 64 bit words
    int triggerSize = (((bitOffset + 8 - 1) / 8 + 8 - 1) / 8) * 8;

    // get diff size in bytes:
    // fedBufferSize - (DAQHeader+TrackHeader+PayloadSize+triggerSize+DAQTrailer)
    int bufferDiff = bufferSize_ - 8 - trackerHeader_.getTrackerHeaderSize() - payloadSize - triggerSize - 8;

    // check if condition data is supposed to be there:
    if (trackerHeader_.getConditionData()) {
      condDataPointer_ = payloadPointer_ + payloadSize + triggerSize;
      // diff must be equal to condition data size
      if (bufferDiff <= 0) {
        std::ostringstream ss;
        ss << "[Phase2Tracker::Phase2TrackerFEDBuffer::" << __func__ << "] "
           << "\n";
        ss << "WARNING: Skipping FED buffer: "
           << "\n";
        ss << "Cause: FED Buffer Size does not match data => missing condition data? : "
           << "\n";
        ss << "Expected Buffer Size " << bufferSize_ << " bytes"
           << "\n";
        ss << "Computed Buffer Size " << bufferSize_ - bufferDiff << " bytes"
           << "\n";
        LogTrace("Phase2TrackerFEDBuffer") << ss.str() << std::endl;
        valid_ = 0;
      }
    }
  }
  /*
  std::map<uint32_t, uint32_t> Phase2TrackerFEDBuffer::conditionData() const {
    std::map<uint32_t, uint32_t> cdata;
    // check if there is condition data
    if (condDataPointer_) {
      const uint8_t* pointer = condDataPointer_;
      const uint8_t* stop = buffer_ + bufferSize_ - 8;
      // first read the size
      uint32_t size = 0;
      // somehow the size is not inverted
      //for (int i=0;i<4;++i) size += *(pointer-4+(i^7)) << (i*8);
      size = *reinterpret_cast<const uint32_t*>(pointer);
      LogTrace("Phase2TrackerFEDBuffer") << "Condition Data size = " << size << std::endl;
      pointer += 8;
      // now the conditions
      while (pointer < stop) {
        // somehow the data is not inverted
        uint32_t data = 0;
        //for (int i = 0, j=3 ; i<4; i++,j--)
        //{ data += (*(pointer+i) << j*8); }
        data = *reinterpret_cast<const uint32_t*>(pointer);
        pointer += 4;

        uint32_t key = 0;
        for (int i = 0, j = 3; i < 4; i++, j--) {
          key += (*(pointer + i) << j * 8);
        }
        pointer += 4;

        cdata[key] = data;
      }
      // final check: cdata size == size
      if (cdata.size() != size) {
        std::ostringstream ss;
        ss << "[Phase2Tracker::Phase2TrackerFEDBuffer::"<<__func__<<"] " << "\n";
        ss << "WARNING: Skipping FED buffer: " << "\n";
        ss << "Cause: FED Buffer Size does not match data => corrupted buffer? : " << "\n";
        ss << "Expected Buffer Size " << bufferSize_ << " bytes" << "\n";
        ss << "Computed Buffer Size " << bufferSize_ - bufferDiff << " bytes" << "\n";
				LogTrace("Phase2TrackerFEDBuffer") << ss.str() << std::endl;
        valid_ = 0;
      }
    } 

  }
  */
  std::map<uint32_t, uint32_t> Phase2TrackerFEDBuffer::conditionData() {
    std::map<uint32_t, uint32_t> cdata;
    // check if there is condition data
    if (condDataPointer_) {
      const uint8_t* pointer = condDataPointer_;
      const uint8_t* stop = buffer_ + bufferSize_ - 8;
      // first read the size
      uint32_t size = 0;
      // somehow the size is not inverted
      //for (int i=0;i<4;++i) size += *(pointer-4+(i^7)) << (i*8);
      size = *reinterpret_cast<const uint32_t*>(pointer);
      LogTrace("Phase2TrackerFEDBuffer") << "Condition Data size = " << size << std::endl;
      pointer += 8;
      // now the conditions
      while (pointer < stop) {
        // somehow the data is not inverted
        uint32_t data = 0;
        //for (int i = 0, j=3 ; i<4; i++,j--)
        //{ data += (*(pointer+i) << j*8); }
        data = *reinterpret_cast<const uint32_t*>(pointer);
        pointer += 4;

        uint32_t key = 0;
        for (int i = 0, j = 3; i < 4; i++, j--) {
          key += (*(pointer + i) << j * 8);
        }
        pointer += 4;

        cdata[key] = data;
      }
      // final check: cdata size == size
      if (cdata.size() != size) {
        std::ostringstream ss;
        ss << "[Phase2Tracker::Phase2TrackerFEDBuffer::" << __func__ << "] "
           << "\n";
        ss << "Number of condition data does not match the announced value"
           << "\n";
        ss << "Expected condition data Size " << size << " entries"
           << "\n";
        ss << "Computed condition data Size " << cdata.size() << " entries"
           << "\n";
        LogTrace("Phase2TrackerFEDBuffer") << ss.str() << std::endl;
        valid_ = 0;
      }
    }
    return cdata;
  }

  FEDReadoutMode Phase2TrackerFEDBuffer::readoutMode() const { return trackerHeader_.getReadoutMode(); }

}  // namespace Phase2Tracker
