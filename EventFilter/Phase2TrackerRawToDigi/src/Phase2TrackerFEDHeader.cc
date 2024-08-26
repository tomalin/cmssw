#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <algorithm>

namespace Phase2Tracker {

  Phase2TrackerFEDHeader::Phase2TrackerFEDHeader(const uint8_t* headerPointer)
      : trackerHeader_(headerPointer), valid_(1) {
    // make a local copy of header (first two words)
    memcpy(headercopy_, headerPointer, 16);
    memcpy(&header_first_word_, headerPointer, 8);
    // decode the Tracker Header and store info
    init();
  }

  void Phase2TrackerFEDHeader::init() {
    dataFormatVersion_ = dataFormatVersion();
    debugMode_ = debugMode();
    // WARNING: eventType must be called before
    // readoutMode, conditionData and dataType
    // as this info is stored in eventType
    eventType_ = eventType();
    readoutMode_ = readoutMode();
    conditionData_ = conditionData();
    dataType_ = dataType();

    glibStatusCode_ = glibStatusCode();
    // numberOfCBC must be called before pointerToData
    numberOfCBC_ = numberOfCBC();
    pointerToData_ = pointerToData();

    LogTrace("Phase2TrackerFEDBuffer") << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "]: \n"
                                       << " Tracker Header contents:\n"
                                       << "  -- Data Format Version : " << uint32_t(dataFormatVersion_) << "\n"
                                       << "  -- Debug Level         : " << debugMode_ << "\n"
                                       << "  -- Operating Mode      : " << readoutMode_ << "\n"
                                       << "  -- Condition Data      : " << (conditionData_ ? "Present" : "Absent")
                                       << "\n"
                                       << "  -- Data Type           : " << (dataType_ ? "Real" : "Fake") << "\n"
                                       << "  -- Glib Stat registers : " << std::hex << std::setw(16) << glibStatusCode_
                                       << "\n"
                                       << "  -- connected CBC       : " << std::dec << numberOfCBC_ << "\n";
  }

  uint8_t Phase2TrackerFEDHeader::dataFormatVersion() {
    uint8_t Version = static_cast<uint8_t>(read_n_at_m_R2L(headercopy_, VERSION_L, VERSION_S));
    if (Version != 2) {
      std::ostringstream ss;
      ss << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "] \n";
      ss << "WARNING: FED has been marked as invalid and will be skipped \n";
      ss << "Cause: Invalid Data Format Version in Tracker Header : ";
      printHex(&header_first_word_, 1, ss);
      LogTrace("Phase2TrackerFEDHeader") << ss.str() << std::endl;
      valid_ = 0;
    }
    return Version;
  }

  void Phase2TrackerFEDHeader::setDataFormatVersion(uint8_t version) {
    write_n_at_m_R2L(headercopy_, VERSION_L, VERSION_S, (uint64_t)version);
  }

  void Phase2TrackerFEDHeader::setDebugMode(READ_MODE mode) {
    write_n_at_m_R2L(headercopy_, HEADER_FORMAT_L, HEADER_FORMAT_S, (uint64_t)mode);
  }

  READ_MODE Phase2TrackerFEDHeader::debugMode() {
    // Read debugMode in Tracker Header
    uint8_t mode = static_cast<uint8_t>(read_n_at_m_R2L(headercopy_, HEADER_FORMAT_L, HEADER_FORMAT_S));

    switch (mode) {  // check if it is one of correct modes
      case SUMMARY:
        return READ_MODE(SUMMARY);
      case FULL_DEBUG:
        return READ_MODE(FULL_DEBUG);
      case CBC_ERROR:
        return READ_MODE(CBC_ERROR);
      default:  // else mark as invalid
        std::ostringstream ss;
        ss << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "] ";
        ss << "WARNING: Skipping FED ";
        ss << "Cause: Invalid Header Format in Tracker Header : ";
        printHex(&header_first_word_, 1, ss);
        LogTrace("Phase2TrackerFEDHeader") << ss.str() << std::endl;
        valid_ = 0;
    }

    return READ_MODE(READ_MODE_INVALID);
  }

  void Phase2TrackerFEDHeader::setEventType(uint8_t event_type) {
    write_n_at_m_R2L(headercopy_, EVENT_TYPE_L, EVENT_TYPE_S, (uint64_t)event_type);
  }

  uint8_t Phase2TrackerFEDHeader::eventType() const {
    return static_cast<uint8_t>(read_n_at_m_R2L(headercopy_, EVENT_TYPE_L, EVENT_TYPE_S));
  }

  // decode eventType_. Read: readoutMode, conditionData and dataType
  FEDReadoutMode Phase2TrackerFEDHeader::readoutMode() {
    // readout mode is first bit of event type
    uint8_t mode = static_cast<uint8_t>(eventType_ >> 2) & 0x3;

    switch (mode) {  // check if it is one of correct modes
      case 2:
        return FEDReadoutMode(READOUT_MODE_PROC_RAW);
      case 1:
        return FEDReadoutMode(READOUT_MODE_ZERO_SUPPRESSED);
      default:  // else create Exception
        std::ostringstream ss;
        ss << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "] ";
        ss << "WARNING: Skipping FED ";
        ss << "Cause: Invalid Readout Mode in Tracker Header : ";
        printHex(&header_first_word_, 1, ss);
        LogTrace("Phase2TrackerFEDHeader") << ss.str() << std::endl;
        valid_ = 0;
        return FEDReadoutMode(READOUT_MODE_INVALID);
    }
  }

  uint8_t Phase2TrackerFEDHeader::conditionData() const { return static_cast<uint8_t>(eventType_ >> 1) & 0x1; }

  uint8_t Phase2TrackerFEDHeader::dataType() const { return static_cast<uint8_t>(eventType_) & 0x1; }

  void Phase2TrackerFEDHeader::setGlibStatusCode(uint64_t status_code) {
    write_n_at_m_R2L(headercopy_, GLIB_STATUS_L, GLIB_STATUS_S, status_code);
  }

  uint64_t Phase2TrackerFEDHeader::glibStatusCode() const {
    return read_n_at_m_R2L(headercopy_, GLIB_STATUS_L, GLIB_STATUS_S);
  }

  std::vector<bool> Phase2TrackerFEDHeader::frontendStatus() const {
    uint8_t fe_status_0 = (uint8_t)(read_n_at_m_R2L(headercopy_, 8, 0));
    uint64_t fe_status_1 = (uint64_t)(read_n_at_m_R2L(headercopy_, 64, 64));
    std::vector<bool> status(72, false);
    for (int i = 0; i < 72; i++) {
      if (i < 64) {
        status[i] = (fe_status_1 >> i) & 0x1;
      } else {
        status[i] = (fe_status_0 >> (i - 64)) & 0x1;
      }
    }
    return status;
  }

  void Phase2TrackerFEDHeader::setFrontendStatus(std::vector<bool> status) {
    uint8_t fe_status_0 = 0x00;
    uint64_t fe_status_1 = 0x00;
    int index = 0;
    std::vector<bool>::iterator sti;
    for (sti = status.begin(); sti < status.end(); sti++) {
      if (*sti) {
        index = (sti - status.begin());
        // IAN. Think this was significant bug. (Compare to read function). Now fixed.
        //if (index < 8) {
        //  fe_status_0 |= 1LL << (index);
        //} else {
        //  fe_status_1 |= 1LL << (index - 8);
        //}
        if (index < 64) {
          fe_status_1 |= 1LL << (index);
        } else {
          fe_status_0 |= 1LL << (index - 64);
        }
      }
    }
    write_n_at_m_R2L(headercopy_, 8, 0, (uint64_t)fe_status_0);
    write_n_at_m_R2L(headercopy_, 64, 64, (uint64_t)fe_status_1);
  }

  void Phase2TrackerFEDHeader::setNumberOfCBC(uint16_t num) {
    write_n_at_m_R2L(headercopy_, CBC_NUMBER_L, CBC_NUMBER_S, (uint64_t)num);
  }

  uint16_t Phase2TrackerFEDHeader::numberOfCBC() const {
    if (debugMode_ != SUMMARY) {
      return static_cast<uint16_t>(read_n_at_m_R2L(headercopy_, CBC_NUMBER_L, CBC_NUMBER_S));
    } else {
      return 0;
    }
  }

  void setCBCStatus() {
    std::ostringstream ss;
    ss << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "] ";
    ss << "CBC status cannot be set : ";
    ss << "Custom header is currently limited to two 64bit words only ";
    throw cms::Exception("Phase2TrackerFEDBuffer") << ss.str();
  }

  std::vector<Phase2TrackerFEDFEDebug> Phase2TrackerFEDHeader::CBCStatus() const {
    // set offset and data to begining
    int offset_bits = 128;
    // number of CBC:
    uint16_t cbc_num = numberOfCBC();
    // get FE status
    std::vector<bool> status = frontendStatus();
    // check that #CBC = CBC_PER_FE_DEBUG x #FE
    int ncbc = CBC_PER_FE_DEBUG;
    int fe_num = std::count(status.begin(), status.end(), true);
    if (cbc_num != fe_num * ncbc) {
      std::ostringstream ss;
      ss << "[Phase2Tracker::Phase2TrackerFEDHeader::" << __func__ << "] ";
      ss << "Number of chips (" << cbc_num << ") should be CBC_PER_FE_DEBUG (default 16) x #FE (" << fe_num << ")";
      throw cms::Exception("Phase2TrackerFEDHeader") << ss.str();
    }
    // store status in list of Phase2TrackerFEDFEDebug
    std::vector<Phase2TrackerFEDFEDebug> all_fe_debug;
    std::vector<bool>::iterator FE_it;
    for (FE_it = status.begin(); FE_it < status.end(); FE_it++) {
      // create empty debug container
      Phase2TrackerFEDFEDebug fe_debug_status(readoutMode_, debugMode_);
      if (*FE_it) {
        fe_debug_status.setOn();
        if (debugMode_ == FULL_DEBUG) {
          if (readoutMode_ == READOUT_MODE_ZERO_SUPPRESSED) {
            for (uint8_t i = 0; i < ncbc; i++) {
              fe_debug_status.setChipDebugStatus(i,
                                                 static_cast<uint32_t>(read_n_at_m_L2R(
                                                     trackerHeader_, CBC_STATUS_SIZE_DEBUG_SPARSIFIED, offset_bits)));
              offset_bits += CBC_STATUS_SIZE_DEBUG_SPARSIFIED;
            }
            fe_debug_status.setFEDebugStatus(
                static_cast<uint32_t>(read_n_at_m_L2R(trackerHeader_, FE_STATUS_SIZE_DEBUG_SPARSIFIED, offset_bits)));
            offset_bits += FE_STATUS_SIZE_DEBUG_SPARSIFIED;
          } else if (readoutMode_ == READOUT_MODE_PROC_RAW) {
            for (uint8_t i = 0; i < ncbc; i++) {
              fe_debug_status.setChipDebugStatus(i,
                                                 static_cast<uint32_t>(read_n_at_m_L2R(
                                                     trackerHeader_, CBC_STATUS_SIZE_DEBUG_UNSPARSIFIED, offset_bits)));
              offset_bits += CBC_STATUS_SIZE_DEBUG_UNSPARSIFIED;
            }
          }
        } else if (debugMode_ == CBC_ERROR) {
          for (uint8_t i = 0; i < ncbc; i++) {
            fe_debug_status.setChipDebugStatus(
                i, static_cast<uint32_t>(read_n_at_m_L2R(trackerHeader_, CBC_STATUS_SIZE_ERROR, offset_bits)));
            offset_bits += CBC_STATUS_SIZE_ERROR;
          }
        }
      }
      all_fe_debug.push_back(fe_debug_status);
    }
    return all_fe_debug;
  }

  // TODO: update this for CBC3 debug format
  const uint8_t* Phase2TrackerFEDHeader::pointerToData() {
    int status_size = 0;
    int cbc_num = numberOfCBC();

    // all sizes in bits here
    if (debugMode_ == FULL_DEBUG) {
      if (readoutMode_ == READOUT_MODE_PROC_RAW) {
        status_size = cbc_num * CBC_STATUS_SIZE_DEBUG_UNSPARSIFIED;
      } else {
        status_size =
            cbc_num * CBC_STATUS_SIZE_DEBUG_SPARSIFIED + 72 * FE_STATUS_SIZE_DEBUG_SPARSIFIED;  // TODO : check this !
      }
    } else if (debugMode_ == CBC_ERROR) {
      status_size = cbc_num * CBC_STATUS_SIZE_ERROR;
    }
    // compute number of additional 64 bit words before payload
    int num_add_words64 = (status_size + 64 - 1) / 64;
    // back to bytes
    trackerHeaderSize_ = (2 + num_add_words64) * 8;
    return &trackerHeader_[trackerHeaderSize_];
  }

}  // namespace Phase2Tracker
