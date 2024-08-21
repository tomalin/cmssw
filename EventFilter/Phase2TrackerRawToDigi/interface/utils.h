#ifndef EventFilter_Phase2TrackerRawToDigi_utils_H  // {
#define EventFilter_Phase2TrackerRawToDigi_utils_H

// common tools
#include <iomanip>
#include <ostream>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cassert>
#include <iostream> // Debug
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

namespace Phase2Tracker {

  // TODO: set this in a common include file.
  // see DataFormats/Phase2TrackerCommon/interface/Constants.h

  // -------------------- FED ids --------------------

  static const uint16_t FED_ID_MIN = static_cast<uint16_t>(FEDNumbering::MINSiStripFEDID);
  static const uint16_t FED_ID_MAX = static_cast<uint16_t>(FEDNumbering::MAXSiStripFEDID);
  static const uint16_t CMS_FED_ID_MAX = static_cast<uint16_t>(FEDNumbering::MAXFEDID);
  static const uint16_t NUMBER_OF_FEDS = static_cast<uint16_t>(FED_ID_MAX - FED_ID_MIN + 1);

  // Assumptions for phase 2

  static const int MAX_FE_PER_FED = 72;
  static const int MAX_CBC_PER_FE = 16;
  static const int STRIPS_PER_CBC = 254;
  static const int PS_ROWS = 120;
  static const int PS_COLS = 32;
  static const int STRIPS_PADDING = 2;
  static const int P_CLUSTER_SIZE_BITS = 18;
  static const int Son2S_CLUSTER_SIZE_BITS = 15;
  static const int SonPS_CLUSTER_SIZE_BITS = 16;
  static const int CBC_STATUS_SIZE_DEBUG_UNSPARSIFIED = 20;
  static const int CBC_STATUS_SIZE_DEBUG_SPARSIFIED = 1;
  static const int FE_STATUS_SIZE_DEBUG_SPARSIFIED = 18;
  static const int CBC_STATUS_SIZE_ERROR = 1;
  static const int STUBS_SIZE_2S = 16;
  static const int STUBS_SIZE_PS = 20;

  // Current dataformat does not allow to know how many CBC there are per FE, set it here
  static const int CBC_PER_FE_DEBUG = 2;

  // definition

  static const uint8_t INVALID = 0xFF;

  // utils

  inline void printNibbleValue(uint8_t value, std::ostream& os) {
    const std::ios_base::fmtflags originalFormatFlags = os.flags();
    os << std::hex << std::setw(1) << value;
    os.flags(originalFormatFlags);
  }

  inline void printHexValue(const uint8_t value, std::ostream& os) {
    const std::ios_base::fmtflags originalFormatFlags = os.flags();
    os << std::hex << std::setfill('0') << std::setw(2);
    os << uint16_t(value);
    os.flags(originalFormatFlags);
  }

  inline void printHexWord(const uint8_t* pointer, const size_t lengthInBytes, std::ostream& os) {
    size_t i = lengthInBytes - 1;
    do {
      printHexValue(pointer[i], os);
      if (i != 0)
        os << " ";
    } while (i-- != 0);
  }

  inline void printHex(const void* pointer, const size_t lengthInBytes, std::ostream& os) {
    const uint8_t* bytePointer = reinterpret_cast<const uint8_t*>(pointer);
    //if there is one 64 bit word or less, print it out
    if (lengthInBytes <= 8) {
      printHexWord(bytePointer, lengthInBytes, os);
    }
    //otherwise, print word numbers etc
    else {
      //header
      os << "word\tbyte\t                       \t\tbyte" << std::endl;
      ;
      const size_t words = lengthInBytes / 8;
      const size_t extraBytes = lengthInBytes - 8 * words;
      //print full words
      for (size_t w = 0; w < words; w++) {
        const size_t startByte = w * 8;
        os << w << '\t' << startByte + 8 << '\t';
        printHexWord(bytePointer + startByte, 8, os);
        os << "\t\t" << startByte << std::endl;
      }
      //print part word, if any
      if (extraBytes) {
        const size_t startByte = words * 8;
        os << words << '\t' << startByte + 8 << '\t';
        //padding
        size_t p = 8;
        while (p-- > extraBytes) {
          os << "00 ";
        }
        printHexWord(bytePointer + startByte, extraBytes, os);
        os << "\t\t" << startByte << std::endl;
      }
      os << std::endl;
    }
  }

  //enum values are values which appear in FED buffer. DO NOT CHANGE!
  enum FEDReadoutMode {
    READOUT_MODE_INVALID = INVALID,
    READOUT_MODE_SCOPE = 0x1,
    READOUT_MODE_VIRGIN_RAW = 0x2,
    READOUT_MODE_PROC_RAW = 0x6,
    READOUT_MODE_ZERO_SUPPRESSED = 0xA,
    READOUT_MODE_ZERO_SUPPRESSED_LITE = 0xC,
    READOUT_MODE_SPY = 0xE
  };

  //to make enums printable
  std::ostream& operator<<(std::ostream& os, const FEDReadoutMode& value);
  inline std::ostream& operator<<(std::ostream& os, const FEDReadoutMode& value) {
    switch (value) {
      case READOUT_MODE_SCOPE:
        os << "Scope mode";
        break;
      case READOUT_MODE_VIRGIN_RAW:
        os << "Virgin raw";
        break;
      case READOUT_MODE_PROC_RAW:
        os << "Processed raw";
        break;
      case READOUT_MODE_ZERO_SUPPRESSED:
        os << "Zero suppressed";
        break;
      case READOUT_MODE_ZERO_SUPPRESSED_LITE:
        os << "Zero suppressed lite";
        break;
      case READOUT_MODE_SPY:
        os << "Spy channel";
        break;
      case READOUT_MODE_INVALID:
        os << "Invalid";
        break;
      default:
        os << "Unrecognized";
        os << " (";
        printHexValue(value, os);
        os << ")";
        break;
    }
    return os;
  }

  // tracker header read modes
  enum READ_MODE { READ_MODE_INVALID = INVALID, SUMMARY = 0, FULL_DEBUG = 1, CBC_ERROR = 2 };

  // module types
  enum DET_TYPE { UNUSED = -1, DET_Son2S = 0, DET_SonPS = 1, DET_PonPS = 2 };

  enum MOD_TYPE { MOD_2S = 0, MOD_PS = 1 };

  enum STACK_LAYER { LAYER_UNUSED = -1, LAYER_INNER = 0, LAYER_OUTER = 1 };

  //to make enums printable
  std::ostream& operator<<(std::ostream& os, const READ_MODE& value);
  inline std::ostream& operator<<(std::ostream& os, const READ_MODE& value) {
    switch (value) {
      case SUMMARY:
        os << "Summary mode";
        break;
      case FULL_DEBUG:
        os << "Full debug mode";
        break;
      case CBC_ERROR:
        os << "CBC error mode";
        break;
      default:
        os << "Unrecognized mode";
        os << " (";
        printHexValue(value, os);
        os << ")";
        break;
    }
    return os;
  }

  // tracker header masks
  enum trackerHeader_m {
    VERSION_M = 0xF000000000000000,
    HEADER_FORMAT_M = 0x0C00000000000000,
    EVENT_TYPE_M = 0x03C0000000000000,
    GLIB_STATUS_M = 0x003FFFFFFFFF0000,
    FRONTEND_STAT_M = 0x000000000000FFFF,
    CBC_NUMBER_M = 0xFFFF000000000000
  };

  // position of first bit
  enum trackerHeader_s {
    VERSION_S = 60,
    HEADER_FORMAT_S = 58,
    EVENT_TYPE_S = 54,
    GLIB_STATUS_S = 24,
    CBC_NUMBER_S = 8,
    FRONTEND_STAT_S = 0
  };

  // number of bits (replaces mask)
  enum trackerheader_l {
    VERSION_L = 4,
    HEADER_FORMAT_L = 2,
    EVENT_TYPE_L = 4,
    GLIB_STATUS_L = 30,
    CBC_NUMBER_L = 16,
    FRONTEND_STAT_L = 0
  };

  // get 64 bits word from data with given offset : only use if at beginning of 64 bits word
  inline uint64_t read64(int offset, const uint8_t* buffer) {
    return *reinterpret_cast<const uint64_t*>(buffer + offset);
  }

  // extract data from a 64 bits word using mask and shift
  inline uint64_t extract64(trackerHeader_m mask, trackerHeader_s shift, uint64_t data) {
    data = (data & mask) >> shift;
    return data;
  }

  // Read "size" bits starting at bit "pos_bit"
  // The bit count pos_bit increases right-to-left inside each line,
  // which is the standard C++ convention.
  // i.e. Top line of buffer has bit 63 at left side and bit 0 at right.
  // Within each word, the MSB is the left-hand one (and on
  // the lower line if the word extends over 2 lines).
  // A word wrapping over 2 lines goes from left-hand-side of line N to
  // right-hand-side of line N+1. 

inline uint64_t read_n_at_m_R2L(const uint8_t* buffer, int size, int pos_bit) {
    // 1) determine which 64 bit word to read
    int iword = pos_bit / 64;
    uint64_t data = *(uint64_t*)(buffer + (iword * 8)); // LSBs
    int start_bit = pos_bit % 64;
    data >>= start_bit;

    // 2) determine if you need to read another
    int bits_on_line = 64 - start_bit;
    bool lineWrap = (size > bits_on_line);
    if (lineWrap) {
      uint64_t data_supp = *(uint64_t*)(buffer + ((iword + 1) * 8)); // MSBs
      data |= data_supp << bits_on_line;
    }

    // 3) mask according to expected size
    if (size < 64) {
      uint64_t mask = (1LL << size) - 1; 
      data &= mask;
    }
    return data;
  }

  // Read "size" bits starting at bit "pos_bit"
  // The bit count pos_bit increases left-to-right inside each line.
  // i.e. Top line of buffer has bit 0 at left side and bit 63 at right.
  // This needed for things like the tracker data buffer,
  // where words are added left to right (e.g words for cluster 0 to
  // left of those for cluster 1).
  // Nonetheless, within each word, the MSB is the left-hand one (and on
  // the upper line if the word extends over 2 lines).
  // A word wrapping over 2 lines goes from right-hand-side of line N to
  // left-hand-side of line N+1. 

inline uint64_t read_n_at_m_L2R(const uint8_t* buffer, int size, int pos_bit) {
    int iword = pos_bit / 64;
    uint64_t data = *(uint64_t*)(buffer + (iword * 8));
    int start_bit = pos_bit % 64;
    int bits_remaining_on_line = 64 - (start_bit + size);
    bool lineWrap = (bits_remaining_on_line < 0);    

    if (not lineWrap) {
      data >>= bits_remaining_on_line;
    } else {
      int bits_on_next_line = - bits_remaining_on_line; 
      data <<= bits_on_next_line; // MSBs
      uint64_t data_supp = *(uint64_t*)(buffer + ((iword + 1) * 8)); // LSBs
      data |= data_supp >> (64 - bits_on_next_line);
    }

    if (size < 64) {
      uint64_t mask = (1LL << size) - 1; 
      data &= mask;
    }
    return data;
  }

  // Write "size " bits of word "data" to "buffer" at bit address "pos_bit".
  // The bit count pos_bit increases right-to-left inside each line,
  // as explained in comment for read_n_at_m_R2L.

inline void write_n_at_m_R2L(uint8_t* buffer, int size, int pos_bit, uint64_t data) {
    //assert (size <= 64 && pos_bit <=64);
    // mask unwanted bits
    uint64_t mask = (size < 64) ? (1LL << size) - 1 : static_cast<uint64_t>(-1);
    data &= mask;
    int iword = pos_bit / 64;
    int start_bit = pos_bit % 64;
    uint64_t curr_data = *(uint64_t*)(buffer + (iword * 8));
    // Set bits to be written to first to zero then to desired value.
    curr_data &= ~(mask << start_bit);
    curr_data |=  (data << start_bit);
    memcpy(buffer + (iword * 8), &curr_data, 8); // LSBs
    
    int bits_on_line = 64 - start_bit;
    bool lineWrap = (size > bits_on_line);

    if (lineWrap) {
      // there are more bits to write
      uint64_t data_supp = *(uint64_t*)(buffer + ((iword + 1) * 8)); 
      data_supp &= ~(mask >> bits_on_line);
      data_supp |= (data >> bits_on_line);
      memcpy(buffer + ((iword + 1) * 8), &data_supp, 8); // MSBs
    }
  }

  // Write "size " bits of word "data" to "buffer" at bit address "pos_bit".
  // The bit count pos_bit increases left-to-right inside each line,
  // as explained in comment for read_n_at_m_L2R.

  inline void write_n_at_m_L2R(uint8_t* buffer, int size, int pos_bit, uint64_t data) {
    //assert (size <= 64 && pos_bit <=64);
    //assert (size <= 64 && pos_bit <=128);
    // assert (size <= 64 && pos_bit <=9999);
    assert (size <= 128 && pos_bit <=9999);

    // mask unwanted bits
    uint64_t mask = (size < 64) ? (1LL << size) - 1 : static_cast<uint64_t>(-1);
    data &= mask;

    int iword = pos_bit / 64;
    int start_bit = pos_bit % 64;
    
    int bits_remaining_on_line = 64 - (start_bit + size);
    bool lineWrap = (bits_remaining_on_line < 0);  
    uint64_t curr_data = *(uint64_t*)(buffer + (iword * 8));

    if (not lineWrap) {
      // Set bits to be written to first to zero then to desired value.
      curr_data &= ~(mask << bits_remaining_on_line);
      curr_data |= (data << bits_remaining_on_line);
      memcpy(buffer + (iword * 8), &curr_data, 8);
    } else {
      int bits_on_next_line = - bits_remaining_on_line;
      curr_data &= ~(mask >> bits_on_next_line);
      curr_data |= (data >> bits_on_next_line);
      memcpy(buffer + (iword * 8), &curr_data, 8); // MSBs

      int pos_bit_next_line = (iword + 1) * 64; // start of line
      write_n_at_m_L2R(buffer, bits_on_next_line, pos_bit_next_line, data); // LSBs
    }
  }

  // Write "size" bits of word "data" to "buffer" at bit address "pos_bit",
  // where assumed size <= 64.
  // boolean indicates if bit count pos_bit increases left to right or visa-versa.
  inline void write_n_at_m(std::vector<uint64_t>& buffer, int size, int pos_bit, uint64_t data, bool L2R = true) {
    assert (size <= 64);
    int iword = pos_bit / 64;
    int toadd = (pos_bit + size - 1 + 64) / 64 - buffer.size();
    if (toadd > 0) {
      buffer.insert(buffer.end(), toadd, (uint64_t)0x00);
    }
    bool lineWrap = (pos_bit % 64 + size > 64);
    uint64_t temp[] = {buffer[iword], 0};
    if (lineWrap) {
      temp[1] = buffer[iword + 1];
    }
    // QUESTION Ian: Why add complexity of casting to uint8?
    uint8_t* tt = (uint8_t*)(temp);
    if (L2R) {
      write_n_at_m_L2R(tt, size, pos_bit % 64, data);
    } else {
      write_n_at_m_R2L(tt, size, pos_bit % 64, data);
    }
    buffer[iword] = *(uint64_t*)(tt);
    if (lineWrap) {
      buffer[iword + 1] = *(uint64_t*)(tt + 8);
    }
  }

  inline void vec_to_array(std::vector<uint64_t> vec, uint8_t* arr) {
    std::vector<uint64_t>::iterator it;
    for (it = vec.begin(); it != vec.end(); it++) {
      memcpy(arr + 8 * (it - vec.begin()), &*it, 8);
    }
  }
}  // namespace Phase2Tracker

#endif  // } end def utils
