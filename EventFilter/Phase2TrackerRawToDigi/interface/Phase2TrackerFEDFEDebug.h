
#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDFEDebug_H // {
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDFEDebug_H
    
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include <stdint.h>
    
namespace Phase2Tracker {

  // tracker FE debug
  class Phase2TrackerFEDFEDebug 
  {
    public: 
      Phase2TrackerFEDFEDebug() : 
        _fedebugstatus(0),
        _chipdebugstatus(),
        _is_on(false),
        _readoutmode(READOUT_MODE_ZERO_SUPPRESSED),
        _debugmode(SUMMARY)
      {
      }
      Phase2TrackerFEDFEDebug(FEDReadoutMode,READ_MODE);
      inline uint32_t getFEDebugStatus() const { return _fedebugstatus; } 
      std::vector<uint16_t> getFEL1ID() const; 
      uint32_t getChipDebugStatus(uint8_t) const;
      uint16_t getChipError(uint8_t) const;
      uint16_t getChipL1ID(uint8_t) const;
      uint16_t getChipPipelineAddress(uint8_t) const;
      inline bool IsOn() const { return _is_on; }
      void setFEDebugStatus(uint32_t);
      void setChipDebugStatus(uint8_t,uint32_t);
      void setReadoutMode(FEDReadoutMode);
      void setDebugMode(READ_MODE);
      inline void setOn() { _is_on = true; }

    private:
      uint32_t _fedebugstatus;
      uint32_t _chipdebugstatus[16];
      bool _is_on;
      FEDReadoutMode _readoutmode;
      READ_MODE _debugmode;
  }; // end of Phase2TrackerFEDFEDebug class

  // Dummy comparator to be able to store in DetSetVector
  inline bool operator<( const Phase2TrackerFEDFEDebug& one, const Phase2TrackerFEDFEDebug& other) { return true; }
}

#endif // } 
