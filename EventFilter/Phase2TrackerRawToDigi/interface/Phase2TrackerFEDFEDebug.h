
#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDFEDebug_H // {
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDFEDebug_H

#include <stdint.h>

namespace Phase2Tracker {

  // tracker FE debug
  class Phase2TrackerFEDFEDebug 
  {
    public: 
      Phase2TrackerFEDFEDebug() : 
        _fedebugstatus(0),
        _chipdebugstatus(),
        _is_on(false)
      {
      }
      inline uint32_t getFEDebugStatus() const { return _fedebugstatus; } 
      inline uint32_t getChipDebugStatus(uint8_t) const;
      inline bool IsOn() const { return _is_on; }
      void setFEDebugStatus(uint32_t);
      void setChipDebugStatus(uint8_t,uint32_t);
      inline void setOn() { _is_on = true; }
    private:
      uint32_t _fedebugstatus;
      uint32_t _chipdebugstatus[16];
      bool _is_on;
  }; // end of Phase2TrackerFEDFEDebug class

}

#endif // } 
