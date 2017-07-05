#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDFEDebug.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Phase2Tracker
{

  inline uint32_t Phase2TrackerFEDFEDebug::getChipDebugStatus(uint8_t chip_id) const
  {
    if (chip_id < 16)
    {
      return _chipdebugstatus[chip_id];
    }
    else
    {
      std::ostringstream ss;
      ss << "[Phase2Tracker::Phase2TrackerFEDHeader::"<<__func__<<"] ";
      ss << "Chip number exceeding maximum";
      throw cms::Exception("Phase2TrackerFEDHeader") << ss.str();
    }
  }

  void Phase2TrackerFEDFEDebug::setFEDebugStatus(uint32_t fe_status)
  {
    _fedebugstatus = fe_status;
  }

  void Phase2TrackerFEDFEDebug::setChipDebugStatus(uint8_t chip_id, uint32_t chip_status)
  {
    _chipdebugstatus[chip_id] = chip_status;
  }

}  // end of Phase2Tracker namespace
