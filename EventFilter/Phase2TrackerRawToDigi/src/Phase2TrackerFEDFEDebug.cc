#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDFEDebug.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Phase2Tracker
{

  Phase2TrackerFEDFEDebug::Phase2TrackerFEDFEDebug(FEDReadoutMode readoutmode, READ_MODE debugmode) :
    _fedebugstatus(0),
    _chipdebugstatus(),
    _is_on(false),
    _readoutmode(readoutmode),
    _debugmode(debugmode)
  {
  }

  uint32_t Phase2TrackerFEDFEDebug::getChipDebugStatus(uint8_t chip_id) const
  {
    if (!(_readoutmode == READOUT_MODE_ZERO_SUPPRESSED) and !(_readoutmode == READOUT_MODE_PROC_RAW))
    {
      std::ostringstream ss;
      ss << "[Phase2Tracker::Phase2TrackerFEDHeader::"<<__func__<<"] ";
      ss << "Invalid Readout mode";
      throw cms::Exception("Phase2TrackerFEDHeader") << ss.str();
    }

    if (chip_id < 16)
    {
      return _chipdebugstatus[chip_id];
    }
    else
    {
      std::ostringstream ss;
      ss << "[Phase2Tracker::Phase2TrackerFEDHeader::"<<__func__<<"] ";
      ss << "Chip id exceeding maximum of 16";
      throw cms::Exception("Phase2TrackerFEDHeader") << ss.str();
    }
    return 0;
  }

  uint16_t Phase2TrackerFEDFEDebug::getChipError(uint8_t chip_id) const 
  {
    if (_readoutmode == READOUT_MODE_ZERO_SUPPRESSED || _debugmode == CBC_ERROR)
    {
      return getChipDebugStatus(chip_id)&0x0001;
    }
    else if (_readoutmode == READOUT_MODE_PROC_RAW)
    {
      return (getChipDebugStatus(chip_id)>>18)&0x0003;
    }
    return 0;
  }

  uint16_t Phase2TrackerFEDFEDebug::getChipL1ID(uint8_t chip_id) const 
  {
    if (_readoutmode == READOUT_MODE_PROC_RAW || _debugmode == FULL_DEBUG)
    {
      return getChipDebugStatus(chip_id)&(0x01FF);
    }
    else 
    {
      #ifdef EDM_ML_DEBUG
      LogTrace("Phase2TrackerFEDFEDebug") << "[Phase2Tracker::Phase2TrackerFEDFEDebug::"<<__func__<<"]: \n";
      LogTrace("Phase2TrackerFEDFEDebug") << "Warning : Chip L1Id only present in Sparsified, Full debug mode" << std::endl;
      #endif
    }
    return 0;
  }

  std::vector<uint16_t> Phase2TrackerFEDFEDebug::getFEL1ID() const 
  {
    std::vector<uint16_t> l1id;
    if (_readoutmode == READOUT_MODE_ZERO_SUPPRESSED && _debugmode == FULL_DEBUG)
    {
        l1id.push_back((uint16_t)(_fedebugstatus&0x01FF));
        l1id.push_back((uint16_t)(_fedebugstatus>>9));
    }
    else 
    {
      #ifdef EDM_ML_DEBUG
      LogTrace("Phase2TrackerFEDFEDebug") << "[Phase2Tracker::Phase2TrackerFEDFEDebug::"<<__func__<<"]: \n";
      LogTrace("Phase2TrackerFEDFEDebug") << "Warning : FE L1Id only present in Unsparsified, Full debug mode" << std::endl;
      #endif
    }
    return l1id;
   
  }

  uint16_t Phase2TrackerFEDFEDebug::getChipPipelineAddress(uint8_t chip_id) const 
  {
    if (_readoutmode == READOUT_MODE_PROC_RAW || _debugmode == FULL_DEBUG)
    {
      return (getChipDebugStatus(chip_id)>>9)&(0x01FF);
    }
    else 
    {
      #ifdef EDM_ML_DEBUG
      LogTrace("Phase2TrackerFEDFEDebug") << "[Phase2Tracker::Phase2TrackerFEDFEDebug::"<<__func__<<"]: \n";
      LogTrace("Phase2TrackerFEDFEDebug") << "Warning : Chip Pipeline Address only present in Unsparsified, Full debug mode" << std::endl;
      #endif
    }
    return 0;
  }

  void Phase2TrackerFEDFEDebug::setFEDebugStatus(uint32_t fe_status)
  {
    _fedebugstatus = fe_status;
  }

  void Phase2TrackerFEDFEDebug::setReadoutMode(FEDReadoutMode readoutmode)
  {
    _readoutmode = readoutmode;
  }

  void Phase2TrackerFEDFEDebug::setDebugMode(READ_MODE debugmode)
  {
    _debugmode = debugmode;
  }

  void Phase2TrackerFEDFEDebug::setChipDebugStatus(uint8_t chip_id, uint32_t chip_status)
  {
    _chipdebugstatus[chip_id] = chip_status;
  }

}  // end of Phase2Tracker namespace
