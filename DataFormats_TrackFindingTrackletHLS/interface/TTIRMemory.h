#ifndef DataFormats_L1TrackTrigger_TTIRMemory_h
#define DataFormats_L1TrackTrigger_TTIRMemory_h

#include "DataFormats/L1TrackTrigger/interface/TTDTC.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "InputRouterTop.h"

#include <vector>

/*! 
 * \class  TTIRMemory
 * \brief  Class to store hardware like structured TTStub Collection used by Track Trigger emulators
 * \author Emyr Clement
 * \date   2020, May
 */
class TTIRMemory {
public:

  typedef std::vector<StubsBarrelPS> IRMemoryCollection;

public:
  TTIRMemory() {}
  TTIRMemory(int numRegions, int numOverlappingRegions, int numDTCsPerRegion);
  ~TTIRMemory() {}


  const StubsBarrelPS& IRMemory(int tfpRegion, int tfpChannel) const;
  void setIRMemory(int tfpRegion, int tfpChannel, const StubsBarrelPS& stubsBarrelPS);


private:
  // DTC Class
  TTDTC TTDTC_;

  // IR Memories
  IRMemoryCollection irMemories_;
};

#endif