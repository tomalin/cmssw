#ifndef DataFormats_L1TrackTrigger_TTIRMemory_h
#define DataFormats_L1TrackTrigger_TTIRMemory_h

#include "L1Trigger/TrackerDTC/interface/Setup.h"
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

  static constexpr int stubWordIRNBits_ = 36;

  typedef std::vector<StubsBarrelPS> IRMemoryCollection;

public:
  TTIRMemory() {}
  TTIRMemory( const trackerDTC::Setup& setup );
  ~TTIRMemory() {}


  const StubsBarrelPS& IRMemory(int dtcId) const;
  void setIRMemory(int dtcId, const StubsBarrelPS& stubsBarrelPS);


private:
  // IR Memories
  IRMemoryCollection irMemories_;

  // TTStubRefs corresponding to stubs in irMemories_ ???
  // std::vector< std::map< std::bitset<stubWordIRNBits_>, TTStubRef > > mapStubWordsToTTStubRefs_;
  // std::vector< std::vector< std::vector< TTStubRef > > > ttStubRefs_;
};

#endif