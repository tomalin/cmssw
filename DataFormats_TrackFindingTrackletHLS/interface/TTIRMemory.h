#ifndef DataFormats_L1TrackTrigger_TTIRMemory_h
#define DataFormats_L1TrackTrigger_TTIRMemory_h

#include "L1Trigger/TrackerDTC/interface/Setup.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "InputRouterTop.h"

#include <vector>
#include <variant>

/*! 
 * \class  TTIRMemory
 * \brief  Class to store output from IR track finding HLS module
 * \author Emyr Clement
 * \date   2020, May
 */
class TTIRMemory {
public:

  // typedef std::vector< std::variant< AllStubMemory<BARRELPS>, AllStubMemory<BARREL2S>, AllStubMemory<DISKPS>, AllStubMemory<DISK2S> > > IRMemories;
  typedef std::vector< InputStubMemory<TRACKER> > IRMemories;
  typedef std::vector<IRMemories> IRMemoryCollection;

  typedef std::vector< std::vector< TTStubRef > > TTStubRefs;
  typedef std::vector< TTStubRefs > TTStubRefsCollection;

public:
  TTIRMemory() {}
  TTIRMemory( const trackerDTC::Setup& setup );
  ~TTIRMemory() {}


  const IRMemories& IRMemory(int streamId) const;
  void setIRMemory(int streamId, const IRMemories& irMemories);

  const TTStubRefs& TTStubs(int streamId ) const;
  void setTTStubs(int streamId, const TTStubRefs& ttStubRefs);

private:
  // IR Memories
  IRMemoryCollection irMemories_;

  TTStubRefsCollection ttStubRefs_;
};

#endif