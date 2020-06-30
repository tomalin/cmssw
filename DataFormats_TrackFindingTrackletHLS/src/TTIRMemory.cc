#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"

using namespace std;
using namespace edm;

TTIRMemory::TTIRMemory( const trackerDTC::Setup& setup ) :
	irMemories_( setup.numRegions() * setup.numOverlappingRegions() * setup.numDTCsPerRegion() ),
	ttStubRefs_( setup.numRegions() * setup.numOverlappingRegions() * setup.numDTCsPerRegion() )
{
}


void TTIRMemory::setIRMemory(int streamId, const IRMemories& irMemories) {

	// TODO : Add checks on arguments, that irMemories_ is large enough
	irMemories_[streamId] = move(irMemories);
}


const TTIRMemory::IRMemories& TTIRMemory::IRMemory(int streamId ) const {

	// TODO : Add checks on arguments, that irMemories_ is large enough
	return irMemories_.at( streamId );
}

void TTIRMemory::setTTStubs(int streamId, const TTStubRefs& ttStubRefs) {

	// TODO : Add checks on arguments, that ttStubRefs_ is large enough
	ttStubRefs_[streamId] = ttStubRefs;
}


const TTIRMemory::TTStubRefs& TTIRMemory::TTStubs(int streamId ) const {

	// TODO : Add checks on arguments, that ttStubRefs_ is large enough
	return ttStubRefs_.at( streamId );
}