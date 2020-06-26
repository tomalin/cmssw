#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"

using namespace std;
using namespace edm;

TTIRMemory::TTIRMemory( const trackerDTC::Setup& setup ) :
	irMemories_( setup.numRegions() * setup.numOverlappingRegions() * setup.numDTCsPerRegion() )
{
}


void TTIRMemory::setIRMemory(int dtcId, const IRMemories& irMemories) {

	// TODO : Add checks on arguments, that irMemories_ is large enough
	irMemories_[dtcId] = move(irMemories);
}


const TTIRMemory::IRMemories& TTIRMemory::IRMemory(int dtcId ) const {

	// TODO : Add checks on arguments, that irMemories_ is large enough
	return irMemories_.at( dtcId );
}