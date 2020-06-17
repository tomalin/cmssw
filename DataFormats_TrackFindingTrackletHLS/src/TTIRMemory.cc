#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"

using namespace std;
using namespace edm;

TTIRMemory::TTIRMemory( const trackerDTC::Setup& setup ) :
	irMemories_( setup.numRegions() * setup.numOverlappingRegions() * setup.numDTCsPerRegion() )
{
}

void TTIRMemory::setIRMemory(int dtcId, const StubsBarrelPS& stubsBarrelPS) {
	// Check arguments, or that the dtc helper function would catch dodgy indices?
	// int dtcId = setup_.dtcId( tfpRegion, tfpChannel );

	// Check irMemories_ is large enough

	irMemories_[dtcId] = move(stubsBarrelPS);
}


const StubsBarrelPS& TTIRMemory::IRMemory(int dtcId ) const {
	// Check arguments, or that the dtc helper function would catch dodgy indices?
	// int dtcId = setup_.dtcId( tfpRegion, tfpChannel );
	
	// Check irMemories_ large enough

	return irMemories_.at( dtcId );
}