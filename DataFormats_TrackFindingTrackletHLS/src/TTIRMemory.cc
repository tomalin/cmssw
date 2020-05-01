#include "L1Trigger/DataFormats_TrackFindingTrackletHLS/interface/TTIRMemory.h"

using namespace std;
using namespace edm;

TTIRMemory::TTIRMemory(int numRegions, int numOverlappingRegions, int numDTCsPerRegion) :
	TTDTC_(numRegions, numOverlappingRegions, numDTCsPerRegion ),
	irMemories_( numRegions * numOverlappingRegions * numDTCsPerRegion )
{
}

void TTIRMemory::setIRMemory(int tfpRegion, int tfpChannel, const StubsBarrelPS& stubsBarrelPS) {
	// Check arguments, or that the dtc helper function would catch dodgy indices?
	int dtcId = TTDTC_.dtcId( tfpRegion, tfpChannel );
	irMemories_[dtcId] = move(stubsBarrelPS);
}


// read one specific stream of TTStubRefs using TFP identifier (region[0-8], channel[0-47])
// tfpRegions aka processing regions are rotated by -0.5 region width w.r.t detector regions
const StubsBarrelPS& TTIRMemory::IRMemory(int tfpRegion, int tfpChannel) const {
	// Check arguments, or that the dtc helper function would catch dodgy indices?
	int dtcId = TTDTC_.dtcId( tfpRegion, tfpChannel );
	return irMemories_.at( dtcId );
}