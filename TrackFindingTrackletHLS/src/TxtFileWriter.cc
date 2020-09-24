#include "L1Trigger/TrackFindingTrackletHLS/interface/TxtFileWriter.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <fstream>
#include <vector>
#include <bitset>

//=== Initialize config

TxtFileWriter::TxtFileWriter(const trackerDTC::Setup& setup,
			     //			     int linkWordNBits,
			     unsigned int linkWord_is2sBit,
			     unsigned int linkWord_nLayersOffset,
			     unsigned int linkWord_layerBarrelEndcapOffset,
			     unsigned int linkWord_layerOrDiskNumberOffset) : 
  // DTC ESProduct
  setup_(&setup), 
  // DTC geometry
  numRegions_(setup.numRegions()), 
  numDTCsPerRegion_(setup.numDTCsPerRegion()),
  // Bit numbers used by writeDtcLinkWords()
  //  linkWordNBits_(linkWordNBits),
  linkWord_is2sBit_(linkWord_is2sBit),
  linkWord_nLayersOffset_(linkWord_nLayersOffset),
  linkWord_layerBarrelEndcapOffset_(linkWord_layerBarrelEndcapOffset),
  linkWord_layerOrDiskNumberOffset_(linkWord_layerOrDiskNumberOffset)
  {}

//=== Write file with layer number & phi range of each DTC.

void TxtFileWriter::writeDtcPhiRange(std::string fileName) const {

  std::ofstream outFile(fileName.c_str());
  assert(outFile.good());

  // Loop over DTCs in this tracker nonant.
  for (int iChan = 0; iChan < numDTCsPerRegion_; iChan++) {
    typedef std::pair<float,float> PhiRange;
    std::map<int, PhiRange> dtcPhiRange;
    int dtcId = iChan; // DTC identifier in nonant 0.

    // Loop over all tracker nonants, taking worst case not all identical.
    for (int iRegion = 0; iRegion < numRegions_; iRegion++) {
      int dtcId_regI = iRegion * numDTCsPerRegion_ + dtcId;
      const std::vector<trackerDTC::SensorModule*>& dtcModules = setup_->dtcModules(dtcId_regI);
      for (const trackerDTC::SensorModule* sm : dtcModules) {
	// Hybrid L1Track currently starts counting endcap layers at 7.
	const int layerOffsetHybrid = 4;
	// Hybrid L1Track currently measures phi w.r.t. lower edge of tracker nonant.
	const float phiOffsetHybrid = M_PI/setup_->numRegions();

	int layer = sm->layerId(); // Barrel = 1-6, Endcap = 11-15;
	if (not sm->barrel()) layer -= layerOffsetHybrid;
	// Inner radius of module.
	float r = sm->barrel()  ?  sm->r()  :  sm->r() - 0.5*sm->numColumns()*sm->pitchCol();
	// phi with respect to tracker nonant centre.
	float phiMin = sm->phi() - 0.5*sm->numRows()*sm->pitchRow() / r;
	float phiMax = sm->phi() + 0.5*sm->numRows()*sm->pitchRow() / r;
	phiMin += phiOffsetHybrid;
	phiMax += phiOffsetHybrid;
	if (dtcPhiRange.find(layer) == dtcPhiRange.end()) {
	  dtcPhiRange[layer] = {phiMin, phiMax};
	} else {
	  dtcPhiRange.at(layer).first  = std::min(phiMin, dtcPhiRange.at(layer).first);
	  dtcPhiRange.at(layer).second = std::max(phiMax, dtcPhiRange.at(layer).second);
	}
      }
    }
    for (const auto& p : dtcPhiRange) {
      outFile<<dtcId<<" "<<this->dtcType(dtcId)<<" "<<p.first<<" "<<p.second.first<<" "<<p.second.second<<std::endl;
    }
  }
}

//=== Write file used by IR HLS with info about layers of each DTC.

void TxtFileWriter::writeDtcLinkWords(std::string fileName) const {

  std::ofstream outFile(fileName.c_str());
  assert(outFile.good());

  // Loop over DTCs in this tracker nonant.
  for (int iChan = 0; iChan < numDTCsPerRegion_; iChan++) {
    typedef std::pair<float,float> PhiRange;
    std::map<int, PhiRange> dtcPhiRange;
    int dtcId = iChan; // DTC identifier in nonant 0.

    std::bitset<linkWordNBits_> linkWord = 0;  // output word in nonant 0
    std::bitset<linkWordNBits_> linkWord_regI = 0;  // output word in nonant i

    // Loop over all tracker nonants, checking all identical.
    for (int iRegion = 0; iRegion < numRegions_; iRegion++) {
      int dtcId_regI = iRegion * numDTCsPerRegion_ + dtcId;

      const std::vector<int>& thisDtcLayerEncodings{ setup_->encodingLayerId( dtcId_regI ) };

      linkWord_regI |= ( thisDtcLayerEncodings.size() << linkWord_nLayersOffset_ );

      bool is2S = not setup_->psModule( dtcId_regI );
      linkWord_regI |= ( is2S << ( linkWord_is2sBit_ ) );

      size_t encodedLayer = 0;
      for ( auto decodedLayer : thisDtcLayerEncodings ) {
        bool isBarrelLayer = decodedLayer < setup_->offsetLayerDisks();
        if ( !isBarrelLayer ) {
          decodedLayer -= setup_->offsetLayerDisks();
        }
        linkWord_regI |= ( ( isBarrelLayer << linkWord_layerBarrelEndcapOffset_ ) ) << ( setup_->hybridNumLayers() * encodedLayer );
        linkWord_regI |= ( decodedLayer << linkWord_layerOrDiskNumberOffset_ ) << ( setup_->hybridNumLayers() * encodedLayer );
        ++encodedLayer;
      }

      if (iRegion == 0) {
        linkWord = linkWord_regI;
	outFile<<dtcId<<" "<<this->dtcType(dtcId)<<std::hex<<" 0x"<<linkWord.to_ulong()<<std::dec<<std::endl;
      }
      if (linkWord_regI != linkWord) throw cms::Exception("LogicError")<<"TxtFileWriter: Link word depends on tracker phi nonant "<<std::hex<<"0x"<<linkWord_regI.to_ulong()<<" 0x"<<linkWord.to_ulong()<<std::dec;
    }
  }
}

// Get string describing DTC type.

std::string TxtFileWriter::dtcType(int dtcId) const {
  std::string dtcType = setup_->psModule(dtcId) ? "PS" : "2S" ;
  if (setup_->slot(dtcId) <= 3) dtcType = "PS10G"; // First 3 slots in ATCA crate used for PS10G.
  bool posZ = setup_->side(dtcId);
  if (not posZ) dtcType = "neg_" + dtcType;
  return dtcType;
}

