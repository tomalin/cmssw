#ifndef L1Trigger_TrackFindingTrackletHLS_interface_TxtFileWriter_h
#define L1Trigger_TrackFindingTrackletHLS_interface_TxtFileWriter_h

#include "L1Trigger/TrackerDTC/interface/Setup.h"

#include <string>

class TxtFileWriter {

public:

  // Initialize config
  TxtFileWriter(const trackerDTC::Setup& setup,
		//		int linkWordNBits,
		unsigned int linkWord_is2sBit,
		unsigned int linkWord_nLayersOffset,
		unsigned int linkWord_layerBarrelEndcapOffset,
		unsigned int linkWord_layerOrDiskNumberOffset); 

  // Write file with layer number & phi range of each DTC.
  void writeDtcPhiRange(std::string fileName) const;

  // Write file used by IR HLS with info about layers of each DTC.
  void writeDtcLinkWords(std::string fileName) const;

private:

  // Get string describing DTC type.
  std::string dtcType(int dtcId) const;

private:

  const trackerDTC::Setup* setup_;

  const int numRegions_;
  const int numDTCsPerRegion_;

  // Annoyingly std::bitset<N> requires constexpr "N", so value of linkWordNBits must be hard-wired
  static constexpr int linkWordNBits_ = 20;
  //int linkWordNBits_;
  const unsigned int linkWord_is2sBit_;
  const unsigned int linkWord_nLayersOffset_;
  const unsigned int linkWord_layerBarrelEndcapOffset_;
  const unsigned int linkWord_layerOrDiskNumberOffset_;
};

#endif
