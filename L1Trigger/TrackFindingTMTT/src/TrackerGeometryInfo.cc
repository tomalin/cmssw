#include "L1Trigger/TrackFindingTMTT/interface/TrackerGeometryInfo.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

using namespace std;

namespace tmtt {

  TrackerGeometryInfo::TrackerGeometryInfo()
      : barrelNTiltedModules_(12), barrelNLayersWithTiltedModules_(3), moduleZoR_(), moduleB_() {}

  void TrackerGeometryInfo::getTiltedModuleInfo(const Settings* settings,
                                                const TrackerTopology* theTrackerTopo,
                                                const TrackerGeometry* theTrackerGeom) {
    for (const GeomDetUnit* gd : theTrackerGeom->detUnits()) {
      float R0 = gd->position().perp();
      float Z0 = std::abs(gd->position().z());
      DetId detId = gd->geographicalId();

      bool barrel = detId.subdetId() == StripSubdetector::TOB || detId.subdetId() == StripSubdetector::TIB;
      enum BarrelModuleType { tiltedMinusZ = 1, tiltedPlusZ = 2, flat = 3 };
      bool tiltedBarrel = barrel && (theTrackerTopo->tobSide(detId) != BarrelModuleType::flat);
      bool isLower = theTrackerTopo->isLower(detId);

      // Calculate and store r, z, B at centre of titled modules
      // Only do this for the "lower" sensor
      if (tiltedBarrel && isLower) {
        double minZOfThisModule = std::abs(gd->surface().zSpan().first);
        double maxZOfThisModule = std::abs(gd->surface().zSpan().second);
        double minROfThisModule = gd->surface().rSpan().first;
        double maxROfThisModule = gd->surface().rSpan().second;

        // Calculate module tilt and B at centre of module
        double moduleTilt = atan(std::abs(maxZOfThisModule - minZOfThisModule) / (maxROfThisModule - minROfThisModule));
        double moduleTheta = atan(R0 / Z0);
        double BCorrection = std::abs(cos(std::abs(moduleTheta) - moduleTilt) / sin(moduleTheta));

        // Only store if this value of B has not been stored already
        bool storeThisBCorrection = true;
        for (unsigned int iCorr = 0; iCorr < moduleB_.size(); ++iCorr) {
          if (std::abs(BCorrection - moduleB_[iCorr]) < 0.0001) {
            storeThisBCorrection = false;
            break;
          }
        }

        if (storeThisBCorrection) {
          moduleZoR_.push_back(Z0 / R0);
          moduleB_.push_back(BCorrection);
        }
      }
    }
  }

}  // namespace tmtt
