///=== Create requested track fitter

#include "L1Trigger/TrackFindingTMTT/interface/TrackFitFactory.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/ChiSquared4ParamsApprox.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFParamsComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/SimpleLR.h"
#ifdef USE_HLS
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFParamsCombCallHLS.h"
#endif
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;

namespace tmtt {

namespace trackFitFactory {

std::unique_ptr<TrackFitGeneric> create(std::string fitterName, const Settings* settings) {
    if (fitterName.compare("ChiSquared4ParamsApprox") == 0) {
      return std::make_unique<ChiSquared4ParamsApprox>(settings, 4);
    } else if (fitterName.compare("KF4ParamsComb") == 0) {
      return std::make_unique<KFParamsComb>(settings, 4, fitterName);
    } else if (fitterName.compare("KF5ParamsComb") == 0) {
	return std::make_unique<KFParamsComb>(settings, 5, fitterName);
    } else if (fitterName.compare("SimpleLR") == 0) {
      return std::make_unique<SimpleLR>(settings);
#ifdef USE_HLS
    } else if (fitterName.compare("KF4ParamsCombHLS") == 0) {
      return std::make_unique<KFParamsCombCallHLS>(settings, 4, fitterName);
    } else if (fitterName.compare("KF5ParamsCombHLS") == 0) {
      return std::make_unique<KFParamsCombCallHLS>(settings, 5, fitterName);
#endif
    } else {
      throw cms::Exception("BadConfig") << "TrackFitFactory: ERROR you requested unknown track fitterName: " << fitterName
                                        << endl;
    }
  }

}  // namespace trackFitFactory

}  // namespace tmtt
