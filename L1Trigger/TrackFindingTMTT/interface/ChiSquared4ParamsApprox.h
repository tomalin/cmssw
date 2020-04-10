#ifndef L1Trigger_TrackFindingTMTT_ChiSquared4ParamsApprox_h
#define L1Trigger_TrackFindingTMTT_ChiSquared4ParamsApprox_h

#include "L1Trigger/TrackFindingTMTT/interface/L1ChiSquared.h"

namespace tmtt {

  class ChiSquared4ParamsApprox : public L1ChiSquared {
  public:
    ChiSquared4ParamsApprox(const Settings* settings, const uint nPar);

  protected:
    TVectorD seed(const L1track3D& l1track3D);
    TVectorD residuals(const TVectorD& x);
    TMatrixD D(const TVectorD& x);
    TMatrixD Vinv();
  };

}  // namespace tmtt

#endif
