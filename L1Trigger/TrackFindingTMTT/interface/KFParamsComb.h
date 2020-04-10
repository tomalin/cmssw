///=== This is the Kalman Combinatorial Filter for 4 & 5 helix parameters track fit algorithm.

#ifndef L1Trigger_TrackFindingTMTT_KFParamsComb_h
#define L1Trigger_TrackFindingTMTT_KFParamsComb_h

#include "L1Trigger/TrackFindingTMTT/interface/L1KalmanComb.h"
#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"

namespace tmtt {

  class KFParamsComb : public L1KalmanComb {
  public:
    KFParamsComb(const Settings* settings, const uint nPar, const std::string& fitterName)
        : L1KalmanComb(settings, nPar, fitterName) {}

    virtual ~KFParamsComb() {}

  protected:
    virtual std::vector<double> getTrackParams(const KalmanState* state) const;
    virtual std::vector<double> getTrackParams_BeamConstr(const KalmanState* state, double& chi2rphi) const;
    virtual std::vector<double> seedx(const L1track3D& l1track3D) const;
    virtual TMatrixD seedP(const L1track3D& l1track3D) const;
    virtual std::vector<double> d(const StubCluster* stubCluster) const;
    virtual TMatrixD H(const StubCluster* stubCluster) const;
    virtual TMatrixD dH(const StubCluster* stubCluster) const;
    virtual TMatrixD F(const StubCluster* stubCluster = 0, const KalmanState* state = 0) const;
    virtual TMatrixD PxxModel(const KalmanState* state, const StubCluster* stubCluster) const;
    virtual TMatrixD PddMeas(const StubCluster* stubCluster, const KalmanState* state) const;
    virtual bool isGoodState(const KalmanState& state) const;
  };

}  // namespace tmtt

#endif
