#ifndef L1Trigger_TrackFindingTMTT_KalmanState_h
#define L1Trigger_TrackFindingTMTT_KalmanState_h

#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFbase.h"
#include <map>

namespace tmtt {

  class KFbase;
  class KalmanState;
  class StubCluster;
  class Settings;

  class KalmanState {
  public:
    KalmanState();
    KalmanState(const Settings* settings,
		const L1track3D &candidate,
                unsigned n_skipped,
                unsigned kLayer_next,
                unsigned layerId,
                const KalmanState *last_state,
                const std::vector<double> &x,
                const TMatrixD &pxx,
                const TMatrixD &K,
                const TMatrixD &dcov,
                const StubCluster *stubcl,
                double chi2rphi,
                double chi2rz);

    KalmanState(const KalmanState &p);
    ~KalmanState() {}

    KalmanState &operator=(const KalmanState &other);

    const Settings* settings() const {return settings_;}
    unsigned nextLayer() const { return kLayerNext_; }
    unsigned layerId() const { return layerId_; }
    unsigned endcapRing() const { return endcapRing_; }
    bool barrel() const { return barrel_; }
    unsigned nSkippedLayers() const { return n_skipped_; }
    // Hit coordinates.
    double r() const { return r_; }
    double z() const { return z_; }
    const KalmanState *last_state() const { return last_state_; }
    // Helix parameters (1/2R, phi relative to sector, z0, tanLambda)
    const std::vector<double> &xa() const { return xa_; }
    // Covariance matrix on helix params.
    const TMatrixD &pxxa() const { return pxxa_; }
    // Kalman Gain matrix
    const TMatrixD &K() const { return K_; }
    // Hit position covariance matrix.
    const TMatrixD &dcov() const { return dcov_; }
    // Hit
    const StubCluster *stubCluster() const { return stubCluster_; }
    double chi2() const { return chi2rphi_ + chi2rz_; }
    double chi2scaled() const { return chi2rphi_ / kalmanChi2RphiScale_ + chi2rz_; }  // Improves electron performance.
    double chi2rphi() const { return chi2rphi_; }
    double chi2rz() const { return chi2rz_; }
    unsigned nStubLayers() const { return n_stubs_; }
    const L1track3D &candidate() const { return l1track3D_; }
    unsigned int hitPattern() const { return hitPattern_; }  // Bit-encoded KF layers the fitted track has stubs in.

    bool good(const TP *tp) const;
    double reducedChi2() const;
    const KalmanState *last_update_state() const;
    std::vector<const Stub *> stubs() const;

    static bool orderChi2(const KalmanState *left, const KalmanState *right);
    static bool orderMinSkipChi2(const KalmanState *left, const KalmanState *right);

    static bool order(const KalmanState *left, const KalmanState *right);
    void setChi2(double chi2rphi, double chi2rz) {
      chi2rphi_ = chi2rphi;
      chi2rz_ = chi2rz;
    }

    // If using HLS, note/get additional output produced by HLS core.
    //void setHLSselect(unsigned int mBinHelix, unsigned int cBinHelix, bool consistent) { mBinHelixHLS_ = mBinHelix; cBinHelixHLS_ = cBinHelix; consistentHLS_ = consistent;}
    //void getHLSselect(unsigned int& mBinHelix, unsigned int& cBinHelix, bool& consistent) const { mBinHelix = mBinHelixHLS_; cBinHelix = cBinHelixHLS_; consistent = consistentHLS_;}

  private:
    const Settings* settings_;
    unsigned kLayerNext_;
    unsigned layerId_;
    unsigned endcapRing_;
    double r_;
    double z_;
    const KalmanState *last_state_;
    std::vector<double> xa_;
    TMatrixD pxxa_;
    TMatrixD K_;
    TMatrixD dcov_;
    const StubCluster *stubCluster_;
    double chi2rphi_;
    double chi2rz_;
    unsigned int kalmanChi2RphiScale_;
    unsigned n_stubs_;
    bool barrel_;
    unsigned n_skipped_;
    L1track3D l1track3D_;
    unsigned int hitPattern_;

    // Additional output from HLS if using it.
    unsigned int mBinHelixHLS_;
    unsigned int cBinHelixHLS_;
    bool consistentHLS_;
  };

}  // namespace tmtt

#endif
