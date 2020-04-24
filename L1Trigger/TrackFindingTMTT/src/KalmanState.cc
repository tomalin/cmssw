#include "L1Trigger/TrackFindingTMTT/interface/KalmanState.h"
#include "L1Trigger/TrackFindingTMTT/interface/StubCluster.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include <TMatrixD.h>

using namespace std;

namespace tmtt {

  KalmanState::KalmanState()
      : settings_(nullptr),
        kLayerNext_(0),
        layerId_(0),
        last_state_(nullptr),
        stubCluster_(nullptr),
        chi2rphi_(0),
        chi2rz_(0),
        barrel_(true),
        n_skipped_(0),
        hitPattern_(0) {}

  KalmanState::KalmanState(const Settings* settings,
			   const L1track3D &candidate,
                           unsigned n_skipped,
                           unsigned kLayer_next,
                           unsigned layerId,
                           const KalmanState *last_state,
                           const std::vector<double> &x,
                           const TMatrixD &pxx,
                           const TMatrixD &K,
                           const TMatrixD &dcov,
                           const StubCluster *stubCluster,
                           double chi2rphi,
                           double chi2rz)
      : settings_(settings),
        kLayerNext_(kLayer_next),
        layerId_(layerId),
        last_state_(last_state),
        xa_(x),
        stubCluster_(stubCluster),
        chi2rphi_(chi2rphi),
        chi2rz_(chi2rz),
        n_skipped_(n_skipped),
        l1track3D_(candidate) {
    pxxa_.Clear();
    pxxa_.ResizeTo(pxx.GetNrows(), pxx.GetNcols());
    pxxa_ = pxx;
    K_.ResizeTo(K.GetNrows(), K.GetNcols());
    K_ = K;
    dcov_.ResizeTo(dcov.GetNrows(), dcov.GetNcols());
    dcov_ = dcov;
    kalmanChi2RphiScale_ = settings_->kalmanChi2RphiScale();

    hitPattern_ = 0;
    if (last_state != nullptr)
      hitPattern_ = last_state->hitPattern();  // Bit encoded list of hit layers
    if (stubCluster != nullptr)
      hitPattern_ |= (1 << (stubCluster->layerKF()));

    r_ = 0.1;
    z_ = 0;
    barrel_ = true;
    endcapRing_ = 0;

    if (stubCluster != nullptr) {
      r_ = stubCluster->r();
      z_ = stubCluster->z();
      barrel_ = stubCluster->barrel();
      endcapRing_ = stubCluster->endcapRing();
    }

    n_stubs_ = kLayerNext_ - n_skipped_;
  }

  KalmanState::KalmanState(const KalmanState &p)
      : settings_(p.settings()),
        kLayerNext_(p.nextLayer()),
        layerId_(p.layerId()),
        endcapRing_(p.endcapRing()),
        r_(p.r()),
        z_(p.z()),
        last_state_(p.last_state()),
        xa_(p.xa()),
        pxxa_(p.pxxa()),
        K_(p.K()),
        dcov_(p.dcov()),
        stubCluster_(p.stubCluster()),
        chi2rphi_(p.chi2rphi()),
        chi2rz_(p.chi2rz()),
        n_stubs_(p.nStubLayers()),
        barrel_(p.barrel()),
        n_skipped_(p.nSkippedLayers()),
        l1track3D_(p.candidate()) {}

  KalmanState &KalmanState::operator=(const KalmanState &other) {
    if (&other == this)
      return *this;

    settings_ = other.settings();
    kLayerNext_ = other.nextLayer();
    layerId_ = other.layerId();
    endcapRing_ = other.endcapRing();
    r_ = other.r();
    z_ = other.z();
    last_state_ = other.last_state();
    xa_ = other.xa();
    pxxa_ = other.pxxa();
    K_ = other.K();
    dcov_ = other.dcov();
    stubCluster_ = other.stubCluster();
    chi2rphi_ = other.chi2rphi();
    chi2rz_ = other.chi2rz();
    n_stubs_ = other.nStubLayers();
    barrel_ = other.barrel();
    n_skipped_ = other.nSkippedLayers();
    l1track3D_ = other.candidate();
    return *this;
  }

  bool KalmanState::good(const TP *tp) const {
    const KalmanState *state = this;
    while (state) {
      const StubCluster *stubCluster = state->stubCluster();
      if (stubCluster) {
        const set<const TP *> &tps = stubCluster->assocTPs();

        if (tps.find(tp) == tps.end())
          return false;
      }
      state = state->last_state();
    }
    return true;
  }

  double KalmanState::reducedChi2() const {
    if (2 * n_stubs_ - xa_.size() > 0)
      return (this->chi2()) / (2 * n_stubs_ - xa_.size());
    else
      return 0;
  }

  const KalmanState *KalmanState::last_update_state() const {
    const KalmanState *state = this;
    while (state) {
      if (state->stubCluster())
        return state;
      state = state->last_state();
    }
    return 0;
  }

  std::vector<const Stub *> KalmanState::stubs() const {
    std::vector<const Stub *> all_stubs;

    const KalmanState *state = this;
    while (state) {
      const StubCluster *stbcl = state->stubCluster();
      if (stbcl) {
        const std::vector<const Stub *>& stubs = stbcl->stubs();
        all_stubs.insert(all_stubs.end(), stubs.begin(), stubs.end());
      }
      state = state->last_state();
    }
    std::reverse(all_stubs.begin(), all_stubs.end());  // Put innermost stub first.
    return all_stubs;
  }

  bool KalmanState::order(const KalmanState *left, const KalmanState *right) {
    return (left->nStubLayers() > right->nStubLayers());
  }

  bool KalmanState::orderMinSkipChi2(const KalmanState *left, const KalmanState *right) {
    return (left->chi2scaled() * (left->nSkippedLayers() + 1) < right->chi2scaled() * (right->nSkippedLayers() + 1));
  }

  bool KalmanState::orderChi2(const KalmanState *left, const KalmanState *right) {
    return (left->chi2scaled() < right->chi2scaled());
  }
}  // namespace tmtt
