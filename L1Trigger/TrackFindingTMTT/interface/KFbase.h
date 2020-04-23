///=== This is the base class for the Kalman Combinatorial Filter track fit algorithm.

#ifndef L1Trigger_TrackFindingTMTT_KFbase_h
#define L1Trigger_TrackFindingTMTT_KFbase_h

#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/TrackFitGeneric.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrack.h"
#include "L1Trigger/TrackFindingTMTT/interface/KalmanState.h"
#include <map>
#include <vector>
#include <fstream>
#include <TString.h>

class TH1F;
class TH2F;

namespace tmtt {

  class TP;
  class KalmanState;
  class StubCluster;

  class KFbase : public TrackFitGeneric {
  public:
    enum PAR_IDS { INV2R, PHI0, T, Z0, D0 };
    enum PAR_IDS_VAR { QOVERPT };
    enum MEAS_IDS { PHI, Z };
    enum OVERLAP_TYPE { TYPE_NORMAL, TYPE_V2, TYPE_NOCLUSTERING, TYPE_TP };

  public:
    KFbase(const Settings *settings, const uint nPar, const std::string &fitterName = "", const uint nMeas = 2);

    virtual ~KFbase() {
      this->resetStates();
      this->deleteStubClusters();
    }

    KFbase(const KFbase &kf) = delete;  // Disable copy & move of this class.
    KFbase(KFbase &&kf) = delete;
    KFbase &operator=(const KFbase &kf) = delete;
    KFbase &operator=(KFbase &&kf) = delete;

    L1fittedTrack fit(const L1track3D &l1track3D);

  protected:
    static std::vector<double> trackParams(const KFbase *p, const KalmanState *state);
    virtual std::vector<double> trackParams(const KalmanState *state) const = 0;

    // Get track params with beam-spot constraint & chi2 (r-phi) after applying it..
    virtual std::vector<double> trackParams_BeamConstr(const KalmanState *state, double &chi2rphi_bcon) const {
      chi2rphi_bcon = 0.0;
      return (this->trackParams(state));  // Returns unconstrained result, unless derived class overrides it.
    }

    double sectorPhi() const {
      float phiCentreSec0 =
          -M_PI / float(settings_->numPhiNonants()) + M_PI / float(settings_->numPhiSectors());
      return 2. * M_PI * float(iCurrentPhiSec_) / float(settings_->numPhiSectors()) + phiCentreSec0;
    }
    //bool kalmanUpdate( const StubCluster *stubCluster, KalmanState &state, KalmanState &new_state, const TP *tpa );
    virtual const KalmanState *kalmanUpdate(
        unsigned skipped, unsigned layer, const StubCluster *stubCluster, const KalmanState &state, const TP *);
    void resetStates();
    void deleteStubClusters();
    const KalmanState *mkState(const L1track3D &candidate,
                               unsigned skipped,
                               unsigned layer,
                               unsigned layerId,
                               const KalmanState *last_state,
                               const std::vector<double> &x,
                               const TMatrixD &pxx,
                               const TMatrixD &K,
                               const TMatrixD &dcov,
                               const StubCluster *stubCluster,
                               double chi2rphi,
                               double chi2rz);

  protected:
    /* Methods */
    std::vector<double> Hx(const TMatrixD &pH, const std::vector<double> &x) const;
    std::vector<double> Fx(const TMatrixD &pF, const std::vector<double> &x) const;
    TMatrixD HxxH(const TMatrixD &pH, const TMatrixD &xx) const;
    void deltaChi2(const TMatrixD &dcov,
                      const std::vector<double> &delta,
                      bool debug,
                      double &delChi2rphi,
                      double &delChi2rz) const;
    TMatrixD GetKalmanMatrix(const TMatrixD &h, const TMatrixD &pxcov, const TMatrixD &dcov) const;
    void GetAdjustedState(const TMatrixD &K,
                          const TMatrixD &pxcov,
                          const std::vector<double> &x,
                          const StubCluster *stubCluster,
                          const std::vector<double> &delta,
                          std::vector<double> &new_x,
                          TMatrixD &new_xcov) const;

    virtual std::vector<double> seedx(const L1track3D &l1track3D) const = 0;
    virtual TMatrixD seedP(const L1track3D &l1track3D) const = 0;
    virtual void barrelToEndcap(double r,
                                const StubCluster *stubCluster,
                                std::vector<double> &x,
                                TMatrixD &cov_x) const {}
    virtual std::vector<double> d(const StubCluster *stubCluster) const = 0;
    virtual TMatrixD H(const StubCluster *stubCluster) const = 0;
    virtual TMatrixD F(const StubCluster *stubCluster = 0, const KalmanState *state = 0) const = 0;
    virtual TMatrixD PddMeas(const StubCluster *stubCluster, const KalmanState *state) const = 0;

    virtual std::vector<double> residual(const StubCluster *stubCluster,
                                         const std::vector<double> &x,
                                         double candQoverPt) const;
    virtual const KalmanState *updateSeedWithStub(const KalmanState &state, const StubCluster *stubCluster) {
      return 0;
    }
    virtual bool isGoodState(const KalmanState &state) const { return true; }

    virtual void calcChi2(const KalmanState &state, double &chi2rphi, double &chi2rz) const;

    virtual unsigned int kalmanLayer(
        unsigned int iEtaReg, unsigned int layerIDreduced, bool barrel, float r, float z) const;
    virtual bool kalmanAmbiguousLayer(unsigned int iEtaReg, unsigned int kfLayer);

    std::vector<const KalmanState *> doKF(const L1track3D &l1track3D,
                                          const std::vector<const StubCluster *> &stubClusters,
                                          const TP *tpa);

    void printTPSummary(std::ostream &os, const TP *tp, bool addReturn = true) const;
    void printTP(std::ostream &os, const TP *tp) const;
    void printStubLayers(std::ostream &os, std::vector<const Stub *> &stubs, unsigned int iEtaReg) const;
    void printStubCluster(std::ostream &os, const StubCluster *stubCluster, bool addReturn = true) const;
    void printStubClusters(std::ostream &os, std::vector<const StubCluster *> &stubClusters) const;
    void printStub(std::ostream &os, const Stub *stub, bool addReturn = true) const;
    void printStubs(std::ostream &os, std::vector<const Stub *> &stubs) const;

    double deltaRphiForClustering(unsigned layerId, unsigned endcapRing);
    double deltaRForClustering(unsigned endcapRing);
    bool isOverlap(const Stub *a, const Stub *b, OVERLAP_TYPE type);

    std::set<unsigned> kalmanDeadLayers(bool &remove2PSCut) const;

    // Function to calculate approximation for tilted barrel modules (aka B) copied from Stub class.
    float approxB(float z, float r) const;

    // Is this HLS code?
    virtual bool isHLS() { return false; };

  protected:
    unsigned nPar_;
    unsigned nMeas_;
    unsigned numEtaRegions_;

    std::vector<KalmanState *> state_list_;
    std::vector<StubCluster *> stbcl_list_;

    unsigned int iCurrentPhiSec_;
    unsigned int iCurrentEtaReg_;
    unsigned int iLastPhiSec_;
    unsigned int iLastEtaReg_;

    unsigned int numUpdateCalls_;

    const TP *tpa_;
  };

}  // namespace tmtt

#endif
