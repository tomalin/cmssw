///=== This is the Kalman Combinatorial Filter for 4 helix parameters track fit algorithm.

#include "L1Trigger/TrackFindingTMTT/interface/KFParamsComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/KalmanState.h"
#include "DataFormats/Math/interface/deltaPhi.h"

using namespace std;

namespace tmtt {

/* Get physical helix params */

  TVectorD KFParamsComb::trackParams(const KalmanState* state) const {
    TVectorD y(nPar_);
    TVectorD x = state->vectorX();
    y[QOVERPT] = 2. * x(INV2R) / settings_->invPtToInvR();
    y[PHI0] = reco::deltaPhi(x(PHI0) + sectorPhi(), 0.);
    y[Z0] = x(Z0);
    y[T] = x(T);
    if (nPar_ == 5) {
      y[D0] = x(D0);
    }
    return y;
  }

  /* If using 5 param helix fit, get track params with beam-spot constraint & track fit chi2 from applying it. */
  /* (N.B. chi2rz unchanged by constraint) */

  TVectorD KFParamsComb::trackParams_BeamConstr(const KalmanState* state, double& chi2rphi) const {
    if (nPar_ == 5) {
      TVectorD y(nPar_);
      TVectorD x = state->vectorX();
      TMatrixD matC = state->matrixC();
      double delChi2rphi = (x(D0) * x(D0)) / matC[D0][D0];
      chi2rphi = state->chi2rphi() + delChi2rphi;
      // Apply beam-spot constraint to helix params in transverse plane only, as most sensitive to it.
      x[INV2R] -= x(D0) * (matC[INV2R][D0] / matC[D0][D0]);
      x[PHI0] -= x(D0) * (matC[PHI0][D0] / matC[D0][D0]);
      x[D0] = 0.0;
      y[QOVERPT] = 2. * x(INV2R) / settings_->invPtToInvR();
      y[PHI0] = reco::deltaPhi(x(PHI0) + sectorPhi(), 0.);
      y[Z0] = x(Z0);
      y[T] = x(T);
      y[D0] = x(D0);
      return y;
    } else {
      return (this->trackParams(state));
    }
  }

/* The Kalman measurement matrix = derivative of helix intercept w.r.t. helix params */
/* Here I always measure phi(r), and z(r) */

  TMatrixD KFParamsComb::matrixH(const Stub* stub) const {
    TMatrixD h(2, nPar_);
    double r = stub->r();
    h(PHI, INV2R) = -r;
    h(PHI, PHI0) = 1;
    if (nPar_ == 5) {
      h(PHI, D0) = -1. / r;
    }
    h(Z, Z0) = 1;
    h(Z, T) = r;
    return h;
  }

  /* Helix state seed  */

  TVectorD KFParamsComb::seedX(const L1track3D& l1track3D) const {
    TVectorD x(nPar_);
    x[INV2R] = settings_->invPtToInvR() * l1track3D.qOverPt() / 2;
    x[PHI0] = reco::deltaPhi(l1track3D.phi0() - sectorPhi(), 0.);
    x[Z0] = l1track3D.z0();
    x[T] = l1track3D.tanLambda();
    if (nPar_ == 5) {
      x[D0] = l1track3D.d0();
    }

    return x;
  }

  /* Helix state seed covariance matrix */

  TMatrixD KFParamsComb::seedC(const L1track3D& l1track3D) const {
    TMatrixD p(nPar_, nPar_);

    double invPtToInv2R = settings_->invPtToInvR() / 2;

    // Assumed track seed (from HT) uncertainty in transverse impact parameter.
    const float d0Sigma = 1.0;

    if (settings_->hybrid()) {
      p(INV2R, INV2R) = 0.0157 * 0.0157 * invPtToInv2R * invPtToInv2R * 4;
      p(PHI0, PHI0) = 0.0051 * 0.0051 * 4;
      p(Z0, Z0) = 5.0 * 5.0;
      p(T, T) = 0.25 * 0.25 * 4;
      // N.B. (z0, tanL, d0) seed uncertainties could be smaller for hybrid, if seeded in PS? -- not tried
      //if (l1track3D.seedPS() > 0) { // Tracklet seed used PS layers
      //  p(Z0,Z0) /= (4.*4.).;
      //  p(T,T) /= (4.*4.);
      // }
      if (nPar_ == 5) {
        p(D0, D0) = d0Sigma * d0Sigma;
      }

    } else {
      // optimised for 18x2 with additional error factor in pt/phi to avoid pulling towards wrong HT params
      p(INV2R, INV2R) = 0.0157 * 0.0157 * invPtToInv2R * invPtToInv2R * 4;  // Base on HT cell size
      p(PHI0, PHI0) = 0.0051 * 0.0051 * 4;                                  // Based on HT cell size.
      p(Z0, Z0) = 5.0 * 5.0;
      p(T, T) = 0.25 * 0.25 * 4;  // IRT: increased by factor 4, as was affecting fit chi2.
      if (nPar_ == 5) {
        p(D0, D0) = d0Sigma * d0Sigma;
      }

      if (settings_->numEtaRegions() <= 12) {
        // Inflate eta errors
        p(T, T) = p(T, T) * 2 * 2;
      }
    }

    return p;
  }

/* The forecast matrix (identity matrix in this KF formulation) */

  TMatrixD KFParamsComb::matrixF(const Stub* stub, const KalmanState* state) const {
    const TMatrixD unitMatrix(TMatrixD::kUnit, TMatrixD(nPar_, nPar_));
    return unitMatrix;
  }

  /* Stub position measurements in (phi,z) */

  TVectorD KFParamsComb::vectorM(const Stub* stub) const {
    TVectorD meas(2);
    meas[PHI] = reco::deltaPhi(stub->phi(), sectorPhi());
    meas[Z] = stub->z();
    return meas;
  }

  // Stub position resolution in (phi,z)

  TMatrixD KFParamsComb::matrixV(const Stub* stub, const KalmanState* state) const {
    double inv2R =
        (settings_->invPtToInvR()) * 0.5 * state->candidate().qOverPt();  // alternatively use state->vectorX()(INV2R)
    double inv2R2 = inv2R * inv2R;

    double tanl = state->vectorX()(T);  // factor of 0.9 improves rejection
    double tanl2 = tanl * tanl;

    TMatrixD p(2, 2);

    double vphi(0);
    double vz(0);
    double vcorr(0);

    // consider error due to integerisation only for z (r in encap) coord when enabled
    double err_digi2(0);
    if (settings_->enableDigitize())
      err_digi2 = 0.15625 * 0.15625 / 12.0;

    double a = stub->sigmaX() * stub->sigmaX();
    double b = stub->sigmaZ() * stub->sigmaZ() + err_digi2;
    double r2 = stub->r() * stub->r();
    double invr2 = 1. / r2;

    // Scattering term scaling as 1/Pt.
    double sigmaScat = settings_->kalmanMultiScattTerm() / (state->candidate().pt());
    double sigmaScat2 = sigmaScat * sigmaScat;

    if (stub->barrel()) {
      vphi = (a * invr2) + sigmaScat2;

      if (stub->tiltedBarrel()) {
        // Convert uncertainty in (r,phi) to (z,phi).
        float scaleTilted = 1.;
        if (settings_->kalmanHOtilted()) {
          if (settings_->useApproxB()) {  // Simple firmware approximation
            scaleTilted = approxB(stub->z(), stub->r());
          } else {  // Exact C++ implementation.
            float tilt = stub->moduleTilt();
            scaleTilted = sin(tilt) + cos(tilt) * tanl;
          }
        }
        float scaleTilted2 = scaleTilted * scaleTilted;
        // This neglects the non-radial strip effect, assumed negligeable for PS.
        vz = b * scaleTilted2;
      } else {
        vz = b;
      }

      if (settings_->kalmanHOdodgy()) {
        // Use original (Dec. 2016) dodgy implementation was this.
        vz = b;
      }

    } else {
      vphi = a * invr2 + sigmaScat2;
      vz = (b * tanl2);

      if (not stub->psModule()) {  // Neglect these terms in PS
        double beta = 0.;
        // Add correlation term related to conversion of stub residuals from (r,phi) to (z,phi).
        if (settings_->kalmanHOprojZcorr() == 2)
          beta += -inv2R;
        // Add alpha correction for non-radial 2S endcap strips..
        if (settings_->kalmanHOalpha() == 2)
          beta += -stub->alpha();  // alpha is 0 except in endcap 2S disks

        double beta2 = beta * beta;
        vphi += b * beta2;
        vcorr = b * (beta * tanl);

        // IRT - for checking efficiency of removing phi-z correlation from projection.
        // "ultimate_off1"
        //vphi  = a * invr2 + b * pow(-stub->alpha(), 2) + b * inv2R2 + sigmaScat2;
        //vcorr = b * ((-stub->alpha()) * tanl);

        // IRT - This higher order correction doesn't significantly improve the track fit performance, so commented out.
        //if (settings_->kalmanHOhelixExp()) {
        //  float dsByDr = 1. + (1./2.)*r2*inv2R2; // Allows for z = z0 + s*tanL, where s is not exactly r due to circle.
        //  vcorr *= dsByDr;
        //  vz *= dsByDr * dsByDr;
        //}

        if (settings_->kalmanHOdodgy()) {
          // Use original (Dec. 2016) dodgy implementation was this.
          vphi = (a * invr2) + (b * inv2R2) + sigmaScat2;
          vcorr = 0.;
          vz = (b * tanl2);
        }
      }
    }

    p(PHI, PHI) = vphi;
    p(Z, Z) = vz;
    p(PHI, Z) = vcorr;
    p(Z, PHI) = vcorr;

    return p;
  }

/* Check if helix state passes cuts */

  bool KFParamsComb::isGoodState(const KalmanState& state) const {
    // Cut values. (Layer 0 entry here is dummy). -- todo : make configurable

    vector<float> z0Cut, ptTolerance, d0Cut, chi2Cut;
    //  Layer   =    0      1      2     3     4      5      6
    ptTolerance = {999., 999., 0.1, 0.1, 0.05, 0.05, 0.05};
    d0Cut = {999., 999., 999., 10., 10., 10., 10.};  // Only used for 5 param fit.
    if (nPar_ == 5) {                                // specific cuts for displaced tracking case.
      //  Layer   =    0      1        2         3         4         5           6
      z0Cut = {
          999., 999., 1.7 * 15., 1.7 * 15., 1.7 * 15., 1.7 * 15., 1.7 * 15.};  // Larger values require digisation change.
      chi2Cut = {999., 999., 10., 30., 80., 120., 160.};                       // Maybe loosen for high d0 ?
    } else {  // specific cuts for prompt tracking case.
      //  Layer   =    0      1      2     3     4      5      6
      z0Cut = {999., 999., 15., 15., 15., 15., 15.};
      chi2Cut = {999., 999., 10., 30., 80., 120., 160.};
    }

    unsigned nStubLayers = state.nStubLayers();
    bool goodState(true);

    TVectorD y = trackParams(&state);
    double qOverPt = y[QOVERPT];
    double pt = std::abs(1 / qOverPt);
    double z0 = std::abs(y[Z0]);

    // state parameter selections

    if (z0 > z0Cut[nStubLayers])
      goodState = false;
    if (pt < settings_->houghMinPt() - ptTolerance[nStubLayers])
      goodState = false;
    if (nPar_ == 5) {
      double d0 = std::abs(state.vectorX()[D0]);
      if (d0 > d0Cut[nStubLayers])
        goodState = false;
    }

    // chi2 selection

    double chi2scaled = state.chi2scaled();  // chi2(r-phi) scaled down to improve electron performance.

    if (settings_->kalmanMultiScattTerm() > 0.0001) {  // Scattering taken into account

      if (chi2scaled > chi2Cut[nStubLayers])
        goodState = false;  // No separate pT selection needed

    } else {  // scattering ignored - HISTORIC

      // N.B. Code below changed by Alexander Morton to allow tracking down to Pt = 2 GeV.
      if (nStubLayers == 2) {
        if (chi2scaled > 15.0)
          goodState = false;  // No separate pT selection needed
      } else if (nStubLayers == 3) {
        if (chi2scaled > 100.0 && pt > 2.7)
          goodState = false;
        if (chi2scaled > 120.0 && pt <= 2.7)
          goodState = false;
      } else if (nStubLayers == 4) {
        if (chi2scaled > 320.0 && pt > 2.7)
          goodState = false;
        if (chi2scaled > 1420.0 && pt <= 2.7)
          goodState = false;
      } else if (nStubLayers == 5) {  // NEEDS TUNING FOR 5 OR 6 LAYERS !!!
        if (chi2scaled > 480.0 && pt > 2.7)
          goodState = false;
        if (chi2scaled > 2130.0 && pt <= 2.7)
          goodState = false;
      } else if (nStubLayers >= 6) {  // NEEDS TUNING FOR 5 OR 6 LAYERS !!!
        if (chi2scaled > 640.0 && pt > 2.7)
          goodState = false;
        if (chi2scaled > 2840.0 && pt <= 2.7)
          goodState = false;
      }
    }

    const bool countUpdateCalls = false;  // Print statement to count calls to Updator.

    if (countUpdateCalls || (settings_->kalmanDebugLevel() >= 2 && tpa_ != nullptr) ||
        (settings_->kalmanDebugLevel() >= 2 && settings_->hybrid())) {
      if (not goodState)
        cout << "State veto:";
      if (goodState)
        cout << "State kept:";
      cout << " nlay=" << nStubLayers << " nskip=" << state.nSkippedLayers() << " chi2_scaled=" << chi2scaled;
      if (tpa_ != nullptr)
        cout << " pt(mc)=" << tpa_->pt();
      cout << " pt=" << pt << " q/pt=" << qOverPt << " tanL=" << y[T] << " z0=" << y[Z0] << " phi0=" << y[PHI0];
      if (nPar_ == 5)
        cout << " d0=" << y[D0];
      cout << " fake" << (tpa_ == nullptr);
      if (tpa_ != nullptr)
        cout << " pt(mc)=" << tpa_->pt();
      cout << endl;
    }

    return goodState;
  }

}  // namespace tmtt
