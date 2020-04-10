#include "L1Trigger/TrackFindingTMTT/interface/ChiSquared4ParamsApprox.h"

namespace tmtt {

  ChiSquared4ParamsApprox::ChiSquared4ParamsApprox(const Settings* settings, const uint nPar)
      : L1ChiSquared(settings, nPar) {
    //parameterStream_ << "4Params_TrackletStyle_MCTruthSeed";
    //configParameters_ = (lsStr_.str());
    largestresid_ = -1.0;
    ilargestresid_ = -1;
  }

  TVectorD ChiSquared4ParamsApprox::seed(const L1track3D& l1track3D) {
    /* Cheat by using MC trutth to initialize helix parameters. Useful to check if conevrgence is the problem */
    TVectorD x(4);
    x[INVR] = getSettings()->invPtToInvR() * l1track3D.qOverPt();
    x[PHI0] = l1track3D.phi0();
    x[T] = l1track3D.tanLambda();
    x[Z0] = l1track3D.z0();
    return x;
  }

  TMatrixD ChiSquared4ParamsApprox::D(const TVectorD& x) {
    TMatrixD D(2 * stubs_.size(), nPar_);  // Empty matrix
    D.Zero();
    int j = 0;
    double rInv = x[INVR];
    double phi0 = x[PHI0];
    double t = x[T];
    for (unsigned i = 0; i < stubs_.size(); i++) {
      double ri = stubs_[i]->r();
      if (stubs_[i]->barrel()) {
        D(j, INVR) = -0.5 * ri * ri;  // Fine for now;
        D(j, PHI0) = ri;              // Fine
        //D(j, T);
        //D(j, Z0);
        j++;
        //D(j, INVR)
        //D(j, PHI0)
        D(j, T) = ri;  // ri; // Fine for now
        D(j, Z0) = 1;   // Fine
        j++;
      } else {
        //here we handle a disk hit
        //first we have the r position
        double phii = stubs_[i]->phi();
        int iphi = stubs_[i]->iphi();

        // N.B. These represent HALF the width and number of strips of sensor.
        double width = stubs_[i]->width() / 2.0;
        double nstrip = stubs_[i]->nstrip() / 2.0;

        double Deltai = width * (iphi - nstrip) / nstrip;  //A bit of a hack...

        if (stubs_[i]->z() > 0.0)
          Deltai = -Deltai;
        double DeltaiOverRi = Deltai / ri;
        double theta0 = (DeltaiOverRi) + 0.67 * (DeltaiOverRi) * (DeltaiOverRi) * (DeltaiOverRi);

        double phi_track = phi0 - 0.5 * rInv * ri;  //Expected phi hit given the track
        //std::cout << phi_track << "/" << phi0 << "/" << rInv << "/" << t << std::endl;

        double tInv = 1 / t;

        D(j, INVR) = -0.167 * ri * ri * ri * rInv;  // Tweaking of constant?
        D(j, PHI0) = 0;                             // Exact
        D(j, T) = -ri * tInv;                    // Fine;
        D(j, Z0) = -1 * tInv;                     // Fine
        j++;
        //second the rphi position
        D(j, INVR) = -0.5 * ri * ri;  // Needs fine tuning, was (phimultiplier*-0.5*(zi-z0)/t+rmultiplier*drdrinv);
        D(j, PHI0) = ri;              // Fine, originally phimultiplier
        D(j, T) = ri * 0.5 * rInv * ri * tInv - ((phi_track - phii) - theta0) * ri * tInv;
        D(j, Z0) = ri * 0.5 * rInv * tInv - ((phi_track - phii) - theta0) * tInv;
        j++;
      }
    }
    return D;
  }

  TMatrixD ChiSquared4ParamsApprox::Vinv() {
    TMatrixD Vinv(2 * stubs_.size(), 2 * stubs_.size());
    for (unsigned i = 0; i < stubs_.size(); i++) {
      if (stubs_[i]->barrel()) {
        Vinv(2 * i, 2 * i) = 1 / stubs_[i]->sigmaX();
        Vinv(2 * i + 1, 2 * i + 1) = 1 / stubs_[i]->sigmaZ();
      } else {
        Vinv(2 * i, 2 * i) = 1 / stubs_[i]->sigmaZ();
        Vinv(2 * i + 1, 2 * i + 1) = 1 / stubs_[i]->sigmaX();
      }
    }
    return Vinv;
  }

  TVectorD ChiSquared4ParamsApprox::residuals(const TVectorD& x) {
    unsigned int n = stubs_.size();

    TVectorD delta(2*n);

    double rInv = x[INVR];
    double phi0 = x[PHI0];
    double t = x[T];
    double z0 = x[Z0];

    double chiSq = 0.0;

    unsigned int j = 0;

    if (getSettings()->debug() == 6)
      std::cout << "Residuals (" << chiSq << ") [" << getSettings()->invPtToInvR() / rInv << "]: ";

    largestresid_ = -1.0;
    ilargestresid_ = -1;

    for (unsigned int i = 0; i < n; i++) {
      double ri = stubs_[i]->r();
      double zi = stubs_[i]->z();
      double phii = stubs_[i]->phi();
      const double sigmax = stubs_[i]->sigmaX();
      const double sigmaz = stubs_[i]->sigmaZ();

      if (stubs_[i]->barrel()) {
        //we are dealing with a barrel stub

        double halfRinvRi = 0.5 * ri * rInv;
        double aSinHalfRinvRi = halfRinvRi + 0.67 * halfRinvRi * halfRinvRi * halfRinvRi;

        double deltaphi = phi0 - aSinHalfRinvRi - phii;
        if (deltaphi > M_PI)
          deltaphi -= 2 * M_PI;
        if (deltaphi < -M_PI)
          deltaphi += 2 * M_PI;
        delta[j++] = (ri * deltaphi) / sigmax;  // TODO this is different from tracklet
        delta[j++] = (z0 + (2.0 / rInv) * t * aSinHalfRinvRi - zi) / sigmaz;
      } else {
        //we are dealing with a disk hit

        double tInv = 1 / t;

        double r_track = (zi - z0) * tInv;
        double phi_track = phi0 - 0.5 * rInv * (zi - z0) * tInv;
        int iphi = stubs_[i]->iphi();

        // N.B. These represent HALF the width and number of strips of sensor.
        double width = stubs_[i]->width() / 2.0;
        double nstrip = stubs_[i]->nstrip() / 2.0;

        double Deltai = width * (iphi - nstrip) / nstrip;  //A bit of a hack...

        if (stubs_[i]->z() > 0.0)
          Deltai = -Deltai;

        double DeltaiOverRi = Deltai / ri;
        double theta0 = (DeltaiOverRi) +
                        0.67 * (DeltaiOverRi) * (DeltaiOverRi) *
                            (DeltaiOverRi);  //+0.125*DeltaiOverRi*DeltaiOverRi*DeltaiOverRi*DeltaiOverRi*DeltaiOverRi;

        double Delta = Deltai - r_track * (theta0 - (phi_track - phii));

        delta[j++] = (r_track - ri) / sigmaz;
        delta[j++] = Delta / sigmax;
      }

      if (getSettings()->debug() == 6)
        std::cout << delta[j - 2] << " " << delta[j - 1] << " ";

      chiSq += delta[j - 2] * delta[j - 2] + delta[j - 1] * delta[j - 1];

      if (fabs(delta[j - 2]) > largestresid_) {
        largestresid_ = fabs(delta[j - 2]);
        ilargestresid_ = i;
      }

      if (fabs(delta[j - 1]) > largestresid_) {
        largestresid_ = fabs(delta[j - 1]);
        ilargestresid_ = i;
      }
      if (getSettings()->debug() == 6)
        std::cout << __FILE__ << ":" << __LINE__ << " - Residuals(): delta[" << j - 2 << "]/delta[" << j - 1
                  << "]: " << delta[j - 2] << "/" << delta[j - 1] << std::endl;
      if (getSettings()->debug() == 6)
        std::cout << __FILE__ << ":" << __LINE__ << " - Residuals(): chisq: " << chiSq << std::endl;
    }

    return delta;
  }

}  // namespace tmtt
