///=== This is the global Linear Regression for 4 helix parameters track fit algorithm.

///=== Written by: Davide Cieri

#include "L1Trigger/TrackFindingTMTT/interface/SimpleLR.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrack.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include <vector>
#include <set>
#include <algorithm>
#include <limits>

using namespace std;

namespace tmtt {

  SimpleLR::SimpleLR(const Settings* settings) : TrackFitGeneric(settings) {
    // Initialize digitization parameters
    phiMult_ = pow(2., getSettings()->phiSBits()) / getSettings()->phiSRange();
    rTMult_ = pow(2., getSettings()->rtBits()) / getSettings()->rtRange();
    zMult_ = pow(2., getSettings()->zBits()) / getSettings()->zRange();
    z0Mult_ = pow(2., getSettings()->slr_z0Bits()) / getSettings()->slr_z0Range();
    phiTMult_ = pow(2., getSettings()->slr_phi0Bits()) / getSettings()->slr_phi0Range();

    qOverPtMult_ = pow(2., getSettings()->slr_oneOver2rBits()) / getSettings()->slr_oneOver2rRange();
    tanLambdaMult_ = pow(2., getSettings()->slr_tanlambdaBits()) / getSettings()->slr_tanlambdaRange();
    chi2Mult_ = pow(2., getSettings()->slr_chisquaredBits()) / getSettings()->slr_chisquaredRange();

    numeratorPtMult_ = rTMult_ * phiMult_;
    numeratorPhiMult_ = rTMult_ * rTMult_ * phiMult_;
    numeratorZ0Mult_ = rTMult_ * rTMult_ * z0Mult_;
    numeratorLambdaMult_ = rTMult_ * z0Mult_;
    denominatorMult_ = rTMult_ * rTMult_;
    resMult_ = rTMult_ * qOverPtMult_;

    digitize_ = getSettings()->digitizeSLR() and getSettings()->enableDigitize();
    dividerBitsHelix_ = getSettings()->dividerBitsHelix();
    dividerBitsHelixZ_ = getSettings()->dividerBitsHelixZ();
    shiftingBitsDenRPhi_ = getSettings()->ShiftingBitsDenRPhi();
    shiftingBitsDenRZ_ = getSettings()->ShiftingBitsDenRZ();

    shiftingBitsPhi_ = getSettings()->ShiftingBitsPhi();
    shiftingBitsz0_ = getSettings()->ShiftingBitsZ0();
    shiftingBitsPt_ = getSettings()->ShiftingBitsPt();
    shiftingBitsLambda_ = getSettings()->ShiftingBitsLambda();

    phiSectorWidth_ = 2. * M_PI / float(getSettings()->numPhiSectors());
    phiNonantWidth_ = 2. * M_PI / float(getSettings()->numPhiNonants());

    chi2cut_ = getSettings()->slr_chi2cut();
    chosenRofPhi_ = getSettings()->chosenRofPhi();
    if (digitize_) chosenRofPhi_ = floor(chosenRofPhi_ * rTMult_) / rTMult_;
  };

  static bool pair_compare(std::pair<const Stub*, float> a, std::pair<const Stub*, float> b) {
    return (a.second < b.second);
  }

  L1fittedTrack SimpleLR::fit(const L1track3D& l1track3D) {
    if (getSettings()->debug() == 6)
      cout << "=============== FITTING TRACK ====================" << endl;

    invPtToDPhi_ = -getSettings()->invPtToDphi();

    double phiCentreSec0 = -0.5 * phiNonantWidth_ + 0.5 * phiSectorWidth_;
    phiSectorCentre_ = phiSectorWidth_ * double(l1track3D.iPhiSec()) - phiCentreSec0;

    if (digitize_)
      phiSectorCentre_ = floor(phiSectorCentre_ * phiTMult_) / phiTMult_;

    // Inizialise track fit parameters
    double qOverPt = 0.;
    double phiT = 0.;
    double phi0 = 0.;
    double z0 = 0.;
    double zT = 0.;
    double tanLambda = 0.;

    // Inizialise Sums
    double SumRPhi = 0.;
    double SumR = 0.;
    double SumPhi = 0.;
    double SumR2 = 0.;
    double SumRZ = 0.;
    double SumZ = 0.;

    unsigned int numStubs = 0;
    // Calc helix parameters on Rphi Plane (STEP 1)
    // This loop calculates the sums needed to calculate the numerators and the denominator to compute the helix parameters in the R-Phi plane (q/pT, phiT)
    for (const Stub* stub : l1track3D.getStubs()) {
      // if((const_cast<Stub*>(stub))->psModule()){
      numStubs++;

      if (digitize_) {
        (const_cast<Stub*>(stub))->digitizeForHTinput(l1track3D.iPhiSec());
        (const_cast<Stub*>(stub))->digitizeForSFinput();
        const DigitalStub digiStub = (const_cast<Stub*>(stub))->digitalStub();

        SumRPhi = SumRPhi + digiStub.rt() * digiStub.phiS();
        SumR = SumR + digiStub.rt();
        SumPhi = SumPhi + digiStub.phiS();
        SumR2 = SumR2 + digiStub.rt() * digiStub.rt();
        if (getSettings()->debug() == 6)
          cout << "phiS " << digiStub.iDigi_PhiS() << " rT " << digiStub.iDigi_Rt() << " z " << digiStub.iDigi_Z()
               << endl;
      } else {
        float phi = 0;
        if (l1track3D.iPhiSec() == 0 and (const_cast<Stub*>(stub))->phi() > 0) {
          phi = (const_cast<Stub*>(stub))->phi() - 2 * M_PI;
        } else if (l1track3D.iPhiSec() == getSettings()->numPhiSectors() and (const_cast<Stub*>(stub))->phi() < 0) {
          phi = (const_cast<Stub*>(stub))->phi() + 2 * M_PI;
        } else {
          phi = (const_cast<Stub*>(stub))->phi();
        }
        SumRPhi = SumRPhi + (const_cast<Stub*>(stub))->r() * phi;
        SumR = SumR + (const_cast<Stub*>(stub))->r();
        SumPhi = SumPhi + phi;
        SumR2 = SumR2 + (const_cast<Stub*>(stub))->r() * (const_cast<Stub*>(stub))->r();
        if (getSettings()->debug() == 6)
          cout << "phi " << phi << " r " << (const_cast<Stub*>(stub))->r() << " z " << (const_cast<Stub*>(stub))->z()
               << endl;
      }
      // }
    }

    double numeratorPt, digiNumeratorPt;
    double denominator, digiDenominator;
    double numeratorPhi, digiNumeratorPhi;
    double reciprocal, digiReciprocal;
    double numeratorZ0, digiNumeratorZ0;
    double numeratorLambda, digiNumeratorLambda;

    digiNumeratorPt = (numStubs * SumRPhi - SumR * SumPhi);
    digiDenominator = (numStubs * SumR2 - SumR * SumR);
    digiNumeratorPhi = (SumR2 * SumPhi - SumR * SumRPhi);

    if (!digitize_) {
      qOverPt = (numStubs * SumRPhi - SumR * SumPhi) / (numStubs * SumR2 - SumR * SumR);
      phi0 = (SumR2 * SumPhi - SumR * SumRPhi) / (numStubs * SumR2 - SumR * SumR);
    } else {
      digiNumeratorPt /= pow(2., shiftingBitsPt_);
      digiNumeratorPt = floor(digiNumeratorPt * numeratorPtMult_);
      numeratorPt = digiNumeratorPt / numeratorPtMult_;

      digiNumeratorPhi /= pow(2., shiftingBitsPhi_);
      digiNumeratorPhi = floor(digiNumeratorPhi * numeratorPhiMult_);
      numeratorPhi = digiNumeratorPhi / numeratorPhiMult_;

      digiDenominator /= pow(2., shiftingBitsDenRPhi_);
      digiDenominator = (floor(digiDenominator * denominatorMult_) + 0.5);
      denominator = digiDenominator / denominatorMult_;
      digiReciprocal = (pow(2., dividerBitsHelix_) - 1) / (denominator);  // To be moved
      digiReciprocal = floor(digiReciprocal / denominatorMult_);
      reciprocal = digiReciprocal * denominatorMult_;

      qOverPt = numeratorPt * reciprocal / pow(2., dividerBitsHelix_ + shiftingBitsDenRPhi_ - shiftingBitsPt_);
      phiT = numeratorPhi * reciprocal / pow(2., dividerBitsHelix_ + shiftingBitsDenRPhi_ - shiftingBitsPhi_);

      qOverPt = floor(qOverPt * qOverPtMult_) / (qOverPtMult_);
      phiT = floor(phiT * phiTMult_) / phiTMult_;
    }

    if (getSettings()->debug() == 6 and digitize_)
      cout << setw(10) << "First Helix parameters: qOverPt = " << qOverPt << " (" << floor(qOverPt * qOverPtMult_)
           << "), phiT = " << phiT << " (" << floor(phiT * phiTMult_) << ") " << endl;

    if (getSettings()->debug() == 6 and !digitize_)
      cout << "First Helix Parameters: qOverPt = " << qOverPt << " phi0 " << phi0 << endl;

    // ================== RESIDUAL CALCULATION ON RPHI ========================
    std::vector<std::pair<const Stub*, double> > vRes;
    unsigned int psStubs = 0;
    for (const Stub* stub : l1track3D.getStubs()) {
      if ((const_cast<Stub*>(stub))->psModule())
        psStubs++;
      double ResPhi;

      if (digitize_) {
        (const_cast<Stub*>(stub))->digitizeForHTinput(l1track3D.iPhiSec());
        (const_cast<Stub*>(stub))->digitizeForSFinput();
        const DigitalStub digiStub = (const_cast<Stub*>(stub))->digitalStub();

        ResPhi =
            digiStub.iDigi_PhiS() * pow(2., shiftingBitsDenRPhi_ - shiftingBitsPt_) -
            floor(phiT * phiTMult_) *
                pow(2., shiftingBitsDenRPhi_ - shiftingBitsPt_ - getSettings()->slr_phi0Bits() + getSettings()->phiSBits()) -
            floor(qOverPt * qOverPtMult_) * digiStub.iDigi_Rt();

        if (getSettings()->debug() == 6) {
          // cout << "floor(phiT*phiTMult_) " << floor(phiT*phiTMult_) << endl;
          // cout << "dsp_PhiSPhiT "<< digiStub.iDigi_PhiS() - floor(phiT*phiTMult_) << " shift_right(dsp_QoverPt_Rt,divider_shift- divider_pt) "<< floor(qOverPt*qOverPtMult_)*digiStub.iDigi_Rt()/(pow(2., ShiftingBits_- shiftingBitsPt_ ) )<< " getSettings()->rtRange() "<< getSettings()->rtRange() << endl;
          cout << "DIGI RESIDUAL " << ResPhi << endl;
        }
        ResPhi = floor(ResPhi) / resMult_;
      }

      else {
        ResPhi = reco::deltaPhi((const_cast<Stub*>(stub))->phi(), phi0 + qOverPt * (const_cast<Stub*>(stub))->r());
      }

      double Res = fabs(ResPhi);
      // if(digitize_) Res = floor(Res*phiMult_)/phiMult_;

      std::pair<const Stub*, double> ResStubPair(stub, Res);
      vRes.push_back(ResStubPair);
      if (getSettings()->debug() == 6) {
        if (const_cast<Stub*>(stub)->assocTP() != nullptr)
          cout << " Stub Residual " << Res << " TP " << const_cast<Stub*>(stub)->assocTP()->index() << endl;
        else
          cout << " Stub Residual " << Res << " TP nullptr" << endl;
      }
    }

    double LargResidual = 9999.;
    // Find largest residuals
    while (vRes.size() > 4 and LargResidual > getSettings()->ResidualCut()) {
      std::vector<std::pair<const Stub*, double> >::iterator maxResIt =
          max_element(vRes.begin(), vRes.end(), pair_compare);
      LargResidual = (*maxResIt).second;
      if (getSettings()->debug() == 6)
        cout << "Largest Residual " << LargResidual << endl;

      if (LargResidual > getSettings()->ResidualCut()) {
        if ((*maxResIt).first->psModule()) {
          if (psStubs > 2) {
            if (getSettings()->debug() == 6)
              cout << "removing PS residual " << (*maxResIt).second << endl;
            vRes.erase(maxResIt);
            psStubs--;
          } else {
            if (getSettings()->debug() == 6)
              cout << "residual " << (*maxResIt).second << " set to -1. " << endl;
            (*maxResIt).second = -1.;
          }
        } else {
          vRes.erase(maxResIt);
          if (getSettings()->debug() == 6)
            cout << "removing residual " << (*maxResIt).second << endl;
        }
      }
    }

    std::vector<const Stub*> fitStubs;
    for (std::pair<const Stub*, double> ResStubPair : vRes) {
      fitStubs.push_back(ResStubPair.first);
    }

    phiT = 0.;
    zT = 0.;

    SumRPhi = 0.;
    SumR = 0.;
    SumPhi = 0.;
    SumR2 = 0.;
    SumRZ = 0.;
    SumZ = 0.;
    double SumR_ps = 0.;
    double SumR2_ps = 0.;

    numStubs = 0;
    psStubs = 0;

    for (const Stub* stub : fitStubs) {
      if ((const_cast<Stub*>(stub))->psModule())
        psStubs++;

      numStubs++;
      if (digitize_) {
        (const_cast<Stub*>(stub))->digitizeForHTinput(l1track3D.iPhiSec());
        // (const_cast<Stub*>(stub))->digitizeForSFinput();
        const DigitalStub digiStub = (const_cast<Stub*>(stub))->digitalStub();
        SumRPhi += digiStub.rt() * digiStub.phiS();
        SumR += digiStub.rt();
        SumPhi += digiStub.phiS();
        SumR2 += digiStub.rt() * digiStub.rt();
        if ((const_cast<Stub*>(stub))->psModule()) {
          SumRZ += digiStub.rt() * digiStub.z();
          SumZ += digiStub.z();
          SumR_ps += digiStub.rt();
          SumR2_ps += digiStub.rt() * digiStub.rt();
        }
        if (getSettings()->debug() == 6) {
          cout << "phiS " << digiStub.iDigi_PhiS() << " rT " << digiStub.iDigi_Rt() << " z " << digiStub.iDigi_Z()
               << endl;
        }
      } else {
        float phi = 0;
        if (l1track3D.iPhiSec() == 0 and (const_cast<Stub*>(stub))->phi() > 0) {
          phi = (const_cast<Stub*>(stub))->phi() - 2 * M_PI;
        } else if (l1track3D.iPhiSec() == getSettings()->numPhiSectors() and (const_cast<Stub*>(stub))->phi() < 0) {
          phi = (const_cast<Stub*>(stub))->phi() + 2 * M_PI;
        } else {
          phi = (const_cast<Stub*>(stub))->phi();
        }

        SumRPhi += (const_cast<Stub*>(stub))->r() * phi;
        SumR += (const_cast<Stub*>(stub))->r();
        SumPhi += phi;
        SumR2 += (const_cast<Stub*>(stub))->r() * (const_cast<Stub*>(stub))->r();
        if ((const_cast<Stub*>(stub))->psModule()) {
          SumRZ += (const_cast<Stub*>(stub))->r() * (const_cast<Stub*>(stub))->z();
          SumZ += (const_cast<Stub*>(stub))->z();
          SumR_ps += (const_cast<Stub*>(stub))->r();
          SumR2_ps += (const_cast<Stub*>(stub))->r() * (const_cast<Stub*>(stub))->r();
        }
        if (getSettings()->debug() == 6)
          cout << "phi " << phi << " r " << (const_cast<Stub*>(stub))->r() << " z " << (const_cast<Stub*>(stub))->z()
               << endl;
      }

      // }
    }

    numeratorZ0 = (SumR2_ps * SumZ - SumR_ps * SumRZ);
    numeratorLambda = (psStubs * SumRZ - SumR_ps * SumZ);
    numeratorPt = (numStubs * SumRPhi - SumR * SumPhi);
    denominator = (numStubs * SumR2 - SumR * SumR);
    double denominatorZ = (psStubs * SumR2_ps - SumR_ps * SumR_ps);
    numeratorPhi = (SumR2 * SumPhi - SumR * SumRPhi);
    double reciprocalZ;
    if (!digitize_) {
      z0 = numeratorZ0 / denominatorZ;
      tanLambda = numeratorLambda / denominatorZ;
      qOverPt = (numStubs * SumRPhi - SumR * SumPhi) / (numStubs * SumR2 - SumR * SumR);
      phi0 = (SumR2 * SumPhi - SumR * SumRPhi) / (numStubs * SumR2 - SumR * SumR);
    } else {
      numeratorPt /= pow(2., shiftingBitsPt_);
      numeratorPt = floor(numeratorPt * numeratorPtMult_) / numeratorPtMult_;

      numeratorPhi /= pow(2., shiftingBitsPhi_);
      numeratorPhi = floor(numeratorPhi * numeratorPhiMult_) / numeratorPhiMult_;

      numeratorLambda /= pow(2., shiftingBitsLambda_);
      numeratorLambda = floor(numeratorLambda * numeratorLambdaMult_) / numeratorLambdaMult_;

      numeratorZ0 /= pow(2., shiftingBitsz0_);
      numeratorZ0 = floor(numeratorZ0 * numeratorZ0Mult_) / numeratorZ0Mult_;

      denominator /= pow(2., shiftingBitsDenRPhi_);
      denominator = (floor(denominator * denominatorMult_) + 0.5) / denominatorMult_;
      reciprocal = (pow(2., dividerBitsHelix_) - 1) / (denominator);
      reciprocal = floor(reciprocal / denominatorMult_) * denominatorMult_;

      denominatorZ /= pow(2., shiftingBitsDenRZ_);
      denominatorZ = (floor(denominatorZ * denominatorMult_) + 0.5) / denominatorMult_;
      reciprocalZ = (pow(2., dividerBitsHelixZ_) - 1) / (denominatorZ);
      reciprocalZ = floor(reciprocalZ / denominatorMult_) * denominatorMult_;

      qOverPt = numeratorPt * reciprocal / pow(2., dividerBitsHelix_ + shiftingBitsDenRPhi_ - shiftingBitsPt_);
      phiT = numeratorPhi * reciprocal / pow(2., dividerBitsHelix_ + shiftingBitsDenRPhi_ - shiftingBitsPhi_);

      tanLambda =
          numeratorLambda * reciprocalZ / pow(2., dividerBitsHelixZ_ + shiftingBitsDenRZ_ - shiftingBitsLambda_);
      zT = numeratorZ0 * reciprocalZ / pow(2., dividerBitsHelixZ_ + shiftingBitsDenRZ_ - shiftingBitsz0_);

      phi0 = phiSectorCentre_ + phiT - qOverPt * getSettings()->chosenRofPhi();
      z0 = zT - tanLambda * getSettings()->chosenRofPhi();

      qOverPt = floor(qOverPt * qOverPtMult_) / qOverPtMult_;
      phiT = floor(phiT * phiTMult_) / phiTMult_;
    }

    // cout << "z0"

    // qOverPt /= -invPtToDPhi_;

    if (getSettings()->debug() == 6 and digitize_) {
      cout << "HT mbin " << int(l1track3D.getCellLocationHT().first) - 16 << " cbin "
           << int(l1track3D.getCellLocationHT().second) - 32 << " iPhi " << l1track3D.iPhiSec() << " iEta "
           << l1track3D.iEtaReg() << endl;
      cout << "Second Helix variables: numeratorPt = " << numeratorPt << ", numeratorPhi = " << numeratorPhi
           << ", numeratorZ0 = " << numeratorZ0 << " numeratorLambda = " << numeratorLambda
           << " denominator =  " << denominator << " reciprocal = " << reciprocal << " denominatorZ =  " << denominatorZ
           << " reciprocalZ = " << reciprocalZ << endl;
      cout << setw(10) << "Final Helix parameters: qOverPt = " << qOverPt << " (" << floor(qOverPt * qOverPtMult_)
           << "), phiT = " << phiT << " (" << floor(phiT * phiTMult_) << "), zT = " << zT << " (" << floor(zT * z0Mult_)
           << "), tanLambda = " << tanLambda << " (" << floor(tanLambda * tanLambdaMult_) << ")"
           << " z0 " << z0 << endl;
    } else if (getSettings()->debug() == 6) {
      cout << setw(10) << "Final Helix parameters: qOverPt = " << qOverPt << ", phi0 = " << phi0 << ", z0 = " << z0
           << ", tanLambda = " << tanLambda << endl;
    }

    double chi2_phi = 0.;
    double chi2_z = 0.;

    for (const Stub* stub : fitStubs) {
      double ResPhi = 0.;
      double ResZ = 0.;
      if (digitize_) {
        (const_cast<Stub*>(stub))->digitizeForHTinput(l1track3D.iPhiSec());
        (const_cast<Stub*>(stub))->digitizeForSFinput();
        const DigitalStub digiStub = (const_cast<Stub*>(stub))->digitalStub();
        ResPhi = digiStub.phiS() - phiT - qOverPt * digiStub.rt();
        ResZ = digiStub.z() - zT - tanLambda * digiStub.rt();
        // ResZ = digiStub.z() - z0 - tanLambda*digiStub.r();
      } else {
        ResPhi = reco::deltaPhi((const_cast<Stub*>(stub))->phi(), phi0 + qOverPt * (const_cast<Stub*>(stub))->r());
        ResZ = (const_cast<Stub*>(stub))->z() - z0 - tanLambda * (const_cast<Stub*>(stub))->r();
      }

      double RPhiSigma = 0.0002;
      float RZSigma = (const_cast<Stub*>(stub))->zErr() + fabs(tanLambda) * (const_cast<Stub*>(stub))->rErr();

      if (not(const_cast<Stub*>(stub))->barrel())
        RPhiSigma = 0.0004;

      if (digitize_) {
        RPhiSigma = floor(RPhiSigma * phiMult_) / phiMult_;
      }

      // if(!(const_cast<Stub*>(stub))->psModule()) RZSigma = 5;

      ResPhi /= RPhiSigma;
      // cout << "zT "<< zT << " tanLambda "<< t<< RZSigma << endl;
      ResZ /= RZSigma;

      chi2_phi += fabs(ResPhi * ResPhi);
      chi2_z += fabs(ResZ * ResZ);
      if (getSettings()->debug() == 6) {
        cout << "Stub ResPhi " << ResPhi * RPhiSigma << " ResSigma " << RPhiSigma << " Res " << ResPhi << " chi2 "
             << chi2_phi << endl;
        cout << "Stub ResZ " << ResZ * RZSigma << " ResSigma " << RZSigma << " Res " << ResZ << " chi2 " << chi2_z
             << endl;
      }
    }
    qOverPt /= invPtToDPhi_;

    bool accepted = false;

    //double chi2 = chi2_phi;  // Ignore r-z residuals due to poor 2S resolution?
    double chi2 = chi2_phi + chi2_z;
    if (digitize_)
      chi2 = floor(chi2 * chi2Mult_) / chi2Mult_;

    // cout << "chi2 "<< chi2 << " phi "<< chi2_phi << " z "<< chi2_z << endl;
    float dof = 2 * fitStubs.size() - 4;
    float chi2dof = chi2 / dof;
    if (chi2 < chi2cut_)
      accepted = true;

    if (getSettings()->debug() == 6)
      cout << "qOverPt " << qOverPt << " phiT " << phiT << endl;

    // This condition can only happen if cfg param TrackFitCheat = True.
    if (fitStubs.size() < 4)
      accepted = false;

    // Kinematic cuts -- NOT YET IN FIRMWARE!!!
    constexpr float tolerance = 0.1;
    if (fabs(qOverPt) > 1. / (getSettings()->houghMinPt() - tolerance))
      accepted = false;
    if (fabs(z0) > 20.)
      accepted = false;

    if (accepted) {
      
      // Create the L1fittedTrack object
      const unsigned int hitPattern = 0;  // FIX: Needs setting
      L1fittedTrack fitTrk(getSettings(), l1track3D, fitStubs, hitPattern, qOverPt, 0., phi0, z0, tanLambda, chi2_phi, chi2_z, 4);

      if (getSettings()->enableDigitize())
	fitTrk.digitizeTrack("SimpleLR");

      if (getSettings()->debug() == 6 and digitize_) {
	cout << "Digitized parameters " << endl;
	cout << "HT mbin " << int(l1track3D.getCellLocationHT().first) - 16 << " cbin "
	     << int(l1track3D.getCellLocationHT().second) - 32 << " iPhi " << l1track3D.iPhiSec() << " iEta "
	     << l1track3D.iEtaReg() << endl;
	cout << setw(10) << "First Helix parameters: qOverPt = " << fitTrk.qOverPt() << " oneOver2r "
	     << fitTrk.digitaltrack().oneOver2r() << " (" << floor(fitTrk.digitaltrack().oneOver2r() * qOverPtMult_)
	     << "), phi0 = " << fitTrk.digitaltrack().phi0() << " (" << fitTrk.digitaltrack().iDigi_phi0rel()
	     << "), zT = " << zT << " (" << floor(zT * z0Mult_) << "), tanLambda = " << tanLambda << " ("
	     << floor(tanLambda * tanLambdaMult_) << ")" << endl;
      }

      if (getSettings()->debug() == 6) {
	cout << "FitTrack helix parameters " << int(fitTrk.getCellLocationFit().first) - 16 << ", "
	     << int(fitTrk.getCellLocationFit().second) - 32 << " HT parameters "
	     << int(fitTrk.getCellLocationHT().first) - 16 << ", " << int(fitTrk.getCellLocationHT().second) - 32 << endl;

	if (fitTrk.getMatchedTP() != nullptr) {
	  cout << "VERY GOOD! " << chi2dof << endl;
	  cout << "TP qOverPt " << fitTrk.getMatchedTP()->qOverPt() << " phi0 " << fitTrk.getMatchedTP()->phi0() << endl;
	  if (!accepted)
	    cout << "BAD CHI2 " << chi2 << " chi2/ndof " << chi2dof << endl;
	} else {
	  cout << "FAKE TRACK!!! " << chi2 << " chi2/ndof " << chi2dof << endl;
	  if (l1track3D.getMatchedTP() != nullptr)
	    cout << "was good" << endl;
	}
	cout << "layers in track " << fitTrk.getNumLayers() << endl;
      }

      return fitTrk;

    } else {

      L1fittedTrack rejectedTrk;
      return rejectedTrk;

    }
  }

}  // namespace tmtt
