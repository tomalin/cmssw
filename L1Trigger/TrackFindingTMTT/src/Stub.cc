#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
#include "L1Trigger/TrackFindingTMTT/interface/PrintL1trk.h"

#include <iostream>

using namespace std;

namespace tmtt {

  //=== Store useful info about the stub (for use with HYBRID code), with hard-wired constants to allow use outside CMSSW.

  Stub::Stub(double phi,
             double r,
             double z,
             double bend,
             unsigned int layerId,
             bool psModule,
             bool barrel,
             unsigned int iphi,
             double alpha,
             const Settings* settings,
             const TrackerTopology* trackerTopology,
             unsigned int ID,
             unsigned int iPhiSec)
      : phi_(phi),
        r_(r),
        z_(z),
        bend_(bend),
        iphi_(iphi),
        alpha_(alpha),
        digitalStub_(std::make_unique<DigitalStub>(settings, r, phi, z, iPhiSec)),
        stubWindowSuggest_(settings),
        psModule_(psModule),
        layerId_(layerId),
        layerIdReduced_(ModuleInfo::calcLayerIdReduced(layerId)),
        endcapRing_(0),
        barrel_(barrel)
  {  //work in progress on better constructor for new hybrid
    if (psModule && barrel) {
      // max z at which non-tilted modules are found in 3 barrel PS layers. (Element 0 not used).
      const vector<float>& zMax = settings->zMaxNonTilted();
      tiltedBarrel_ = (std::abs(z) > zMax[layerId]);
    } else {
      tiltedBarrel_ = false;
    }
    if (psModule) {
      stripPitch_ = settings->psPixelPitch();
      stripLength_ = settings->psPixelLength();
      nStrips_ = settings->psNStrips();
    } else {
      stripPitch_ = settings->ssStripPitch();
      stripLength_ = settings->ssStripLength();
      nStrips_ = settings->ssNStrips();
    }
    index_in_vStubs_ = ID;  // A unique ID to label the stub.
  }

  //=== Store useful info about stub (for use with TMTT tracking).

  Stub::Stub(const TTStubRef& ttStubRef,
             unsigned int index_in_vStubs,
             const Settings* settings,
             const TrackerGeometry* trackerGeometry,
             const TrackerTopology* trackerTopology,
	     const ModuleInfo* moduleInfo,
	     const StubKiller* stubKiller)
      : ttStubRef_(ttStubRef),
        settings_(settings),
        index_in_vStubs_(index_in_vStubs),
        assocTP_(nullptr),  // Initialize in case job is using no MC truth info.
        digitizedForGPinput_(false),  // Has stub has been digitized for GP input?
        digitizedForHTinput_(false),  // Has stub been digitized for HT input?
        digitizedForSForTFinput_(""),  // Has stub been digitized for SF/TF input?
        digitizeWarningsOn_(true),
        moduleInfo_(moduleInfo),  // Info about tracker module containing stub 
        stubWindowSuggest_(settings, trackerTopology),  // TMTT recommended stub window sizes.
        degradeBend_(trackerTopology),                   // Used to degrade stub bend information.
  // Module related variables (need to be stored for Hybrid)
    psModule_(moduleInfo->psModule()),
    layerId_(moduleInfo->layerId()),
    layerIdReduced_(moduleInfo->layerIdReduced()),
    endcapRing_(moduleInfo->endcapRing()),
    barrel_(moduleInfo->barrel()),
    tiltedBarrel_(moduleInfo->tiltedBarrel()),
    stripPitch_(moduleInfo->stripPitch()),
    stripLength_(moduleInfo->stripLength()),
    nStrips_(moduleInfo->nStrips())
  {
    // Get coordinates of stub.
    const TTStub<Ref_Phase2TrackerDigi_>* ttStubP = ttStubRef_.get();

    const PixelGeomDetUnit* specDet = moduleInfo_->specDet();
    const PixelTopology* specTopol = moduleInfo_->specTopol();
    MeasurementPoint measurementPoint = ttStubRef_->clusterRef(0)->findAverageLocalCoordinatesCentered();
    LocalPoint clustlp = specTopol->localPosition(measurementPoint);
    GlobalPoint pos = specDet->surface().toGlobal(clustlp);

    phi_ = pos.phi();
    r_ = pos.perp();
    z_ = pos.z();

    if (r_ < settings_->trackerInnerRadius() || r_ > settings_->trackerOuterRadius() ||
        std::abs(z_) > settings_->trackerHalfLength()) {
      throw cms::Exception("BadConfig") << "Stub: Stub found outside assumed tracker volume. Please update tracker "
                                           "dimensions specified in Settings.h!"
                                        << " r=" << r_ << " z=" << z_ << " id=" << moduleInfo_->detId().subdetId();
    }

    // Get the coordinates of the two clusters that make up this stub, measured in units of strip pitch, and measured
    // in the local frame of the sensor. They have a granularity  of 0.5*pitch.
    for (unsigned int iClus = 0; iClus <= 1; iClus++) {  // Loop over two clusters in stub.
      localU_cluster_[iClus] = ttStubP->clusterRef(iClus)->findAverageLocalCoordinatesCentered().x();
      localV_cluster_[iClus] = ttStubP->clusterRef(iClus)->findAverageLocalCoordinatesCentered().y();
    }

    // Get location of stub in module in units of strip number (or pixel number along finest granularity axis).
    // Range from 0 to (nStrips - 1) inclusive.
    // N.B. Since iphi is integer, this degrades the granularity by a factor 2. This seems silly, but track fit wants it.
    iphi_ = localU_cluster_[0];  // granularity 1*strip (unclear why we want to degrade it ...)

    // Determine alpha correction for non-radial strips in endcap 2S modules.
    // (If true hit at larger r than stub r by deltaR, then stub phi needs correcting by +alpha*deltaR).
    alpha_ = 0.;
    if ((not barrel()) && (not psModule())) {
      float fracPosInModule = (float(iphi_) - 0.5 * float(nStrips())) / float(nStrips());
      float phiRelToModule = sensorWidth() * fracPosInModule / r_;
      if (z_ < 0)
        phiRelToModule *= -1;
      if (moduleInfo_->outerModuleAtSmallerR())
        phiRelToModule *= -1;  // Module flipped.
      // If true hit at larger r than stub r by deltaR, then stub phi needs correcting by +alpha*deltaR.
      alpha_ = -phiRelToModule / r_;
    }

    // Calculate variables giving ratio of track intercept angle to stub bend.
    this->calcDphiOverBend();

    // Get stub bend that is available in front-end electronics, where bend is displacement between
    // two hits in stubs in units of strip pitch.
    bendInFrontend_ = ttStubRef_->bendFE();
    if ((not barrel()) && pos.z() > 0)
      bendInFrontend_ *= -1;
    if (barrel())
      bendInFrontend_ *= -1;

    // Get stub bend that is available in off-detector electronics, allowing for degredation of
    // bend resolution due to bit encoding by FE chip if required.
    bool rejectStub = false;  // indicates if bend is outside window assumed in DegradeBend.h
    numMergedBend_ = 1;       // Number of bend values merged into single degraded one.
    if (settings->degradeBendRes() == 2) {
      float degradedBend;  // degraded bend
      this->degradeResolution(
          bendInFrontend_, degradedBend, rejectStub, numMergedBend_);  // sets value of last 3 arguments.
      bend_ = degradedBend;
    } else if (settings->degradeBendRes() == 1) {
      bend_ = ttStubRef_->bendBE();  // Degraded bend from official CMS recipe.
      if ((not barrel()) && pos.z() > 0)
        bend_ *= -1;
      if (barrel())
        bend_ *= -1;
    } else {
      bend_ = bendInFrontend_;
    }

    // Fill frontendPass_ flag, indicating if frontend readout electronics will output this stub.
    this->setFrontend(rejectStub, stubKiller);

    // Calculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub.
    this->calcQoverPtrange();

    // Initialize class used to produce digital version of stub, with original stub parameters pre-digitization.
    digitalStub_ = std::make_unique<DigitalStub>(settings_, 
		      phi_,
                      r_,
                      z_,
                      min_qOverPt_bin_,
                      max_qOverPt_bin_,
                      layerId(),
                      moduleInfo_->layerIdReduced(),
                      bend_,
                      stripPitch(),
                      moduleInfo_->sensorSpacing(),
                      barrel(),
                      tiltedBarrel(),
		      psModule());

    // Update recommended stub window sizes that TMTT recommends that CMS should use in FE electronics.
    if (settings_->printStubWindows())
      stubWindowSuggest_.process(this);

    // Initialize truth info to false in case job is using no MC truth info.
    for (unsigned int iClus = 0; iClus <= 1; iClus++) {
      assocTPofCluster_[iClus] = nullptr;
    }
  }

  //=== Calculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub.

  void Stub::calcQoverPtrange() {
    // First determine bin range along q/Pt axis of HT array
    const int nbinsPt =
        (int)settings_->houghNbinsPt();  // Use "int" as nasty things happen if multiply "int" and "unsigned int".
    const int min_array_bin = 0;
    const int max_array_bin = nbinsPt - 1;
    // Now calculate range of q/Pt bins allowed by bend filter.
    float qOverPtMin = this->qOverPtOverBend() * (this->bend() - this->bendRes());
    float qOverPtMax = this->qOverPtOverBend() * (this->bend() + this->bendRes());
    int houghNbinsPt = settings_->houghNbinsPt();
    const float houghMaxInvPt = 1. / settings_->houghMinPt();
    float qOverPtBinSize = (2. * houghMaxInvPt) / houghNbinsPt;
    if (settings_->shape() == 2 || settings_->shape() == 1 || settings_->shape() == 3)  // Non-square HT cells.
      qOverPtBinSize = 2. * houghMaxInvPt / (houghNbinsPt - 1);
    // Convert to bin number along q/Pt axis of HT array.
    // N.B. For square HT cells, setting "tmp = -0.5" causeas cell to be accepted if q/Pt at its centre is consistent
    // with the stub bend. Instead using "tmp = 0.0" accepts cells if q/Pt at any point in cell is consistent with bend.
    // So if you use change from -0.5 to 0.0, you have to tighten the bend cut (by ~0.05) to get similar performance.
    // Decision to set tmp = 0.0 taken in softare & GP firmware on 9th August 2016.
    //float tmp = ( settings_->shape() == 2 || settings_->shape() == 1 || settings_->shape() == 3 ) ? 1. : -0.5;

    float tmp = (settings_->shape() == 2 || settings_->shape() == 1 || settings_->shape() == 3) ? 1. : 0.;
    int min_bin = std::floor(-tmp + (qOverPtMin + houghMaxInvPt) / qOverPtBinSize);
    int max_bin = std::floor(tmp + (qOverPtMax + houghMaxInvPt) / qOverPtBinSize);

    // Limit it to range of HT array.
    min_bin = max(min_bin, min_array_bin);
    max_bin = min(max_bin, max_array_bin);
    // If min_bin > max_bin at this stage, it means that the Pt estimated from the bend is below the cutoff for track-finding.
    // Keep min_bin > max_bin, so such stubs can be rejected, but set both variables to values inside the HT bin range.
    if (min_bin > max_bin) {
      min_bin = max_array_bin;
      max_bin = min_array_bin;
    }
    min_qOverPt_bin_ = (unsigned int)min_bin;
    max_qOverPt_bin_ = (unsigned int)max_bin;
  }

  //=== Digitize stub for input to Geographic Processor, with digitized phi coord. measured relative to closest phi sector.
  //=== (This approximation is valid if their are an integer number of digitisation bins inside each phi nonant).
  //=== However, you should also call digitizeForHTinput() before accessing digitized stub data, even if you only care about that going into GP! Otherwise, you will not identify stubs assigned to more than one nonant.

  void Stub::digitizeForGPinput(unsigned int iPhiSec) {
    if (settings_->enableDigitize()) {
      // Save CPU by not redoing digitization if stub was already digitized for this phi sector.
      if (!(digitizedForGPinput_ && digitalStub_->iGetNonant(iPhiSec) == digitalStub_->iDigi_Nonant())) {
        // Digitize
        digitalStub_->makeGPinput(iPhiSec);

        // Replace stub coordinates with those degraded by digitization process.
        phi_ = digitalStub_->phi();
        r_ = digitalStub_->r();
        z_ = digitalStub_->z();
        bend_ = digitalStub_->bend();

        // If the Stub class contains any data members that are not input to the GP, but are derived from variables that
        // are, then be sure to update these here too, unless Stub.h uses the check*() functions to declare them invalid.

        // Update variables giving ratio of track intercept angle to stub bend.
        this->calcDphiOverBend();

        // Note that stub has been digitized for GP input
        digitizedForGPinput_ = true;
      }
      digitizedForHTinput_ = false;
    }
  }

  //=== Digitize stub for input to Hough transform, with digitized phi coord. measured relative to specified phi sector.

  void Stub::digitizeForHTinput(unsigned int iPhiSec) {
    if (settings_->enableDigitize()) {
      // Save CPU by not redoing digitization if stub was already digitized for this phi sector.
      if (!(digitizedForHTinput_ && iPhiSec == digitalStub_->iDigi_PhiSec())) {
        // Call digitization for GP in case not already done. (Needed for variables that are common to GP & HT).
        this->digitizeForGPinput(iPhiSec);

        // Digitize
        digitalStub_->makeHTinput(iPhiSec);

        // Since GP and HT use same digitisation in r and z, don't bother updating their values.
        // (Actually, the phi digitisation boundaries also match, except for systolic array, so could skip updating phi too).

        // Replace stub coordinates and bend with those degraded by digitization process. (Don't bother with r & z, as already done by GP digitisation).
        phi_ = digitalStub_->phi();

        // Recalculate bin range along q/Pt axis of r-phi Hough transform array
        // consistent with bend of this stub, since it depends on r & z which have now been digitized.
        // (This recalculation should really be done in DigitalStub::makeHTinput(), but too lazy to move it there ...).
        this->calcQoverPtrange();

        // If the Stub class contains any data members that are not input to the HT, but are derived from variables that
        // are, then be sure to update these here too, unless Stub.h uses the check*() functions to declare them invalid.
        // - currently none.

        // Note that stub has been digitized.
        digitizedForHTinput_ = true;
      }
    }
  }

  //=== Digitize stub for input to r-z Seed Filter or Track Fitter.
  //=== Argument is "SeedFilter" or name of Track Fitter.

  void Stub::digitizeForSForTFinput(string SForTF) {
    if (settings_->enableDigitize()) {
      if (digitizedForSForTFinput_ != SForTF) {
        // Digitize variables specific to seed filter or track fittr if not already done.
        digitalStub_->makeSForTFinput(SForTF);

        // Must replace stub r coordinate, as seed filter & fitters work with digitized r instead of digitized rT.
        r_ = digitalStub_->r();
        // And KF may also redigitize z.
        z_ = digitalStub_->z();

        digitizedForSForTFinput_ = SForTF;
      }
    }
  }

  //=== Digitize stub for input to r-z Seed Filter.

  void Stub::digitizeForDRinput(unsigned int stubId) {
    if (settings_->enableDigitize()) {
      // Digitize variables specific to seed filter if not already done.
      digitalStub_->makeDRinput(stubId);
      // digitizedForDRinput_ = true;
    }
  }

  //===  Restore stub to pre-digitized state. i.e. Undo what function digitize() did.

  void Stub::reset_digitize() {
    if (settings_->enableDigitize()) {
      // Save CPU by not undoing digitization if stub was not already digitized.
      if (digitizedForGPinput_ || digitizedForHTinput_) {
        // Replace stub coordinates and bend with original coordinates stored prior to any digitization.
        phi_ = digitalStub_->orig_phi();
        r_ = digitalStub_->orig_r();
        z_ = digitalStub_->orig_z();
        bend_ = digitalStub_->orig_bend();

        // Note that stub is (no longer) digitized.
        digitizedForGPinput_ = false;
        digitizedForHTinput_ = false;
        digitizedForSForTFinput_ = "";

        // If the Stub class contains any data members that are not input to the GP or HT, but are derived from
        // variables that are, then be sure to update these here too.

        // Update variables giving ratio of track intercept angle to stub bend.
        this->calcDphiOverBend();
      }
    }
  }

  //=== Degrade assumed stub bend resolution.
  //=== Also return boolean indicating if stub bend was outside assumed window, so stub should be rejected
  //=== and return an integer indicating how many values of bend are merged into this single one.

  void Stub::degradeResolution(float bend, float& degradedBend, bool& reject, unsigned int& num) const {
    // If TMTT code is tightening official CMS FE stub window cuts, then calculate TMTT stub windows.
    float windowFE;
    if (settings_->killLowPtStubs()) {
      // Window size corresponding to Pt cut used for tracking.
      float invPtMax = 1. / (settings_->houghMinPt());
      windowFE = invPtMax / std::abs(this->qOverPtOverBend());
      // Increase half-indow size to allow for resolution in bend.
      windowFE += this->bendResInFrontend();
    } else {
      windowFE = 99999.;  // TMTT is not tightening windows.
    }

    static std::atomic<bool> firstErr = true;
    degradeBend_.degrade(bend, psModule(), moduleInfo_->detId(), windowFE, degradedBend, reject, num);
  }

  //=== Set flag indicating if stub will be output by front-end readout electronics
  //=== (where we can reconfigure the stub window size and rapidity cut).
  //=== Argument indicates if stub bend was outside window size encoded in DegradeBend.h
  //=== Note that this should run on quantities as available inside front-end chip, which are not
  //=== degraded by loss of bits or digitisation.

  void Stub::setFrontend(bool rejectStub, const StubKiller* stubKiller) {
    frontendPass_ = true;              // Did stub pass cuts applied in front-end chip
    stubFailedDegradeWindow_ = false;  // Did it only fail cuts corresponding to windows encoded in DegradeBend.h?
    // Don't use stubs at large eta, since it is impossible to form L1 tracks from them, so they only contribute to combinatorics.
    if (std::abs(this->eta()) > settings_->maxStubEta())
      frontendPass_ = false;
    // Don't use stubs whose Pt is significantly below the Pt cut used in the L1 tracking, allowing for uncertainty in q/Pt due to stub bend resolution.
    if (settings_->killLowPtStubs()) {
      const float qOverPtCut = 1. / settings_->houghMinPt();
      // Apply this cut in the front-end electronics.
      if (std::abs(this->bendInFrontend()) - this->bendResInFrontend() > qOverPtCut / this->qOverPtOverBend())
        frontendPass_ = false;
      // Reapply the same cut using the degraded bend information available in the off-detector electronics.
      // The reason is  that the bend degredation can move the Pt below the Pt cut, making the stub useless to the off-detector electronics.
      if (std::abs(this->bend()) - this->bendRes() > qOverPtCut / this->qOverPtOverBend())
        frontendPass_ = false;
    }
    // Don't use stubs whose bend is outside the window encoded into DegradeBend.h
    if (rejectStub) {
      if (frontendPass_)
        stubFailedDegradeWindow_ = true;
      frontendPass_ = false;
    }

    // Emulate stubs in dead tracker regions..
    StubKiller::KillOptions killScenario = static_cast<StubKiller::KillOptions>(settings_->killScenario()); 
    if (killScenario != StubKiller::KillOptions::none) {
      bool kill = stubKiller->killStub(ttStubRef_.get());
      if (kill)
        frontendPass_ = false;
    }
  }

  //=== Function to calculate approximation for dphiOverBendCorrection aka B
  double Stub::approxB() {
    if (tiltedBarrel()) {
      return settings_->bApprox_gradient() * std::abs(z_) / r_ + settings_->bApprox_intercept();
    } else {
      return barrel() ? 1 : std::abs(z_) / r_;
    }
  }

    //=== Calculate variables giving ratio of track intercept angle to stub bend.

  void Stub::calcDphiOverBend() {
   // Uses stub (r,z) instead of module (r,z). Logically correct but has negligable effect on results.
    dphiOverBendCorrection_ = std::abs(cos(this->theta() - moduleInfo_->moduleTilt()) / sin(this->theta()));
    dphiOverBendCorrection_approx_ = approxB();
    if (settings_->useApproxB()) {
      dphiOverBend_ = moduleInfo_->pitchOverSep() * dphiOverBendCorrection_approx_;
    } else {
      dphiOverBend_ = moduleInfo_->pitchOverSep() * dphiOverBendCorrection_;
    }

    } 

  //=== Note which tracking particle(s), if any, produced this stub.
  //=== The 1st argument is a map relating TrackingParticles to TP.

  void Stub::fillTruth(const map<edm::Ptr<TrackingParticle>, const TP*>& translateTP,
                       const edm::Handle<TTStubAssMap>& mcTruthTTStubHandle,
                       const edm::Handle<TTClusterAssMap>& mcTruthTTClusterHandle) {
    //--- Fill assocTP_ info. If both clusters in this stub were produced by the same single tracking particle, find out which one it was.

    bool genuine = mcTruthTTStubHandle->isGenuine(ttStubRef_);  // Same TP contributed to both clusters?
    assocTP_ = nullptr;

    // Require same TP contributed to both clusters.
    if (genuine) {
      edm::Ptr<TrackingParticle> tpPtr = mcTruthTTStubHandle->findTrackingParticlePtr(ttStubRef_);
      if (translateTP.find(tpPtr) != translateTP.end()) {
        assocTP_ = translateTP.at(tpPtr);
        // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
      }
    }

    // Fill assocTPs_ info.

    if (settings_->stubMatchStrict()) {
      // We consider only stubs in which this TP contributed to both clusters.
      if (assocTP_ != nullptr)
        assocTPs_.insert(assocTP_);

    } else {
      // We consider stubs in which this TP contributed to either cluster.

      for (unsigned int iClus = 0; iClus <= 1; iClus++) {  // Loop over both clusters that make up stub.
        const TTClusterRef& ttClusterRef = ttStubRef_->clusterRef(iClus);

        // Now identify all TP's contributing to either cluster in stub.
        vector<edm::Ptr<TrackingParticle> > vecTpPtr = mcTruthTTClusterHandle->findTrackingParticlePtrs(ttClusterRef);

        for (edm::Ptr<TrackingParticle> tpPtr : vecTpPtr) {
          if (translateTP.find(tpPtr) != translateTP.end()) {
            assocTPs_.insert(translateTP.at(tpPtr));
            // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
          }
        }
      }
    }

    //--- Also note which tracking particles produced the two clusters that make up the stub

    for (unsigned int iClus = 0; iClus <= 1; iClus++) {  // Loop over both clusters that make up stub.
      const TTClusterRef& ttClusterRef = ttStubRef_->clusterRef(iClus);

      bool genuineCluster = mcTruthTTClusterHandle->isGenuine(ttClusterRef);  // Only 1 TP made cluster?
      assocTPofCluster_[iClus] = nullptr;

      // Only consider clusters produced by just one TP.
      if (genuineCluster) {
        edm::Ptr<TrackingParticle> tpPtr = mcTruthTTClusterHandle->findTrackingParticlePtr(ttClusterRef);

        if (translateTP.find(tpPtr) != translateTP.end()) {
          assocTPofCluster_[iClus] = translateTP.at(tpPtr);
          // N.B. Since not all tracking particles are stored in InputData::vTPs_, sometimes no match will be found.
        }
      }
    }
  }

  //=== Estimated phi angle at which track intercepts a given radius rad, based on stub bend info. Also estimate uncertainty on this angle due to endcap 2S module strip length.
  //=== N.B. This is identical to Stub::beta() if rad=0.

  pair<float, float> Stub::trkPhiAtR(float rad) const {
    float rStubMax = r_ + sigmaR();  // Uncertainty in radial stub coordinate due to strip length.
    float rStubMin = r_ - sigmaR();
    float trkPhi1 = (phi_ + dphi() * (1. - rad / rStubMin));
    float trkPhi2 = (phi_ + dphi() * (1. - rad / rStubMax));
    float trkPhi = 0.5 * (trkPhi1 + trkPhi2);
    float errTrkPhi = 0.5 * std::abs(trkPhi1 - trkPhi2);
    return pair<float, float>(trkPhi, errTrkPhi);
  }

  //=== Note if stub is a crazy distance from the tracking particle trajectory that produced it.
  //=== If so, it was probably produced by a delta ray.

  bool Stub::crazyStub() const {
    bool crazy;
    if (assocTP_ == nullptr) {
      crazy = false;  // Stub is fake, but this is not crazy. It happens ...
    } else {
      // Stub was produced by TP. Check it lies not too far from TP trajectory.
      crazy = std::abs(reco::deltaPhi(phi_, assocTP_->trkPhiAtStub(this))) > settings_->crazyStubCut();
    }
    return crazy;
  }
}  // namespace tmtt
