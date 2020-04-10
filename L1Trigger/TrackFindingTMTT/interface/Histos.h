#ifndef L1Trigger_TrackFindingTMTT_Histos_h
#define L1Trigger_TrackFindingTMTT_Histos_h

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/TrackerGeometryInfo.h"

#include "boost/numeric/ublas/matrix.hpp"
using boost::numeric::ublas::matrix;

#include <vector>
#include <map>
#include <string>

class TH1F;
class TH2F;
class TH2Poly;
class TF1;
class TProfile;
class TGraphAsymmErrors;
class TGraph;
class TEfficiency;

namespace tmtt {

  class InputData;
  class TP;
  class Sector;
  class HTrphi;
  class Get3Dtracks;
  class L1fittedTrack;
  class L1fittedTrk4and5;

  class Histos {
  public:
    // Store cfg parameters.
    Histos(const Settings* settings);

    virtual ~Histos() {}

    // Book & fill all histograms.
    virtual void book();
    virtual void fill(const InputData& inputData,
                      const matrix<Sector>& mSectors,
                      const matrix<HTrphi>& mHtPhis,
                      const matrix<Get3Dtracks> mGet3Dtrks,
                      const std::map<std::string, std::vector<L1fittedTrack>>& fittedTracks);

    // Print tracking performance summary & make tracking efficiency histograms.
    virtual void endJobAnalysis();

    // Determine "B" parameter, used in GP firmware to allow for tilted modules.
    virtual void trackerGeometryAnalysis(const TrackerGeometryInfo trackerGeometryInfo);

    // Did user request output histograms via the TFileService in their cfg?
    virtual bool available() const { return fs_.isAvailable(); }

    // Should histograms be produced?
    virtual bool enabled() const { return (settings_->enableHistos() && available()); }

  protected:
    // Book histograms for specific topics.
    virtual TFileDirectory bookInputData();
    virtual TFileDirectory bookEtaPhiSectors();
    virtual TFileDirectory bookRphiHT();
    virtual TFileDirectory bookRZfilters();
    virtual TFileDirectory bookStudyBusyEvents();
    virtual TFileDirectory bookTrackCands(const std::string& tName);
    virtual std::map<std::string, TFileDirectory> bookTrackFitting();

    // Fill histograms for specific topics.
    virtual void fillInputData(const InputData& inputData);
    virtual void fillEtaPhiSectors(const InputData& inputData, const matrix<Sector>& mSectors);
    virtual void fillRphiHT(const matrix<HTrphi>& mHtRphis);
    virtual void fillRZfilters(const matrix<Get3Dtracks>& mGet3Dtrks);
    virtual void fillStudyBusyEvents(const InputData& inputData,
                                     const matrix<Sector>& mSectors,
                                     const matrix<HTrphi>& mHtRphis,
                                     const matrix<Get3Dtracks>& mGet3Dtrks);
    virtual void fillTrackCands(const InputData& inputData, const std::vector<L1track3D>& tracks, std::string tName);
    virtual void fillTrackFitting(const InputData& inputData,
                                  const std::map<std::string, std::vector<L1fittedTrack>>& fittedTracks);

    // Produce plots of tracking efficiency after HZ or after r-z track filter (run at end of job)
    virtual TFileDirectory plotTrackEfficiency(const std::string& tName);
    // Produce plots of tracking efficiency after track fit (run at end of job).
    virtual TFileDirectory plotTrackEffAfterFit(const std::string& fitName);

    // For Hybrid tracking
    // Produce plots of tracklet seed finding efficiency before track reco
    virtual void plotTrackletSeedEfficiency(){};
    // Produce plots of hybrid duplicate removal efficiency after track reco
    virtual void plotHybridDupRemovalEfficiency(){};

    virtual void makeEfficiencyPlot(
        TFileDirectory& inputDir, TEfficiency* outputEfficiency, TH1F* pass, TH1F* all, TString name, TString title);

    // Print summary of track-finding performance after track pattern reco.
    virtual void printTrackPerformance(const std::string& tName);

    // Print summary of track-finding performance after helix fit for given track fitter.
    virtual void printFitTrackPerformance(const std::string& fitName);

    // For Hybrid tracking
    // Print summary of seed finding and extrapolation performance during track pattern reco.
    virtual void printTrackletSeedFindingPerformance(){};

    // Print summary of duplicate removal performance after track pattern reco.
    virtual void printHybridDupRemovalPerformance(){};

    // Understand why not all tracking particles were reconstructed.
    // Returns list of tracking particles that were not reconstructed and an integer indicating why.
    // Only considers TP used for algorithmic efficiency measurement.
    virtual std::map<const TP*, std::string> diagnoseTracking(const std::vector<TP>& allTPs,
                                                    const std::vector<L1track3D>& tracks,
                                                    bool withRZfilter) const;

  protected:
    // Configuration parameters.
    const Settings* settings_;
    unsigned int genMinStubLayers_;
    unsigned int numPhiSectors_;
    unsigned int numEtaRegions_;
    float houghMinPt_;
    unsigned int houghNbinsPt_;
    unsigned int houghNbinsPhi_;
    float chosenRofZ_;
    std::vector<std::string> trackFitters_;
    std::vector<std::string> useRZfilter_;
    bool ranRZfilter_;
    bool resPlotOpt_;

    edm::Service<TFileService> fs_;

    bool oldSumW2opt_;

    // Histograms of input data.
    TH1F* hisNumEvents_;
    TProfile* profNumStubs_;
    TH1F* hisStubsVsEta_;
    TH1F* hisStubsVsR_;
    TH2F* hisStubsVsRVsZ_;
    TH2F* hisStubsModuleVsRVsZ_;
    TH2F* hisStubsModuleTiltVsZ_;
    TH2F* hisStubsdPhiCorrectionVsZ_;
    TH2F* hisStubsVsRVsPhi_;
    TH2F* hisStubsModuleVsRVsPhi_;

    TH2F* hisStubsVsRVsZ_outerModuleAtSmallerR_;
    TH2F* hisStubsVsRVsPhi_outerModuleAtSmallerR_;

    TProfile* profNumTPs_;
    TH1F* hisNumStubsPerTP_;
    TH1F* hisNumPSStubsPerTP_;
    TH1F* hisNum2SStubsPerTP_;
    TH1F* hisNumLayersPerTP_;
    TH1F* hisNumPSLayersPerTP_;
    TH1F* hisNum2SLayersPerTP_;

    TH1F* hisNumLayersPerTP_lowPt_;
    TH1F* hisNumPSLayersPerTP_lowPt_;
    TH1F* hisNum2SLayersPerTP_lowPt_;

    TH1F* hisNumLayersPerTP_mediumPt_;
    TH1F* hisNumPSLayersPerTP_mediumPt_;
    TH1F* hisNum2SLayersPerTP_mediumPt_;

    TH1F* hisNumLayersPerTP_highPt_;
    TH1F* hisNumPSLayersPerTP_highPt_;
    TH1F* hisNum2SLayersPerTP_highPt_;

    TH1F* hisNumLayersPerTP_muons_;
    TH1F* hisNumPSLayersPerTP_muons_;
    TH1F* hisNum2SLayersPerTP_muons_;

    TH1F* hisNumLayersPerTP_electrons_;
    TH1F* hisNumPSLayersPerTP_electrons_;
    TH1F* hisNum2SLayersPerTP_electrons_;

    TH1F* hisNumLayersPerTP_pions_;
    TH1F* hisNumPSLayersPerTP_pions_;
    TH1F* hisNum2SLayersPerTP_pions_;

    TProfile* hisStubKillFE_;
    TProfile* hisStubIneffiVsInvPt_;
    TProfile* hisStubIneffiVsEta_;
    TProfile* hisStubKillDegradeBend_;
    TH1F* hisPtStub_;
    TH1F* hisPtResStub_;
    TH1F* hisBendFilterPower_;
    TH1F* hisDelPhiStub_;
    TH1F* hisDelPhiResStub_;
    TH1F* hisDelPhiResStub_tilted_;
    TH1F* hisDelPhiResStub_notTilted_;
    TH1F* hisBendStub_;
    TH1F* hisBendResStub_;
    TH1F* hisNumMergedBend_;
    TH2F* hisBendVsLayerOrRingPS_;
    TH2F* hisBendVsLayerOrRing2S_;
    TH2F* hisBendFEVsLayerOrRingPS_;
    TH2F* hisBendFEVsLayerOrRing2S_;
    TH1F* hisPhiStubVsPhiTP_;
    TH1F* hisPhiStubVsPhi0TP_;
    TH1F* hisPhi0StubVsPhi0TP_;
    TH1F* hisPhi0StubVsPhi0TPres_;
    TH1F* hisPhiStubVsPhi65TP_;
    TH1F* hisPhi65StubVsPhi65TP_;
    TH1F* hisPhi65StubVsPhi65TPres_;
    TH1F* hisPitchOverSep_;
    TH1F* hisRhoParameter_;
    TH2F* hisAlphaCheck_;
    TH1F* hisFracStubsSharingClus0_;
    TH1F* hisFracStubsSharingClus1_;

    // Histograms of B
    TH1F* hisStubB_;
    TH1F* hisStubBApproxDiff_tilted_;
    TGraph* graphBVsZoverR_;

    // Histograms checking that (eta,phi) sector definition is good.
    TH1F* hisFracStubsInSec_;
    TH1F* hisFracStubsInEtaSec_;
    TH1F* hisFracStubsInPhiSec_;
    TH1F* hisNumSecsPerStub_;
    TH1F* hisNumEtaSecsPerStub_;
    TH1F* hisNumPhiSecsPerStub_;
    TH1F* hisNumStubsPerSec_;
    TProfile* profNumStubsPerEtaSec_;
    TH2F* hisLayerIDvsEtaSec_;
    TH2F* hisLayerIDreducedvsEtaSec_;

    // Histograms checking filling of r-phi HT array.
    TH2Poly* hisArrayHT_;
    TF1* hisStubHT_;
    TH1F* hisIncStubsPerHT_;
    TH1F* hisExcStubsPerHT_;
    TH2F* hisNumStubsInCellVsEta_;
    TH1F* hisStubsOnRphiTracksPerHT_;
    TH1F* hisHTstubsPerTrack_;
    TH1F* hisHTmBin_;
    TH1F* hisHTcBin_;

    // Histograms about r-z track filters (or other filters applied after r-phi HT array).
    TH1F* hisNumZtrkSeedCombinations_;
    TH1F* hisNumSeedCombinations_;
    TH1F* hisNumGoodSeedCombinations_;
    TH1F* hisCorrelationZTrk_;

    // Histograms for studying freak, large events with too many stubs.
    TH1F* hisNumBusySecsInPerEvent_;
    TH1F* hisNumBusySecsOutPerEvent_;
    TProfile* profFracBusyInVsEtaReg_;
    TProfile* profFracBusyOutVsEtaReg_;
    TProfile* profFracStubsKilledVsEtaReg_;
    TProfile* profFracTracksKilledVsEtaReg_;
    TProfile* profFracTracksKilledVsInvPt_;
    TProfile* profFracTPKilledVsEta_;
    TProfile* profFracTPKilledVsInvPt_;
    TH1F* hisNumTPkilledBusySec_;
    std::map<std::string, TH1F*> hisNumInputStubs_;
    std::map<std::string, TH1F*> hisQoverPtInputStubs_;
    std::map<std::string, TH1F*> hisNumOutputStubs_;
    std::map<std::string, TH1F*> hisNumTracks_;
    std::map<std::string, TH1F*> hisNumStubsPerTrack_;
    std::map<std::string, TH1F*> hisTrackQoverPt_;
    std::map<std::string, TH1F*> hisTrackPurity_;
    std::map<std::string, TH1F*> hisNumTPphysics_;
    std::map<std::string, TH1F*> hisNumTPpileup_;
    std::map<std::string, TH1F*> hisSumPtTPphysics_;
    std::map<std::string, TH1F*> hisSumPtTPpileup_;

    // Histograms studying 3D track candidates found by Hough Transform or r-z Track Filter.
    std::map<std::string, TProfile*> profNumTrackCands_;
    std::map<std::string, TProfile*> profNumTracksVsEta_;
    std::map<std::string, TH1F*> hisNumTracksVsQoverPt_;
    std::map<std::string, TH1F*> hisNumTrksPerSect_;
    std::map<std::string, TH1F*> hisNumTrksPerNon_;
    std::map<std::string, TProfile*> profStubsOnTracks_;
    std::map<std::string, TProfile*> profStubsOnTracksVsEta_;
    std::map<std::string, TH1F*> hisStubsOnTracksPerSect_;
    std::map<std::string, TH1F*> hisStubsOnTracksPerNon_;
    std::map<std::string, TH1F*> hisUniqueStubsOnTrksPerSect_;
    std::map<std::string, TH1F*> hisUniqueStubsOnTrksPerNon_;
    std::map<std::string, TH1F*> hisStubsPerTrack_;
    std::map<std::string, TH1F*> hisLayersPerTrack_;
    std::map<std::string, TH1F*> hisPSLayersPerTrack_;
    std::map<std::string, TH1F*> hisLayersPerTrueTrack_;
    std::map<std::string, TH1F*> hisPSLayersPerTrueTrack_;

    std::map<std::string, TH1F*> hisNumStubsPerLink_;
    std::map<std::string, TH2F*> hisNumStubsVsLink_;
    std::map<std::string, TProfile*> profMeanStubsPerLink_;
    std::map<std::string, TH1F*> hisNumTrksPerLink_;
    std::map<std::string, TH2F*> hisNumTrksVsLink_;
    std::map<std::string, TProfile*> profMeanTrksPerLink_;

    std::map<std::string, TProfile*> profExcessStubsPerTrackVsPt_;
    std::map<std::string, TH1F*> hisFracMatchStubsOnTracks_;
    std::map<std::string, TH1F*> hisDeltaBendTrue_;
    std::map<std::string, TH1F*> hisDeltaBendFake_;
    std::map<std::string, TProfile*> profFracTrueStubsVsLayer_;
    std::map<std::string, TProfile*> profDupTracksVsEta_;
    std::map<std::string, TProfile*> profDupTracksVsInvPt_;
    //map<std::string, TH2F*> hisWrongSignStubRZ_pBend_;
    //map<std::string, TH2F*> hisWrongSignStubRZ_nBend_;

    // Histos of track params after HT.
    std::map<std::string, TH1F*> hisQoverPt_;
    std::map<std::string, TH1F*> hisPhi0_;
    std::map<std::string, TH1F*> hisEta_;
    std::map<std::string, TH1F*> hisZ0_;

    // Histograms of track parameter resolution after HT transform.
    std::map<std::string, TH1F*> hisQoverPtRes_;
    std::map<std::string, TH1F*> hisPhi0Res_;
    std::map<std::string, TH1F*> hisEtaRes_;
    std::map<std::string, TH1F*> hisZ0Res_;

    std::map<std::string, TH2F*> hisRecoVsTrueQinvPt_;
    std::map<std::string, TH2F*> hisRecoVsTruePhi0_;
    std::map<std::string, TH2F*> hisRecoVsTrueD0_;
    std::map<std::string, TH2F*> hisRecoVsTrueZ0_;
    std::map<std::string, TH2F*> hisRecoVsTrueEta_;

    // Diagnosis of failed tracking.
    std::map<std::string, TH1F*> hisRecoFailureReason_;
    std::map<std::string, TH1F*> hisRecoFailureLayer_;

    std::map<std::string, TH1F*> hisNumStubsOnLayer_;

    // Histos used for denominator of tracking efficiency plots.
    TH1F* hisTPinvptForEff_;
    TH1F* hisTPptForEff_;
    TH1F* hisTPetaForEff_;
    TH1F* hisTPphiForEff_;
    TH1F* hisTPd0ForEff_;
    TH1F* hisTPz0ForEff_;
    //
    TH1F* hisTPinvptForAlgEff_;
    TH1F* hisTPptForAlgEff_;
    TH1F* hisTPetaForAlgEff_;
    TH1F* hisTPphiForAlgEff_;
    TH1F* hisTPd0ForAlgEff_;
    TH1F* hisTPz0ForAlgEff_;
    //
    TH1F* hisTPphisecForAlgEff_;
    TH1F* hisTPetasecForAlgEff_;
    TH1F* hisTPinvptForAlgEff_inJetPtG30_;
    TH1F* hisTPinvptForAlgEff_inJetPtG100_;
    TH1F* hisTPinvptForAlgEff_inJetPtG200_;

    // Histograms used to make efficiency plots with 3D track candidates prior to fit.
    std::map<std::string, TH1F*> hisRecoTPinvptForEff_;
    std::map<std::string, TH1F*> hisRecoTPptForEff_;
    std::map<std::string, TH1F*> hisRecoTPetaForEff_;
    std::map<std::string, TH1F*> hisRecoTPphiForEff_;
    //
    std::map<std::string, TH1F*> hisPerfRecoTPinvptForEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPptForEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPetaForEff_;
    //
    std::map<std::string, TH1F*> hisRecoTPd0ForEff_;
    std::map<std::string, TH1F*> hisRecoTPz0ForEff_;
    //
    std::map<std::string, TH1F*> hisRecoTPinvptForAlgEff_;
    std::map<std::string, TH1F*> hisRecoTPptForAlgEff_;
    std::map<std::string, TH1F*> hisRecoTPetaForAlgEff_;
    std::map<std::string, TH1F*> hisRecoTPphiForAlgEff_;
    //
    std::map<std::string, TH1F*> hisPerfRecoTPinvptForAlgEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPptForAlgEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPetaForAlgEff_;
    //
    std::map<std::string, TH1F*> hisRecoTPd0ForAlgEff_;
    std::map<std::string, TH1F*> hisRecoTPz0ForAlgEff_;
    //
    std::map<std::string, TH1F*> hisRecoTPphisecForAlgEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPphisecForAlgEff_;
    std::map<std::string, TH1F*> hisRecoTPetasecForAlgEff_;
    std::map<std::string, TH1F*> hisPerfRecoTPetasecForAlgEff_;

    std::map<std::string, TH1F*> hisRecoTPinvptForAlgEff_inJetPtG30_;
    std::map<std::string, TH1F*> hisRecoTPinvptForAlgEff_inJetPtG100_;
    std::map<std::string, TH1F*> hisRecoTPinvptForAlgEff_inJetPtG200_;

    // Histograms for track fitting evaluation, where std::map index specifies name of track fitting algorithm used.

    std::map<std::string, TProfile*> profNumFitTracks_;
    std::map<std::string, TH1F*> hisNumFitTrks_;
    std::map<std::string, TH1F*> hisNumFitTrksPerNon_;
    std::map<std::string, TH1F*> hisNumFitTrksPerSect_;

    std::map<std::string, TH1F*> hisStubsPerFitTrack_;
    std::map<std::string, TProfile*> profStubsOnFitTracks_;

    std::map<std::string, TH1F*> hisFitQinvPtMatched_;
    std::map<std::string, TH1F*> hisFitPhi0Matched_;
    std::map<std::string, TH1F*> hisFitD0Matched_;
    std::map<std::string, TH1F*> hisFitZ0Matched_;
    std::map<std::string, TH1F*> hisFitEtaMatched_;

    std::map<std::string, TH1F*> hisFitQinvPtUnmatched_;
    std::map<std::string, TH1F*> hisFitPhi0Unmatched_;
    std::map<std::string, TH1F*> hisFitD0Unmatched_;
    std::map<std::string, TH1F*> hisFitZ0Unmatched_;
    std::map<std::string, TH1F*> hisFitEtaUnmatched_;

    std::map<std::string, TH1F*> hisKalmanNumUpdateCalls_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay0Matched_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay1Matched_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay2Matched_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay0Unmatched_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay1Unmatched_;
    std::map<std::string, TH1F*> hisKalmanChi2DofSkipLay2Unmatched_;

    std::map<std::string, TH1F*> hisFitChi2Matched_;
    std::map<std::string, TH1F*> hisFitChi2DofMatched_;
    std::map<std::string, TH1F*> hisFitChi2DofRphiMatched_;
    std::map<std::string, TH1F*> hisFitChi2DofRzMatched_;
    std::map<std::string, TH1F*> hisFitBeamChi2Matched_;
    std::map<std::string, TH1F*> hisFitBeamChi2DofMatched_;
    std::map<std::string, TProfile*> profFitChi2VsEtaMatched_;
    std::map<std::string, TProfile*> profFitChi2DofVsEtaMatched_;
    std::map<std::string, TProfile*> profFitChi2VsInvPtMatched_;
    std::map<std::string, TProfile*> profFitChi2DofVsInvPtMatched_;
    std::map<std::string, TProfile*> profFitChi2VsTrueD0Matched_;
    std::map<std::string, TProfile*> profFitChi2DofVsTrueD0Matched_;
    std::map<std::string, TH1F*> hisFitChi2PerfMatched_;
    std::map<std::string, TH1F*> hisFitChi2DofPerfMatched_;

    std::map<std::string, TH1F*> hisFitChi2Unmatched_;
    std::map<std::string, TH1F*> hisFitChi2DofUnmatched_;
    std::map<std::string, TH1F*> hisFitChi2DofRphiUnmatched_;
    std::map<std::string, TH1F*> hisFitChi2DofRzUnmatched_;
    std::map<std::string, TH1F*> hisFitBeamChi2Unmatched_;
    std::map<std::string, TH1F*> hisFitBeamChi2DofUnmatched_;
    std::map<std::string, TProfile*> profFitChi2VsEtaUnmatched_;
    std::map<std::string, TProfile*> profFitChi2DofVsEtaUnmatched_;
    std::map<std::string, TProfile*> profFitChi2VsInvPtUnmatched_;
    std::map<std::string, TProfile*> profFitChi2DofVsInvPtUnmatched_;

    std::map<std::string, TProfile*> profFitChi2VsPurity_;
    std::map<std::string, TProfile*> profFitChi2DofVsPurity_;

    std::map<std::string, TH1F*> hisDeltaPhitruePSbarrel_;
    std::map<std::string, TH1F*> hisDeltaRorZtruePSbarrel_;
    std::map<std::string, TH1F*> hisDeltaPhitrue2Sbarrel_;
    std::map<std::string, TH1F*> hisDeltaRorZtrue2Sbarrel_;
    std::map<std::string, TH1F*> hisDeltaPhitruePSendcap_;
    std::map<std::string, TH1F*> hisDeltaRorZtruePSendcap_;
    std::map<std::string, TH1F*> hisDeltaPhitrue2Sendcap_;
    std::map<std::string, TH1F*> hisDeltaRorZtrue2Sendcap_;
    std::map<std::string, TH1F*> hisDeltaPhifakePSbarrel_;
    std::map<std::string, TH1F*> hisDeltaRorZfakePSbarrel_;
    std::map<std::string, TH1F*> hisDeltaPhifake2Sbarrel_;
    std::map<std::string, TH1F*> hisDeltaRorZfake2Sbarrel_;
    std::map<std::string, TH1F*> hisDeltaPhifakePSendcap_;
    std::map<std::string, TH1F*> hisDeltaRorZfakePSendcap_;
    std::map<std::string, TH1F*> hisDeltaPhifake2Sendcap_;
    std::map<std::string, TH1F*> hisDeltaRorZfake2Sendcap_;
    std::map<std::string, TProfile*> profRecalcRphiChi2VsEtaTrue1_;
    std::map<std::string, TProfile*> profRecalcRzChi2VsEtaTrue1_;
    std::map<std::string, TProfile*> profRecalcChi2VsEtaTrue1_;
    std::map<std::string, TProfile*> profRecalcChi2VsEtaTrue2_;
    std::map<std::string, TProfile*> profNsigmaPhivsInvPt_;
    std::map<std::string, TProfile*> profNsigmaPhivsR_;
    std::map<std::string, TProfile*> profNsigmaPhivsTanl_;

    std::map<std::string, TH2F*> hisFitVsSeedQinvPtMatched_;
    std::map<std::string, TH2F*> hisFitVsSeedPhi0Matched_;
    std::map<std::string, TH2F*> hisFitVsSeedD0Matched_;
    std::map<std::string, TH2F*> hisFitVsSeedZ0Matched_;
    std::map<std::string, TH2F*> hisFitVsSeedEtaMatched_;

    std::map<std::string, TH2F*> hisFitVsSeedQinvPtUnmatched_;
    std::map<std::string, TH2F*> hisFitVsSeedPhi0Unmatched_;
    std::map<std::string, TH2F*> hisFitVsSeedD0Unmatched_;
    std::map<std::string, TH2F*> hisFitVsSeedZ0Unmatched_;
    std::map<std::string, TH2F*> hisFitVsSeedEtaUnmatched_;

    std::map<std::string, TH2F*> hisNumStubsVsPurityMatched_;
    std::map<std::string, TProfile*> profFitFracTrueStubsVsLayerMatched_;
    std::map<std::string, TProfile*> profFitFracTrueStubsVsEtaMatched_;

    std::map<std::string, TH2F*> hisFitVsTrueQinvPt_;
    std::map<std::string, TH2F*> hisFitVsTruePhi0_;
    std::map<std::string, TH2F*> hisFitVsTrueD0_;
    std::map<std::string, TH2F*> hisFitVsTrueZ0_;
    std::map<std::string, TH2F*> hisFitVsTrueEta_;

    std::map<std::string, TH1F*> hisFitQinvPtRes_;
    std::map<std::string, TH1F*> hisFitPhi0Res_;
    std::map<std::string, TH1F*> hisFitD0Res_;
    std::map<std::string, TH1F*> hisFitZ0Res_;
    std::map<std::string, TH1F*> hisFitEtaRes_;

    std::map<std::string, TProfile*> hisQoverPtResVsTrueEta_;
    std::map<std::string, TProfile*> hisPhi0ResVsTrueEta_;
    std::map<std::string, TProfile*> hisEtaResVsTrueEta_;
    std::map<std::string, TProfile*> hisZ0ResVsTrueEta_;
    std::map<std::string, TProfile*> hisD0ResVsTrueEta_;

    std::map<std::string, TProfile*> hisQoverPtResVsTrueInvPt_;
    std::map<std::string, TProfile*> hisPhi0ResVsTrueInvPt_;
    std::map<std::string, TProfile*> hisEtaResVsTrueInvPt_;
    std::map<std::string, TProfile*> hisZ0ResVsTrueInvPt_;
    std::map<std::string, TProfile*> hisD0ResVsTrueInvPt_;

    std::map<std::string, TProfile*> hisQoverPtResBeamVsTrueEta_;
    std::map<std::string, TProfile*> hisPhi0ResBeamVsTrueEta_;
    std::map<std::string, TProfile*> hisQoverPtResBeamVsTrueInvPt_;
    std::map<std::string, TProfile*> hisPhi0ResBeamVsTrueInvPt_;

    std::map<std::string, TH2F*> hisFitEfficiencyVsChi2Dof_;
    std::map<std::string, TH2F*> hisNumStubsVsChi2Dof_;
    std::map<std::string, TH2F*> hisNumLayersVsChi2Dof_;
    std::map<std::string, TH2F*> hisAvgNumStubsPerLayerVsChi2Dof_;

    std::map<std::string, TProfile*> profDupFitTrksVsEta_;
    std::map<std::string, TProfile*> profDupFitTrksVsInvPt_;

    // Histograms used for efficiency plots made with fitted tracks.
    std::map<std::string, TH1F*> hisFitTPinvptForEff_;
    std::map<std::string, TH1F*> hisFitTPptForEff_;
    std::map<std::string, TH1F*> hisFitTPetaForEff_;
    std::map<std::string, TH1F*> hisFitTPphiForEff_;
    std::map<std::string, TH1F*> hisPerfFitTPinvptForEff_;
    std::map<std::string, TH1F*> hisPerfFitTPptForEff_;
    std::map<std::string, TH1F*> hisPerfFitTPetaForEff_;
    std::map<std::string, TH1F*> hisFitTPd0ForEff_;
    std::map<std::string, TH1F*> hisFitTPz0ForEff_;
    std::map<std::string, TH1F*> hisFitTPinvptForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPptForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPetaForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPphiForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPptForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPetaForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPd0ForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPz0ForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPphisecForAlgEff_;
    std::map<std::string, TH1F*> hisFitTPetasecForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPphisecForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPetasecForAlgEff_;
    std::map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_inJetPtG30_;
    std::map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_inJetPtG100_;
    std::map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_inJetPtG200_;

    // Histograms of tracking efficiency & fake rate after Hough transform or after r-z track filter.
    std::map<std::string, TEfficiency*> teffEffVsInvPt_;
    std::map<std::string, TEfficiency*> teffEffVsPt_;
    std::map<std::string, TEfficiency*> teffEffVsEta_;
    std::map<std::string, TEfficiency*> teffEffVsPhi_;
    //
    std::map<std::string, TEfficiency*> teffPerfEffVsInvPt_;
    std::map<std::string, TEfficiency*> teffPerfEffVsPt_;
    std::map<std::string, TEfficiency*> teffPerfEffVsEta_;
    //
    std::map<std::string, TEfficiency*> teffEffVsD0_;
    std::map<std::string, TEfficiency*> teffEffVsZ0_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffVsInvPt_;
    std::map<std::string, TEfficiency*> teffAlgEffVsPt_;
    std::map<std::string, TEfficiency*> teffAlgEffVsEta_;
    std::map<std::string, TEfficiency*> teffAlgEffVsPhi_;
    std::map<std::string, TEfficiency*> teffAlgEffVsInvPt_inJetPtG30_;
    std::map<std::string, TEfficiency*> teffAlgEffVsInvPt_inJetPtG100_;
    std::map<std::string, TEfficiency*> teffAlgEffVsInvPt_inJetPtG200_;
    //
    std::map<std::string, TEfficiency*> teffPerfAlgEffVsInvPt_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffVsPt_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffVsEta_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffVsD0_;
    std::map<std::string, TEfficiency*> teffAlgEffVsZ0_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffVsPhiSec_;
    std::map<std::string, TEfficiency*> teffAlgEffVsEtaSec_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffVsPhiSec_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffVsEtaSec_;

    // Histograms of tracking efficiency & fake rate after Hough transform based on tracks after the track fit.
    std::map<std::string, TEfficiency*> teffEffFitVsInvPt_;
    std::map<std::string, TEfficiency*> teffEffFitVsPt_;
    std::map<std::string, TEfficiency*> teffEffFitVsEta_;
    std::map<std::string, TEfficiency*> teffEffFitVsPhi_;
    //
    std::map<std::string, TEfficiency*> teffPerfEffFitVsInvPt_;
    std::map<std::string, TEfficiency*> teffPerfEffFitVsPt_;
    std::map<std::string, TEfficiency*> teffPerfEffFitVsEta_;
    //
    std::map<std::string, TEfficiency*> teffEffFitVsD0_;
    std::map<std::string, TEfficiency*> teffEffFitVsZ0_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffFitVsInvPt_;
    std::map<std::string, TEfficiency*> teffAlgEffFitVsPt_;
    std::map<std::string, TEfficiency*> teffAlgEffFitVsEta_;
    std::map<std::string, TEfficiency*> teffAlgEffFitVsPhi_;
    //
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsInvPt_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsPt_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsEta_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsInvPt_inJetPtG30_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsInvPt_inJetPtG100_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsInvPt_inJetPtG200_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffFitVsD0_;
    std::map<std::string, TEfficiency*> teffAlgEffFitVsZ0_;
    //
    std::map<std::string, TEfficiency*> teffAlgEffFitVsPhiSec_;
    std::map<std::string, TEfficiency*> teffAlgEffFitVsEtaSec_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsPhiSec_;
    std::map<std::string, TEfficiency*> teffPerfAlgEffFitVsEtaSec_;

    bool plotFirst_;

    // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
    std::map<std::string, unsigned int> numFitAlgEff_;
    std::map<std::string, unsigned int> numFitPerfAlgEff_;

    // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
    std::map<std::string, unsigned int> numFitAlgEffPass_;
    std::map<std::string, unsigned int> numFitPerfAlgEffPass_;

    // Range in r of each barrel layer.
    std::map<unsigned int, float> mapBarrelLayerMinR_;
    std::map<unsigned int, float> mapBarrelLayerMaxR_;
    // Range in z of each endcap wheel.
    std::map<unsigned int, float> mapEndcapWheelMinZ_;
    std::map<unsigned int, float> mapEndcapWheelMaxZ_;

    // Range in (r,z) of each module type.
    std::map<unsigned int, float> mapModuleTypeMinR_;
    std::map<unsigned int, float> mapModuleTypeMaxR_;
    std::map<unsigned int, float> mapModuleTypeMinZ_;
    std::map<unsigned int, float> mapModuleTypeMaxZ_;
    // Extra std::maps for wierd barrel layers 1-2 & endcap wheels 3-5.
    std::map<unsigned int, float> mapExtraAModuleTypeMinR_;
    std::map<unsigned int, float> mapExtraAModuleTypeMaxR_;
    std::map<unsigned int, float> mapExtraAModuleTypeMinZ_;
    std::map<unsigned int, float> mapExtraAModuleTypeMaxZ_;
    std::map<unsigned int, float> mapExtraBModuleTypeMinR_;
    std::map<unsigned int, float> mapExtraBModuleTypeMaxR_;
    std::map<unsigned int, float> mapExtraBModuleTypeMinZ_;
    std::map<unsigned int, float> mapExtraBModuleTypeMaxZ_;
    std::map<unsigned int, float> mapExtraCModuleTypeMinR_;
    std::map<unsigned int, float> mapExtraCModuleTypeMaxR_;
    std::map<unsigned int, float> mapExtraCModuleTypeMinZ_;
    std::map<unsigned int, float> mapExtraCModuleTypeMaxZ_;
    std::map<unsigned int, float> mapExtraDModuleTypeMinR_;
    std::map<unsigned int, float> mapExtraDModuleTypeMaxR_;
    std::map<unsigned int, float> mapExtraDModuleTypeMinZ_;
    std::map<unsigned int, float> mapExtraDModuleTypeMaxZ_;

    bool bApproxMistake_;
  };

}  // namespace tmtt
#endif
