#ifndef L1Trigger_TrackFindingTMTT_InputData_h
#define L1Trigger_TrackFindingTMTT_InputData_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
#include "L1Trigger/TrackFindingTMTT/interface/TrackerModule.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <list>

namespace tmtt {

  class Settings;

  //=== Unpacks stub & tracking particle (truth) data into user-friendlier format in Stub & TP classes.
  //=== Also makes B-field available to Settings class.

  class InputData {
  public:
    InputData(const edm::Event& iEvent,
              const edm::EventSetup& iSetup,
              Settings* settings,
              const TrackerGeometry* trackerGeometry,
              const TrackerTopology* trackerTopology,
	      const std::list<TrackerModule>& listTrackerModule,
              const edm::EDGetTokenT<TrackingParticleCollection> tpToken,
              const edm::EDGetTokenT<TTStubDetSetVec> stubToken,
              const edm::EDGetTokenT<TTStubAssMap> stubTruthToken,
              const edm::EDGetTokenT<TTClusterAssMap> clusterTruthToken,
              const edm::EDGetTokenT<reco::GenJetCollection> genJetToken);

    // Info about each tracker module
    const std::list<TrackerModule>& trackerModules() const { return trackerModules_; };

    // Get tracking particles
    const std::list<TP>& getTPs() const { return vTPs_; }
    // Get stubs that would be output by the front-end readout electronics
    const std::list<const Stub*>& stubs() const { return vStubs_; }

    //--- of minor importance ...

    // Get number of stubs prior to applying tighted front-end readout electronics cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
    const std::list<Stub>& allStubs() const { return vAllStubs_; }

  private:
    bool enableMCtruth_;  // Notes if job will use MC truth info.

    std::list<TrackerModule> trackerModules_; // Info about each tracker module.

    std::list<TP> vTPs_;             // tracking particles
    std::list<const Stub*> vStubs_;  // stubs that would be output by the front-end readout electronics.

    //--- Used for a few minor studies ...

    // all stubs, even those that would fail any tightened front-end readout electronic cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
    std::list<Stub> vAllStubs_;
  };

}  // namespace tmtt
#endif
