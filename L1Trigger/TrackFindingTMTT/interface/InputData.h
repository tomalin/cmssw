#ifndef L1Trigger_TrackFindingTMTT_InputData_h
#define L1Trigger_TrackFindingTMTT_InputData_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <vector>

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
              const edm::EDGetTokenT<TrackingParticleCollection> tpToken,
              const edm::EDGetTokenT<TTStubDetSetVec> stubToken,
              const edm::EDGetTokenT<TTStubAssMap> stubTruthToken,
              const edm::EDGetTokenT<TTClusterAssMap> clusterTruthToken,
              const edm::EDGetTokenT<reco::GenJetCollection> genJetToken);

    // Get tracking particles
    const std::vector<TP>& getTPs() const { return vTPs_; }
    // Get stubs that would be output by the front-end readout electronics
    const std::vector<const Stub*>& getStubs() const { return vStubs_; }

    //--- of minor importance ...

    // Get number of stubs prior to applying tighted front-end readout electronics cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
    const std::vector<Stub>& getAllStubs() const { return vAllStubs_; }

  private:
    // const edm::EDGetTokenT<TrackingParticleCollection> inputTag;

    // Can optionally be used to sort stubs by bend.
    struct SortStubsInBend {
      inline bool operator()(const Stub* stub1, const Stub* stub2) {
        return (fabs(stub1->bend()) < fabs(stub2->bend()));
      }
    };

  private:
    bool enableMCtruth_;  // Notes if job will use MC truth info.

    std::vector<TP> vTPs_;             // tracking particles
    std::vector<const Stub*> vStubs_;  // stubs that would be output by the front-end readout electronics.

    //--- of minor importance ...

    std::vector<Stub>
        vAllStubs_;  // all stubs, even those that would fail any tightened front-end readout electronic cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
  };

}  // namespace tmtt
#endif
