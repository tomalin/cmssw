#ifndef L1Trigger_TrackFindingTMTT_TMTrackProducer_h
#define L1Trigger_TrackFindingTMTT_TMTrackProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/TrackerGeometryInfo.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <vector>
#include <map>
#include <string>

namespace tmtt {

  class Settings;
  class Histos;
  class TrackFitGeneric;

  class TMTrackProducer : public edm::EDProducer {
  public:
    explicit TMTrackProducer(const edm::ParameterSet &);
    ~TMTrackProducer() {}

  private:
    typedef std::vector<TTTrack<Ref_Phase2TrackerDigi_> > TTTrackCollection;

    virtual void beginRun(const edm::Run &, const edm::EventSetup &);
    virtual void produce(edm::Event &, const edm::EventSetup &);
    virtual void endJob();

  private:
    // ES tokens
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
    edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
    edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
    // ED tokens
    edm::EDGetTokenT<TTStubDetSetVec> stubToken_;
    edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
    edm::EDGetTokenT<TTStubAssMap> stubTruthToken_;
    edm::EDGetTokenT<TTClusterAssMap> clusterTruthToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;

    const TrackerGeometry *trackerGeometry_;
    const TrackerTopology *trackerTopology_;

    // Configuration parameters
    Settings *settings_;
    std::vector<std::string> trackFitters_;
    std::vector<std::string> useRZfilter_;
    bool runRZfilter_;

    Histos *hists_;
    std::map<std::string, TrackFitGeneric *> fitterWorkerMap_;

    TrackerGeometryInfo trackerGeometryInfo_;
  };

}  // namespace tmtt

#endif
