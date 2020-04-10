#ifndef L1Trigger_TrackFindingTMTT_DigiProducer_h
#define L1Trigger_TrackFindingTMTT_DigiProducer_h

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
#include "DataFormats/Common/interface/DetSetVector.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/Histos.h"
#include "L1Trigger/TrackFindingTMTT/interface/TrackFitGeneric.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Demonstrator/DataFormats/interface/DigiKF4Track.hpp"
#include "Demonstrator/DataFormats/interface/DigiHTStub.hpp"
#include "Demonstrator/DataFormats/interface/DigiHTMiniStub.hpp"
#include "Demonstrator/DataFormats/interface/DigiDTCStub.hpp"

#include <vector>
#include <map>
#include <string>

namespace demo {

  /*class Settings;
class Histos;
class TrackFitGeneric;*/

  class DigiProducer : public edm::EDProducer {
  public:
    explicit DigiProducer(const edm::ParameterSet &);
    ~DigiProducer() {}

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
    const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
    const edm::EDGetTokenT<tmtt::TTStubDetSetVec> stubToken_;
    const edm::EDGetTokenT<tmtt::TTStubAssMap> stubTruthToken_;
    const edm::EDGetTokenT<tmtt::TTClusterAssMap> clusterTruthToken_;
    const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;

    const TrackerGeometry* trackerGeometry_; 
    const TrackerTopology* trackerTopology_;

    // Configuration parameters
    tmtt::Settings *settings_;
    std::vector<std::string> trackFitters_;
    std::vector<std::string> useRZfilter_;
    bool runRZfilter_;

    tmtt::Histos *hists_;
   std::map<std::string, tmtt::TrackFitGeneric *> fitterWorkerMap_;

    tmtt::TrackerGeometryInfo trackerGeometryInfo_;
  };

}  // namespace demo

#endif /* __DEMONSTRATOR_PRODUCER_DIGIPRODUCER_HPP__ */
