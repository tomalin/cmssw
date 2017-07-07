#ifndef EventFilter_Phase2TrackerRawToStub_Phase2TrackerStubProducer_H
#define EventFilter_Phase2TrackerRawToStub_Phase2TrackerStubProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>

namespace Phase2Tracker {

  class Phase2TrackerStubProducer : public edm::EDProducer
  {
  public:
    Phase2TrackerStubProducer( const edm::ParameterSet& pset );
    ~Phase2TrackerStubProducer();
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

  private:
    unsigned int runNumber_;
    edm::EDGetTokenT<FEDRawDataCollection> token_;
    const Phase2TrackerCabling * cabling_;
    std::map< int, std::pair<int,int> > stackMap_;
    uint32_t cacheId_;
    DetIdCollection detids_;
  };
}
#endif // EventFilter_Phase2TrackerRawToStub_Phase2TrackerStubProducer_H
