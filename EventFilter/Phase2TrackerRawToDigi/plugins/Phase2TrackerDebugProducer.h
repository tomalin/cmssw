#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerDebugProducer_H
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerDebugProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>

namespace Phase2Tracker {

  class Phase2TrackerDebugProducer : public edm::EDProducer
  {
  public:
    Phase2TrackerDebugProducer( const edm::ParameterSet& pset );
    ~Phase2TrackerDebugProducer();
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

  private:
    edm::EDGetTokenT<FEDRawDataCollection> token_;
    const Phase2TrackerCabling * cabling_;
    std::map< int, std::pair<int,int> > stackMap_;
    DetIdCollection detids_;
  };
}
#endif // EventFilter_Phase2TrackerRawToDigi_Phase2TrackerDebugProducer_H
