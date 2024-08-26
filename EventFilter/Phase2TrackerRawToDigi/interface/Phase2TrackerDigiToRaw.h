#ifndef  EventFilter_Phase2TrackerRawToDigi_DigiToRaw_H 
#define EventFilter_Phase2TrackerRawToDigi_DigiToRaw_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerStackedClus.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

namespace Phase2Tracker {

  class Phase2TrackerDigiToRaw {
  public:

    typedef edmNew::DetSetVector<Phase2TrackerCluster1D> DetSetVecClus;
    typedef edmNew::DetSet<Phase2TrackerCluster1D> DetSetClus;

    Phase2TrackerDigiToRaw() {}
    Phase2TrackerDigiToRaw(const Phase2TrackerCabling*,
                           const TrackerGeometry* tGeom,
                           const TrackerTopology* tTopo,
                           std::map<int, std::pair<int, int> > stackMap,
                           edm::Handle<DetSetVecClus>,
                           int);
    ~Phase2TrackerDigiToRaw() {}
    // loop on FEDs to create buffers
    void buildFEDBuffers(std::unique_ptr<FEDRawDataCollection>&);

  private:
    // builds a single FED buffer
    std::vector<uint64_t> makeBuffer(std::vector<DetSetVecClus::const_iterator>);
    // write FE Header to buffer
    void writeFeHeaderSparsified(std::vector<uint64_t>&, uint64_t&, int, int, int);
    // determine if a P or S cluster should be written
    void writeCluster(std::vector<uint64_t>& buffer, uint64_t&, const StackedClus&);
    // write S cluster to buffer
    void writeSCluster(std::vector<uint64_t>& buffer, uint64_t&, const StackedClus&, bool);
    // write P cluster to buffer
    void writePCluster(std::vector<uint64_t>& buffer, uint64_t&, const StackedClus&);

  private:
    // data you get from outside
    const Phase2TrackerCabling* cabling_;
    const TrackerTopology* tTopo_;
    const TrackerGeometry* tGeom_;
    std::map<int, std::pair<int, int> > stackMap_;
    edm::Handle<DetSetVecClus > clusshandle_;
    int mode_;
    // headers
    FEDDAQHeader FedDaqHeader_;
    FEDDAQTrailer FedDaqTrailer_;
    Phase2TrackerFEDHeader FedHeader_;
  };
}  // namespace Phase2Tracker

#endif
