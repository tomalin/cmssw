#include "EventFilter/Phase2TrackerRawToDigi/plugins/Phase2TrackerDigiProducer.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/src/fed_header.h"
#include "DataFormats/FEDRawData/src/fed_trailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDRawChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDZSChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ext/algorithm>

using namespace std;

namespace Phase2Tracker {

  Phase2TrackerStubProducer::Phase2TrackerStubProducer( const edm::ParameterSet& pset ) :
    runNumber_(0),
    cabling_(0),
    cacheId_(0)
  {
    // define product
    produces< edm::DetSetVector<TTStub> >("Unsparsified");
    token_ = consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"));
  }
  
  Phase2TrackerStubProducer::~Phase2TrackerStubProducer()
  {
  }
  
  void Phase2TrackerStubProducer::beginJob( )
  {
  }
  
  void Phase2TrackerStubProducer::beginRun( edm::Run const& run, edm::EventSetup const& es)
  {
    // fetch cabling from event setup
    edm::ESHandle<Phase2TrackerCabling> c;
    es.get<Phase2TrackerCablingRcd>().get( c );
    cabling_ = c.product();

    // FIXME: build map of stacks to compensate for missing trackertopology methods
    edm::ESHandle<TrackerTopology> tTopoHandle;
    es.get<TrackerTopologyRcd>().get(tTopoHandle);
    const TrackerTopology* tTopo = tTopoHandle.product();

    edm::ESHandle< TrackerGeometry > tGeomHandle;
    es.get< TrackerDigiGeometryRecord >().get( tGeomHandle );
    const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

    for (auto iu = theTrackerGeom->detUnits().begin(); iu != theTrackerGeom->detUnits().end(); ++iu) {
      unsigned int detId_raw = (*iu)->geographicalId().rawId();
      DetId detId = DetId(detId_raw);
      if (detId.det() == DetId::Detector::Tracker) {
          // build map of upper and lower for each module
          if ( tTopo->isLower(detId) != 0 ) {
              stackMap_[tTopo->stack(detId)].first = detId;
          }
          if ( tTopo->isUpper(detId) != 0 ) {
              stackMap_[tTopo->stack(detId)].second = detId;
          }
      }
    } // end loop on detunits
  }
  
  void Phase2TrackerStubProducer::endJob()
  {
  }
  
  void Phase2TrackerStubProducer::produce( edm::Event& event, const edm::EventSetup& es)
  {
    std::unique_ptr<edmNew::DetSetVector<TTStub>> stubs( new edmNew::DetSetVector<TTStub>() ); 

    // Retrieve FEDRawData collection
    edm::Handle<FEDRawDataCollection> buffers;
    event.getByToken( token_, buffers );

    // Analyze strip tracker FED buffers in data
    std::vector<int> feds = cabling_->listFeds();
    std::vector<int>::iterator fedIndex;
    for(fedIndex = feds.begin(); fedIndex != feds.end(); ++fedIndex)
    {
      const FEDRawData& fed = buffers->FEDData(*fedIndex);
      if(fed.size()==0) continue;
	  // construct buffer
	  Phase2Tracker::Phase2TrackerFEDBuffer buffer(fed.data(),fed.size());
      // Skip FED if buffer is not a valid tracker FEDBuffer
      if(buffer.isValid() == 0) 
      { 
        LogTrace("Phase2TrackerStubProducer") << "[Phase2Tracker::Phase2TrackerStubProducer::"<<__func__<<"]: \n";
        LogTrace("Phase2TrackerStubProducer") << "Skipping invalid buffer for FED nr " << *fedIndex << endl;
        continue; 
      }
      // loop on channels
      // TODO : check loops consistency
      // TODO : use proper TTStub constructor
      int ichan = 0;
      for ( int ife = 0; ife < MAX_FE_PER_FED; ife++ )
      {
        // get fedid from cabling
        const Phase2TrackerModule mod = cabling_->findFedCh(std::make_pair(*fedIndex, ife));
        uint32_t detid = mod.getDetid();
        // container for this module's digis
        std::vector<TTStub> stubsTop;
        std::vector<TTStub> stubsBottom;
        const Phase2TrackerFEDChannel& channel = buffer.channel(ichan);
        if(channel.length() > 0)
        {
          // create appropriate unpacker
          if (channel.dettype() == DET_Son2S) 
          {
            Phase2TrackerFEDStubSon2SChannelUnpacker unpacker = Phase2TrackerFEDStubSon2SChannelUnpacker(channel);
            while (unpacker.hasData())
            {
              if (unpacker.rawX()%2) 
              {
                stubsTop.push_back(TTStub());
              }
              else 
              {
                stubsBottom.push_back(TTStub());
              }
              unpacker++;
            }
          } 
          else if (channel.dettype() == DET_SonPS)
          {
            Phase2TrackerFEDStubSonPSChannelUnpacker unpacker = Phase2TrackerFEDStubSonPSChannelUnpacker(channel);
            while (unpacker.hasData())
            {
              stubsTop.push_back(TTStub());
              unpacker++;
            }
          }
          else if (channel.dettype() == DET_PonPS)
          {
            Phase2TrackerFEDStubPonPSChannelUnpacker unpacker = Phase2TrackerFEDStubPonPSChannelUnpacker(channel);
            while (unpacker.hasData())
            {
              stubsBottom.push_back(TTStub());
              unpacker++;
            }
          }
        } // end reading stubs channel
        ichan++;
        if(detid > 0)
        {
          std::vector<TTStub>::iterator it;
          {
            // outer detid is defined as inner detid + 1 or module detid + 2
            edmNew::DetSetVector<TTStub>::FastFiller spct(*stubs, stackMap_[detid].second);
            for(it=stubsTop.begin();it!=stubsTop.end();it++)
            {
              spct.push_back(*it);
            }
          }
          {
            edmNew::DetSetVector<TTStub>::FastFiller spcb(*stubs, stackMap_[detid].first);
            for(it=stubsBottom.begin();it!=stubsBottom.end();it++)
            {
              spcb.push_back(*it);
            }
          }
        } // end if detid > 0 (?)
      } // end loop on FE
    } // en loop on FED   
    // store stubs
    event.put(std::move(stubs), "Stubs" );
  } 
}
