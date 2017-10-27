#include "EventFilter/Phase2TrackerRawToDigi/plugins/Phase2TrackerDebugProducer.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/src/fed_header.h"
#include "DataFormats/FEDRawData/src/fed_trailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDFEDebug.h"
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

  Phase2TrackerDebugProducer::Phase2TrackerDebugProducer( const edm::ParameterSet& pset ) :
    cabling_(0)
  {
    // define product
    produces< edmNew::DetSetVector<Phase2TrackerFEDFEDebug> >("Debugs");
    token_ = consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"));
  }
  
  Phase2TrackerDebugProducer::~Phase2TrackerDebugProducer()
  {
  }
  
  void Phase2TrackerDebugProducer::beginJob( )
  {
  }
  
  void Phase2TrackerDebugProducer::beginRun( edm::Run const& run, edm::EventSetup const& es)
  {
    // fetch cabling from event setup
    edm::ESHandle<Phase2TrackerCabling> c;
    es.get<Phase2TrackerCablingRcd>().get( c );
    cabling_ = c.product();
  }
  
  void Phase2TrackerDebugProducer::endJob()
  {
  }
  
  void Phase2TrackerDebugProducer::produce( edm::Event& event, const edm::EventSetup& es)
  {
    std::unique_ptr<edmNew::DetSetVector<Phase2TrackerFEDFEDebug>> debugs( new edmNew::DetSetVector<Phase2TrackerFEDFEDebug>() ); 

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
        LogTrace("Phase2TrackerDebugProducer") << "[Phase2Tracker::Phase2TrackerDebugProducer::"<<__func__<<"]: \n";
        LogTrace("Phase2TrackerDebugProducer") << "Skipping invalid buffer for FED nr " << *fedIndex << endl;
        continue; 
      }
      std::vector<Phase2TrackerFEDFEDebug> all_fed_debugs = buffer.trackerHeader().CBCStatus();
      std::vector<bool> fe_status = buffer.trackerHeader().frontendStatus();
      // loop on FE
      for ( int ife = 0; ife < MAX_FE_PER_FED; ife++ )
      {
        if ( fe_status[ife] )
        {
          // get fedid from cabling
          const Phase2TrackerModule mod = cabling_->findFedCh(std::make_pair(*fedIndex, ife));
          uint32_t detid = mod.getDetid();
          // fill one debug using fastfiller
          edmNew::DetSetVector<Phase2TrackerFEDFEDebug>::FastFiller spct(*debugs, detid);
          spct.push_back(all_fed_debugs[ife]);
        }
      } // end loop on FE
    } // en loop on FED   
    // store debugs
    event.put(std::move(debugs), "Debugs" );
  } 
}
