#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerCommissioningDigi.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace Phase2Tracker {

  class Phase2TrackerCommissioningDigiProducer : public edm::global::EDProducer<> {
  public:
    /// constructor
    Phase2TrackerCommissioningDigiProducer(const edm::ParameterSet& pset);
    /// default constructor
    ~Phase2TrackerCommissioningDigiProducer() override = default;
    void produce(edm::StreamID, edm::Event& ev, const edm::EventSetup& es) const override;

  private:
    edm::EDGetTokenT<FEDRawDataCollection> token_;
  };
}  // namespace Phase2Tracker

#include "FWCore/Framework/interface/MakerMacros.h"
typedef Phase2Tracker::Phase2TrackerCommissioningDigiProducer Phase2TrackerCommissioningDigiProducer;
DEFINE_FWK_MODULE(Phase2TrackerCommissioningDigiProducer);

using namespace std;

namespace Phase2Tracker {

  typedef std::vector< std::vector<Phase2TrackerCommissioningDigi> > condata_map;  

  Phase2Tracker::Phase2TrackerCommissioningDigiProducer::Phase2TrackerCommissioningDigiProducer( const edm::ParameterSet& pset )
  {
    produces<condata_map>("ConditionData");
    token_ = consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"));
  }
  
  Phase2Tracker::Phase2TrackerCommissioningDigiProducer::~Phase2TrackerCommissioningDigiProducer()
  {
  }
  
  void Phase2Tracker::Phase2TrackerCommissioningDigiProducer::produce( edm::Event& event, const edm::EventSetup& es)
  {
    // Retrieve FEDRawData collection
    edm::Handle<FEDRawDataCollection> buffers;
    event.getByToken( token_, buffers );

     // fill collection
     std::auto_ptr<condata_map> cdigis( new condata_map ); /* switch to unique_ptr in CMSSW 7 */

     size_t fedIndex;
     for( fedIndex = Phase2Tracker::FED_ID_MIN; fedIndex < Phase2Tracker::CMS_FED_ID_MAX; ++fedIndex )
     {
       // reading
       const FEDRawData& fed = buffers->FEDData(fedIndex);
       if(fed.size()==0) continue;
       Phase2Tracker::Phase2TrackerFEDBuffer buffer(fed.data(),fed.size());
       std::map<uint32_t,uint32_t> cond_data = buffer.conditionData();

       // DEBUG
       LogTrace("Phase2TrackerCommissioningDigiProducer") << "--- Condition data debug ---" << std::endl;
       for(auto it = cond_data.begin(); it != cond_data.end(); it++)
       {
         LogTrace("Phase2TrackerCommissioningDigiProducer") << std::hex << "key: " << it->first
                                                            << std::hex << " value: " << it->second << " (hex) "
                                                            << std::dec               << it->second << " (dec) " << std::endl;
       }
       LogTrace("Phase2TrackerCommissioningDigiProducer") << "----------------------------" << std::endl;

       // storage
       std::vector<Phase2TrackerCommissioningDigi> fed_com;
       for(auto it = cond_data.begin(); it != cond_data.end(); it++)
       {
         fed_com.push_back(Phase2TrackerCommissioningDigi(it->first,it->second)); 
       }
       cdigis->push_back(fed_com);
     }
     event.put(cdigis, "ConditionData");
  }
} // end of Phase2Tracker namespace
