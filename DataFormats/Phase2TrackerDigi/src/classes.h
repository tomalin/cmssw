#ifndef PHASE2TRACKERDIGI_CLASSES_H
#define PHASE2TRACKERDIGI_CLASSES_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include <vector>

namespace {
  struct dictionary {
    edm::Wrapper<Phase2TrackerDigi> zs0;
    edm::Wrapper< std::vector<Phase2TrackerDigi>  > zs1;
    edm::Wrapper< edm::DetSet<Phase2TrackerDigi> > zs2;
    edm::Wrapper< std::vector<edm::DetSet<Phase2TrackerDigi> > > zs3;
    edm::Wrapper< edm::DetSetVector<Phase2TrackerDigi> > zs4;
  };
}

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerStub.h"
namespace {
  struct dictionary2 {
    edm::Wrapper<Phase2TrackerStub> zs0;
    edm::Wrapper< std::vector<Phase2TrackerStub>  > zs1;
    edm::Wrapper< edm::DetSet<Phase2TrackerStub> > zs2;
    edm::Wrapper< std::vector<edm::DetSet<Phase2TrackerStub> > > zs3;
    edm::Wrapper< edm::DetSetVector<Phase2TrackerStub> > zs4;
  };
}

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerCommissioningDigi.h"
namespace {
  struct dictionary4 {
    edm::Wrapper<Phase2TrackerCommissioningDigi > pcom0;
    edm::Wrapper<edm::DetSet<Phase2TrackerCommissioningDigi> > pcom1;
    edm::Wrapper<std::vector< std::vector<Phase2TrackerCommissioningDigi> > > pcom2;
  };
}

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerHeaderDigi.h"
namespace {
  struct dictionary5 {
    edm::Wrapper<Phase2TrackerHeaderDigi > pcom0;
    edm::Wrapper<std::vector<Phase2TrackerHeaderDigi> > pcom1;
  };
}

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDFEDebug.h"
namespace {
  struct dictionary6 {
    edm::Wrapper<Phase2Tracker::Phase2TrackerFEDFEDebug> zs0;
    edm::Wrapper< std::vector<Phase2Tracker::Phase2TrackerFEDFEDebug>  > zs1;
    edm::Wrapper< edm::DetSet<Phase2Tracker::Phase2TrackerFEDFEDebug> > zs2;
    edm::Wrapper< std::vector<edm::DetSet<Phase2Tracker::Phase2TrackerFEDFEDebug> > > zs3;
    edm::Wrapper< edm::DetSetVector<Phase2Tracker::Phase2TrackerFEDFEDebug> > zs4;
  };
}

#endif // PHASE2TRACKERDIGI_CLASSES_H
