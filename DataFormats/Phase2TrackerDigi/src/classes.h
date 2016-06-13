#ifndef PHASE2TRACKERDIGI_CLASSES_H
#define PHASE2TRACKERDIGI_CLASSES_H

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include <vector>

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerCommissioningDigi.h"

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerHeaderDigi.h"
namespace {
  struct dictionary5 {
    edm::Wrapper<Phase2TrackerHeaderDigi > pcom0;
    edm::Wrapper<std::vector<Phase2TrackerHeaderDigi> > pcom1;
  };
}

#endif // PHASE2TRACKERDIGI_CLASSES_H
