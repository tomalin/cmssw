#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"

namespace Phase2Tracker {
  class StackedDigi {
  public:
    // layer is 0/1 for lower/upper sensor in module (lower is p-type for PS).
    // moduletype is 0/1 for 2S/PS module.
    StackedDigi() {}
    StackedDigi(const Phase2TrackerCluster1D *, STACK_LAYER layer, int moduletype);
    StackedDigi(int digix, int digiy, int sizex, STACK_LAYER layer, int moduletype);
    StackedDigi(int digix, int digiy, int sizex, STACK_LAYER layer, int moduletype, int threshold);
    ~StackedDigi() {}
    bool operator<(StackedDigi) const;
    inline STACK_LAYER getLayer() const { return layer_; }
    inline int getModuleType() const { return moduletype_; }
    inline int getRawX() const { return rawx_; }
    inline int getRawY() const { return rawy_; }
    inline int getDigiX() const { return digix_; }
    inline int getDigiY() const { return digiy_; }
    inline int getSizeX() const { return sizex_; }
    inline int getSide() const { return side_; } // Which end of module: 0 or 1.
    inline int getThreshold() const { return threshold_; }
    // get side and type, to map to concentrators (0 = P-left, 1 = P-right, 2 = S-left, 3 = S-right)
    inline int getSideType() const { return side_ + 2 * (1 - moduletype_) + 2 * layer_ * moduletype_; }
    inline int getChipId() const { return chipid_; }
    inline int getStripsX() const { return (moduletype_ == 1) ? PS_ROWS : (STRIPS_PER_CBC / 2); }
    void setPosSizeX(int, int);
    bool shouldSplit() const;
    std::vector<StackedDigi> splitDigi() const;

  private:
    void calcchipid();
    int digix_, digiy_, sizex_;
    int threshold_;
    STACK_LAYER layer_;
    int moduletype_;
    int rawx_, rawy_;
    int side_;
    int chipid_;
  };

  StackedDigi::StackedDigi(const Phase2TrackerCluster1D *digi, STACK_LAYER layer, int moduletype)
      : digix_(digi->firstStrip()),
        digiy_(digi->column()),
        sizex_(digi->size()),
        threshold_(digi->threshold()),
        layer_(layer),
        moduletype_(moduletype) {
    calcchipid();
  }

  StackedDigi::StackedDigi(int digix, int digiy, int sizex, STACK_LAYER layer, int moduletype)
      : digix_(digix), digiy_(digiy), sizex_(sizex), threshold_(0), layer_(layer), moduletype_(moduletype) {
    calcchipid();
  }

  StackedDigi::StackedDigi(int digix, int digiy, int sizex, STACK_LAYER layer, int moduletype, int threshold)
      : digix_(digix), digiy_(digiy), sizex_(sizex), threshold_(threshold), layer_(layer), moduletype_(moduletype) {
    calcchipid();
  }

  bool StackedDigi::shouldSplit() const {
    if (getSizeX() > 8 or getDigiX() % getStripsX() + getSizeX() > getStripsX())
      return true;
    return false;
  }

  std::vector<StackedDigi> StackedDigi::splitDigi() const {
    std::vector<StackedDigi> parts;
    if (shouldSplit()) {
      int pos = getDigiX();
      int end = pos + getSizeX();
      int isize;
      while (pos < end) {
        int nextchip = (pos / getStripsX() + 1) * getStripsX();
        isize = std::min(std::min(8, end - pos), nextchip - pos);
        StackedDigi ndig(*this);
        ndig.setPosSizeX(pos, isize);
        parts.push_back(ndig);
        pos += isize;
      }
#ifdef EDM_ML_DEBUG
      std::ostringstream ss;
      ss << "--- Split digi at " << getDigiX() << ", raw: " << getRawX() << ", length: " << getSizeX() << std::endl;
      for (auto id = parts.begin(); id < parts.end(); id++) {
        ss << " -- " << id->getDigiX() << ", raw:  " << id->getRawX() << " " << id->getSizeX() << std::endl;
      }
      ss << "--- End of split ---" << std::endl;
      LogTrace("Phase2TrackerDigiProducer") << ss.str();
      ss.clear();
      ss.str("");
#endif
    } else {
      parts.push_back(*this);
    }
    return parts;
  }

  void StackedDigi::setPosSizeX(int pos, int size) {
    digix_ = pos;
    sizex_ = size;
    calcchipid();
  }

  void StackedDigi::calcchipid() {
    int x = digix_;
    int y = digiy_;
    side_ = 0;
    if (moduletype_ == 1) {
      if (layer_ == LAYER_INNER) {
        // PonPS
        chipid_ = x / PS_ROWS;
        // which side ?
        if (y >= PS_COLS / 2) {
          chipid_ += 8;
          side_ = 1;
          rawy_ = y - PS_COLS / 2;
        } else {
          rawy_ = y;
        }
        rawx_ = x % PS_ROWS;
      } else {
        // SonPS
        chipid_ = x / PS_ROWS;
        if (y > 0) {
          chipid_ += 8;
          side_ = 1;
        }
        rawx_ = x % PS_ROWS;
        rawy_ = y;
      }
    } else {
      x *= 2;
      if (layer_ == LAYER_OUTER) {
        x += 1;
      }
      chipid_ = x / STRIPS_PER_CBC;
      // which side ?
      if (y > 0) {
        chipid_ += 8;
        side_ = 1;
      }
      rawx_ = x % STRIPS_PER_CBC;
      rawy_ = y;
    }
  }

// Used for sorting of StackedDigis
bool StackedDigi::operator<(const StackedDigi d2) const {
    // Ensures lower sensor before upper one. (lower is p-type for PS).
    if (layer_ < d2.getLayer()) {
      return true;
    }
    if (layer_ > d2.getLayer()) {
      return false;
    }
    // Ensures if two digis have same layer, then ordered by left-end / right-end.
    if (chipid_ < d2.getChipId()) {
      return true;
    }
    if (chipid_ == d2.getChipId() and rawx_ < d2.getRawX()) {
      return true;
    }
    if (chipid_ == d2.getChipId() and rawx_ == d2.getRawX() and rawy_ < d2.getRawY()) {
      return true;
    }
    return false;
  }
}  // namespace Phase2Tracker
