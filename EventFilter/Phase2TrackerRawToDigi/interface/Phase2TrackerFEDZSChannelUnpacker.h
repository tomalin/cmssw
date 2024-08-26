#ifndef EventFilter_Phase2TrackerRawToDigi_FEDZSChannelUnpacker_H
#define EventFilter_Phase2TrackerRawToDigi_FEDZSChannelUnpacker_H

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include <cstdint>

namespace Phase2Tracker {

  class Phase2TrackerFEDZSChannelUnpacker {
  public:
    Phase2TrackerFEDZSChannelUnpacker(const Phase2TrackerFEDChannel& channel);
    virtual ~Phase2TrackerFEDZSChannelUnpacker() = default;
    bool hasData() const { return (clustersLeft_ > 0); }
    // go to next clusters, merge adjacent clusters (default method)
    Phase2TrackerFEDZSChannelUnpacker& operator++();
    Phase2TrackerFEDZSChannelUnpacker& operator++(int);
    inline int clusterX() const { return clusterx_; }
    inline int clusterY() const { return clustery_; }
    inline int clusterSize() const { return clustersize_; }
    virtual int unMergedX() const = 0;
    virtual int unMergedY() const = 0;
    virtual int unMergedSize() const = 0;
    // merge next clusters if adjacent
    void Merge();

  protected:
    // go to next cluster without merging adjacent clusters
    Phase2TrackerFEDZSChannelUnpacker& advance();
    // raw methods (bitwise operations)
    virtual uint8_t chipId() const = 0;
    virtual uint8_t rawX() const = 0;
    virtual uint8_t rawSize() const = 0;
    virtual bool gluedToNextCluster() const = 0;
    const uint8_t* data_;
    uint16_t currentOffset_;  // caution : this is in bits, not bytes
    int clustersLeft_;
    uint8_t clusterdatasize_;
    int clusterx_;
    int clustery_;
    int clustersize_;
  };

  class Phase2TrackerFEDZSSon2SChannelUnpacker : public Phase2TrackerFEDZSChannelUnpacker {
  public:
    Phase2TrackerFEDZSSon2SChannelUnpacker(const Phase2TrackerFEDChannel& channel)
        : Phase2TrackerFEDZSChannelUnpacker(channel) {
      clusterdatasize_ = Son2S_CLUSTER_SIZE_BITS;
      clustersLeft_ = channel.length() * 8 / clusterdatasize_;
    }
    ~Phase2TrackerFEDZSSon2SChannelUnpacker() = default;
    inline uint8_t chipId() const { return (uint8_t)read_n_at_m_L2R(data_, 4, currentOffset_); }
    inline uint8_t rawX() const { return (uint8_t)read_n_at_m_L2R(data_, 8, currentOffset_ + 4); }
    inline uint8_t rawSize() const { return (uint8_t)read_n_at_m_L2R(data_, 3, currentOffset_ + 12) + 1; }
    inline int Plane() const { return rawX() % 2; }
    int unMergedX() const;
    int unMergedY() const;
    inline int unMergedSize() const { return rawSize(); }
    bool gluedToNextCluster() const;
    Phase2TrackerFEDZSSon2SChannelUnpacker next() const;
  };

  class Phase2TrackerFEDZSSonPSChannelUnpacker : public Phase2TrackerFEDZSChannelUnpacker {
  public:
    Phase2TrackerFEDZSSonPSChannelUnpacker(const Phase2TrackerFEDChannel& channel)
        : Phase2TrackerFEDZSChannelUnpacker(channel) {
      clusterdatasize_ = SonPS_CLUSTER_SIZE_BITS;
      clustersLeft_ = channel.length() * 8 / clusterdatasize_;
    }
    ~Phase2TrackerFEDZSSonPSChannelUnpacker() = default;
    inline uint8_t chipId() const { return (uint8_t)read_n_at_m_L2R(data_, 4, currentOffset_); }
    inline uint8_t rawX() const { return (uint8_t)read_n_at_m_L2R(data_, 8, currentOffset_ + 4); }
    inline uint8_t rawSize() const { return (uint8_t)read_n_at_m_L2R(data_, 3, currentOffset_ + 12) + 1; }
    inline uint8_t threshold() const { return (uint8_t)read_n_at_m_L2R(data_, 1, currentOffset_ + 15); }
    int unMergedX() const;
    int unMergedY() const;
    inline int unMergedSize() const { return rawSize(); }
    bool gluedToNextCluster() const;
    Phase2TrackerFEDZSSonPSChannelUnpacker next() const;
  };

  class Phase2TrackerFEDZSPonPSChannelUnpacker : public Phase2TrackerFEDZSChannelUnpacker {
  public:
    Phase2TrackerFEDZSPonPSChannelUnpacker(const Phase2TrackerFEDChannel& channel)
        : Phase2TrackerFEDZSChannelUnpacker(channel) {
      clusterdatasize_ = P_CLUSTER_SIZE_BITS;
      clustersLeft_ = channel.length() * 8 / clusterdatasize_;
    }
    ~Phase2TrackerFEDZSPonPSChannelUnpacker() = default;
    inline uint8_t chipId() const { return (uint8_t)read_n_at_m_L2R(data_, 4, currentOffset_); }
    inline uint8_t rawX() const { return (uint8_t)read_n_at_m_L2R(data_, 7, currentOffset_ + 4); }
    inline uint8_t rawY() const { return (uint8_t)read_n_at_m_L2R(data_, 4, currentOffset_ + 11); }
    inline uint8_t rawSize() const { return (uint8_t)read_n_at_m_L2R(data_, 3, currentOffset_ + 15) + 1; }
    int unMergedX() const;
    int unMergedY() const;
    inline int unMergedSize() const { return rawSize(); }
    bool gluedToNextCluster() const;
    Phase2TrackerFEDZSPonPSChannelUnpacker next() const;
  };
}  // namespace Phase2Tracker
#endif  // } end def EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDZSChannelUnpacker_H
