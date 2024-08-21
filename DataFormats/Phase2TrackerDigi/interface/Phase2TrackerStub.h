#ifndef DataFormats_Phase2TrackerDigi_Phase2TrackerStub_H
#define DataFormats_Phase2TrackerDigi_Phase2TrackerStub_H

#include <cstdint>

class Phase2TrackerStub {
public:
  Phase2TrackerStub() : posx_(0), bendcode_(0), posy_(0) {}
  Phase2TrackerStub(uint32_t posx, uint8_t bendcode, uint8_t posy = 0)
      : posx_(posx), bendcode_(bendcode), posy_(posy) {}
  ~Phase2TrackerStub() {}
  uint32_t getPositionX() const { return posx_; }
  uint8_t getPositionY() const { return posy_; }
  uint8_t getBendCode() const { return bendcode_; }

private:
  uint32_t posx_;
  uint8_t bendcode_;
  uint8_t posy_;
};

inline bool operator<(const Phase2TrackerStub& one, const Phase2TrackerStub& other) {
  return one.getPositionX() < other.getPositionX();
}

#endif
