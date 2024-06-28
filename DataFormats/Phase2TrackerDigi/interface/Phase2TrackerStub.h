#ifndef DataFormats_Phase2TrackerDigi_Phase2TrackerStub_H
#define DataFormats_Phase2TrackerDigi_Phase2TrackerStub_H

#include <cstdint>

// stub bend codes to actual binary values
static const float STUB_BEND_MAP[16] = {
    0,     // -0.5,0,0.5
    1,     // 1
    2,     // 1.5,2,2.5
    3.5,   // 3,3.5,4
    4.5,   // 4.5, 5
    5.5,   // 5.5
    6,     // 6
    6.5,   // 6.5, 7
    15,    // NOT USED
    -6.5,  // -6.5, -7
    -6,    // -6
    -5.5,  // -5.5
    -4.5,  // -4.5, -5
    -3.5,  // -3,-3.5,-4
    -2,    // -1.5,-2,-2.5
    -1,    // -1
};

class Phase2TrackerStub {
public:
  Phase2TrackerStub() : posx_(0), bendcode_(0), posy_(0) {}
  Phase2TrackerStub(uint32_t posx, uint8_t bendcode, uint8_t posy = 0)
      : posx_(posx), bendcode_(bendcode), posy_(posy) {}
  ~Phase2TrackerStub() {}
  uint32_t getPositionX() const { return posx_; }
  uint8_t getPositionY() const { return posy_; }
  uint8_t getBendCode() const { return bendcode_; }
  float getBend() const { return STUB_BEND_MAP[bendcode_]; }

private:
  uint32_t posx_;
  uint8_t bendcode_;
  uint8_t posy_;
};

inline bool operator<(const Phase2TrackerStub& one, const Phase2TrackerStub& other) {
  return one.getPositionX() < other.getPositionX();
}

#endif
