#ifndef CondFormats_SiStripObjects_Phase2TrackerCabling_H
#define CondFormats_SiStripObjects_Phase2TrackerCabling_H

#include "CondFormats/Serialization/interface/Serializable.h"
#include <CondFormats/SiStripObjects/interface/Phase2TrackerModule.h>
#include <vector>
#include <algorithm>

class Phase2TrackerCabling {
  typedef std::vector<Phase2TrackerModule> Store;

public:
  typedef Store::const_iterator key;
  typedef std::vector<key> Cabling;

  // Constructor taking FED channel connection objects as input.
  Phase2TrackerCabling(const std::vector<Phase2TrackerModule>& cons);

  // Copy ocnstructor
  Phase2TrackerCabling(const Phase2TrackerCabling& src);

  // Default constructor
  Phase2TrackerCabling() {}

  // Default destructor
  virtual ~Phase2TrackerCabling() {}

  // Initialize the internal maps
  void initializeCabling();

  // get the list of modules
  const Store& connections() const { return connections_; }

  // get ordered collections
  const Cabling& orderedConnections(int) const;

  // find a connection for a given fed channel
  const Phase2TrackerModule& findFedCh(std::pair<unsigned int, unsigned int> fedch) const;

  // find a connection for a given detid
  const Phase2TrackerModule& findDetid(uint32_t detid) const;

  // find a connection for a given gbtid
  const Phase2TrackerModule& findGbtid(uint32_t gbtid) const;

  // return all the modules connected to a given cooling line
  Phase2TrackerCabling filterByCoolingLine(uint32_t coolingLine) const;

  // return all the modules connected to a given HV group
  Phase2TrackerCabling filterByPowerGroup(uint32_t powerGroup) const;

  // return all fedids
  std::vector<int> listFeds() const;

  // Says which of all 72 inputs of specified DTC are connected to a module.
  std::vector<bool> connectedInputs(unsigned int fedid) const;

  // print a summary of the content
  std::string summaryDescription() const;

  // print the details of the content
  std::string description(bool compact = false) const;

private:
  // Dummy to represent unconnected channel.
  Phase2TrackerModule dummyModule_ COND_TRANSIENT;

  // the connections
  Store connections_;

  // indices for fast searches
  Cabling fedCabling_ COND_TRANSIENT;
  Cabling gbtCabling_ COND_TRANSIENT;
  Cabling detCabling_ COND_TRANSIENT;

private:
  // sorting functions
  static bool chOrdering(key a, key b);
  static bool chComp(key a, std::pair<unsigned int, unsigned int> b);
  static bool fedeq(key a, key b);
  static bool detidOrdering(key a, key b);
  static bool detidComp(key a, uint32_t b);
  static bool gbtidOrdering(key a, key b);
  static bool gbtidComp(key a, uint32_t b);
  static bool coolingOrdering(const Phase2TrackerModule& a, const Phase2TrackerModule& b);
  static bool coolingComp(const Phase2TrackerModule& a, uint32_t b);
  static bool cooleq(const Phase2TrackerModule& a, const Phase2TrackerModule& b);
  static bool powerOrdering(const Phase2TrackerModule& a, const Phase2TrackerModule& b);
  static bool powerComp(const Phase2TrackerModule& a, uint32_t b);
  static bool poweq(const Phase2TrackerModule& a, const Phase2TrackerModule& b);

  COND_SERIALIZABLE;
};

#endif
