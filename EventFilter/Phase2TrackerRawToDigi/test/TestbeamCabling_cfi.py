import FWCore.ParameterSet.Config as cms

# set cabling by hand (typically for testbeam)
Phase2TrackerCabling = cms.ESSource("Phase2TrackerCablingCfgESSource",
    modules = cms.VPSet(
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50000), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(0), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50004), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(1), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50008), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(2), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50012), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(3), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50016), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(4), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50020), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(5), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(6), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(7), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(8), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(9), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(10), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(11), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(12), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(13), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(14), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(50024), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(15), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 ),
                 cms.PSet( # Phase2 tracker module connection
                   moduleType=cms.string("2S"),
                   detid=cms.uint32(51028), 
                   gbtid=cms.uint32(0), 
                   fedid=cms.uint32(51), 
                   fedch=cms.uint32(16), 
                   powerGroup=cms.uint32(0), 
                   coolingLoop=cms.uint32(0)
                 )
              )
)
