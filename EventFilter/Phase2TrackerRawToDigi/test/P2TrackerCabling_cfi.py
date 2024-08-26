import FWCore.ParameterSet.Config as cms


def make_vpset_fromfile(filename):
    psets = []
    with open(filename, 'r') as f:
        fedch = 0
        for line in f:
            line = line.split()
            detid  = int(line[0])
            fedid  = int(line[1])
            fedch = int(line[2])
            modtype = str(line[3])
            psets.append(cms.PSet(
                detid=cms.uint32(detid),
                gbtid=cms.uint32(0),
                fedid=cms.uint32(fedid),
                fedch=cms.uint32(fedch),
                moduleType = cms.string(modtype),
                powerGroup=cms.uint32(0),
                coolingLoop=cms.uint32(0))
                )

        """
        while fedch != 72:
            detid = 0 
            psets.append(cms.PSet(
                moduleType=cms.string("dummy"),
                detid=cms.uint32(detid),
                gbtid=cms.uint32(0),
                fedid=cms.uint32(fedid),
                fedch=cms.uint32(fedch),
                powerGroup=cms.uint32(0),
                coolingLoop=cms.uint32(0))
                )
            fedch += 1
        """
        return psets

# Original file from 8 years ago
#my_psets = make_vpset_fromfile('detids_phase2.txt')
# New file generated from TkLayout by parseTkLayoutCabling.py
my_psets = make_vpset_fromfile('tkLayoutCabling_parsed.txt')

Phase2TrackerCabling = cms.ESSource("Phase2TrackerCablingCfgESSource", modules = cms.VPSet( *my_psets ))
