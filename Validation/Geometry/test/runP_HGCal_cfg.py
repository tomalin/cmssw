# In order to produce everything that you need in one go, use the command:
#
# for t in {'BeamPipe','Tracker','PixBar','PixFwdMinus','PixFwdPlus','TIB','TOB','TIDB','TIDF','TEC','TkStrct','InnerServices'}; do cmsRun runP_Tracker_cfg.py geom=XYZ label=$t >& /dev/null &; done


import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys, re

process = cms.Process("PROD")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# The default geometry is Extended2023D28. If a different geoemtry
# is needed, the appropriate flag has to be passed at command line,
# e.g.: cmsRun runP_HGCal_cfg.py geom="XYZ"

# The default component to be monitored is the HGCal. If other
# components need to be studied, they must be supplied, one at a time,
# at the command line, e.g.: cmsRun runP_HGCal_cfg.py
# label="XYZ"

from Validation.Geometry.plot_hgcal_utils import _LABELS2COMPS

_ALLOWED_LABELS = _LABELS2COMPS.keys()

options = VarParsing('analysis')
options.register('geom',             #name
                 'Extended2023D28',      #default value
                 VarParsing.multiplicity.singleton,   # kind of options
                 VarParsing.varType.string,           # type of option
                 "Select the geometry to be studied"  # help message
                )

options.register('label',         #name
                 'HGCal',              #default value
                 VarParsing.multiplicity.singleton,   # kind of options
                 VarParsing.varType.string,           # type of option
                 "Select the label to be used to create output files. Default to HGCal. If multiple components are selected, it defaults to the join of all components, with '_' as separator."  # help message
                )

options.setDefault('inputFiles', ['file:single_neutrino_random.root'])

options.parseArguments()
# Option validation

if options.label not in _ALLOWED_LABELS:
    print "\n*** Error, '%s' not registered as a valid components to monitor." % options.label
    print "Allowed components:", _ALLOWED_LABELS
    print
    raise RuntimeError("Unknown label")

_components = _LABELS2COMPS[options.label]

# Load geometry either from the Database of from files
process.load("Configuration.Geometry.Geometry%s_cff" % options.geom)

#
#Magnetic Field
#
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# Output of events, etc...
#
# Explicit note : since some histos/tree might be dumped directly,
#                 better NOT use PoolOutputModule !
# Detector simulation (Geant4-based)
#
process.load("SimG4Core.Application.g4SimHits_cfi")

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.g4SimHits.initialSeed = 9876

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

process.p1 = cms.Path(process.g4SimHits)
process.g4SimHits.StackingAction.TrackNeutrino = cms.bool(True)
process.g4SimHits.UseMagneticField = False
process.g4SimHits.Physics.type = 'SimG4Core/Physics/DummyPhysics'
process.g4SimHits.Physics.DummyEMPhysics = True
process.g4SimHits.Physics.CutsPerRegion = False
process.g4SimHits.Watchers = cms.VPSet(cms.PSet(
    type = cms.string('MaterialBudgetAction'),
    MaterialBudgetAction = cms.PSet(
        HistosFile = cms.string('matbdg_%s.root' % options.label),
        AllStepsToTree = cms.bool(True),
        HistogramList = cms.string('HGCal'),
        SelectedVolumes = cms.vstring(_components),
        TreeFile = cms.string('None'), ## is NOT requested

        StopAfterProcess = cms.string('None'),
#        TextFile = cms.string("matbdg_HGCal.txt")
        TextFile = cms.string('None')
    )
))