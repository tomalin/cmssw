# CMSSW_TrackFinding_HLSFramework

### Setup instructions
```
cmsrel CMSSW_11_1_0_pre6
cd CMSSW_11_1_0_pre6/src/
cmsenv
git-cms-init

mkdir L1Trigger
cd L1Trigger/
git clone ssh://git@gitlab.cern.ch:7999/ejclemen/cmssw_trackfinding_hlsframework.git .
cd ../

# Get latest DTC branch, not yet in release
# This is getting a branch called "emyr" from Thomas' repo
git cms-addpkg L1Trigger/TrackerDTC
git cms-addpkg DataFormats/L1TrackTrigger
git remote add tschuh git@github.com:tschuh/cmssw.git
git fetch tschuh emyr
git checkout -b tschuh_emyr tschuh/emyr

scram b -j 8

cd L1Trigger/TrackFindingTrackletHLS/test
cmsRun config.py Events=1
```
The producer currently prints, for each DTC, the link word needed by the IR, and the 39 bit word for each stub.
