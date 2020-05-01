# CMSSW_TrackFinding_HLSFramework

### Setup instructions
```
export TOPDIR=$PWD
cmsrel CMSSW_11_1_0_pre6
cd $CMSSW_BASW/src/
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

### Setting up HLS repo

Currently assuming you have set up CMSSW as in previous instructions, and you are in top level dir (i.e. directory containing CMSSW_11_1_0_pre6)
```
cd $TOPDIR
git clone git@github.com:sseifeln/firmware-hls.git
cd firmware-hls
export HLS_BASE=$PWD
git checkout -b Sarah_dev origin/Dev

cd $CMSSW_BASE/src/
cmsenv
export hlsIncludeDir=`scram tool info hls | grep 'HLS_BASE' | cut -d '=' -f 2`/include

cd $HLS_BASE/project
c++ -std=c++11 -c -Wall  -fpic -I/software/CAD/Xilinx/2019.1/Vivado/2019.1/include/ -I$HLS_BASE/firmware-hls/TrackletAlgorithm/ ../TrackletAlgorithm/InputRouterTop.cpp
c++ -shared -o libTFIR.so InputRouterTop.o 
export HLS_LIB_DIR=$PWD
```

### Going back to CMSSW
```
cd $CMSSW_BASE/src
cd L1Trigger/TrackFindingTrackletHLS
sed -i 's|.*environment name="TFIR_BASE".*|     <environment name="TFIR_BASE" default="'$HLS_BASE'"/>|g' myTrackFindingLib.xml
scram setup myTrackFindingLib.xml
cd $CMSSW_BASE/src
scram b -j 8
```