# CMSSW_TrackFinding_HLSFramework

### Setup instructions
```
export TOPDIR=$PWD
export CMSSW_VERSION=CMSSW_11_2_0_pre1
cmsrel $CMSSW_VERSION
cd $CMSSW_VERSION/src/
cmsenv
git-cms-init

mkdir L1Trigger
cd L1Trigger/
git clone ssh://git@gitlab.cern.ch:7999/ejclemen/cmssw_trackfinding_hlsframework.git .
git checkout -b AddHLSIR origin/AddHLSIR
cd ../

git cms-addpkg L1Trigger/TrackerDTC
git cms-addpkg DataFormats/L1TrackTrigger
git remote add emyr git@github.com:EmyrClement/cmssw.git
git fetch emyr DTC_For_HLS
git checkout -b DTC_For_HLS emyr/DTC_For_HLS
```

### Setting up HLS repo

Currently assuming you have set up CMSSW as in previous instructions, and you are in top level dir (i.e. directory containing $CMSSW_VERSION)
```
cd $TOPDIR
git clone git@github.com:EmyrClement/firmware-hls.git
cd firmware-hls
export HLS_BASE=$PWD
git checkout -b IR_final origin/IR_final
sed -i 's|#include "hls_math.h"|//#include "hls_math.h"|g' TrackletAlgorithm/InputRouterTop.h

cd $HLS_BASE/emData
./clean.sh
./download.sh
cd ../
cp -r emData $CMSSW_BASE/src/L1Trigger/TrackFindingTrackletHLS/plugins/

cd $CMSSW_BASE/src/
cmsenv
export hlsIncludeDir=`scram tool info hls | grep 'HLS_BASE' | cut -d '=' -f 2`/include

cd $HLS_BASE/project
c++ -std=c++11 -c -Wall  -fpic -I$hlsIncludeDir -I$HLS_BASE/firmware-hls/TrackletAlgorithm/ ../TrackletAlgorithm/InputRouterTop.cc
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
cd L1Trigger/TrackFindingTrackletHLS/
cmsRun config.py Events=10
```

### All commands for quick copy&paste
```
export TOPDIR=$PWD
export CMSSW_VERSION=CMSSW_11_2_0_pre1
cmsrel $CMSSW_VERSION
cd $CMSSW_VERSION/src/
cmsenv
git-cms-init
mkdir L1Trigger
cd L1Trigger/
git clone ssh://git@gitlab.cern.ch:7999/ejclemen/cmssw_trackfinding_hlsframework.git .
git checkout -b AddHLSIR origin/AddHLSIR
cd ../
git cms-addpkg L1Trigger/TrackerDTC
git cms-addpkg DataFormats/L1TrackTrigger
git remote add emyr git@github.com:EmyrClement/cmssw.git
git fetch emyr DTC_For_HLS
git checkout -b DTC_For_HLS emyr/DTC_For_HLS
cd $TOPDIR
git clone git@github.com:EmyrClement/firmware-hls.git
cd firmware-hls
export HLS_BASE=$PWD
git checkout -b IR_final origin/IR_final
sed -i 's|#include "hls_math.h"|//#include "hls_math.h"|g' TrackletAlgorithm/InputRouterTop.h
cd $HLS_BASE/emData
./clean.sh
./download.sh
cd ../
cp -r emData $CMSSW_BASE/src/L1Trigger/TrackFindingTrackletHLS/plugins/
cd $CMSSW_BASE/src/
cmsenv
export hlsIncludeDir=`scram tool info hls | grep 'HLS_BASE' | cut -d '=' -f 2`/include
cd $HLS_BASE/project
c++ -std=c++11 -c -Wall  -fpic -I$hlsIncludeDir -I$HLS_BASE/firmware-hls/TrackletAlgorithm/ ../TrackletAlgorithm/InputRouterTop.ccc++ -shared -o libTFIR.so InputRouterTop.o 
export HLS_LIB_DIR=$PWD
cd $CMSSW_BASE/src
cd L1Trigger/TrackFindingTrackletHLS
sed -i 's|.*environment name="TFIR_BASE".*|     <environment name="TFIR_BASE" default="'$HLS_BASE'"/>|g' myTrackFindingLib.xml
scram setup myTrackFindingLib.xml
cd $CMSSW_BASE/src
scram b -j 8
cd L1Trigger/TrackFindingTrackletHLS/
cmsRun config.py Events=10
```


### Setup when not starting from scratch
Assume you have done cmsenv and are in CMSSW_X_Y_Z/src
```
export TOPDIR=$PWD/../../
cd $TOPDIR/firmware-hls
export HLS_BASE=$PWD
cd $CMSSW_BASE/src/
cmsenv
export hlsIncludeDir=`scram tool info hls | grep 'HLS_BASE' | cut -d '=' -f 2`/include
cd $HLS_BASE/project
export HLS_LIB_DIR=$PWD
cd $CMSSW_BASE/src
```