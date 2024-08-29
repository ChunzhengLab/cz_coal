
### CZCOAL
- A coalescence model for studying baryon azimuthal correlations and the background of CVE
- CMake project (CMake version > 3.2.0)
- C++20 is recommended and adapted to C++11 on PC farm. Root package required to read the root file.
- Classes:
  - Namespace: par
  - Base: Event class, Particle class (inherited by Parton class and Hadrons class), DistanceFun
  - Core: Coalescence class, CalculateObvs class
  - Process: EventReader class, EventWriter class
- Included utility tools
  - Random event generator: rand_parton_gen, with 2 run mode: kPureRandom and  kAMPTAnchored
  - Coalescence center comparator
  - Automatic path getter and job submission script generator create_conf.sh and sort_paths.py
 
Usage:
```bash
./bin/Coalescence config.conf
```
Configuration example:

```
isDebug=false
isLocalDraw=false
isWriteEvents=true
isCalculateObvs=true
isRemoveHFQuarks=true
isEnableMassVarify=true
isEnableQuarkMoveOn=false
isBalanceQuarkNumber=true
eventType=kAMPT
r_bm=0.5
flavourBreakTolerance=0.001
coalescenceAlgorithm=kFromParton
inputFile=../test_data/zpc-1.root
outputFile=../test_data/dataCoalHadrons.root
obvsFile=../test_data/obvsCoalHadrons.root
```







