
### CZCOAL -- a coalescence model for studying baryon azimuthal correlations and the background of CVE
- CMake project (Cmake version > 3.2)
- C++20 is recommended and adapted to C++11 on PC farm. Root package required to read the root file.
- Classes:
  - Namespace: par
  - Base: Event class, Particle class (inherited by Parton class and Hadrons class), DistanceFun
  - Core: Coalescence class, CalculateObvs class
  - Process: EventReader class, EventWriter class
- Included Utility Tools
  - Random event generator: rand_parton_gen, with 2 run mode: kPureRandom and  kAMPTAnchored
  - Coalescence center comparator
  - Automatic path getter and job submission script generator create_conf.sh and sort_paths.py
