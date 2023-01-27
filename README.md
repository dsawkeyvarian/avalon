# VirtualLinac
Monte Carlo simulation code for various Varian treatment heads. Based on [GEANT4](https://geant4.web.cern.ch) Monte Carlo code. Based on the original code and modeling work by Daren Sawkey.

## Features
* Simulation of radiation physics with Geant4.
* Geometry and materials for multiple treatment heads
    - TrueBeam with Millennium 120 or HD multileaf collimator
    - Halcyon with RDS multileaf collimator
* Multileaf collimator geometry is from accurate CAD models tested against measurements.
* Multiple validated beams are available (more coming).
    - TrueBeam: 6X-FFF, 6X
    - Halcyon: 6X-FFF
* Custom particle sources (photons, electrons, protons) with or without treatment head geometry.
* Phantoms
    - Water box
    - Slab model, different materials for each slab
    - Voxel box, where each voxel can have different material
* Record dose and particles (phase space) for further analysis
* Static any dynamic plans
* Optional python tools for analyzing dose and particles.

## Quick start
Install GEANT4 and the physics lists following the instructions <https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/>.
Compile with cmake. Depending on how you installed Geant4 you might need to provide the install location with CMAKE_PREFIX_PATH.
```
mkdir <build_folder>
cd <build_folder>
cmake <cloned_repository_path> -DCMAKE_PREFIX_PATH=<G4_installation_path>
cmake --build . --config Release
```
See options by running `linac.exe --help`.

## Input
VirtualLinac accepts JSON file describing the simulation as an input. See examples/ folder for input .jsonc files (.jsonc is .json with comments). VirtualLinac command line interface also accepts second JSON file and some command line argument for ease-of-use, see `linac.exe --help`.

### Coordinate system
VirtualLinac has a global coordinate system used for inputs and outputs (dose, phase-space). 
* Origin is at the isocenter. 
* +Z-direction is from the isocenter to the target.
* Looking towards the machine from the front: +X is from left to right, and +Y is into the monitor.
Outputs are (like 3D dose grids) are ordered so that Z is the fastest running index then Y and X.

## Testing
VirtualLinac comes with a test suite. It tests the functionality but do not validate physical accuracy. Tests use [doctest](https://github.com/doctest/doctest) C++ testing framework. There are two types of tests: unit tests of VirtuaLinac parts and system tests, which run the simulation and inspect the output.
* Run unit tests with `test.exe --test-suite-exclude="System"`.
* Each system test initializes Geant4 and therefore each test needs to be run in a separate process. List all available system tests `test.exe --test-suite="System" -ltc` and run each tests separately with `test.exe --test-case="<test_name> -- --out <output_folder>"` where output foldary can be temporary one but is required for VirtualLinac output files.
All tests require that VirtualLinac root folder is supplied either through environment variable "VIRTUALLINAC_ROOT" or by extra argument to tests `test.exe --test-suite-exclude="System" -- --root <root_folder>`. VirtualLinac root folder is the repository root where stl / gdml / example folders are located).
