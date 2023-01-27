# Foreword
If you contribute back to the main repo, please make sure to read through README.md and README_developer.md and make sure they are up to date with your changes. Make sure that all configs in example folder are working. The main objective of the virtuallinac is to provide 
* Treatment head modelling 
* PhaseSpace recording 
* Dose measurement 
* Flux measurement
If there is need for a new goal then it should be carefully evaluated if the changes should be kept separate or integrated to the main branch.

# Big picture
The virtuallinac has the following major parts
* Treatment head (Truebeam, Halcyon). Currently not refactored. Has the monitor chamber sensitive detector.
* Paralled worlds which contain the detectors, DoseGrid, FluxGrid and PhaseSpacePlane.
* LinacRunAction, LinacRunAction. These are major parts for collecting data during run (hits on sensitive detectors). Each LinacRun corresponds to one thread where each thread has its own instance of sensitive detectors (see ConstructSD method for UserDetector and ParallelWorlds). Then LinacRun::Merge is used to collect the data from multiple threads into one LinacRun (main thread) and data of that LinacRun is written to files in LinacRunAction::EndOfRunAction (on main thread). 
* LinacActionInitialization is a "dummy" class which the only purpose is to "install" our own classes into the Geant4 process.
* Main point which allows threading is that the sensitive detectors are created and stored in ParallelWorld::ConstructSD (thread-local method) and then they are fetched in LinacRun (which is also thread-local).
* Biggest annoyance is that for example logic of PhaseSpacePlane has to be splitten up to multiple places. E.g. PhaseSpacewWriter handles the hits, ParallelWorld handles construction of geometry and sensitive detector, LinacRun collects and sums the particle data (by accessing PhaseSpaceWriter) and then RunAction writes it to file. So keeping this sane is pretty difficult.
* Visualization works with Qt library
* Millenium 120 MLC modelled with CADmesh and STL files from drawings
* Halcyon MLC checked with visualisation and fixed accordingly to drawings

GEANT4 users might wonder that macros are not used to configure the program. Currently the treatment head is constructed according to the configuration at the beginning and then it is considered to be complete and not touched anymore. This is to make it more simple and straightforward.

# Configuration
Configuration of the virtual linac is currently done with .toml file instead of GEANT4 macros.
Only criteria for the configuration file format is that there are ready libraries for reading/writing the format in Python/C++. JSON config files would have been the other choices but the choice of toml is because for plain dictionaries toml is more simple and comments are possible. Some documentation of .toml format (which is pretty simple):

* TOML-format: https://toml.io/en/
* TOML library used to read in the configuration: https://marzer.github.io/tomlplusplus/


# Misc notes
* Uncertainty in monitor chamber dose output is uncertainty due to different amount of monitor dose in different threads. It is calculated with running mean (and variance). For additional info see: <https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm>
* If you want to create research version where you do one thing then consider creating new branch for that version alone so that the main version is not too complicated. You can also ditch the configuration format and just create the configuration in C++ code.

# TODO
* Refactor treatment heads.
* Check geometries: flattening filter, target

# Known issues
* Minor issue: WaterPhantom (with SSD > 96cm) will overlap with gantry mother volume in TrueBeam. We get around this issue by defining phantom in parallel world using layeredMass feature.
* Minor issue: Visualization doesn't work with Windows server computer