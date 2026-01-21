# LeptonInjector

LeptonInjector is a group of modules used to create events in IceCube. This code represents a standalone version of the original LeptonInjector which has been trimmed of the proprietary Icetray dependencies. It is currently fully functional and compatible with LeptonWeighter. 

To use it, you will

    1. Prepare a Injector object (or list of Injector objects).

    2. Use one injector object, along with several generation parameters, to create a Controller object. These were called MultiLeptonInjector in the original code. 

    3. Add more injectors to the controller object using the add injector function. Verify that the controller is properly configured.
    
    4. Specify the full paths and names of the destination output file and LeptonInjector Configuration (LIC) file.

    5. Call Controller.Execute(). This will run the simulation. 

For an example of this in action, see $root/resources/example/inject_muons.py

To learn about the LIC files and weighting, see https://github.com/IceCubeOpenSource/LeptonWeighter

# Installation

## Quick Start: pip install (Recommended)

The easiest way to install LeptonInjector is via pip:

```bash
pip install .
```

This will automatically:
- Build the C++ library with pybind11 bindings
- Install all Python dependencies
- Set up the package for immediate use

### Requirements for pip install

- Python >= 3.8
- A C++ compiler with C++17 support
- CMake >= 3.15
- The following libraries must be installed on your system:
  - HDF5 C libraries
  - Photospline (https://github.com/IceCubeOpenSource/photospline)
  - Boost (headers only)

On **macOS** with Homebrew:
```bash
brew install hdf5 boost cmake
# Install photospline separately following its documentation
```

On **Ubuntu/Debian**:
```bash
sudo apt-get install libhdf5-dev libboost-all-dev cmake
# Install photospline separately following its documentation
```

### Optional dependencies

For HDF5 file handling in Python:
```bash
pip install ".[hdf5]"
```

## Python Usage

Once installed, you can import and use LeptonInjector directly:

```python
import LeptonInjector as LI

# Create an injector for muon neutrino charged-current interactions
injector = LI.Injector(
    NEvents=1000,
    FinalType1=LI.Particle.ParticleType.MuMinus,
    FinalType2=LI.Particle.ParticleType.Hadrons,
    DoublyDifferentialCrossSectionFile="path/to/differential_xs.fits",
    TotalCrossSectionFile="path/to/total_xs.fits",
    Ranged=True
)

# Create a controller to manage the simulation
controller = LI.Controller(
    injectors=injector,
    minimum_energy=1e2,  # GeV
    maximum_energy=1e6,  # GeV
    spectral_index=2.0,
    minimum_azimuth=0,
    maximum_azimuth=2*LI.Constants.pi,
    minimum_zenith=0,
    maximum_zenith=LI.Constants.pi
)

# Configure output files
controller.NameOutfile("output.h5")
controller.NameLicFile("config.lic")

# Run the simulation
controller.Execute()
```

### Earth Model Services

LeptonInjector includes Earth model services for density calculations:

```python
import LeptonInjector as LI

# Create an Earth model service with PREM
ems = LI.EarthModelService(
    "mymodel",
    tablepath="/path/to/earthparams/",
    earthmodels=["PREM_mmc.dat"],
    materialmodels=["Standard.dat"],
    icecapname="SimpleIceCap",
    icecapangle=20.0 * LI.Constants.degrees,
    detectordepth=1948.0 * LI.Constants.m
)

# Get density at a position (detector-centered coordinates)
pos = LI.LI_Position(0, 0, -1000)  # 1km below detector center
density = ems.GetDensityInCGS(pos)
print(f"Density: {density} g/cm^3")

# Convert to Earth-centered coordinates
earth_pos = ems.GetEarthCoordPosFromDetCoordPos(pos)

# Get the medium type at that position
medium = ems.GetMedium(earth_pos)
print(f"Medium: {LI.EarthModelService.GetMediumTypeString(medium)}")

# Calculate column depth between two points
pos1 = LI.LI_Position(0, 0, 0)
pos2 = LI.LI_Position(0, 0, -10000)
column_depth = ems.GetColumnDepthInCGS(pos1, pos2)
print(f"Column depth: {column_depth} g/cm^2")

# Get lepton range
direction = LI.LI_Direction(LI.Constants.pi, 0)  # straight down
energy = 1e5  # 100 TeV
muon_range = ems.GetLeptonRangeInMeterFrom(pos, direction, energy, isTau=False)
print(f"Muon range: {muon_range} m")
```

### Available Earth Model Enums

```python
# Medium types
LI.MediumType.INNERCORE
LI.MediumType.OUTERCORE
LI.MediumType.MANTLE
LI.MediumType.ROCK
LI.MediumType.ICE
LI.MediumType.WATER
LI.MediumType.AIR
LI.MediumType.VACUUM

# Ice cap types
LI.IceCapType.NOICE
LI.IceCapType.ICESHEET
LI.IceCapType.SIMPLEICECAP

# Integration types for density integration
LI.IntegType.PATH    # Along a path
LI.IntegType.RADIUS  # Projected on radial direction
LI.IntegType.CIRCLE  # 2*pi*r weighted
LI.IntegType.SPHERE  # 4*pi*r^2 weighted (volume mass)

# Lepton range calculation options
LI.LeptonRangeOption.DEFAULT  # Dima's fitting function
LI.LeptonRangeOption.LEGACY   # Legacy NuGen equation
LI.LeptonRangeOption.NUSIM    # Gary's fitting function
```

### Utility Functions

```python
# Geometry calculations
impact, t, closest_pos = LI.GetImpactParameter(position, direction)
n_intersections, start_pos, end_pos = LI.GetIntersectionsWithSphere(pos, dir, radius)
n_intersections, enter_dist, exit_dist = LI.GetDistsToIntersectionsWithSphere(pos, dir, radius)

# Lepton range in meter-water-equivalent
range_mwe = LI.GetLeptonRange(energy, isTau=False)

# Unit conversions
mwe = LI.ColumnDepthCGStoMWE(column_depth_cgs)
cgs = LI.MWEtoColumnDepthCGS(range_mwe)
```

---

# Dependencies (for manual installation)

All of the dependencies are already installed on the CVMFS environments on the IceCube Cobalt testbeds.

For local installations, you need the following:

* A C++ compiler with C++17 support.

* The `HDF5` C libraries. Read more about it here: https://portal.hdfgroup.org/display/support. These libraries are, of course, used to save the data files.

* It also requires Photospline to create and to read cross sections. Read more about it, and its installation at https://github.com/IceCubeOpenSource/photospline. Note that Photospline has dependencies that you will need that are not listed here.

* LeptonInjector requires Photospline's `SuiteSparse` capabilities, whose dependencies are available here http://faculty.cse.tamu.edu/davis/suitesparse.html

For building py-bindings:

* Python >= 3.8

* `pybind11` (used for the pip-installable Python bindings)

* `BOOST` headers (for some internal functionality)


# Included Dependencies

These are not ostensibly a part of LeptonInjector, but are included for its functionality. They were developed by the IceCube Collaboration and modified slightly to use the LeptonInjector datatypes instead of the IceCube proprietary ones. 

* I3CrossSections: provies the tools for sampling DIS and GR cross sections. 

* Earthmodel Services: provides the PREM for column depth calculations. 

# Manual Compilation with CMake (Alternative)

If you prefer to build manually or need more control over the build process, you can use CMake directly.

First, clone the repository:

```bash
git clone git@github.com:IceCubeOpenSource/LeptonInjector.git
cd LeptonInjector
```

Create build and install directories:

```bash
mkdir build install
cd build
```

Configure with CMake:

```bash
cmake -DCMAKE_INSTALL_PREFIX=../install ..
```

Build and install:

```bash
make -j4 && make install
```

Set up environment variables (add to your `.bashrc` or `.bash_profile`):

```bash
# Allow Python to find the installed library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/your/install/path/

# Allow the EarthModel to find Earth parameter files
export EARTH_PARAMS=/path/to/LeptonInjector/resources/earthparams/
```

**Note:** When using pip install, these environment variables are not needed as the package is installed into your Python environment and resources are bundled with the package.

# Structure
The code base is divided into several files. 
* Constants: a header defining various constants. 
* Controller: implements a class for managing the simulation
* DataWriter: writes event properties and MCTrees to an HDF5 file
* EventProps: implements a few structures used to write events in the hdf5 file. 
* h5write: may be renamed soon. This will be used to write the configurations onto a file
* LeptonInjector (the file): defines the Injector objects described above in addition to several configuration objects and event generators 
* Particle: simple implementation of particles. Includes a big enum. 
* Random: object for random number sampling.

# Cross Sections
For generating events you will need fits files of splines specifying the cross sections (total and differential cross sections). These should be made with photospline. 

# Making Contributions
If you would like to make contributions to this project, please create a branch off of the `master` branch and name it something following the template: `$YourLastName/$YourSubProject`. 
Work on this branch until you have made the changes you wished to see and your branch is _stable._ 
Then, pull from master, and create a pull request to merge your branch back into master. 

# Detailed Author Contributions and Citation

The LeptonInjector and LeptonWeighter modules were motivated by the high-energy light sterile neutrino search performed by B. Jones and C. Argüelles. C. Weaver wrote the first implementation of LeptonInjector using the IceCube internal software framework, icetray, and wrote the specifications for LeptonWeighter. In doing so, he also significantly enhanced the functionality of IceCube's Earth-model service. These weighting specifications were turned into code by C. Argüelles in LeptonWeighter. B. Jones performed the first detailed Monte Carlo comparisons that showed that this code had similar performance to the standard IceCube neutrino generator at the time for throughgoing muon neutrinos.

It was realized that these codes could have use beyond IceCube and could benefit the broader neutrino community. The codes were copied from IceCube internal subversion repositories to this GitHub repository; unfortunately, the code commit history was not preserved in this process. Thus the current commits do not represent the contributions from the original authors, particularly from the initial work by C. Weaver and C. Argüelles. 

The transition to this public version of the code has been spearheaded by A. Schneider and B. Smithers, with significant input and contributions from C. Weaver and C. Argüelles. B. Smithers isolated the components of the code needed to make the code public, edited the examples, and improved the interface of the code. A. Schneider contributed to improving the weighting algorithm, particularly to making it work for volume mode cascades, as well as in writing the general weighting formalism that enables joint weighting of volume and range mode.

This project also received contributions and suggestions from internal IceCube reviewers and the collaboration as a whole. Please cite this work as:

LeptonInjector and LeptonWeighter: A neutrino event generator and weighter for neutrino observatories
IceCube Collaboration
https://arxiv.org/abs/2012.10449

## CRediT

**Austin Schneider**: Software, Validation, Writing - Original Draft, Writing - Review & Editing;
**Benjamin Jones**: Conceptualization, Validation;
**Benjamin Smithers**: Software, Validation, Writing - Original Draft, Visualization, Writing - Review & Editing;
**Carlos Argüelles**: Conceptualization, Software, Writing - Original Draft, Writing - Review & Editing, Supervision;
**Chris Weaver**: Methodology, Software, Writing - Review & Editing
