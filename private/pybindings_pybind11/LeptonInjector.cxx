/*
 * pybind11 bindings for LeptonInjector
 *
 * This provides the same interface as the Boost.Python bindings but using
 * pybind11, which is header-only and easier to package for pip installation.
 */

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Random.h>
#include <LeptonInjector/Constants.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace LeptonInjector {
    struct LIConstants {};
}

PYBIND11_MODULE(LeptonInjector, m) {
    using namespace LeptonInjector;

    m.doc() = "LeptonInjector - Generates neutrino events for large volume Cherenkov telescopes";

    // Note: std::pair<double, double> is automatically converted by pybind11/stl.h
    // to Python tuples, so DoublePair is not needed. If explicit access is needed,
    // users can use Python tuples directly: (first, second)

    // Random number generator
    py::class_<LI_random>(m, "RNG")
        .def(py::init<unsigned int>(), "seed"_a = 1)
        .def("Uniform", &LI_random::Uniform);

    // LI_Position
    py::class_<LI_Position>(m, "LI_Position")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def(py::init<const LI_Position&>())
        .def(py::init<std::array<double, n_dimensions>>())
        .def("at", &LI_Position::at)
        .def("Magnitude", &LI_Position::Magnitude)
        .def("GetX", &LI_Position::GetX)
        .def("GetY", &LI_Position::GetY)
        .def("GetZ", &LI_Position::GetZ)
        .def("SetX", &LI_Position::SetX)
        .def("SetY", &LI_Position::SetY)
        .def("SetZ", &LI_Position::SetZ)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self += py::self)
        .def(py::self == py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self * LI_Direction())
        .def(py::self * py::self);

    // LI_Direction
    py::class_<LI_Direction>(m, "LI_Direction")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def(py::init<std::array<double, 2>>())
        .def(py::init<std::pair<double, double>>())
        .def(py::init<const LI_Direction&>())
        .def(py::init<const LI_Position&>())
        .def("GetX", &LI_Direction::GetX)
        .def("GetY", &LI_Direction::GetY)
        .def("GetZ", &LI_Direction::GetZ)
        .def_readwrite("zenith", &LI_Direction::zenith)
        .def_readwrite("azimuth", &LI_Direction::azimuth)
        .def(-py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self * LI_Position());

    // Coordinate functions
    m.def("RotateX", &RotateX);
    m.def("RotateY", &RotateY);
    m.def("RotateZ", &RotateZ);
    m.def("rotateRelative", &rotateRelative);
    m.def("computeCylinderIntersections", &computeCylinderIntersections);

    // Controller
    py::class_<Controller>(m, "Controller")
        .def(py::init<Injector, double, double, double, double, double, double, double, double, double, double, double>(),
             "injectors"_a, "minimum_energy"_a, "maximum_energy"_a, "spectral_index"_a,
             "minimum_azimuth"_a, "maximum_azimuth"_a, "minimum_zenith"_a, "maximum_zenith"_a,
             "injection_radius"_a = 1200., "endcap_length"_a = 1200.,
             "cylinder_radius"_a = 1200., "cylinder_height"_a = 1200.)
        .def("Execute", &Controller::Execute)
        .def("SetEarthModel", &Controller::setEarthModel)
        .def("AddInjector", &Controller::AddInjector)
        .def("NameOutfile", &Controller::NameOutfile)
        .def("NameLicFile", &Controller::NameLicFile)
        .def("Overwrite", &Controller::Overwrite)
        .def("setSeed", &Controller::setSeed);

    // Particle class and enum
    py::class_<Particle> particle(m, "Particle");
    particle.def(py::init<>())
        .def_readwrite("type", &Particle::type)
        .def_readwrite("energy", &Particle::energy)
        .def_readwrite("direction", &Particle::direction)
        .def_property("position",
            [](Particle& p) {
                return py::array_t<double>({3}, {sizeof(double)}, p.position, py::cast(p));
            },
            [](Particle& p, py::array_t<double> arr) {
                auto buf = arr.request();
                if (buf.size != 3) throw std::runtime_error("Position must have 3 elements");
                double* ptr = static_cast<double*>(buf.ptr);
                for (int i = 0; i < 3; i++) p.position[i] = ptr[i];
            })
        .def("GetMass", &Particle::GetMass)
        .def("HasMass", &Particle::HasMass)
        .def("GetTypeString", &Particle::GetTypeString);

    // ParticleType enum
    py::enum_<Particle::ParticleType>(particle, "ParticleType")
        .value("unknown", Particle::unknown)
        .value("Gamma", Particle::Gamma)
        .value("EPlus", Particle::EPlus)
        .value("EMinus", Particle::EMinus)
        .value("MuPlus", Particle::MuPlus)
        .value("MuMinus", Particle::MuMinus)
        .value("Pi0", Particle::Pi0)
        .value("PiPlus", Particle::PiPlus)
        .value("PiMinus", Particle::PiMinus)
        .value("K0_Long", Particle::K0_Long)
        .value("KPlus", Particle::KPlus)
        .value("KMinus", Particle::KMinus)
        .value("Neutron", Particle::Neutron)
        .value("PPlus", Particle::PPlus)
        .value("PMinus", Particle::PMinus)
        .value("K0_Short", Particle::K0_Short)
        .value("Eta", Particle::Eta)
        .value("Lambda", Particle::Lambda)
        .value("SigmaPlus", Particle::SigmaPlus)
        .value("Sigma0", Particle::Sigma0)
        .value("SigmaMinus", Particle::SigmaMinus)
        .value("Xi0", Particle::Xi0)
        .value("XiMinus", Particle::XiMinus)
        .value("OmegaMinus", Particle::OmegaMinus)
        .value("NeutronBar", Particle::NeutronBar)
        .value("LambdaBar", Particle::LambdaBar)
        .value("SigmaMinusBar", Particle::SigmaMinusBar)
        .value("Sigma0Bar", Particle::Sigma0Bar)
        .value("SigmaPlusBar", Particle::SigmaPlusBar)
        .value("Xi0Bar", Particle::Xi0Bar)
        .value("XiPlusBar", Particle::XiPlusBar)
        .value("OmegaPlusBar", Particle::OmegaPlusBar)
        .value("DPlus", Particle::DPlus)
        .value("DMinus", Particle::DMinus)
        .value("D0", Particle::D0)
        .value("D0Bar", Particle::D0Bar)
        .value("DsPlus", Particle::DsPlus)
        .value("DsMinusBar", Particle::DsMinusBar)
        .value("LambdacPlus", Particle::LambdacPlus)
        .value("WPlus", Particle::WPlus)
        .value("WMinus", Particle::WMinus)
        .value("Z0", Particle::Z0)
        .value("NuE", Particle::NuE)
        .value("NuEBar", Particle::NuEBar)
        .value("NuMu", Particle::NuMu)
        .value("NuMuBar", Particle::NuMuBar)
        .value("TauPlus", Particle::TauPlus)
        .value("TauMinus", Particle::TauMinus)
        .value("NuTau", Particle::NuTau)
        .value("NuTauBar", Particle::NuTauBar)
        .value("H2Nucleus", Particle::H2Nucleus)
        .value("He3Nucleus", Particle::He3Nucleus)
        .value("He4Nucleus", Particle::He4Nucleus)
        .value("Li6Nucleus", Particle::Li6Nucleus)
        .value("Li7Nucleus", Particle::Li7Nucleus)
        .value("Be9Nucleus", Particle::Be9Nucleus)
        .value("B10Nucleus", Particle::B10Nucleus)
        .value("B11Nucleus", Particle::B11Nucleus)
        .value("C12Nucleus", Particle::C12Nucleus)
        .value("C13Nucleus", Particle::C13Nucleus)
        .value("N14Nucleus", Particle::N14Nucleus)
        .value("N15Nucleus", Particle::N15Nucleus)
        .value("O16Nucleus", Particle::O16Nucleus)
        .value("O17Nucleus", Particle::O17Nucleus)
        .value("O18Nucleus", Particle::O18Nucleus)
        .value("F19Nucleus", Particle::F19Nucleus)
        .value("Ne20Nucleus", Particle::Ne20Nucleus)
        .value("Ne21Nucleus", Particle::Ne21Nucleus)
        .value("Ne22Nucleus", Particle::Ne22Nucleus)
        .value("Na23Nucleus", Particle::Na23Nucleus)
        .value("Mg24Nucleus", Particle::Mg24Nucleus)
        .value("Mg25Nucleus", Particle::Mg25Nucleus)
        .value("Mg26Nucleus", Particle::Mg26Nucleus)
        .value("Al26Nucleus", Particle::Al26Nucleus)
        .value("Al27Nucleus", Particle::Al27Nucleus)
        .value("Si28Nucleus", Particle::Si28Nucleus)
        .value("Si29Nucleus", Particle::Si29Nucleus)
        .value("Si30Nucleus", Particle::Si30Nucleus)
        .value("Si31Nucleus", Particle::Si31Nucleus)
        .value("Si32Nucleus", Particle::Si32Nucleus)
        .value("P31Nucleus", Particle::P31Nucleus)
        .value("P32Nucleus", Particle::P32Nucleus)
        .value("P33Nucleus", Particle::P33Nucleus)
        .value("S32Nucleus", Particle::S32Nucleus)
        .value("S33Nucleus", Particle::S33Nucleus)
        .value("S34Nucleus", Particle::S34Nucleus)
        .value("S35Nucleus", Particle::S35Nucleus)
        .value("S36Nucleus", Particle::S36Nucleus)
        .value("Cl35Nucleus", Particle::Cl35Nucleus)
        .value("Cl36Nucleus", Particle::Cl36Nucleus)
        .value("Cl37Nucleus", Particle::Cl37Nucleus)
        .value("Ar36Nucleus", Particle::Ar36Nucleus)
        .value("Ar37Nucleus", Particle::Ar37Nucleus)
        .value("Ar38Nucleus", Particle::Ar38Nucleus)
        .value("Ar39Nucleus", Particle::Ar39Nucleus)
        .value("Ar40Nucleus", Particle::Ar40Nucleus)
        .value("Ar41Nucleus", Particle::Ar41Nucleus)
        .value("Ar42Nucleus", Particle::Ar42Nucleus)
        .value("K39Nucleus", Particle::K39Nucleus)
        .value("K40Nucleus", Particle::K40Nucleus)
        .value("K41Nucleus", Particle::K41Nucleus)
        .value("Ca40Nucleus", Particle::Ca40Nucleus)
        .value("Ca41Nucleus", Particle::Ca41Nucleus)
        .value("Ca42Nucleus", Particle::Ca42Nucleus)
        .value("Ca43Nucleus", Particle::Ca43Nucleus)
        .value("Ca44Nucleus", Particle::Ca44Nucleus)
        .value("Ca45Nucleus", Particle::Ca45Nucleus)
        .value("Ca46Nucleus", Particle::Ca46Nucleus)
        .value("Ca47Nucleus", Particle::Ca47Nucleus)
        .value("Ca48Nucleus", Particle::Ca48Nucleus)
        .value("Sc44Nucleus", Particle::Sc44Nucleus)
        .value("Sc45Nucleus", Particle::Sc45Nucleus)
        .value("Sc46Nucleus", Particle::Sc46Nucleus)
        .value("Sc47Nucleus", Particle::Sc47Nucleus)
        .value("Sc48Nucleus", Particle::Sc48Nucleus)
        .value("Ti44Nucleus", Particle::Ti44Nucleus)
        .value("Ti45Nucleus", Particle::Ti45Nucleus)
        .value("Ti46Nucleus", Particle::Ti46Nucleus)
        .value("Ti47Nucleus", Particle::Ti47Nucleus)
        .value("Ti48Nucleus", Particle::Ti48Nucleus)
        .value("Ti49Nucleus", Particle::Ti49Nucleus)
        .value("Ti50Nucleus", Particle::Ti50Nucleus)
        .value("V48Nucleus", Particle::V48Nucleus)
        .value("V49Nucleus", Particle::V49Nucleus)
        .value("V50Nucleus", Particle::V50Nucleus)
        .value("V51Nucleus", Particle::V51Nucleus)
        .value("Cr50Nucleus", Particle::Cr50Nucleus)
        .value("Cr51Nucleus", Particle::Cr51Nucleus)
        .value("Cr52Nucleus", Particle::Cr52Nucleus)
        .value("Cr53Nucleus", Particle::Cr53Nucleus)
        .value("Cr54Nucleus", Particle::Cr54Nucleus)
        .value("Mn52Nucleus", Particle::Mn52Nucleus)
        .value("Mn53Nucleus", Particle::Mn53Nucleus)
        .value("Mn54Nucleus", Particle::Mn54Nucleus)
        .value("Mn55Nucleus", Particle::Mn55Nucleus)
        .value("Fe54Nucleus", Particle::Fe54Nucleus)
        .value("Fe55Nucleus", Particle::Fe55Nucleus)
        .value("Fe56Nucleus", Particle::Fe56Nucleus)
        .value("Fe57Nucleus", Particle::Fe57Nucleus)
        .value("Fe58Nucleus", Particle::Fe58Nucleus)
        .value("Qball", Particle::Qball)
        .value("CherenkovPhoton", Particle::CherenkovPhoton)
        .value("Nu", Particle::Nu)
        .value("Monopole", Particle::Monopole)
        .value("Brems", Particle::Brems)
        .value("DeltaE", Particle::DeltaE)
        .value("PairProd", Particle::PairProd)
        .value("NuclInt", Particle::NuclInt)
        .value("MuPair", Particle::MuPair)
        .value("Hadrons", Particle::Hadrons)
        .value("ContinuousEnergyLoss", Particle::ContinuousEnergyLoss)
        .value("FiberLaser", Particle::FiberLaser)
        .value("N2Laser", Particle::N2Laser)
        .value("YAGLaser", Particle::YAGLaser)
        .value("STauPlus", Particle::STauPlus)
        .value("STauMinus", Particle::STauMinus)
        .value("SMPPlus", Particle::SMPPlus)
        .value("SMPMinus", Particle::SMPMinus)
        .export_values();

    // Particle-related functions
    m.def("isLepton", &isLepton);
    m.def("isCharged", &isCharged);
    m.def("particleName", &particleName);
    m.def("particleMass", &particleMass);
    m.def("kineticEnergy", &kineticEnergy);
    m.def("particleSpeed", &particleSpeed);
    m.def("decideShape", &decideShape);
    m.def("deduceInitialType", &deduceInitialType);
    m.def("getInteraction", &getInteraction);

    // Injector
    py::class_<Injector, std::shared_ptr<Injector>>(m, "Injector")
        .def(py::init<unsigned int, Particle::ParticleType, Particle::ParticleType, std::string, std::string, bool>(),
             "NEvents"_a, "FinalType1"_a, "FinalType2"_a,
             "DoublyDifferentialCrossSectionFile"_a, "TotalCrossSectionFile"_a, "Ranged"_a)
        .def_readwrite("events", &Injector::events)
        .def_readwrite("finalType1", &Injector::finalType1)
        .def_readwrite("finalType2", &Injector::finalType2)
        .def_readwrite("crossSectionPath", &Injector::crossSectionPath)
        .def_readwrite("totalCrossSectionPath", &Injector::totalCrossSectionPath)
        .def_readwrite("ranged", &Injector::ranged);

    // Constants
    py::class_<LIConstants> constants(m, "Constants");
    constants.attr("pi") = Constants::pi;
    constants.attr("tau") = Constants::tau;
    constants.attr("degrees") = Constants::degrees;
    constants.attr("deg") = Constants::deg;
    constants.attr("radian") = Constants::radian;
    constants.attr("m") = Constants::m;
    constants.attr("meter") = Constants::meter;
    constants.attr("cm") = Constants::cm;
    constants.attr("centimeter") = Constants::centimeter;
    constants.attr("second") = Constants::second;
    constants.attr("c") = Constants::c;
    constants.attr("protonMass") = Constants::protonMass;
    constants.attr("neutronMass") = Constants::neutronMass;
    constants.attr("isoscalarMass") = Constants::isoscalarMass;
    constants.attr("electronMass") = Constants::electronMass;
    constants.attr("muonMass") = Constants::muonMass;
    constants.attr("tauMass") = Constants::tauMass;
    constants.attr("s") = Constants::s;
    constants.attr("tauLifeTime") = Constants::tauLifeTime;
    constants.attr("MuonLifeTime") = Constants::MuonLifeTime;
    constants.attr("wMass") = Constants::wMass;
    constants.attr("wWidth") = Constants::wWidth;
    constants.attr("zMass") = Constants::zMass;
    constants.attr("WBranchE") = Constants::WBranchE;
    constants.attr("WBranchMuon") = Constants::WBranchMuon;
    constants.attr("WBranchTau") = Constants::WBranchTau;
    constants.attr("WBranchHadronic") = Constants::WBranchHadronic;
    constants.attr("nuEMass") = Constants::nuEMass;
    constants.attr("nuMuMass") = Constants::nuMuMass;
    constants.attr("nuTauMass") = Constants::nuTauMass;
    constants.attr("GeV") = Constants::GeV;
    constants.attr("EeV") = Constants::EeV;
    constants.attr("PeV") = Constants::PeV;
    constants.attr("TeV") = Constants::TeV;
    constants.attr("MeV") = Constants::MeV;
    constants.attr("keV") = Constants::keV;
    constants.attr("eV") = Constants::eV;
    constants.attr("Joule") = Constants::Joule;
    constants.attr("FermiConstant") = Constants::FermiConstant;
    constants.attr("avogadro") = Constants::avogadro;
    constants.attr("thetaWeinberg") = Constants::thetaWeinberg;
    constants.attr("gravConstant") = Constants::gravConstant;
    constants.attr("fineStructure") = Constants::fineStructure;
}
