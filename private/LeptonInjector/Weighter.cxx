#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/Weighter.h"

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"

#include <rk/rk.hh>

namespace LeptonInjector {

//---------------
// class LeptonWeighter
//---------------

namespace {
    template <class InIt>
    typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
        typedef typename std::iterator_traits<InIt>::value_type real;
        real sum = real(0);
        real running_error = real(0);
        real temp;
        real difference;

        for (; begin != end; ++begin) {
            difference = *begin;
            difference -= running_error;
            temp = sum;
            temp += difference;
            running_error = temp;
            running_error -= sum;
            running_error -= difference;
            sum = std::move(temp);
        }
        return sum;
    }

    template<typename T>
    T accumulate(std::initializer_list<T> list) {
        return accumulate(list.begin(), list.end());
    }

    double one_minus_exp_of_negative(double x) {
        if(x < 1e-1) {
            return std::exp(std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0);
        } else {
            return 1.0 - std::exp(-x);
        }
    }
}


double LeptonWeighter::InteractionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const {
    std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
    return InteractionProbability(bounds, record);
}

double LeptonWeighter::InteractionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<LeptonInjector::Particle::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            total_xs += xs->TotalCrossSection(record);
        }
        total_cross_sections.push_back(total_xs);
    }
    std::vector<double> particle_counts = earth_model->GetParticleColumnDepth(intersections, bounds.first, bounds.second, targets);
    std::vector<double> interaction_probabilities;
    interaction_probabilities.reserve(particle_counts.size());
    for(unsigned int i=0; i<particle_counts.size(); ++i) {
        interaction_probabilities.push_back(particle_counts[i] * total_cross_sections[i]);
    }
    double exponent = accumulate(interaction_probabilities.begin(), interaction_probabilities.end());
    return one_minus_exp_of_negative(exponent);
}

double LeptonWeighter::UnnormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LeptonInjector::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            total_xs += xs->TotalCrossSection(record);
        }
        total_cross_sections.push_back(total_xs);
    }

    std::vector<double> particle_counts = earth_model->GetParticleColumnDepth(intersections, bounds.first, interaction_vertex, targets);
    std::vector<double> interaction_probabilities; interaction_probabilities.reserve(n_targets);
    for(unsigned int i=0; i<n_targets; ++i) {
        interaction_probabilities.push_back(particle_counts[i] * total_cross_sections[i]);
    }
    double exponent = accumulate(interaction_probabilities.begin(), interaction_probabilities.end());
    double density = std::exp(-exponent);

    std::vector<double> particle_densities = earth_model->GetParticleDensity(intersections, interaction_vertex, targets.begin(), targets.end());
    std::vector<double> jacobians; jacobians.reserve(n_targets);
    for(unsigned int i=0; i<n_targets; ++ i) {
        jacobians.push_back(particle_densities[i] * total_cross_sections[i]);
    }
    double jacobian = accumulate(jacobians.begin(), jacobians.end()) * 100; // cm^-1 --> m^-1

    return density * jacobian;
}

double LeptonWeighter::NormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    double norm = InteractionProbability(bounds, record);
    return UnnormalizedPositionProbability(bounds, record) / norm;
}

void LeptonWeighter::Initialize() {
    
}

LeptonWeighter::LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<CrossSectionCollection> cross_sections, std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions)
    : injectors(injectors)
    , earth_model(earth_model)
    , cross_sections(cross_sections)
    , physical_distributions(physical_distributions)
{}

double LeptonWeighter::EventWeight(InteractionRecord const & record) const {
    // The weight is given by
    //  w = (\sum_i p_gen^i / p_phys^i)^-1

    // The generation probabilities are different between each injector.
    // Most of the physical probabilities are common between all injectors.
    // The physical interaction probability and physical position distribution
    //  depend on the position boundaries of the injection
    //  and thus are different for each injection.
    // Thus the weighting can be given by
    //  w = p_physCommon / (\sum_i p_gen^i / (p_physPos^i * p_physInt^i))

    // However, the normalization of the physical position distribution is identical to the interaction probability.
    // Thus, the two will cancel out and we are left with only the unnormalized position probability
    //  w = p_physCommon / (\sum_i p_gen^i / p_physPosNonNorm^i)

    // The ratio between physical and generation probabilities that differ between injectors
    std::vector<double> gen_over_phys;
    gen_over_phys.reserve(injectors.size());

    // From each injector we need the generation probability and the unnormalized position probability (interaction probability * position probability)
    for(auto injector : injectors) {
        double generation_probability = injector->GenerationProbability(record);
        std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
        double physical_probability = UnnormalizedPositionProbability(bounds, record);

        gen_over_phys.push_back(generation_probability / physical_probability);
    }

    // The denominator is the sum over the ratios for each injector
    double injection_specific_factors = accumulate(gen_over_phys.begin(), gen_over_phys.end());

    // One physical probability density is computed for each distribution, independent of the injectors
    double common_physical_probability = 1.0;
    for(auto physical_distribution : physical_distributions) {
        common_physical_probability *= physical_distribution->GenerationProbability(earth_model, cross_sections, record);
    }

    return common_physical_probability / injection_specific_factors;
}

} // namespace LeptonInjector

