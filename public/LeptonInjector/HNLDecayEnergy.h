#ifndef HNLDECAYENERGY_H_INCLUDED    // To make sure you don't declare the function more than once by including the header multiple times.
#define HNLDECAYENERGY_H_INCLUDED

#include <iostream>
#include "math.h"
#include "string.h"
#include <typeinfo> 
#include <vector>

double TwoBodyRestE1(double m_HNL, double m1, double m2);
double GammaFromE(double E, double m);
double GammaFromBeta(double Beta);
double Beta(double E, double m);
std::vector<double> Lorentz(double E, std::vector<double> p, double boost[3], double gamma);
std::vector<double> RotateVector(std::vector<double> p, double alpha, int dir);
double FullWidth(double m_HNL);
int DecayChannel(double rand, double m_HNL);
void TwoBodyDecay(std::vector<double> &decay_products, double pdgId_daughter, double m_daughter, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double zenith_daughter, double azimuth_daughter);
void ThreeBodyDecay(std::vector<double> &decay_products, double pdgId_1, double pdgId_2, double m_1, double m_2, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_energy_1, double rand_energy_2, double azimuth_p1p2, double azimuth_rotZ, double alpha_rotX);
std::vector<double> GetDaughterVector(bool isParticle, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_channel, double rand_energy_1, double rand_energy_2, double azimuth_p1p2, double azimuth_rotZ, double alpha_rotX);

#endif
