#include "LeptonInjector/DecayWidths.h"
#include "LeptonInjector/HNLDecayEnergy.h"
#include <time.h>

double TwoBodyRestE1(double m_HNL, double m1, double m2){return (pow(m_HNL,2)+pow(m1,2)-pow(m2,2))/(2*m_HNL);}

double GammaFromE(double E, double m){return E/m;}

double GammaFromBeta(double Beta){return 1./sqrt(1-pow(Beta,2));}

double Beta(double E, double m){return sqrt(1-pow(m/E,2));}

std::vector<double> Lorentz(double E, std::vector<double> p, double boost[3], double gamma){
	double c = 299792458;

	std::vector<double> p_prime;

	double boost_normsq = pow(boost[0],2)+pow(boost[1],2)+pow(boost[2],2);

	p_prime.push_back(gamma*boost[0]*E/c + (1+(gamma-1)*pow(boost[0],2)/boost_normsq)*p.at(0) + ((gamma-1)*boost[0]*boost[1]/boost_normsq)*p.at(1) + ((gamma-1)*boost[0]*boost[2]/boost_normsq)*p.at(2));
	p_prime.push_back(gamma*boost[1]*E/c + ((gamma-1)*boost[1]*boost[0]/boost_normsq)*p.at(0) + (1+(gamma-1)*pow(boost[1],2)/boost_normsq)*p.at(1) + ((gamma-1)*boost[1]*boost[2]/boost_normsq)*p.at(2));
	p_prime.push_back(gamma*boost[2]*E/c + ((gamma-1)*boost[2]*boost[0]/boost_normsq)*p.at(0) + ((gamma-1)*boost[2]*boost[1]/boost_normsq)*p.at(1) + (1+(gamma-1)*pow(boost[2],2)/boost_normsq)*p.at(2));

	return p_prime;
}

std::vector<double> RotateVector(std::vector<double> p, double alpha, int dir){
	std::vector<double> p_rot;

	if(dir==0){
		p_rot.push_back(p.at(0));
		p_rot.push_back(cos(alpha)*p.at(1) - sin(alpha)*p.at(2));
		p_rot.push_back(sin(alpha)*p.at(1) + cos(alpha)*p.at(2));
	}
	else if(dir==2){
		p_rot.push_back(cos(alpha)*p.at(0) - sin(alpha)*p.at(1));
		p_rot.push_back(sin(alpha)*p.at(0) + cos(alpha)*p.at(1));
		p_rot.push_back(p.at(2));
	}

	return p_rot;
}

double FullWidth(double m_HNL){
    // Declare global physical constants
    double G_F = 1.17E-05; // Fermi constant in GeV^-2
    double sin2thetaW = 0.23122; // sin^2(weak angle)

    //return 2*(gamma_3nu(G_F,m_HNL)+                   // 3nu
    return 2*(gamma_nu_ll(G_F,sin2thetaW,511E-06,m_HNL)+    // nu_ee
              gamma_pi0_nu(G_F,0.130,0.135,m_HNL)+            // pi0_nu
              gamma_nu_ll(G_F,sin2thetaW,0.1056,m_HNL)+    // nu_mumu
              gamma_pi0_nu(G_F,0.156,0.548,m_HNL)+            // eta_nu
              gamma_rho0_nu(0.102,G_F,0.775,m_HNL)+      // rho0_nu
              gamma_pi0_nu(G_F,-0.0585,0.958,m_HNL)+          // etaprime_nu
              gamma_nu1_l1_l2(G_F,511E-06,1.776,m_HNL)+      // nue_e_tau
              gamma_nu1_l1_l2(G_F,0.1056,1.776,m_HNL)+      // numu_mu_tau
              gamma_H_l(G_F,0.97427,0.130,0.140,1.776,m_HNL)); // pi_tau
}

int DecayChannel(double rand, double m_HNL){
	// Declare global physical constants
    double G_F = 1.17E-05; // Fermi constant in GeV^-2
    double sin2thetaW = 0.23122; // sin^2(weak angle)

	//double x1 = 2*gamma_3nu(G_F,m_HNL)/FullWidth(m_HNL);
	double x1 = 2*gamma_nu_ll(G_F,sin2thetaW,511E-06,m_HNL)/FullWidth(m_HNL);
	double x2 = 2*gamma_pi0_nu(G_F,0.130,0.135,m_HNL)/FullWidth(m_HNL);
	double x3 = 2*gamma_nu_ll(G_F,sin2thetaW,0.1056,m_HNL)/FullWidth(m_HNL);
	double x4 = 2*gamma_pi0_nu(G_F,0.156,0.548,m_HNL)/FullWidth(m_HNL);
	double x5 = 2*gamma_rho0_nu(0.102,G_F,0.775,m_HNL)/FullWidth(m_HNL);
	double x6 = 2*gamma_pi0_nu(G_F,-0.0585,0.958,m_HNL)/FullWidth(m_HNL);
	double x7 = 2*gamma_nu1_l1_l2(G_F,511E-06,1.776,m_HNL)/FullWidth(m_HNL);
	double x8 = 2*gamma_nu1_l1_l2(G_F,0.1056,1.776,m_HNL)/FullWidth(m_HNL);

	double x[8] = {x1,
                   x1+x2,
                   x1+x2+x3,
                   x1+x2+x3+x4,
                   x1+x2+x3+x4+x5,
                   x1+x2+x3+x4+x5+x6,
                   x1+x2+x3+x4+x5+x6+x7,
                   x1+x2+x3+x4+x5+x6+x7+x8};

	for(int i=0; i<8; i++){
		if(rand<x[i]) return i+1;
	}
	return 9;
}

// Pushes back the following elements to the vector: the pdgid, the energy and the direction of the daughter produced in the HNL 2-body decay 
void TwoBodyDecay(std::vector<double> &decay_products, double pdgId_daughter, double m_daughter, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double zenith_daughter, double azimuth_daughter){
	double c = 299792458;
    // Compute boost
    double v_HNL[3] = {c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*cos(azimuth_HNL),c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*sin(azimuth_HNL),c*Beta(E_HNL,m_HNL)*cos(zenith_HNL)};
    // Compute rest frame variables
	double E_daughter = TwoBodyRestE1(m_HNL,m_daughter,0.0);
	double beta_daughter = Beta(E_daughter,m_daughter);
	double gamma_daughter = GammaFromBeta(beta_daughter);
    double p_daughter_x = beta_daughter*gamma_daughter*m_daughter*sin(zenith_daughter)*cos(azimuth_daughter);
    double p_daughter_y = beta_daughter*gamma_daughter*m_daughter*sin(zenith_daughter)*sin(azimuth_daughter);
    double p_daughter_z = beta_daughter*gamma_daughter*m_daughter*cos(zenith_daughter);
	std::vector<double> p_daughter{p_daughter_x,p_daughter_y,p_daughter_z};
    // Lorentz-transform the momentum components
	std::vector<double> p_daughter_prime;
	p_daughter_prime = Lorentz(E_daughter,p_daughter,v_HNL,GammaFromE(E_HNL,m_HNL));
	// Compute daughter energy and angles
	double E_daughter_prime = sqrt(pow(p_daughter_prime.at(0),2)+pow(p_daughter_prime.at(1),2)+pow(p_daughter_prime.at(2),2)+pow(m_daughter,2));
	double zenith_daughter_prime = acos(p_daughter_prime.at(2)/sqrt(pow(p_daughter_prime.at(0),2)+pow(p_daughter_prime.at(1),2)+pow(p_daughter_prime.at(2),2)));
	double azimuth_daughter_prime = atan(p_daughter_prime.at(1)/p_daughter_prime.at(0));
	if(p_daughter_prime.at(0)>0 && p_daughter_prime.at(1)<0) azimuth_daughter_prime = 2*M_PI + azimuth_daughter_prime;
	else if(p_daughter_prime.at(0)<0 && p_daughter_prime.at(1)>0) azimuth_daughter_prime = M_PI + azimuth_daughter_prime;
	else if(p_daughter_prime.at(0)<0 && p_daughter_prime.at(1)<0) azimuth_daughter_prime = M_PI + azimuth_daughter_prime;
    // Fill vector
	decay_products.push_back(pdgId_daughter);
	decay_products.push_back(m_daughter);
	decay_products.push_back(E_daughter_prime);
	decay_products.push_back(zenith_daughter_prime);
	decay_products.push_back(azimuth_daughter_prime);
}

// Pushes back the following elements to the vector: the pdgid, the energy and the direction of the two visible daughters produced in the HNL 3-body decay 
void ThreeBodyDecay(std::vector<double> &decay_products, double pdgId_1, double pdgId_2, double m_1, double m_2, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_energy_1, double rand_energy_2, double azimuth_p1p2, double azimuth_rotZ, double alpha_rotX){
	double c = 299792458;
    // Compute boost
    double v_HNL[3] = {c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*cos(azimuth_HNL),c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*sin(azimuth_HNL),c*Beta(E_HNL,m_HNL)*cos(zenith_HNL)};
	// Compute energy-momentum of two visible daughters
	double E_1 = 0;
	double E_2 = 0;
	double p_1_norm = 0;
	double p_2_norm = 0;
	double p_3_norm = 0;
	srand (100000*m_HNL);
	double in_sqrt = -1;
	while(in_sqrt < 0){
		rand_energy_1 = ((double) rand() / (RAND_MAX));
		E_1 = m_1 + rand_energy_1*(((pow(m_HNL,2)+pow(m_1,2)-pow(m_2,2))/(2*m_HNL))-m_1);
		double m23_sq = pow(m_HNL,2)-(2*E_1*m_HNL)+pow(m_1,2);
		double lambda = pow(m23_sq,2)+pow(m_2,4)-(2*m23_sq*pow(m_2,2));
		rand_energy_2 = ((double) rand() / (RAND_MAX));
    	E_2 = 1/(2*m23_sq)*((m_HNL-E_1)*(m23_sq+pow(m_2,2))+(rand_energy_2*sqrt((pow(E_1,2)-pow(m_1,2))*lambda)));
		p_1_norm = sqrt(pow(E_1,2)-pow(m_1,2));
		p_2_norm = sqrt(pow(E_2,2)-pow(m_2,2));
		p_3_norm = m_HNL - E_1 - E_2;
		in_sqrt = 4*pow(p_2_norm/p_3_norm,2) - 4*(p_1_norm/p_3_norm) - 3;
		double cos_phi2_minus = (-1-sqrt(in_sqrt))/(2*(p_2_norm/p_3_norm));
		double cos_phi2_plus = (-1+sqrt(in_sqrt))/(2*(p_2_norm/p_3_norm));
		azimuth_p1p2 = acos(cos_phi2_plus);
	}
    // Compute rest frame momenta in Z = 0 plane
	std::vector<double> p_1_Z0{p_1_norm,0,0};
	std::vector<double> p_2_Z0{p_2_norm*cos(azimuth_p1p2),p_2_norm*sin(azimuth_p1p2),0};
	// Rotate around Z axis (XY symmetry)
	std::vector<double> p_1_rotZ, p_2_rotZ;
	p_1_rotZ = RotateVector(p_1_Z0,azimuth_rotZ,2);
	p_2_rotZ = RotateVector(p_2_Z0,azimuth_rotZ,2);
	// Rotate around X axis to get final rest-frame momenta
	std::vector<double> p_1, p_2;
	p_1 = RotateVector(p_1_rotZ,alpha_rotX,0);
	p_2 = RotateVector(p_2_rotZ,alpha_rotX,0);
    // Lorentz-transform the momentum components
	std::vector<double> p_1_prime, p_2_prime;
	p_1_prime = Lorentz(E_1,p_1,v_HNL,GammaFromE(E_HNL,m_HNL));
	p_2_prime = Lorentz(E_2,p_2,v_HNL,GammaFromE(E_HNL,m_HNL));
	// Compute daughter_1 energy and angles
	double E_1_prime = sqrt(pow(p_1_prime.at(0),2)+pow(p_1_prime.at(1),2)+pow(p_1_prime.at(2),2)+pow(m_1,2));
	double zenith_1_prime = acos(p_1_prime.at(2)/sqrt(pow(p_1_prime.at(0),2)+pow(p_1_prime.at(1),2)+pow(p_1_prime.at(2),2)));
	double azimuth_1_prime = atan(p_1_prime.at(1)/p_1_prime.at(0));
	if(p_1_prime.at(0)>0 && p_1_prime.at(1)<0) azimuth_1_prime = 2*M_PI + azimuth_1_prime;
	else if(p_1_prime.at(0)<0 && p_1_prime.at(1)>0) azimuth_1_prime = M_PI + azimuth_1_prime;
	else if(p_1_prime.at(0)<0 && p_1_prime.at(1)<0) azimuth_1_prime = M_PI + azimuth_1_prime;
	// Compute daughter_2 energy and angles
	double E_2_prime = sqrt(pow(p_2_prime.at(0),2)+pow(p_2_prime.at(1),2)+pow(p_2_prime.at(2),2)+pow(m_2,2));
	double zenith_2_prime = acos(p_2_prime.at(2)/sqrt(pow(p_2_prime.at(0),2)+pow(p_2_prime.at(1),2)+pow(p_2_prime.at(2),2)));
	double azimuth_2_prime = atan(p_2_prime.at(1)/p_2_prime.at(0));
	if(p_2_prime.at(0)>0 && p_2_prime.at(1)<0) azimuth_2_prime = 2*M_PI + azimuth_2_prime;
	else if(p_2_prime.at(0)<0 && p_2_prime.at(1)>0) azimuth_2_prime = M_PI + azimuth_2_prime;
	else if(p_2_prime.at(0)<0 && p_2_prime.at(1)<0) azimuth_2_prime = M_PI + azimuth_2_prime;
    // Fill vector
	decay_products.push_back(pdgId_1);
	decay_products.push_back(m_1);
	decay_products.push_back(E_1_prime);
	decay_products.push_back(zenith_1_prime);
	decay_products.push_back(azimuth_1_prime);
	decay_products.push_back(pdgId_2);
	decay_products.push_back(m_2);
	decay_products.push_back(E_2_prime);
	decay_products.push_back(zenith_2_prime);
	decay_products.push_back(azimuth_2_prime);
}

std::vector<double> GetDaughterVector(bool isParticle, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_channel, double rand_energy_1, double rand_energy_2, double azimuth_p1p2, double azimuth_rotZ, double alpha_rotX){
	// Declare vector
	std::vector<double> decay_products;
	// Particle or antiparticle
	int isParticle_int = 1;
	if(!isParticle) isParticle_int = -1;
	// Determine the decay channel
	int channel = DecayChannel(rand_channel,m_HNL);
	if(channel==1) ThreeBodyDecay(decay_products,11,-11,511E-06,511E-06,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,rand_energy_1,rand_energy_2,azimuth_p1p2,azimuth_rotZ,alpha_rotX);
	else if(channel==2) TwoBodyDecay(decay_products,111,0.135,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
	else if(channel==3) ThreeBodyDecay(decay_products,13,-13,0.1056,0.1056,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,rand_energy_1,rand_energy_2,azimuth_p1p2,azimuth_rotZ,alpha_rotX);
	else if(channel==4) TwoBodyDecay(decay_products,221,0.548,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
	else if(channel==5) TwoBodyDecay(decay_products,113,0.775,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
	else if(channel==6) TwoBodyDecay(decay_products,331,0.958,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
	else if(channel==7) ThreeBodyDecay(decay_products,-1*isParticle_int*11,isParticle_int*15,511E-06,1.776,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,rand_energy_1,rand_energy_2,azimuth_p1p2,azimuth_rotZ,alpha_rotX);
	else if(channel==8) ThreeBodyDecay(decay_products,-1*isParticle_int*13,isParticle_int*15,0.1056,1.776,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,rand_energy_1,rand_energy_2,azimuth_p1p2,azimuth_rotZ,alpha_rotX);
	else if(channel==9){
		TwoBodyDecay(decay_products,isParticle_int*211,0.140,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
		TwoBodyDecay(decay_products,isParticle_int*15,1.776,m_HNL,E_HNL,zenith_HNL,azimuth_HNL,azimuth_p1p2/2,azimuth_rotZ);
	}
	// Return the vector
	return decay_products;
}
