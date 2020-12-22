#include "LeptonInjector/DecayWidths.h"

double gamma_pi0_nu(double G_F, double f, double m_pi0, double m_HNL){
	if(m_HNL < m_pi0) return 0;
	double mass_ratio = pow(m_pi0/m_HNL,2);
	double factor = pow(G_F,2)*pow(f,2)*pow(m_HNL,3)/(32*M_PI);
	return factor*pow(1-mass_ratio,2);
}

double gamma_H_l(double G_F, double V, double f, double m_H, double m_l, double m_HNL){
	if(m_HNL < (m_H+m_l)) return 0;
	double mass_ratio1 = pow(m_l/m_HNL,2);
	double mass_ratio2 = pow(m_H/m_HNL,2);
	double mass_diff_minus = pow((m_H-m_l)/m_HNL,2);
	double mass_diff_plus = pow((m_H+m_l)/m_HNL,2);
	double factor = pow(G_F,2)*pow(V,2)*pow(f,2)*pow(m_HNL,3)/(16*M_PI);
	return factor*(pow(1-mass_ratio1,2)-mass_ratio2*(1+mass_ratio1))*sqrt((1-mass_diff_minus)*(1-mass_diff_minus));
}

double gamma_rho_l(double g_p, double G_F, double V, double m_rho, double m_l, double m_HNL){
	if(m_HNL < (m_rho+m_l)) return 0;
	double mass_ratio1 = pow(m_l/m_HNL,2);
	double mass_ratio2 = pow(m_rho/m_HNL,2);
	double mass_diff_minus = pow((m_rho-m_l)/m_HNL,2);
	double mass_diff_plus = pow((m_rho+m_l)/m_HNL,2);
	double mass_diff_minus_2 = (pow(m_l,2)-2*pow(m_rho,2))/pow(m_HNL,2);
	double factor = pow(g_p/m_rho,2)*pow(G_F,2)*pow(V,2)*pow(m_HNL,3)/(8*M_PI);
	return factor*(pow(1-mass_ratio1,2)-mass_ratio2*(1+mass_diff_minus_2))*sqrt((1-mass_diff_minus)*(1-mass_diff_minus));
}

double gamma_rho0_nu(double g_p, double G_F, double m_rho0, double m_HNL){
	if(m_HNL < m_rho0) return 0;
	double mass_ratio = pow(m_rho0/m_HNL,2);
	double factor = pow(g_p/m_rho0,2)*pow(G_F,2)*pow(m_HNL,3)/(16*M_PI);
	return factor*(1+2*mass_ratio)*pow(1-mass_ratio,2);
}

double gamma_3nu(double G_F, double m_HNL){
	return pow(G_F,2)*pow(m_HNL,5)/(192*pow(M_PI,3));
}

double gamma_nu1_l1_l2(double G_F, double m_l1, double m_l2, double m_HNL){
	if(m_HNL < (m_l1+m_l2)) return 0;
	double factor = pow(G_F,2)*pow(m_HNL,5)/(192*pow(M_PI,3));
	double mass_ratio = std::max(m_l1,m_l2)/m_HNL;
	double x_factor = 1-(8*pow(mass_ratio,2))+(8*pow(mass_ratio,6))-pow(mass_ratio,8)-(12*pow(mass_ratio,4)*log(pow(mass_ratio,2)));
	return factor*x_factor;
}

double gamma_nu_ll(double G_F, double sin2thetaW, double m_l, double m_HNL){
	if(m_HNL < 2*m_l) return 0;
	double x = pow(m_l/m_HNL,2);
	double L = 0;
	double num = 0, den = 0;
	if(x<0.25){
		num = 1-3*x-(1-x)*sqrt(1-4*x);
		den = x*(1+sqrt(1-4*x));
	}
	if(den!=0 && num/den>0) L = log((1-3*x-(1-x)*sqrt(1-4*x))/(x*(1+sqrt(1-4*x))));
	double C1 = 0.25*(1-4*sin2thetaW+8*pow(sin2thetaW,2));
	double C2 = 0.5*sin2thetaW*(2*sin2thetaW-1);
	double C3 = 0.25*(1+4*sin2thetaW+8*pow(sin2thetaW,2));
	double C4 = 0.5*sin2thetaW*(2*sin2thetaW+1);
	double factor = pow(G_F,2)*pow(m_HNL,5)/(192*pow(M_PI,3));
	double gamma_same_flavor = (C3*((1-14*x-2*pow(x,2)-12*pow(x,3))*sqrt(1-4*x)+12*pow(x,2)*(pow(x,2)-1)*L)+(4*C4*(x*(2+10*x-12*pow(x,2))*sqrt(1-4*x)+6*pow(x,2)*(1-2*x+2*pow(x,2))*L)));
	double gamma_opposite_flavor = (C1*((1-14*x-2*pow(x,2)-12*pow(x,3))*sqrt(1-4*x)+12*pow(x,2)*(pow(x,2)-1)*L)+(4*C2*(x*(2+10*x-12*pow(x,2))*sqrt(1-4*x)+6*pow(x,2)*(1-2*x+2*pow(x,2))*L)));
	return factor*(gamma_same_flavor+2*gamma_opposite_flavor);
}
