#ifndef TOOLS_HEADER_INCLUDED
#define TOOLS_HEADER_INCLUDED

#include <cmath>

#include <LeptonInjector.h>
#include <EarthModelService.h>
#include <Particle.h>

#include <cmath>
#include <map> //used by a test registry
#include <string>
#include <sstream>
#include <memory> // shared_ptr

// ---------------TEST DEFINITIONS BORROWED FROM PHOTOSPLINE-----------------------------------
// https://github.com/IceCubeOpenSource/photospline
// ^^  Published under the BSD 2-Clause "Simplified License"  ^^


void emit_error(const std::string& file, size_t line,
				const std::string& criterion, const std::string& message="");

#define ENSURE(cond,...) \
	do{ \
		if(!(cond)) \
			emit_error(__FILE__,__LINE__,#cond,##__VA_ARGS__); \
	}while(0)

namespace{
	template<typename T1, typename T2>
	std::string express_comparison(const std::string& e1, const T1& v1,
	                               const std::string& e2, const T2& v2){
		std::ostringstream ss;
		ss.precision(16);
		ss << v1 << " (" << e1 << ") != " << v2 << " (" << e2 << ")";
		return(ss.str());
	}
	template<typename T>
	std::string express_comparison(const std::string& e1, const T& v1,
	                               const std::string& e2, const T& v2,
	                               const T& tolerance=0){
		std::ostringstream ss;
		ss.precision(16);
		ss << v1 << " (" << e1 << ") != " << v2 << " (" << e2 << ")";
		if(tolerance!=0)
			ss << " to within " << tolerance;
		return(ss.str());
	}

	bool ensure_equal_impl(float v1, float v2){
		return(v1==v2 || (std::isnan(v1) && std::isnan(v2)));
	}
	bool ensure_equal_impl(double v1, double v2){
		return(v1==v2 || (std::isnan(v1) && std::isnan(v2)));
	}
	template<typename T1, typename T2>
	bool ensure_equal_impl(const T1& v1, const T2& v2){
		return(v1==v2);
	}
}

#define ENSURE_EQUAL(first,second,...) \
	do{ \
		if(!ensure_equal_impl(first,second)) \
			emit_error(__FILE__,__LINE__, \
			  express_comparison(#first,first,#second,second),##__VA_ARGS__); \
	}while(0)

#define ENSURE_DISTANCE(first,second,tolerance,...) \
do{ \
	if(!(std::abs((first)-(second))<(tolerance))) \
		emit_error(__FILE__,__LINE__, \
		  express_comparison(#first,first,#second,second,tolerance),##__VA_ARGS__); \
}while(0)

#define FAIL(...) \
	emit_error(__FILE__,__LINE__,"FAIL",##__VA_ARGS__)

std::map<std::string,void(*)()>& test_registry();

struct register_test{
	register_test(const std::string& test_name, void(*test)()){
		test_registry().insert(std::make_pair(test_name,test));
	}
};

#define TEST(name) \
	void test_func ## name (); \
	static register_test register_ ## name (#name,&test_func ## name); \
	void test_func ## name ()


// --------------------------------------------------------------

extern const std::string earthModelName;
extern const std::string defaultCrosssectionPath;
extern const std::string defaultTotalCrosssectionPath;

extern std::shared_ptr<LeptonInjector::LI_random> random_machine;
extern std::shared_ptr<earthmodel::EarthModelService> earth;
extern std::shared_ptr<LeptonInjector::MinimalInjectionConfiguration> minimal_ranged;
extern std::shared_ptr<LeptonInjector::MinimalInjectionConfiguration> minimal_volume;


//Welford's on-line variance algorithm, extended to higher moments, following
//implmentation by John D. Cook: http://www.johndcook.com/blog/skewness_kurtosis/
class MomentAccumulator{
public:
	MomentAccumulator():
	n(0),m1(0),m2(0),m3(0),m4(0){}
	
	void Insert(double x){
		const double oldN=n++;
		
		const double d = x-m1;
		const double don = d/n;
		const double don2=don*don;
		const double t = d*don*oldN;
		m1+=don;
		m4+=t*don2*(double(n)*n - 3.*n + 3) + 6*don2*m2 - 4*don*m3;
		m3+=(t*(n-2) - 3*m2)*don;
		m2+=t;
	}
	
	unsigned long NumDataValues() const{
		return n;
	}
	
	double Mean() const{
		return(m1);
	}
	
	double Variance() const{
		return((n>1) ? m2/(n-1) : 0.0);
	}
	
	double StandardDeviation() const{
		return(sqrt(Variance()));
	}
	
	double Skewness() const{
		return(sqrt(double(n)) * m3/ pow(m2, 1.5));
	}
    double Kurtosis() const{
		return n*m4 / (m2*m2) - 3.0;
	}
	
private:
	unsigned long n;
	double m1, m2, m3, m4;
};

double predictPowerLawMoment(double index, double a, double b, unsigned int moment);

void testPowerLawness(double powerlawIndex, double min, double max, unsigned int count,
					  const MomentAccumulator& moments, const std::string& file, unsigned line);


#endif 