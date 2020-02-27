#include <fstream>

#include <LeptonInjector.h>
#include <BasicInjectionConfiguration.h>
#include <Controller.h>
#include <Particle.h>

#include "tools.h"

namespace LeptonInjector{ 

//None of the default values in the *Configuration objects should set off
//errors in the module
TEST(1_sane_defaults_in_config){
	std::cout << "creating minimal config" << std::endl;
	MinimalInjectionConfiguration config;

	Controller inj(config);
}

//attempting to set nonsensical parameter values should be rejected
TEST(2_reject_invalid_params){
	MinimalInjectionConfiguration config;

	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 0, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi);
		FAIL("configuring with a non-positive minimum energy should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 0, 2, 0, 2*Constants::pi, 0, Constants::pi);
		FAIL("configuring with a non-positive maximum energy should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e2, 2, 0, 2*Constants::pi, 0, Constants::pi);
		FAIL("configuring with a maximum energy below the minimum energy should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, -1, 2*Constants::pi, 0, Constants::pi);
		FAIL("configuring with a negative minimum azimuth should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 7, 0, Constants::pi);
		FAIL("configuring with a maximum azimuth greater than 2 pi should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 2., 1., 0, Constants::pi);
		FAIL("configuring with a maximum azimuth less than the minimum azimuth should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, -1, Constants::pi);
		FAIL("configuring with a negative minimum zenith should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, 2*Constants::pi);
		FAIL("configuring with a maximum zenith greater than pi should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 2, 1);
		FAIL("configuring with a maximum zenith less than the minimum zenith should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	
	try{
		MinimalInjectionConfiguration config;
		config.crossSectionPath = "";
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi);
		inj.Configure();
		FAIL("configuring without a cross section should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi,-200);
		FAIL("configuring with a negative injection radius should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi,1200,-200);
		FAIL("configuring with a negative endcap length should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi, 1200,1200,-200);
		FAIL("configuring with a negative cylinder radius should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
	
	try{
		MinimalInjectionConfiguration config;
		Controller inj(config, 1e3, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi,1200,1200,1200,-200);
		FAIL("configuring with a negative cylinder height should be rejected");
	}catch(std::runtime_error& e){/*squash*/}
}

TEST(3_number_of_events){
	const unsigned int nEvents=100;
	const std::string filename="Multi_NEvents_test.i3";
	
	
	MinimalInjectionConfiguration mic1(nEvents,Particle::MuMinus,Particle::Hadrons,defaultCrosssectionPath,defaultTotalCrosssectionPath,true); //CC
	MinimalInjectionConfiguration mic2(nEvents,Particle::NuMu,Particle::Hadrons,defaultCrosssectionPath,defaultTotalCrosssectionPath,false); //NC

	Controller this_test( mic1, 1e2, 1e6, 2, 0, 2*Constants::pi, 0, Constants::pi );
	this_test.AddInjector( mic2 );
	this_test.NameOutfile(filename);
	this_test.NameLicFile("Multi_NEvents_test.lic");
	this_test.Execute();

	std::ifstream resultsFile(filename.c_str());
	ENSURE(resultsFile.good()); //the file should exist and be readable

	//If the test was successful, clean up
	unlink(filename.c_str());
}

} // LeptonInjector Namespace