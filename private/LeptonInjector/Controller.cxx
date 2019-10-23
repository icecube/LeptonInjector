#include <LeptonInjector/Controller.h>

// This is the "controller" of the whole thing
// It is a reimplementation of the multileptoninjector 

namespace LeptonInjector {

    Controller::Controller(){
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = 1200*Constants::m;
        this->endcapLength = 1200*Constants::m;
        this->cylinderRadius = 1200*Constants::m;
        this->cylinderHeight = 1200*Constants::m;

    }

    // constructor if the user only provides one injector and no parameters
    Controller::Controller(MinimalInjectionConfiguration configs_received){
        this->configs.push_back( conconfigs_receivedfigs );
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = 1200*Constants::m;
        this->endcapLength = 1200*Constants::m;
        this->cylinderRadius = 1200*Constants::m;
        this->cylinderHeight = 1200*Constants::m;
    }

    //constructor if the user provides a list of injectors and no parameters 
    Controller::Controller(	std::vector<MinimalInjectionConfiguration> configs_received ){
        this->configs = configs_received;
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = 1200*Constants::m;
        this->endcapLength = 1200*Constants::m;
        this->cylinderRadius = 1200*Constants::m;
        this->cylinderHeight = 1200*Constants::m;
    }

    Controller::Controller(std::vector<MinimalInjectionConfiguration> configs_received, double minimumEnergy
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith,
            double injectionRadius=1200*Constants::m, double endcapLength=1200*Constants::m,
            double cylinderRadius=1200*Constants::m, double cylinderHeight= 1200*Constants::m){

        this->configs = configs_received;
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = injectionRadius;
        this->endcapLength    = endcapLength;
        this->cylinderRadius  = cylinderRadius;
        this->cylinderHeight  = cylinderHeight;

    }

    Controller::~Controller(){
    }

    void Controller::AddInjector( MinimalInjectionConfiguration configs_received ){
        this->configs.push_back( configs_received );
    }

    void Controller::Execute(){
        // setup the injectors! 

        bool hasRanged=false, hasVolume=false;
		for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=this->configs.begin(), end=this->configs.end(); genSet!=end; genSet++){
			hasRanged |= genSet->ranged;
			hasVolume |= !genSet->ranged;
		}

        LI_random random(this->seed);

        // sanity check! 
        if (this->energyMinium <= 0 ){ throw "minimum energy must be positive"; }
        if (this->energyMaximum <= 0 ){ throw "maximum energy must be positive"; }
        if (this->energyMaximum < this->energyMinimum ){ throw "Max energy must be greater or equal to minimum energy"; }
        if (this->azimuthMinimum < 0 ){ throw "minimum azimuth must be positive"; }
        if (this->azimuthMaximum > 2*Constants::pi ){ throw "maximum azimuth must be less than 2pi"; }
        if (this->azimuthMinimum < this->azimuthMaximum ){ throw "Max azimuth must be greater or equal to min."; }
        if (this->zenithMinimum < 0.0 ){ throw "minimum zenith must be positive"; }
        if (this->zenithMaximum > Constants::pi ){ throw "maximum zenith must be less than or equal to pi"; }
        if (this->zenithMinimum >this->zenithMaximum){throw "Max zenith must be greater or equal to min.";}

        // first, construct the template injector configuration objects
        this->rangedConfig.energyMinimum = this->minimumEnergy; 
        this->rangedConfig.energyMaximum = this->maximumEnergy; 
        this->rangedConfig.powerlawIndex = this->powerlawIndex; 
        this->rangedConfig.azimuthMinimum = this->minimumAzimuth; 
        this->rangedConfig.azimuthMaximum = this->maximumAzimuth; 
        this->rangedConfig.zenithMinimum = this->minimumZenith; 
        this->rangedConfig.zenithMaximum = this->maximumZenith; 

        this->volumeConfig.energyMinimum = this->minimumEnergy; 
        this->volumeConfig.energyMaximum = this->maximumEnergy; 
        this->volumeConfig.powerlawIndex = this->powerlawIndex; 
        this->volumeConfig.azimuthMinimum = this->minimumAzimuth; 
        this->volumeConfig.azimuthMaximum = this->maximumAzimuth; 
        this->volumeConfig.zenithMinimum = this->minimumZenith; 
        this->volumeConfig.zenithMaximum = this->maximumZenith; 

        //  SETUP EARTHMODEL

        if(hasRanged){
            this->rangedConfig.injectionRadius = this->injectionRadius; 
            this->rangedConfig.endcapLength = this->endcapLength; 
            // set pointer to earthmodel -- GetParameter("EarthModel",earthModelName);
            
            if(this->rangedConfig.injectionRadius<0){throw": InjectionRadius must be non-negative"; }
            if(this->rangedConfig.endcapLength<0){ throw ": EndcapLength must be non-negative"; }
            /*(
            earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
            
            if(!earthModel)
                log_fatal_stream(GetName() << ": an Earth model service is required");
            
            innerContext.Put(earthModel,earthModelName); */
        }
        
        //get the properties for volume injectors
        if(hasVolume){
            this->rangedConfig.cylinderRadius = this->cylinderRadius; 
            this->rangedConfig.cylinderHeight = this->cylinderHeight;
            
            if(this->rangedConfig.cylinderRadius<0){throw": InjectionRadius must be non-negative"; }
            if(this->rangedConfig.cylinderHeight<0){ throw ": EndcapLength must be non-negative"; }

        }

        //construct all generators
        unsigned int i=0;
        for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=this->configs.begin(), end=this->configs.end(); genSet!=end; genSet++){
//            log_debug_stream("Configuring injector " << i << ":");
            LeptonInjectorBase* generator=NULL;
            try{
                if(genSet->ranged){
                    //log_debug_stream(" this is a ranged injector");
                    generator=new RangedLeptonInjector(this->rangedConfig);
                    generator->GetConfiguration().Set("EarthModel",boost::python::object(earthModelName));
                }
                else{ //volume
                    //log_debug_stream(" this is a volume injector");
                    generator=new VolumeLeptonInjector(this->volumeConfig);
                }
                
                //set properties not shared with other injectors, or which are not part of the config object
                generator->GetConfiguration().Set("NEvents",boost::python::object(genSet->events));
                generator->GetConfiguration().Set("FinalType1",boost::python::object(genSet->finalType1));
                generator->GetConfiguration().Set("FinalType2",boost::python::object(genSet->finalType2));
                generator->GetConfiguration().Set("RandomService",boost::python::object(randomServiceName));
                generator->GetConfiguration().Set("DoublyDifferentialCrossSectionFile",boost::python::object(genSet->crossSectionPath));
                generator->GetConfiguration().Set("TotalCrossSectionFile",boost::python::object(genSet->totalCrossSectionPath));
                generator->GetConfiguration().Set("SuspendOnCompletion",boost::python::object(false));
                
                generator->SetName(GetName()+"_Generator_"+boost::lexical_cast<std::string>(i++));
                generator->Configure();
            }catch(...){
                delete generator;
                throw;
            }
            generators.push_back(generator);
        }
        
        // set up the reference to the first generator
        this->ActiveGenerator = *(generators.front());

    }

	Controller::Generate(void) {

	}


}