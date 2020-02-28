#include "LeptonInjector.h"
#include "EventProps.h"

#include <cassert>
#include <fstream>

// namespace constants = boost::math::constants;

namespace LeptonInjector{

	MinimalInjectionConfiguration::MinimalInjectionConfiguration():
		events(1), 
		finalType1(Particle::EMinus),
		finalType2(Particle::Hadrons),
		crossSectionPath(""),
		totalCrossSectionPath(""),
		ranged(true){}

	LeptonInjectorBase::LeptonInjectorBase(){
		eventsGenerated = 0;
		wroteConfigFrame = false;
		suspendOnCompletion = true;
		// also do nothing with this...
		//config = BasicInjectionConfiguration();
	}

	
	void LeptonInjectorBase::Configure(MinimalInjectionConfiguration basic){// , std::shared_ptr<LI_random> pass){
		

		//std::cout<<"configured the random thing" << std::endl;
		config.events = basic.events;
		config.finalType1 = basic.finalType1;
		config.finalType2 = basic.finalType2;
		
		//std::cout << "wrote more" << std::endl;
		
		if(this->config.events==0){
			std::cout<< "there's no point in running this if you don't generate at least one event" <<std::endl;
			throw std::runtime_error("");
		}
		if(this->config.energyMinimum<=0){
			std::cout << ": minimum energy must be positive" << std::endl;
			throw std::runtime_error("");
		}
		if(this->config.energyMaximum<=0){
			std::cout<<": maximum energy must be positive" << std::endl;
			throw std::runtime_error("");
		}
		if(this->config.energyMaximum<this->config.energyMinimum){
			std::cout << ": maximum energy must be greater than or equal to minimum energy" << std::endl;
			throw std::runtime_error("");
		}
		if(this->config.azimuthMinimum<0.0){
			std::cout << ": minimum azimuth angle must be greater than or equal to zero" << std::endl;
			throw std::runtime_error("");
		}
		if(this->config.azimuthMaximum>2*Constants::pi){
			std::cout << ": maximum azimuth angle must be less than or equal to 2 pi" <<std::endl;
			throw std::runtime_error("");
		}
		if(this->config.azimuthMinimum>this->config.azimuthMaximum){
			std::cout<<": minimum azimuth angle must be less than or equal to maximum azimuth angle" <<std::endl;
			throw std::runtime_error("");
		}
		if(this->config.zenithMinimum<0.0){
			std::cout<<": minimum zenith angle must be greater than or equal to zero" << std::endl;
			throw std::runtime_error("");
		}
		if(this->config.zenithMaximum>Constants::pi){
			std::cout<< ": maximum zenith angle must be less than or equal to pi" <<std::endl;
			throw std::runtime_error("");
		}
		if(this->config.zenithMinimum>this->config.zenithMaximum){
			std::cout<<": minimum zenith angle must be less than or equal to maximum zenith angle"<<std::endl;
			throw std::runtime_error("");
		}
		try{
			Particle::ParticleType initialType=deduceInitialType(this->config.finalType1,this->config.finalType2);
		}catch(std::runtime_error& re){
			std::cout<<"Something bad happened while deducing the Initial particle type" << std::endl;
			throw std::runtime_error("");
		}
		this->initialType = initialType;
		// write the pointer to the RNG
		if(!random){
			std::cout<<"A random service is required"<<std::endl;
			throw std::runtime_error("");
		}
		if(basic.crossSectionPath.empty()){
			std::cout<<": DoublyDifferentialCrossSectionFile must be specified"<<std::endl;
			throw std::runtime_error("");
		}
		else if(basic.totalCrossSectionPath.empty()){
			std::cout<<": TotalCrossSectionFile must be specified"<<std::endl;
			throw std::runtime_error("");
		}else{
			crossSection.load(basic.crossSectionPath,basic.totalCrossSectionPath);
		}
		this->crossSection.insert_blobs( this->config );
	}

	BasicInjectionConfiguration& LeptonInjectorBase::getConfig(void){
		return( this->config );
	}

	

	void LeptonInjectorBase::Print_Configuration(){
		std::cout << "    Events   : " << this->config.events << std::endl;
		std::cout << "    Emin     : " << this->config.energyMinimum << std::endl;
		std::cout << "    Emax     : " << this->config.energyMaximum << std::endl;
		std::cout << "    PowerLaw : " << this->config.powerlawIndex << std::endl;
		std::cout << "    AzMin    : " << this->config.azimuthMinimum << std::endl;
		std::cout << "    AzMax    : " << this->config.azimuthMaximum << std::endl;
		std::cout << "    ZenMin   : " << this->config.zenithMinimum << std::endl;
		std::cout << "    ZenMax   : " << this->config.zenithMaximum << std::endl;
		std::cout << "    FT1      : " << static_cast<int32_t>(this->config.finalType1) << std::endl;
		std::cout << "    FT2      : " << static_cast<int32_t>(this->config.finalType2) << std::endl;
	}

	
	LI_Position LeptonInjectorBase::SampleFromDisk(double radius, double zenith, double azimuth){
		//choose a random point on a disk laying in the xy plane
		double t=this->random->Uniform(0,2*Constants::pi);
		double u=this->random->Uniform()+this->random->Uniform();
		double r=(u>1.?2.-u:u)*radius;
		LI_Position  pos(r*cos(t) ,r*sin(t), 0.0);
		//now rotate to make the disc perpendicular to the requested normal vector
		pos = RotateY(pos, zenith );
		pos = RotateZ(pos, azimuth);
		return(pos);
	}
	
	double LeptonInjectorBase::SampleEnergy(){
		if(config.energyMinimum==config.energyMaximum)
			return(config.energyMinimum); //return the only allowed energy
			
		if(config.powerlawIndex==1.0) //sample uniformly in log space
			return(pow(10.0,this->random->Uniform(log10(config.energyMinimum),log10(config.energyMaximum))));
		else{
			double u=this->random->Uniform();
			double energyP=(1-u)*pow(config.energyMinimum,1-config.powerlawIndex) + u*pow(config.energyMaximum,1-config.powerlawIndex);
			return(pow(energyP,1/(1-config.powerlawIndex)));
		}
	}
	
    // this function returns a pair of angles
	std::pair<double,double> LeptonInjectorBase::computeFinalStateAngles(double E_total, double x, double y){
		const double M_N = crossSection.GetTargetMass();
		double theta1=0, theta2=0;
		
		//first particle is a lepton, which covers CC, NC, and leptonic GR
		if(isLepton(config.finalType1)){
			double m1=Particle(config.finalType1).GetMass();
			double E1 = (1 - y) * E_total;
			double cos_theta1, kE1;
			
			if(!isLepton(config.finalType2)){ //CC and NC have Hadrons as particle 2
				//squared kinetic energy of final state particle 1:
				double kE1sq=E1*E1 - m1*m1;
				if(kE1sq<=0){
					std::cout<< "Negative kinetic energy detected in `computeFinalStateAngles` " <<std::endl;
					throw std::runtime_error("");
                }
				cos_theta1=(E1 - x*y*M_N - m1*m1/(2*E_total))/sqrt(kE1sq);
				kE1=sqrt(kE1sq);
			}
			else{ //leptonic GR
				double m_e = Constants::electronMass;
				
				if(E1<=0){
					std::cout << "Bjorken Y > 1?" << std::endl;
					throw std::runtime_error("");
				}
                
				
				cos_theta1=1 - (m_e*m_e + 2*m_e*E_total - m1*m1)/(2*E_total*E1);
				kE1=E1;
			}
			
			if(cos_theta1<-1){
                // commented out until new logger is implemented
//				log_warn_stream("cos(theta) underflow (" << cos_theta1 << "); rounding up to -1"
//					"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y << ")");
				cos_theta1=-1;
			}
			else if(cos_theta1>1){
				//tell the user if the difference was large enough to plausibly not be just round-off
                
                // need new logger 
//				if((cos_theta1-1)>1e-3)
//					log_warn_stream("cos(theta) overflow (" << cos_theta1 << "); rounding down to 1"
//						"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y <<")");
				cos_theta1=1;
			}
			
			theta1=acos(cos_theta1);
			
			//longitudinal component of final state particle 2 momentum:
			double p_long=E_total-kE1*cos_theta1;
			//transverse component of final state particle 2 momentum:
			double p_trans=kE1*sin(theta1);
			theta2=atan(p_trans/p_long);
		}
		//otherwise we have hadronic GR, so both final state masses are unknown
		//and there isn't much we can do, so leave everything colinear
		// angle of particle 1 from initial dir
		return(std::make_pair(theta1,theta2));
	}
	
	// take a direction, deflect that direction by a distance /zenith/
	//		rotate the new direction around the initial direction by /azimuth/
	// So the zenith and azimuth are only what their names would suggest in the coordinate system where 
	//		/base/ is the \hat{z} axis 
	
	void LeptonInjectorBase::FillTree(LI_Position vertex, LI_Direction dir, double energy, BasicEventProperties& properties, std::array<h5Particle,3>& particle_tree){
		const LICrossSection::finalStateRecord& fs=crossSection.sampleFinalState(energy,config.finalType1,this->random);
		
		std::pair<double,double> relativeZeniths=computeFinalStateAngles(energy,fs.x,fs.y);
		double azimuth1=this->random->Uniform(0,2*Constants::pi);
		double azimuth2=azimuth1+(azimuth1<Constants::pi ? 1 : -1)*Constants::pi ;
		
		(particle_tree)[0]=  h5Particle(true,
					static_cast<int32_t>(this->initialType),
					vertex,
					dir,
					energy
		);

		//Make the first final state particle
		(particle_tree)[1] = h5Particle( false, 
					static_cast<int32_t>(config.finalType1),
					vertex,
					rotateRelative(dir,relativeZeniths.first,azimuth1),
					kineticEnergy(config.finalType1,(1-fs.y)*energy)
		);

		(particle_tree)[2] = h5Particle(false,
					static_cast<int32_t>(config.finalType2),
					vertex,
					rotateRelative(dir,relativeZeniths.second,azimuth2),
					kineticEnergy(config.finalType2,fs.y*energy)
		);

		properties.totalEnergy=energy;
		properties.zenith=dir.zenith;
		properties.azimuth=dir.azimuth;
		properties.finalStateX=fs.x;
		properties.finalStateY=fs.y;
		properties.finalType1= static_cast<int32_t>(config.finalType1);
		properties.finalType2= static_cast<int32_t>(config.finalType2);
		properties.initialType=static_cast<int32_t>(deduceInitialType(config.finalType1, config.finalType2));
		properties.interaction=getInteraction(config.finalType1, config.finalType2);

	}
	
	//-----------------------
	//Ranged injection module
	
	RangedLeptonInjector::RangedLeptonInjector(){
	}
	
	RangedLeptonInjector::RangedLeptonInjector( BasicInjectionConfiguration config_, std::shared_ptr<earthmodel::EarthModelService> earth_, std::shared_ptr<LI_random> random_){
		config = config_;
		eventsGenerated = 0;
		wroteConfigFrame = false;
		suspendOnCompletion = true;
		random = random_;
		earthModel = earth_;
		if(config.injectionRadius<0){
			std::cout << ": InjectionRadius must be non-negative" <<std::endl;
			throw std::runtime_error("");
		}
		if(config.endcapLength<0){
			std::cout << ": EndcapLength must be non-negative" <<std::endl;
			throw std::runtime_error("");
		}
		if(!earthModel){
			std::cout<< ": an Earth model service is required" << std::endl;
			throw std::runtime_error("");
		}
	}
	

	
	
	bool RangedLeptonInjector::Generate(){
		
		//Choose an energy
		double energy=SampleEnergy();
		
		//Pick a direction on the sphere
		LI_Direction dir(acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
						random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
		//log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
		
		//decide the point of closest approach
		LI_Position pca=SampleFromDisk(config.injectionRadius,dir.zenith,dir.azimuth);
		//log_trace_stream("pca=(" << pca.GetX() << ',' << pca.GetY() << ',' << pca.GetZ() << ')');
		
		//Figure out where we want the vertex
		//Add up the column depth for the range of a muon at this energy with the
		//column depth for the fixed endcaps to ensure that the whole detector is
		//covered
		using namespace earthmodel::EarthModelCalculator;
		bool isTau = (this->config.finalType1==Particle::ParticleType::TauMinus) || (this->config.finalType1==Particle::ParticleType::TauPlus); 

		bool use_electron_density = getInteraction(this->config.finalType1, this->config.finalType2 ) == 2;
		double totalColumnDepth=MWEtoColumnDepthCGS(GetLeptonRange(energy, isTau=isTau))
		+earthModel->GetColumnDepthInCGS(pca-config.endcapLength*dir,pca+config.endcapLength*dir, use_electron_density);
		//See whether that much column depth actually exists along the chosen path
		{
			double maxDist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth, use_electron_density)-config.endcapLength;
			double actualColumnDepth=earthModel->GetColumnDepthInCGS(pca+config.endcapLength*dir,pca-maxDist*dir, use_electron_density);
			if(actualColumnDepth<(totalColumnDepth-1)){ //if actually smaller, clip as needed, but for tiny differences we don't care
				//log_debug_stream("Wanted column depth of " << totalColumnDepth << " but found only " << actualColumnDepth << " g/cm^2");
				totalColumnDepth=actualColumnDepth;
			}
		}
		//Choose how much of the total column depth this event should have to traverse
		double traversedColumnDepth=totalColumnDepth*random->Uniform();
		//endcapLength is subtracted so that dist==0 corresponds to pca
		double dist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth-traversedColumnDepth, use_electron_density)-config.endcapLength;
		
		{ //ensure that the point we picked is inside the atmosphere
			LI_Position atmoEntry, atmoExit;
			int isect=GetIntersectionsWithSphere(earthModel->GetEarthCoordPosFromDetCoordPos(pca),
												 earthModel->GetEarthCoordDirFromDetCoordDir(dir),
												 earthModel->GetAtmoRadius(),atmoEntry,atmoExit);
			if(isect<2)
				throw("PCA not inside atmosphere");
			atmoEntry=earthModel->GetDetCoordPosFromEarthCoordPos(atmoEntry);
			double atmoDist=(pca-atmoEntry).Magnitude();
			if(std::abs(dist-atmoDist)<100.0)
				dist=std::min(dist,atmoDist);
		}
		LI_Position vertex=pca-(dist*dir);
		
		//assemble the MCTree

		RangedEventProperties properties;
		std::array<h5Particle, 3> particle_tree;

		properties.impactParameter=(pca-LI_Position(0,0,0)).Magnitude();
		properties.totalColumnDepth=totalColumnDepth;
		FillTree(vertex,dir,energy, properties, particle_tree);

		//set subclass properties
		

		//update event count and check for completion
		eventsGenerated++;
		writer_link->WriteEvent( properties, particle_tree[0] , particle_tree[1], particle_tree[2]);

		return(  (eventsGenerated < this->config.events)  );

	}
	
	
	//-----------------------
	//Volume injection module
	
	VolumeLeptonInjector::VolumeLeptonInjector(){
	}
	
	VolumeLeptonInjector::VolumeLeptonInjector(BasicInjectionConfiguration config_, std::shared_ptr<LI_random> random_){
		eventsGenerated = 0;
		wroteConfigFrame = false;
		suspendOnCompletion = true;
		random = random_;
		config = config_;
		if(config.cylinderRadius<0){
			std::cout<< ": CylinderRadius must be non-negative" << std::endl;
			throw std::runtime_error("");
		}
		if(config.cylinderHeight<0){
			std::cout << ": CylinderHeight must be non-negative" << std::endl;
			throw std::runtime_error("");
		}
	}
	
	
	
	bool VolumeLeptonInjector::Generate(){

		//Choose an energy
		double energy=SampleEnergy();
		
		//Pick a direction on the sphere
		LI_Direction dir(acos(this->random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
						this->random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
		//log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
		
		//Pick a position in the xy-plane
		LI_Position vertex=SampleFromDisk(config.cylinderRadius);
		//Add on the vertical component
		vertex.SetZ(this->random->Uniform(-config.cylinderHeight/2,config.cylinderHeight/2));
		//log_trace_stream("vtx=(" << vertex.GetX() << ',' << vertex.GetY() << ',' << vertex.GetZ() << ')');
		
		//assemble the MCTree
		VolumeEventProperties properties;
		std::array<h5Particle, 3> particle_tree;
		
		//set subclass properties
		properties.radius=vertex.Magnitude();
		properties.z=vertex.GetZ();

		// write hdf5 file! 

		FillTree(vertex,dir,energy, properties, particle_tree);
		
		


		//package up output and send it
		writer_link->WriteEvent( properties, particle_tree[0] , particle_tree[1], particle_tree[2]);
		//update event count and check for completion
		eventsGenerated++;

		return(  (eventsGenerated < config.events)  );

	}
	
	bool operator == (const MinimalInjectionConfiguration& one , const MinimalInjectionConfiguration& two){
		return( one.crossSectionPath == two.crossSectionPath
			&& one.events == two.events
			&& one.finalType1 == two.finalType2
			&& one.finalType2 == two.finalType2
			&& one.totalCrossSectionPath == two.totalCrossSectionPath 
			&& one.ranged == two.ranged 
		);
	}

	
	
} //namespace LeptonInjector
