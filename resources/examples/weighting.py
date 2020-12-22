import LeptonWeighter as LW
import h5py as h5
import numpy as np
import os

"""
A Basic Weighting example 
This calculates the weight of each event in the standalone LeptonInjector example script. 
Ben Smithers
benjamin.smithers@mavs.uta.edu
To convert this to one working with i3-LeptonInjector, you will need to modify the event loop and the get_weight function. 
"""

#                         GeV      unitless        GeV
flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
#liveTime   = 3.1536e7 #s
liveTime   = 1 #s

# Create generator
#    if there were multiple LIC files, you would instead make a list of Generators
net_generation = LW.MakeGeneratorsFromLICFile("config.lic")

# Create flux, load cross sections 
# This cross section object takes four differential cross sections (dS/dEdxdy) 
#            Neutrino CC-DIS xs
#       Anti-Neutrino CC-DIS xs
#            Neutrino NC-DIS xs
#       Anti-Neutrino NC-DIS xs
#cross_section_location = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/"
cross_section_location = os.path.join( os.path.dirname(__file__), '..' )
xs = LW.CrossSectionFromSpline(
                    cross_section_location+"/dsdxdy-hnl-N-cc-HERAPDF15NLO_EIG_central.fits",
                    cross_section_location+"/dsdxdy-hnlbar-N-cc-HERAPDF15NLO_EIG_central.fits",
                    cross_section_location+"/dsdxdy-hnl-N-nc-HERAPDF15NLO_EIG_central.fits",
                    cross_section_location+"/dsdxdy-hnlbar-N-nc-HERAPDF15NLO_EIG_central.fits")
flux = LW.PowerLawFlux( flux_params['constant'] , flux_params['index'] , flux_params['scale'] )
#Flux = "/data/user/dvannerom/HNL/Fluxes/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000_2006_2014.hdf5"
#flux = LW.nuSQUIDSAtmFlux(Flux)

# build weighter
weight_event = LW.Weighter( flux, xs, net_generation )

# Write utility function to get the weight 
def get_weight( props ):
	"""
	Accepts the properties list of an event and returns the weight
	"""
	LWevent = LW.Event()
	LWevent.energy = props[0]
	LWevent.zenith = props[1]
	LWevent.azimuth = props[2]
	
	LWevent.interaction_x = props[3]
	LWevent.interaction_y = props[4]
	LWevent.final_state_particle_0 = LW.ParticleType( props[5] )
	LWevent.final_state_particle_1 = LW.ParticleType( props[6] )
	LWevent.primary_type = LW.ParticleType( props[7] )
	LWevent.radius = props[8]
	# Ranged = props[9], Volume = props[10]
	LWevent.total_column_depth = props[10]
	LWevent.x = 0
	LWevent.y = 0
	LWevent.z = 0

	weight = weight_event(LWevent)
	
	# this would alert us that something bad is happening 
	if weight==np.nan:
	    raise ValueError("Bad Weight!")
	
	return( weight*liveTime )


# load data
data_file = h5.File('data_output.h5', 'a')
injector_list = data_file.keys()
weight = []

for injector in injector_list:
	print(injector)
	for event in range(len( data_file[injector]['properties'] )):
		weight.append(get_weight( data_file[injector]['properties'][event]))
		#print("Event Weight: {}".format( get_weight( data_file[injector]['properties'][event]) ))

data_file.create_dataset('weight', data=weight)
data_file.close()
