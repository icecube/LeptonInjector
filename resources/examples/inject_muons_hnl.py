# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu

# this example script ...
#   + imports the LeptonInjector libraries 
#   + creates two injectors, and their operator 
#   + tells the operator to execute the process 

import LeptonInjector as LI
from math import pi
import os 

# use this if you are on the Cobalt testbed 
#xs_folder = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0"

# for now, just use the example cross sections that come pre-installed. 
# this looks in the parent folder to this file's containing folder 
xs_folder = os.path.join( os.path.dirname(__file__), '..' )

n_events        = 5000
is_ranged       = True

# hnl events
diff_xs_hnl     = xs_folder + "/dsdxdy-hnl-N-nc-HERAPDF15NLO_EIG_central.fits"
total_xs_hnl    = xs_folder + "/sigma-hnl-N-nc-HERAPDF15NLO_EIG_central.fits"
final_1_hnl     = LI.Particle.HNL
final_2_hnl     = LI.Particle.Hadrons
injector_hnl    = LI.Injector( n_events , final_1_hnl, final_2_hnl, diff_xs_hnl, total_xs_hnl, is_ranged)
# hnlbar events
diff_xs_hnlbar  = xs_folder + "/dsdxdy-hnlbar-N-nc-HERAPDF15NLO_EIG_central.fits"
total_xs_hnlbar = xs_folder + "/sigma-hnlbar-N-nc-HERAPDF15NLO_EIG_central.fits"
final_1_hnlbar  = LI.Particle.HNLBar
final_2_hnlbar  = LI.Particle.Hadrons
injector_hnlbar = LI.Injector( n_events , final_1_hnlbar, final_2_hnlbar, diff_xs_hnlbar, total_xs_hnlbar, is_ranged)


deg = pi/180.

# define some defaults 
minE        = 1.     # [GeV]
maxE        = 10000.   # [GeV]
gamma       = 2. 
minZenith   = 80.*deg
maxZenith   = 180.*deg
minAzimuth  = 0.*deg
maxAzimuth  = 360.*deg

# construct the controller 
controller  = LI.Controller(injector_hnl, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)  
controller.AddInjector(injector_hnlbar)

# specify the output
controller.NameOutfile("./data_output.h5")
controller.NameLicFile("./config.lic")

# run the simulation
controller.Execute()
