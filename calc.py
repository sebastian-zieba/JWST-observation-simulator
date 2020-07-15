import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import numpy as np
import os

#Parameters
magK = 7.853 #https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf
Ts =6290 #https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf
metal = -0.24#https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf
logg = 4.29#https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf
Rp_earth = 3.95 #https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf
Rs_sun = 1.18 #https://iopscience.iop.org/article/10.3847/1538-3881/aa6e01/pdf


exo_dict = jdi.load_exo_dict()

exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
exo_dict['observation']['sat_unit'] = '%'
exo_dict['observation']['noccultations'] = 1 #number of transits
exo_dict['observation']['R'] = 50          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                             #without having to redo the calcualtion
exo_dict['observation']['baseline_unit'] = 'frac'  #Defines how you specify out of transit observing time
                                                    #'frac' : fraction of time in transit versus out = in/out
                                                    #'total' : total observing time (seconds)
exo_dict['observation']['baseline'] = 1 #in accordance with what was specified above (total observing time)

exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                             #to a wavelength dependent noise floor solution (units are ppm)
exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
exo_dict['star']['mag'] = magK             #magnitude of the system
exo_dict['star']['ref_wave'] = 2.22         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
exo_dict['star']['temp'] = Ts            #in K
exo_dict['star']['metal'] = metal            # as log Fe/H
exo_dict['star']['logg'] = logg              #log surface gravity cgs

exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
exo_dict['planet']['transit_duration'] = 4.693*60.0*60.0   #transit duration
exo_dict['planet']['td_unit'] = 's'
exo_dict['planet']['radius'] = Rp_earth
exo_dict['planet']['r_unit'] = 'R_earth'            #Any unit of distance in accordance with astropy.units can be added here
exo_dict['star']['radius'] = Rs_sun
exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit

#ORRRRR....
#if you wanted to instead to secondary transit at constant temperature
#exo_dict['planet']['f_unit'] = 'fp/f*'
#exo_dict['planet']['temp'] = T_day(Ts, ars, 0)[0]

#jdi.print_instruments()

result = jdi.run_pandexo(exo_dict,['NIRSpec G395H'])

n_lam = 2/(exo_dict['planet']['transit_duration']*result['FinalSpectrum']['error_w_floor']**2)

# in 1 sec
print('one sec:', np.sqrt(sum(n_lam))/(sum(n_lam))*1e6)
# in 1 min #"it only takes a minute girl"
print('one min:', np.sqrt(sum(n_lam)*60)/(sum(n_lam)*60)*1e6)
# in 30 minute
print('30 mins:', np.sqrt(sum(n_lam)*60*30)/(sum(n_lam)*60*30)*1e6)