import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO

import os

from astropy.io import ascii
from datetime import datetime
from shutil import copyfile

import yaml


with open('./config/params.yaml', 'r') as file:
    params = yaml.safe_load(file)

for k,v in params.items():
    print(k, type(v))

class AncillaryData:
    """
    doc
    """

    def __init__(self, params):
        self.magK = params['magK']
        self.Ts = params['Ts']
        self.metal = params['metal']
        self.logg = params['logg']
        self.Rp_earth = params['Rp_earth']
        self.Rs_sun = params['Rs_sun']
        self.D = params['D']

        self.noccultations = params['noccultations']
        self.R = params['R']
        self.baseline = params['baseline']
        self.output = params['output']
        self.path = params['path']

        self.instrument = params['instrument']


ancil = AncillaryData(params)


### SAVE OUTPUT IF WISHED

if ancil.output == True:
    dirname = "runs_dir/" + datetime.strftime(datetime.now(), '%m_%d_%H_%M')
    if not os.path.exists(dirname): os.makedirs(dirname)

    resultfile = open(dirname+'/results.txt', 'w')

    copyfile("./config/params.txt", dirname+"/params.txt")        #stores obs_par.txt


### PLOT SPECTRUM

datafile_nometal = 'trans_spect_hd106315c_LKrescale_m0.5_co1.0nc_f0.1.txt'
datafile_metal = 'trans_spect_hd106315c_LKrescale_m2.5_co1.0nc.txt'

print(ancil.path)

wvl_nometal, depth_nometal = np.loadtxt(ancil.path + datafile_nometal, skiprows=1).T
wvl_metal, depth_metal = np.loadtxt(ancil.path + datafile_metal, skiprows=1).T

print(np.median(depth_nometal))
print(np.median(depth_metal))

fig, ax = plt.subplots(1,1,figsize=(11,5))

ax.errorbar(wvl_nometal, depth_nometal*1e6, yerr = [67.14]*len(depth_nometal))

ax.set_xlim(0, 21)
#ax.set_xlim(0.7, 21)
ax.set_ylim(1e3, 3e3)

ax.set_xlabel('wavelength (microns)')
ax.set_ylabel('transit depth (ppm)')

#ax.set_xscale('log')

plt.savefig('test.png')


### RUN PANDEXO

exo_dict = jdi.load_exo_dict()

exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
exo_dict['observation']['sat_unit'] = '%'
exo_dict['observation']['noccultations'] = ancil.noccultations #number of transits
exo_dict['observation']['R'] = ancil.R          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                             #without having to redo the calcualtion
exo_dict['observation']['baseline_unit'] = 'frac'  #Defines how you specify out of transit observing time
                                                    #'frac' : fraction of time in transit versus out = in/out
                                                    #'total' : total observing time (seconds)
exo_dict['observation']['baseline'] = ancil.baseline #in accordance with what was specified above (total observing time)

exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                             #to a wavelength dependent noise floor solution (units are ppm)
exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
exo_dict['star']['mag'] = ancil.magK             #magnitude of the system
exo_dict['star']['ref_wave'] = 2.22         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
exo_dict['star']['temp'] = ancil.Ts            #in K
exo_dict['star']['metal'] = ancil.metal            # as log Fe/H
exo_dict['star']['logg'] = ancil.logg             #log surface gravity cgs

exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
exo_dict['planet']['transit_duration'] = ancil.D*60.0*60.0   #transit duration
exo_dict['planet']['td_unit'] = 's'
exo_dict['planet']['radius'] = ancil.Rp_earth
exo_dict['planet']['r_unit'] = 'R_earth'            #Any unit of distance in accordance with astropy.units can be added here
exo_dict['star']['radius'] = ancil.Rs_sun
exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit

#ORRRRR....
#if you wanted to instead to secondary transit at constant temperature
#exo_dict['planet']['f_unit'] = 'fp/f*'
#exo_dict['planet']['temp'] = T_day(Ts, ars, 0)[0]

#jdi.print_instruments()

result = jdi.run_pandexo(exo_dict,[ancil.instrument])
#result = jdi.run_pandexo(exo_dict,['NIRSpec PRISM'])

n_lam = 2/(exo_dict['planet']['transit_duration']*result['FinalSpectrum']['error_w_floor']**2)

# in 1 sec
print('one sec:', np.sqrt(sum(n_lam))/(sum(n_lam))*1e6)
# in 1 min #"it only takes a minute girl"
print('one min:', np.sqrt(sum(n_lam)*60)/(sum(n_lam)*60)*1e6)
# in 30 minute
print('30 mins:', np.sqrt(sum(n_lam)*60*30)/(sum(n_lam)*60*30)*1e6)

if ancil.output == True:
    for k, v in result['timing'].items():
        print("{:<40} {:<25}".format(k, v), file=resultfile)

    print('\n', file=resultfile)

    for k, v in result['warning'].items():
        print("{:<25} {:<25}".format(k, v), file=resultfile)

    print('\n', file=resultfile)

    print('rms one sec:', np.sqrt(sum(n_lam)) / (sum(n_lam)) * 1e6, file=resultfile)
    print('rms one min:', np.sqrt(sum(n_lam) * 60) / (sum(n_lam) * 60) * 1e6, file=resultfile)
    print('rms 30 mins:', np.sqrt(sum(n_lam) * 60 * 30) / (sum(n_lam) * 60 * 30) * 1e6, file=resultfile)

if ancil.output == True:
    resultfile.close()