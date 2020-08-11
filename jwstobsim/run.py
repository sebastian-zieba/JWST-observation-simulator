import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO

import os

from jwstobsim.utils import *


#from astropy.io import ascii
from datetime import datetime
from shutil import copyfile

import yaml

import sys
#__all__ = ['AncillaryData']
#lib_dir = '/home/zieba/Desktop/Projects/JWST-observation-simulator/sim/'
#sys.path.insert(0,lib_dir)
sys.path.insert(0, os.path.abspath('..'))





### BINNING FUNCTION ###
# Note that the binning is constant
# One would have to modify that function if the binsize should be a function of wavelength

###MY OLD RELAYABLE BINNING FUNCTION
#def bins(x, y, y_err, n_bins):
#    indexes = np.array_split(np.arange(len(x)), n_bins)
#    #print((indexes[0]).size)
#    binned_x = np.array([np.mean(x[a]) for a in indexes])
#    binned_y = np.array([np.mean(y[a]) for a in indexes])
#    binned_y_err = np.array([np.sqrt(sum(y_err[a]**2)) / np.size(a) for a in indexes])
#    return (binned_x, binned_y, binned_y_err)



### IMPORT CONFIG FILE ###

yaml_path = './config/params.yaml'
#yaml_path = '/home/zieba/Desktop/Projects/JWST-observation-simulator/jwstobsim/config/params.yaml'

with open(yaml_path, 'r') as file:
    params = yaml.safe_load(file)



#for k,v in params.items():
#    print(k, type(v))



ancil = AncillaryData(params)

### SAVE OUTPUT IF WISHED ###

if ancil.output == True:
    dirname = "runs_dir/" + datetime.strftime(datetime.now(), '%b_%d_%H:%M:%S')
    if not os.path.exists(dirname): os.makedirs(dirname)

    resultfile = open(dirname+'/results.txt', 'w')

    copyfile(yaml_path, dirname+"/params.yaml")        #stores a copy of params.yaml in runs_dir



### GET THROUGHPUT DATA ###

try:
    thru_dict = jdi.get_thruput(ancil.instrument)
except:
    print(jdi.print_instruments())
    raise

wvl_tp = thru_dict['wave']
pce_tp = thru_dict['pce']

mask_pce=[]
for i in pce_tp:
    if i/max(pce_tp)>0.05: #only take data where photon conversion effiency is above 5% max
        mask_pce.append(True)
    else:
        mask_pce.append(False)

wvl_range = [wvl_tp[mask_pce][0], wvl_tp[mask_pce][-1]] #range in which instrument operates above 5%max

print('Wavelength range of the instrument: ', [float('{:.3f}'.format(num)) for num in wvl_range])

if ancil.path_to_model != False:
    wvl_model, flux_model = np.loadtxt(ancil.path_to_model).T
    print('Wavelength range of the model: ', [min(wvl_model), max(wvl_model)])

    wvl_range = [max([wvl_tp[mask_pce][0], min(wvl_model)]), min([wvl_tp[mask_pce][-1], max(wvl_model)])]

    print('Final Wavelength range: ', [float('{:.3f}'.format(num)) for num in wvl_range])

### RUN PANDEXO

exo_dict = jdi.load_exo_dict()

### OBSERVATION INPUTS ###

exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
exo_dict['observation']['sat_unit'] = '%'
exo_dict['observation']['noccultations'] = ancil.noccultations #number of transits
exo_dict['observation']['R'] = None#ancil.R          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                             #without having to redo the calcualtion
exo_dict['observation']['baseline_unit'] = 'frac'  #Defines how you specify out of transit observing time
                                                    #'frac' : fraction of time in transit versus out = in/out
                                                    #'total' : total observing time (seconds)
exo_dict['observation']['baseline'] = ancil.baseline #in accordance with what was specified above (total observing time)

exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                             #to a wavelength dependent noise floor solution (units are ppm)

### HOST STAR INPUTS ###


exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
exo_dict['star']['mag'] = ancil.magK             #magnitude of the system
exo_dict['star']['ref_wave'] = 2.22         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
exo_dict['star']['temp'] = ancil.Ts            #in K
exo_dict['star']['metal'] = ancil.metal            # as log Fe/H
exo_dict['star']['logg'] = ancil.logg             #log surface gravity cgs


### EXOPLANET INPUTS ###

### IF NO MODEL, USE CONSTANT

if ancil.path_to_model == False:
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

### IF PATH HAS MODEL, USE IT

else:
    exo_dict['planet']['type'] = 'user'  # tells pandexo you are uploading your own spectrum
    exo_dict['planet']['exopath'] = ancil.path_to_model
    exo_dict['planet']['w_unit'] = ancil.w_unit  # other options include "um","nm" ,"Angs", "sec" (for phase curves)
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'  # other options are 'fp/f*'
    exo_dict['planet']['transit_duration'] = ancil.D*60.0*60.0  # transit duration
    exo_dict['planet']['td_unit'] = 's'  # Any unit of time in accordance with astropy.units can be added


### RUN IT ###

result = jdi.run_pandexo(exo_dict,[ancil.instrument])


# ONLY USE THE PARTS WHERE THE INSTUMENT OPERATES WELL
mask = [wvl_range[0] < wvli < wvl_range[1] for wvli in result['FinalSpectrum']['wave']]


# PLOT IT
fig, ax = plt.subplots(1,1,figsize=(11,5))

ax.axvline(wvl_range[0], c='r', ls='--', label='wvl range')
ax.axvline(wvl_range[1], c='r', ls='--')

ax.set_xlabel('wavelength (microns)')
ax.set_ylabel('tansit depth (rprs**2)')


### CALC NUMBER OF BINS ###

n_bins = (wvl_range[1]-wvl_range[0])/(wvl_range[1]+wvl_range[0]) * 2 * ancil.R
n_bins = int(round(n_bins)) #rounding will lead to a resolution which is a bit off the wanted value (could be improved)
print('Number of bins = ', n_bins)
print('Final resolution = ', '{:.3f}'.format( n_bins/2 * (wvl_range[1]+wvl_range[0])/(wvl_range[1]-wvl_range[0]) ))

a,b,c=bins_new(result['FinalSpectrum']['wave'][mask], result['FinalSpectrum']['spectrum'][mask], result['FinalSpectrum']['error_w_floor'][mask], n_bins)
#dont judge me on the naming

ind = np.isnan(c)
a, b, c = a[~ind], b[~ind], c[~ind]

b_rand = b + c * np.random.normal(0,1,len(b))

ax.errorbar(a, b_rand, yerr=c, fmt='.', ls='', label='spectrum')

print('median error bar size: ', '{:.3e}'.format( np.median(c)) )


# PLOT TRANSMISSION AND MAKE LIMITS BETTER

xlim0, xlim1 = ax.get_xlim()

ylim0, ylim1 = ax.get_ylim()
ylim_range = ylim1 - ylim0
ylim_scaler=0.4
ax.set_ylim(ylim0 - ylim_scaler * ylim_range, ylim1 + ylim_scaler*1.2 * ylim_range)

ylim0, ylim1 = ax.get_ylim()
pce_tp_rescale = pce_tp * (ylim1-ylim0)/(max(pce_tp)*3-0) + ylim0
ax.plot(wvl_tp, pce_tp_rescale, ls='-.', label='bandpass (arb. units)')

ax.plot(wvl_model, flux_model, c='k', alpha=0.2, label='model')
ax.scatter(a, b, marker='s', facecolors='none', edgecolors='r', label='mean bin')

ax.set_xlim(xlim0,xlim1)

plt.title('{0}, #Transits = {1}, R = {2:.1f}'.format(ancil.instrument, ancil.noccultations, n_bins/2 * (wvl_range[1]+wvl_range[0])/(wvl_range[1]-wvl_range[0])))

plt.legend(loc=1)

plt.savefig(dirname + '/spectrum.png', dpi=200)
plt.show()

### SAVE SPECTRUM ###

np.savetxt(dirname + "/spectrum.txt", list(zip(a,b,b_rand,c)))

n_lam = 2/(exo_dict['planet']['transit_duration']*result['FinalSpectrum']['error_w_floor']**2)

# in 1 sec
print('RMS one sec:', '{:.3f}'.format(np.sqrt(sum(n_lam))/(sum(n_lam))*1e6))
# in 1 min #"it only takes a minute girl"
print('RMS one min:', '{:.3f}'.format(np.sqrt(sum(n_lam)*60)/(sum(n_lam)*60)*1e6))
# in 30 minute
print('RMS 30 mins:', '{:.3f}'.format(np.sqrt(sum(n_lam)*60*30)/(sum(n_lam)*60*30)*1e6))


if ancil.output == True:
    for k, v in result['timing'].items():
        print("{:<40} {:<25}".format(k, v), file=resultfile)

    print('\n', file=resultfile)

    for k, v in result['warning'].items():
        print("{:<25} {:<25}".format(k, v), file=resultfile)

    print('\n', file=resultfile)

    print('Number of bins = ', n_bins, file=resultfile)
    print('Final resolution = ',
          '{:.3f}'.format(n_bins / 2 * (wvl_range[1] + wvl_range[0]) / (wvl_range[1] - wvl_range[0])), file=resultfile)
    print('median error bar size: ', '{:.3e}'.format(np.median(c)), file=resultfile)

    print('\n', file=resultfile)

    print('rms one sec:', '{:.3f}'.format(np.sqrt(sum(n_lam)) / (sum(n_lam)) * 1e6), file=resultfile)
    print('rms one min:', '{:.3f}'.format(np.sqrt(sum(n_lam) * 60) / (sum(n_lam) * 60) * 1e6), file=resultfile)
    print('rms 30 mins:', '{:.3f}'.format(np.sqrt(sum(n_lam) * 60 * 30) / (sum(n_lam) * 60 * 30) * 1e6), file=resultfile)

if ancil.output == True:
    resultfile.close()
