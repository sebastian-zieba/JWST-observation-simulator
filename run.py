import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import numpy as np
import os

from astropy.io import ascii
from datetime import datetime
from shutil import copyfile

def str_to_num(str):
    'Return value of numeric literal string or ValueError exception'

    # Handle '0'
    if str == '0': return 0

    # Int/Float/Complex
    try:
        return int(str)
    except ValueError:
        pass
    try:
        return float(str)
    except ValueError:
        pass
    return str

def convert_to_bool(str):
    if str=="True": return True
    elif str=="False": return False
    else: return "String not equal to True or False"

class AncillaryData:
    """
    doc
    """

    def __init__(self, params):
        self.magK = str_to_num(params['magK'])
        self.Ts = str_to_num(params['Ts'])
        self.metal = str_to_num(params['metal'])
        self.logg = str_to_num(params['logg'])
        self.Rp_earth = str_to_num(params['Rp_earth'])
        self.Rs_sun = str_to_num(params['Rs_sun'])
        self.D = str_to_num(params['D'])
        self.instrument = 'NIRSpec ' + params['NIRSpec']
        self.noccultations = str_to_num(params['noccultations'])
        self.R = str_to_num(params['R'])
        self.baseline = str_to_num(params['baseline'])
        self.output = convert_to_bool(params['output'])


def make_dict(table):
    return {x['parameter']: x['value'] for x in table}


params = make_dict(ascii.read("config/params.txt", Reader=ascii.CommentedHeader))
ancil = AncillaryData(params)





### SAVE OUTPUT IF WISHED

if ancil.output == True:
    dirname = "runs_dir/" + datetime.strftime(datetime.now(), '%m_%d_%H_%M')
    if not os.path.exists(dirname): os.makedirs(dirname)

    resultfile = open(dirname+'/results.txt', 'w')

    copyfile("./config/params.txt", dirname+"/params.txt")        #stores obs_par.txt



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