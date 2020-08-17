#import numpy as np

class AncillaryData:
    """
	.. note::
		- Units for the orbital period and ephemeris can be anything as long as they are consistent (e.g. both in days).
		- The orbital path is calculated based on `t0` for primary transits and `t_secondary` for secondary eclipses.

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
        self.path_to_model = params['path_to_model']
        self.w_unit = params['w_unit']

        self.instrument = params['instrument']


def bins_new(x, y, y_err, n_bins):
    """
    Calculate maximum error for transit light curve calculation.

    :param plot: If ``True``, plots the error in the light curve model as a function of separation of centers.
    :type plot: bool
    :return: Truncation error (parts per million)
    :rtype: float
    """
    binned_x, binned_y, binned_y_err = np.zeros(n_bins), np.zeros(n_bins), np.zeros(n_bins)
    xmin = min(x)
    xmax = max(x)
    stepsize = (xmax-xmin)/n_bins
    for i in range(n_bins):
        ind = (x > xmin + i*stepsize) & (x < xmin + (i+1)*stepsize)
        binned_x[i] = np.mean(x[ind])
        binned_y[i] = np.mean(y[ind])
        binned_y_err[i] = np.sqrt(sum(y_err[ind]**2)) / sum(ind)
    return (binned_x, binned_y, binned_y_err)