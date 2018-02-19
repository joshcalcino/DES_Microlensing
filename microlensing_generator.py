import numpy as np
import astropy.constants as const
from utils.loaddata import load_csv_data
import matplotlib.pyplot as plt


class GenerateMicrolensingEvent(object):
    """This class will generate a random microlensing event for every item in a dictionary of lightcurves.
     All of the parameters are randomly generated, over a default or specified range."""

    def __init__(self, u_min_min=0., u_min_max=1.2, t_max_min=0., t_max_max=1500.):
        self.u_min_min = u_min_min
        self.u_min_max = u_min_max
        self.t_max_min = t_max_min
        self.t_max_max = t_max_max

    def get_t_max(self):  # time of maximum approach, in years since the start of DES.. for now
        t_max = np.random.uniform(low=self.t_max_min, high=self.t_max_max)
        return t_max

    def get_u_min(self):
        return np.random.uniform(self.u_min_min, self.u_min_max)

    def get_u(self, t, t_hat, u_min, t_max):
        u = np.sqrt(u_min ** 2 + ((t - t_max) / t_hat) ** 2)
        return u

    def get_amplification(self, t, t_hat, u_min, t_max):
        u = self.get_u(t, t_hat, u_min, t_max)
        return (u**2 + 2.)/(u*np.sqrt(u**2 + 4))

    def get_delta_mag(self, t, t_hat, u_min, t_max):  # change in the magnitude of the star due to the lensing
        u = self.get_u(t, t_hat, u_min, t_max)
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        delta_mag = 2.5 * np.log10(A)
        return delta_mag

    def mag_to_flux(self, mag):
        m_vega = -21.10
        flux = 10**((m_vega - mag)/2.5)
        return flux

    def flux_to_mag(self, flux):
        m_vega = -21.10
        mag = -2.5*np.log10(flux) + m_vega
        return mag

    def generate_microlensing_catalogue(self, data, t_hat):
        """ Return a catalogue of lightcurves with microlensing events implanted into them.

        data must be a dictionary where each key is an object, and each key points to a dictionary that contains the
        observational data. e.g.
                                    -- obj = '1234'
                                    -- data[obj] = {'g_mag_obs': [22., 22.1, 22.4, 22.1, 22.],
                                                    'g_obs_time': [1., 21, 43, 65, 83],
                                                    'g_mag_err': [ 0.1, 0.1, 0.1, 0.1, 0.1],
                                                    .....,
                                    }

        Thus, each object has a dictionary of all of the observational parameters, such as the time of the observation,
        the observed magnitude, and the error on the observed magnitude. This is done for each band:
                                                                                        ['g', 'r', 'i', 'z', 'Y']
        """

        microlensing_dict = {}

        des_colour = ['g', 'r', 'i', 'z', 'Y']

        for key in data:
            obj = data[key]

            # generate the microlensing parameters
            tmp_t_max = self.get_t_max()
            tmp_u_min = self.get_u_min()

            tmp_dict = {}

            for colour in des_colour:
                if colour + '_obs_time' in obj:
                    obs_time = obj[colour + '_obs_time']
                    obs_mag_err = obj[colour + '_mag_err']
                    obs_mag = obj[colour + '_mag_obs']

                    obs_time -= 56255.09246  # the time of the first observation

                    obs_flux = self.mag_to_flux(obs_mag)

                    A = self.get_amplification(obs_time, t_hat, tmp_u_min, tmp_t_max)

                    flux_sim = obs_flux + (A - 1.)*np.median(obs_flux)

                    mag_sim = self.flux_to_mag(flux_sim)

                    tmp_dict[colour + '_obs_time'] = obs_time
                    tmp_dict[colour + '_mag_err'] = obs_mag_err
                    tmp_dict[colour + '_mag_obs'] = mag_sim

                else:
                    pass

            tmp_dict['u_min'] = tmp_u_min
            tmp_dict['t_max'] = tmp_t_max

            t_max_rnd = np.round(tmp_t_max, 4)
            u_min_rnd = np.round(tmp_t_max, 4)

            # microlensing_dict['tmax'+str(t_max_rnd)+'_umin'+str(u_min_rnd)] = tmp_dict
            microlensing_dict[key] = tmp_dict

        return microlensing_dict


data = load_csv_data('data/high_nepochs_objects.csv')

gme = GenerateMicrolensingEvent()

ml_dicts = gme.generate_microlensing_catalogue(data, 300.)

key_list = ml_dicts.keys()

test_dict = ml_dicts[key_list[300]]
test_dict2 = data[key_list[300]]

des_filter_arr = ['g', 'r', 'i', 'z', 'Y']

# for des_filter in des_filter_arr:

des_filter = 'r'

t = np.linspace(0, 1200, 1200)


for i in range(0, 1000):
    test_dict = ml_dicts[key_list[i]]
    test_dict2 = data[key_list[i]]

    plt.errorbar(test_dict[des_filter + '_obs_time'], test_dict[des_filter + '_mag_obs'],
                 yerr=test_dict[des_filter + '_mag_err'], label=des_filter, fmt='o')

    plt.errorbar(test_dict2[des_filter + '_obs_time'], test_dict2[des_filter + '_mag_obs'],
                 yerr=test_dict2[des_filter + '_mag_err'], label=des_filter, fmt='o')

    median = np.median(test_dict2[des_filter + '_mag_obs'])

    u_min, t_max = test_dict['u_min'], test_dict['t_max']

    lc = median - gme.get_delta_mag(t, 300., u_min, t_max)

    plt.plot(t, lc)

    plt.gca().invert_yaxis()  # So that decreasing y-axis is up, consistent with smaller magnitudes being brighter

    plt.legend()
    plt.show()

# t_hat = np.logspace(0, 3, 100)
#
# ml_dicts = []
#
# for i in range(0, len(t_hat)):
#     ml_dicts.append(gme.generate_microlensing_catalogue(data, t_hat[i]))


