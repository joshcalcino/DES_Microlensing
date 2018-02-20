import numpy as np


class DetectMicrolensingEvents(object):
    """ Detect microlensing events which are stored in a dictionary. The input dictionary must be of the following
    form:

    Data must be a dictionary where each key is an object, and each key points to a dictionary that contains the
    observational data. e.g.
                                -- obj = '1234'
                                -- data[obj] = {'g_mag_obs': [22., 22.1, 22.4, 22.1, 22.],
                                                'g_obs_time': [1., 21, 43, 65, 83],
                                                'g_mag_err': [ 0.1, 0.1, 0.1, 0.1, 0.1],
                                                .....,
                                }

    Maybe make it such that only one object is passed into the functions at a time, rather than an entire dict of
    objects.

    I should make this such that we do not require all of the 5 bands to have loads of observations in each.. But the
    observations that do exist must be consistent with a microlensing event (i.e. achromatic) for the object to pass.
    """

    def __init__(self):
        self.bands = ['g', 'r', 'i', 'z', 'Y']
        self.variability_threshold = 3.  # the threshold to class something as variable
        self.size_threshold = 3.

    def detect_events(self, large_data):
        """ Run through the microlensing pipeline and try and fine microlensing events. Record how many objects pass
        each cut.

        large_data --- a dictionary with object IDs as keys, and the object ID keys point to dictionaries that
        contain the observations.

        """

        total_num_objs = len(large_data)

        num_pass_first_cut = 0
        num_pass_second_cut = 0
        num_pass_third_cut = 0
        num_pass_fourth_cut = 0

        for key in large_data:
            """ Iterate over every single object and see if it passes the microlensing cuts. """

            object = large_data[key]

            if self.variability_cut(object):
                num_pass_first_cut += 1
            else:
                pass

            if self.red_cut(object):
                num_pass_second_cut += 1
            else:
                pass

            if self.achromaticity_cut(object):
                num_pass_third_cut += 1
            else:
                pass

        return num_pass_fourth_cut

    def variability_cut(self, data):
        """ Very basic cut for now.. just see if the sum of the errors is much smaller than the variance of the
        observations.. """

        variable_bands = 0
        non_variable_bands = 0

        for band in self.bands:
            obs_key = band + '_mag_obs'
            obs_err_key = band + '_mag_err'

            if obs_key in data:  # check if the key exists..
                if len(data[obs_key]) > self.size_threshold:
                    if np.std(data[obs_key]) > self.variability_threshold * np.median(data[obs_err_key]):
                        variable_bands += 1
                    else:
                        non_variable_bands += 1
                else:
                    pass
            else:
                pass

        return variable_bands > non_variable_bands

    def vary_only_once_cut(self, data):
        """ For a single lens and source, the source should only appear to vary once.. """

        pass

    def achromaticity_cut(self, data):
        pass

    def red_cut(self, data):
        """ Get rid of red stars, since they're far more likely to be variabile stars than others.. """

        if 'r_mag_obs' in data and 'g_mag_obs' in data:
            avg_r = np.average(data['r_mag_obs'])
            avg_g = np.average(data['g_mag_obs'])

            med_r = np.median(data['r_mag_obs'])
            med_g = np.median(data['g_mag_obs'])

            if avg_g - avg_r > 0.9 and med_g - med_r > 0.9:
                return False
            else:
                return True
