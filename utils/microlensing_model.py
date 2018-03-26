import numpy as np


""" This module contains the equations to plot the Paczynski curve, and a few other useful functions for
    computing and handling light-curves. """


def get_amplification(t, t_hat, u_min, t_max):
    """ The amplification caused by a point lens, on a point source """
    u = get_u(t, t_hat, u_min, t_max)
    return (u ** 2 + 2.) / (u * np.sqrt(u ** 2 + 4))


def get_u(t, t_hat, u_min, t_max):
    u = np.sqrt(u_min ** 2 + ((t - t_max) / t_hat) ** 2)
    return u


def get_delta_mag(t, t_hat, u_min, t_max):  # change in the magnitude of the star due to the lensing
    # u = mm.get_u(t, t_hat, u_min, t_max)
    A = get_amplification(t, t_hat, u_min, t_max)
    delta_mag = 2.5 * np.log10(A)
    return delta_mag


def mag_to_flux(mag):
    m_vega = -21.10
    flux = 10**((m_vega - mag)/2.5)
    return flux


def flux_to_mag(flux):
    m_vega = -21.10
    mag = -2.5*np.log10(flux) + m_vega
    return mag


def microlensing_fit_function(t, t_hat, u_min, t_max, *f_bands):

    print(t)
    print(t[0, :])

    A = get_amplification(t[0, :], t_hat, u_min, t_max)

    tmp_band = t[1, 0]
    mag_arr = []

    i = 0

    for i in range(0, len(f_bands)):
        band = t[1, i]
        if band == tmp_band:
            mag_arr = np.append(mag_arr, f_bands[i])
        else:
            tmp_band = band
            mag_arr = np.append(mag_arr, f_bands[i])
    return A + mag_arr

