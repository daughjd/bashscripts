# fakeps.py

from __future__ import division

__doc__ = """Generate point-like simulation given diffuse simulation."""

import numpy as np


def _get (o, key):
    """Get member key from object o.

    First subscript (o[key]) and then member (getattr (o, key)) access will
    be tried.

    """
    try:
        return o[key]
    except:
        try:
            return getattr (o, key)
        except:
            raise ValueError ('could not find "{0}" key in {1}'.format (
                key, o))


def delta_angle (fit1, fit2, azimuth1=None, azimuth2=None):
    """Return the space angle between ``fit1`` and ``fit2``.

    :type   fit1: float or mapping
    :param  fit1: Zenith 1, or a mapping containing array values with keys
        'zenith' and 'azimuth'.

    :type   fit2: float or mapping
    :param  fit2: Zenith 2, or a mapping containing array values with keys
        'zenith' and 'azimuth'.

    :type   azimuth1: float
    :param  azimuth1: Azimuth 1.

    :type   azimuth2: float
    :param  azimuth2: Azimuth 2.

    :return: ``numpy.ndarray`` of delta_angle values.

    """
    sin, cos, arccos = np.sin, np.cos, np.arccos
    if azimuth1 is None:
        assert azimuth2 is None
        z1, a1 = misc._get (fit1, 'zenith'), misc._get (fit1, 'azimuth')
        z2, a2 = misc._get (fit2, 'zenith'), misc._get (fit2, 'azimuth')
    else:
        z1, a1 = fit1, azimuth1
        z2, a2 = fit2, azimuth2

    return arccos (
            sin (z1) * cos (a1)  *  sin (z2) * cos (a2)
            + sin (z1) * sin (a1)  *  sin (z2) * sin (a2)
            + cos (z1) * cos (z2))


def ps_events (point, diffsim,
        mode='circle',
        frac_of_sphere=.05):
    """Get the idx (array of bools) for representative events.

    :type   point: mapping
    :param  point: Object with 'zenith' and 'azimuth' keys and float values.

    :type   diffsim: mapping
    :param  diffsim: Object with 'true_zenith' and 'true_azimuth' keys and
        ndarray values.

    :type   mode: str
    :param  mode: Either 'circle' (nearby events) or 'band' (in zenith)

    :type   frac_of_sphere: float
    :param  frac_of_sphere: Determines width of circle or zenith band.

    :return: ndarray of bools.

    Note: where a mapping is required, either subscript (object[key]) or
    attribute (getattr (object, key)) is allowed.

    """
    if mode == 'circle':
        dpsi = delta_angle (
                _get (point, 'zenith'), _get (diffsim, 'true_zenith'),
                _get (point, 'azimuth'), _get (diffsim, 'true_azimuth'))
        dpsi_thresh = np.arccos (1 - 2 * frac_of_sphere)
        idx = dpsi <= dpsi_thresh
        return idx
    elif mode == 'band':
        dcz = frac_of_sphere * 2
        cz = np.cos (_get (point, 'zenith'))
        diff_cz = np.cos (_get (diffsim, 'true_zenith'))
        diff_cz_min = diff_cz.min ()
        diff_cz_max = diff_cz.max ()
        # if cz - .5 * dcz < diff_cz_min:
        #     cz_min, cz_max = (diff_cz_min, diff_cz_min + dcz)
        if cz + .5 * dcz > diff_cz_max:
            cz_min, cz_max = (diff_cz_max - dcz, diff_cz_max)
        else:
            cz_min, cz_max = cz - .5 * dcz, cz + .5 * dcz
        idx = (cz_min <= diff_cz) * (diff_cz < cz_max)
        return idx
    else:
        raise NotImplementedError ('mode="{0}" not implemented'.format (mode))

def ps_weights (diff_fluence, point, diffsim,
        mode='circle',
        frac_of_sphere=.05):
    """Get the weights of representative events.

    :type   diff_fluence: function
    :param  diff_fluence: Differential fluence as a function of energy.

    All other arguments are passed thru to ps_events. For this function,
    diffsim must also contain 'true_energy' and 'oneweight' arrays, as well as
    an 'n_gen'->float item.

    :return: ndarray of floats.

    """
    idx = ps_events (point, diffsim,
            mode=mode, frac_of_sphere=frac_of_sphere)
    oneweight = _get (diffsim, 'oneweight')
    try:
        n_gen = _get (diffsim, 'n_gen')
    except:
        n_gen = diffsim.meta.n_gen
    weights = oneweight[idx] \
            * diff_fluence (_get (diffsim, 'true_energy')[idx]) \
            / (4 * np.pi * frac_of_sphere * n_gen)
    return weights

def corrected_frac_of_sphere (frac_of_sphere, point, horizon,
        mode='circle'):
    """Get the corrected fraction of sphere if there is a horizon.

    :type   frac_of_sphere: float
    :param  frac_of_sphere: The original fraction of sphere

    :type   point: mapping
    :param  point: Object with 'zenith' and 'azimuth' keys and float values.

    :type   horizon: float
    :param  horizon: The horizon in IceCube coordinates in radians.

    :type   mode: str
    :param  mode: Either 'circle' (nearby events) or 'band' (in zenith)

    :return: float: the fraction of the sphere actually drawn from.

    """

    if mode == 'circle':
        delta = np.arccos (1 - 2 * frac_of_sphere)
        alpha = _get (point, 'zenith') - horizon
        if alpha < 0:
            raise NotImplementedError (
                    'not implemented for Southern sky (zenith < horizon)')
        if alpha > delta:
            return frac_of_sphere
        frac_remain = 1 + 1 / (np.pi * delta**2) \
                * (alpha * np.sqrt (delta**2 - alpha**2)
                        - delta**2 * np.arccos (alpha / delta))
        return frac_remain * frac_of_sphere
    else:
        raise NotImplementedError (
                'not implemented for mode "{0}"'.format (mode))
