# grblist.py


from __future__ import division, print_function

__doc__ = """Hold information about a list of GRBs."""

import copy
import datetime
from itertools import izip
import numpy as np
pi = np.pi
import re

from icecube.umdtools.vars_class import Vars
import datetime

class GRB (object):

    """Specifics of a single GRB."""

    def __init__ (self, name, t_start, t_end, ra, dec, sigma,
                  run_number=0, gbm_pos=False,z=None):
        """Construct a GRB.

        A new :class:`grbdb.GRB` instance is created with the given properties.
        The GRB zenith and azimuth in IceCube coordinates are calculated using
        icecube.coordinate_service.

        :type   name: str
        :param  name: The name.

        :type   t_start: datetime.datetime
        :param  t_start: The start time.

        :type   t_end: datetime.datetime
        :param  t_end: The end time.

        :type   ra: float
        :param  ra: The right ascension (radians).

        :type   dec: float
        :param  dec: The declination (radians).

        :type   sigma: float
        :param  dec: The angular error estimate (radians).

        :type   run_number: int
        :param  run_number: The IceCube run during which this GRB took place.

        """

        from icecube.coordinate_service import calendar_date_2_mjd
        from icecube.coordinate_service import equa_2_local_zenith__inv
        from icecube.coordinate_service import equa_2_local_azimuth__inv
        self.name = name
        self.t_start = t_start
        self.t_end = t_end
        self.ra = ra
        self.dec = dec
        self.sigma = sigma
        self.run_number = run_number
        MJD = calendar_date_2_mjd (
                t_start.year, t_start.month, t_start.day,
                t_start.hour, t_start.minute, t_start.second)
        self.zenith = equa_2_local_zenith__inv (ra, dec, MJD, 2000)
        self.azimuth = equa_2_local_azimuth__inv (ra, dec, MJD, 2000)
        dt = self.t_end - self.t_start
        self.t_100 = dt.microseconds / 1e6 + dt.seconds + dt.days * 86400.
        self.gbm_pos = gbm_pos
	self.z = z
        # self.t100_flnc = t100_flnc
        # self.flnc = flnc
        # self.alpha = alpha
        # self.beta = beta
        # self.epeak = epeak
        # self.emin = emin
        # self.emax = emax
        # self.z = z


    def __repr__ (self):
        props = [
            'name={0}'.format (self.name),
            't_start={0}'.format (self.t_start),
            't_end={0}'.format (self.t_end),
            't_100={0:.2f}'.format (self.t_100),
            'ra={0:.2f}'.format (self.ra),
            'dec={0:.2f}'.format (self.dec),
            'zenith={0:.2f}'.format (self.zenith),
            'azimuth={0:.2f}'.format (self.azimuth),
            'sigma={0:.2e}'.format (self.sigma),
            'run_number={0}'.format (self.run_number or 'unknown'),
            'gbm_pos={0}'.format (self.gbm_pos),
            # 't100_flnc={0}'.format (self.t100_flnc),
            # 'flnc={0}'.format (self.flnc),
            # 'alpha={0}'.format (self.alpha),
            # 'beta={0}'.format (self.beta),
            # 'epeak={0}'.format (self.epeak),
            # 'emin={0}'.format (self.emin),
            # 'emax={0}'.format (self.emax),
            'z={0}'.format (self.z),
            ]
        out = '{0}.{1}({2})'.format (
                self.__class__.__module__,
                self.__class__.__name__,
                ','.join (props)
                )
        return out

    def __cmp__ (self, other):
        return cmp (self.name, other.name) or cmp (self.t_start, other.t_start)


class GRBColumn (object):

    """Column of data for a :class:`GRBList`."""

    def __init__ (self, grblist, array=False):
        self._grblist = grblist
        self._array = array
        # props = ['name', 't_start', 't_end', 't_100',
        #         'ra', 'dec', 'zenith', 'azimuth', 'sigma', 'run_number',
        #         'gbm_pos', 't100_flnc', 'flnc', 'alpha', 'beta', 
        #          'epeak', 'emin', 'emax', 'z']
        props = ['name', 't_start', 't_end', 't_100',
                'ra', 'dec', 'zenith', 'azimuth', 'sigma', 'run_number',
                'gbm_pos','z']

        for prop in props:
            setattr (self, prop, None)

    def __getattribute__ (self, name):
        if name[0] == '_':
            return object.__getattribute__ (self, name)
        else:
            return self[name]

    def __getitem__ (self, name):
        x = [getattr (grb, name) for grb in self._grblist]
        if self._array:
            return np.array (x)
        else:
            return x


class GRBList (object):

    """A list of :class:`GRB` s.

    The list is kept sorted at all times.
    """

    def __init__ (self, name='', grbs=[]):
        """Construct a new GRBList.

        :type   name: str
        :param  name: The name of this GRBList.

        :type   grbs: list
        :param  grbs: Some bursts with which to populate the new GRBList.

        """
        self.name = name
        self.grbs = []
        for grb in grbs:
            self.add (grb)
        self.lists = GRBColumn (self, False)
        self.arrays = GRBColumn (self, True)

    def __len__ (self):
        return len (self.grbs)

    def __iter__ (self):
        return iter (self.grbs)

    def __getitem__ (self, i):
        if isinstance (i, str):
            index = self.lists.name.index (i)
        else:
            index = i
        return self.grbs[index]

    def add (self, grb):
        """Add a GRB to the sample."""
        self.grbs.append (grb)
        self.grbs.sort ()

    def remove (self, grb):
        """Remove a GRB from the sample.

        :type   grb: :class:`GRB` or str
        :param  grb: A GRB object, or the name of the GRB to remove.

        """
        if isinstance (grb, GRB):
            index = self.grbs.index (grb)
        else:
            index = self.lists.name.index (grb)
            grb = self.grbs[index]
        self.grbs.remove (grb)

    def apply_cut (self, i):
        """Only keep GRBs where the corresponding element in ``i`` is True.

        :type   i: numpy.ndarray of bools
        :param  i: The indexing array.

        """
        keep_grbs = list (np.array (self.grbs)[i])
        self.grbs = []
        for grb in keep_grbs:
            self.add (grb)


    def get_copy (self, name='', i=None):
        """Make a copy GRBList.

        :type   name: str
        :param  name: The name of this GRBList.

        :type   i: numpy.ndarray of bools
        :param  i: The indexing array -- which elements to keep.

        :return: :class:`GRBList`

        """
        out = copy.deepcopy (self)
        if i is not None:
            out.apply_cut (i)
        return out


def grb_list_from_grb_web (
        list_name='',
        t_min=None, t_max=None,
        zenith_min=None, zenith_max=None,
        grb_name_blacklist=[]):
    """Get a GRBList from the GRB-web database.

    :type   list_name: str
    :param  name: The name of this GRB list.

    :type   t_min: datetime.datetime
    :param  t_min: Earliest GRB time to include.

    :type   t_max: datetime.datetime
    :param  t_max: Latest GRB time to include.

    :type   zenith_min: float
    :param  zenith_min: Minimum zenith to include (radians).

    :type   zenith_max: float
    :param  zenith_max: Maximum zenith to include (radians).

    :type   grb_name_blacklist: sequence
    :param  grb_name_blacklist: Sequence of names of GRBs to exclude.

    """
    # imported here because MySQL is only needed for import, not general use
    try:
        import MySQLdb as pymsql
    except:
        import pymysql
        MySQLdb = pymysql

    # first, get run start and stop times from IceCube live
    #db = MySQLdb.connect (
    #        host='cygnus.icecube.wisc.edu',
    #        user='icecube', passwd='skua', db='live')
    db = pymysql.Connect (
            host='cygnus.icecube.wisc.edu',
            user='icecube', passwd='skua', db='live')
    cursor = db.cursor ()
    cursor.execute ("""SELECT runNumber,tStart,tStop FROM livedata_run;""")
    data = cursor.fetchall ()
    i3_run_number = np.array ([d[0] for d in data if None not in d[1:]])
    i3_t_start = np.array ([d[1] for d in data if None not in d[1:]])
    i3_t_end = np.array ([d[2] for d in data if None not in d[1:]])


    def get_run_number (grb_start, grb_end):
        for i in xrange (len (i3_run_number)):
            if i3_t_start[i] is None or i3_t_end[i] is None:
                continue
            if i3_t_start[i] < grb_start and grb_end < i3_t_end[i]:
                return i3_run_number[i]
        return 0


    # now build GRBList from grbweb database
    grbs = GRBList (name=list_name)
    #db = MySQLdb.connect (
    #        host='dbs3.icecube.wisc.edu',
    #        user='grbweb-ro', passwd='Toshavmai', db='grbweb')
    db = pymysql.Connect (
            host='dbs3.icecube.wisc.edu',
            user='grbweb-ro', passwd='Toshavmai', db='grbweb')
    cursor = db.cursor ()
    db_cols = {}
    # colnames = (
    #         'GRBNAME', 'RA', 'DECL', 'ERR',
    #         'UTTIMET1', 'T1', 'T2', 'POS_TEXT',
    #         'FLUENCENU', 'EPSILON1', 'EPSILON2',
    #         'ALPHANU', 'BETANU', 'GAMMANU',
    #         'RUNNUMBER')
    # colnames = (
    #         'GRBNAME', 'RA', 'DECL', 'ERR',
    #         'UTTIMET1', 'T1', 'T2', 'POS_TEXT',
    #         'TFLUENCE', 'FLUENCE', 'ALPHA', 'BETA',
    #         'EPEAK', 'EMIN', 'EMAX', 'Z',
    #         'FLUENCENU', 'EPSILON1', 'EPSILON2',
    #         'ALPHANU', 'BETANU', 'GAMMANU',
    #         'RUNNUMBER')
    colnames = (
            'GRBNAME', 'RA', 'DECL', 'ERR',
            'UTTIMET1', 'T1', 'T2', 'POS_TEXT',
            'FLUENCENU', 'EPSILON1', 'EPSILON2',
            'ALPHANU', 'BETANU', 'GAMMANU',
            'RUNNUMBER')
    n_grb = -1
    for colname in colnames:
        cursor.execute (""" SELECT {0} FROM summary """.format (colname))
        data = cursor.fetchall ()
        col = np.array ([item[0] for item in data])
        col[col == '-'] = 'nan'
        db_cols[colname] = col
        if n_grb < 0:
            n_grb = len (col)
        else:
            assert (len (col) == n_grb)
    db_cols['RUNNUMBER'][db_cols['RUNNUMBER']=='nan'] = '0'
    for i in xrange (n_grb):
        name = db_cols['GRBNAME'][i]
        if name in grb_name_blacklist:
            continue
        ra = float (db_cols['RA'][i]) / 180. * pi
        dec = float (db_cols['DECL'][i]) / 180. * pi
        if np.isnan (ra) or np.isnan (dec):
            continue
        t_trigger = db_cols['UTTIMET1'][i]
        dt = float (db_cols['T1'][i])
        if t_trigger is None or np.isnan (dt):
            continue
        t_start = t_trigger + datetime.timedelta (
                microseconds=dt * 1e6)
        dt = float (db_cols['T2'][i])
        if t_trigger is None or np.isnan (dt):
            continue
        t_end = t_trigger + datetime.timedelta (
                microseconds=dt * 1e6)
        sigma = float (db_cols['ERR'][i]) / 180. * pi
        gbm_pos = ('GBM' in db_cols['POS_TEXT'][i])

        run_number = int (db_cols['RUNNUMBER'][i])
        if t_min is not None and t_start < t_min:
            continue
        if t_max is not None and t_end > t_max:
            continue

        # t100_flnc = float (db_cols['TFLUENCE'][i])
        # flnc = float (db_cols['FLUENCE'][i])
        # alpha = float (db_cols['ALPHA'][i])
        # beta = float (db_cols['BETA'][i])
        # epeak = float (db_cols['EPEAK'][i])
        # emin = float (db_cols['EMIN'][i])
        # emax = float (db_cols['EMAX'][i])
        # z = float (db_cols['Z'][i])

        # grb = GRB (name=name,
        #            t_start=t_start, t_end=t_end, ra=ra, dec=dec, sigma=sigma,
        #            run_number=run_number, gbm_pos=gbm_pos,
        #            t100_flnc=t100_fln, flnc=flnc, alpha=alpha, beta=beta, 
        #            epeak=epeak, emin=emin, emax=emax, z=z)
        grb = GRB (name=name,
                   t_start=t_start, t_end=t_end, ra=ra, dec=dec, sigma=sigma,
                   run_number=run_number, gbm_pos=gbm_pos)
        if zenith_min is not None and grb.zenith < zenith_min:
            continue
        if zenith_max is not None and grb.zenith > zenith_max:
            continue

        f = float (db_cols['FLUENCENU'][i])
        e1 = float (db_cols['EPSILON1'][i]) * 1e6
        e2 = float (db_cols['EPSILON2'][i]) * 1e6
        alpha = -float (db_cols['ALPHANU'][i])
        beta = -float (db_cols['BETANU'][i])
        gamma = -float (db_cols['GAMMANU'][i])

        grb.Guetta = Vars ()
        grb.Guetta.f = f
        grb.Guetta.e1 = e1
        grb.Guetta.e2 = e2
        grb.Guetta.alpha = alpha
        grb.Guetta.beta = beta
        grb.Guetta.gamma = gamma

        #@np.vectorize
        #def diff_fluence (E):
        #    if E < e1:
        #        return f * e1**-alpha * E**alpha
        #    elif e1 <= E < e2:
        #        return f * e1**-beta * E**beta
        #    else:
        #        return f * e1**-beta \
        #                * e2**(beta-gamma) \
        #                * E**gamma

        #grb.Guetta.diff_fluence = diff_fluence

        grbs.add (grb)
    return grbs

class Source:
    """ This class holds all the sources information.  Individual sources can be
    accessed through indexing or iteration.  The initial layout of this class
    originated in grbllh with some parts modified to better meet my current
    needs.
    """
    # FIXME: change zenith and azimuth to ra/dec
    def __init__ (self, ra, dec, sigma, duration, t, name):
        """Initialize Sources with arrays of event properties.

        :type   zenith: ndarray
        :param  zenith: Per-GRB zenith (radians).

        :type   azimuth: ndarray
        :param  azimuth: Per-GRB azimuth (radians).

        :type   sigma: ndarray
        :param  sigma: Per-GRB angular uncertainty (radians).

        :type   duration: ndarray
        :param  duration: Per-GRB start time (seconds).

        :type   t: ndarray
        :param  t: Per-GRB start time (datetime.datetime objects).
        """

        self.ra = ra.copy ()
        self.dec = dec.copy ()
        self.sigma = sigma.copy ()
        self.duration = duration.copy ()
        self.t = t.copy ()
        self.name = name.copy ()

    # The following methods should be fairly self explainitory
    def __iter__ (self):
        for i in xrange (len(self.ra)):
            yield Source (
                    self.ra[i], self.dec[i], self.sigma[i],
                    self.duration[i], self.t[i], self.name[i])

    def __len__ (self):
        return self.ra.size

    def __getitem__(self,grbid):
	if type(grbid) == int:
            i = grbid
        else:
            i = int(np.argwhere (self.name == grbid))

        return Source ( self.ra[i], self.dec[i], self.sigma[i],
                        self.duration[i], self.t[i], self.name[i])


def get_grbs_from_grbweb_file(grb_file_name='lowen_trans_grblist.txt',list_name='IC86IIGRBS',grb_name_blacklist=[],t_min=None,t_max=None,bad_run_list=[]):
    """ Load GRBs from grbweb text file. This should be done once and the output
    stored as a pickle file where the data is accessed 
    """
    grbs = GRBList (name=list_name)
    db_cols = {}

    infile = np.loadtxt(grb_file_name,dtype='string',delimiter=',')

    labels = infile[0]
    labels_old = labels.copy()

    labels = []
    for label in labels_old:
        labels.append(label.replace(' ',''))

    indata = infile[1:]

    colnames = (
            'GRBNAME', 'RA', 'DECL', 'ERR',
            'UTTIMET1', 'T1', 'T2', 'POS_TEXT',
            'FLUENCENU', 'EPSILON1', 'EPSILON2',
            'ALPHANU', 'BETANU', 'GAMMANU',
            'RUNNUMBER','Z')
    n_grb = -1

    for indy,col in enumerate(indata.transpose()):
	col[col == '-'] = 'nan'
	db_cols[labels[indy]] = col
	if n_grb < 0:
		n_grb = len(col)
	else:
		assert (len(col) == n_grb)
    

    db_cols['RunNumber'][db_cols['RunNumber']=='nan'] = '0'
    for i in xrange (n_grb):
        name = db_cols['Name'][i]
        if name in grb_name_blacklist:
            continue
        ra = float (db_cols['RA'][i]) / 180. * pi
        dec = float (db_cols['Decl'][i]) / 180. * pi
        if np.isnan (ra) or np.isnan (dec):
            continue
        t_trigger = db_cols['UT(Trigger)'][i]
	t_trigger_date = db_cols['Date'][i]
        dt = float (db_cols['T1'][i])
        if t_trigger is None or np.isnan (dt):
            continue
        fulltime = datetime.datetime(int(t_trigger_date[0:4]),int(t_trigger_date[5:7]),int(t_trigger_date[8:10]),int(t_trigger[0:2]),int(t_trigger[3:5]),int(t_trigger[6:8]))
        t_start = fulltime + datetime.timedelta (
                microseconds=dt * 1e6)
        dt = float (db_cols['T2'][i])
        if t_trigger is None or np.isnan (dt):
            continue
        t_end = fulltime + datetime.timedelta (
                microseconds=dt * 1e6)
        sigma = float (db_cols['ERR'][i]) / 180. * pi
	z = float(db_cols['z'][i])
        #gbm_pos = ('GBM' in db_cols['POS_TEXT'][i])
	gbm_pos = bool(True)

        run_number = int (db_cols['RunNumber'][i])
        if t_min is not None and t_start < t_min:
            continue
        if t_max is not None and t_end > t_max:
            continue

        grb = GRB (name=name,
                   t_start=t_start, t_end=t_end, ra=ra, dec=dec, sigma=sigma,
                   run_number=run_number, gbm_pos=gbm_pos, z=z)
        f = float (db_cols['Fluence'][i])
        e1 = float (db_cols['epsilon1'][i]) * 1e6
        e2 = float (db_cols['epsilon2'][i]) * 1e6
        alpha = -float (db_cols['AlphaNu'][i])
        beta = -float (db_cols['BetaNu'][i])
        gamma = -float (db_cols['GammaNu'][i])

        grb.Guetta = Vars ()
        grb.Guetta.f = f
        grb.Guetta.e1 = e1
        grb.Guetta.e2 = e2
        grb.Guetta.alpha = alpha
        grb.Guetta.beta = beta
        grb.Guetta.gamma = gamma

        #@np.vectorize
        #def diff_fluence (E):
        #    if E < e1:
        #        return f * e1**-alpha * E**alpha
        #    elif e1 <= E < e2:
        #        return f * e1**-beta * E**beta
        #    else:
        #        return f * e1**-beta \
        #                * e2**(beta-gamma) \
        #                * E**gamma

        #grb.Guetta.diff_fluence = diff_fluence

        grbs.add (grb)

    return grbs

def convert_date(in_date,in_time):
    """ Convert date from GRBweb file to datetime object """
    out_date = datetime.datetime(2000 + int(in_date[0:2]),
                                 int(in_date[2:4]),
                                 int(in_date[4:6]),
                                 int(in_time[0:2]),
                                 int(in_time[3:5]),
                                 int(in_time[6:8]))
    return out_date


