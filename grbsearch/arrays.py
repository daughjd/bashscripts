# cutvars.py


__doc__ = """Provide analysis-level objects holding equal-length arrays."""

import copy

from vars_class import Vars

import numpy as np

class Arrays (Vars):

    def __init__ (self, initializer={}):
        Vars.__init__ (self, initializer)
        self.meta = Vars ()

    def __len__ (self):
        for key in self:
            if key != 'meta':
                return len (self[key])

    def apply_cut (self, i):
        for k, a in self._dict.iteritems ():
            if k != 'meta':
                self[k] = a[i]

    def get_copy (self, i=None):
        out = self.__class__ ()
        for k, a in self._dict.iteritems ():
            if k == 'meta' and not a._dict:
                out[k] = Vars ()
            elif k == 'meta' or i is None:
                out[k] = copy.deepcopy (a)
            else:
                out[k] = a[i]
        return out


def combined_arrays (arrays):
    """Combine Arrays objects.

    :type   arrays: list of :class:`Arrays`
    :param  arrays: The Arrays objects.

    It is up to the user to ensure that all the Arrays objects are compatible!
    
    """
    a1 = arrays[0]
    out = a1.__class__ ()
    for k in a1._dict.iterkeys ():
        if k == 'meta':
            continue
        out[k] = np.copy (np.concatenate (tuple (a[k] for a in arrays)))
    out.meta = Vars ()
    return out
