#!/usr/bin/env python
# setup_sample_validator.py

from __future__ import print_function

__doc__ = """Setup a pybdt.validate.Validator for the sample BDT."""

from pybdt.util import save
from pybdt.validate import Validator

print ('Creating Validator...')
v = Validator ('output/fiveparam_final_hes.bdt')

def add_data (name, label):
    filename = 'hesdata/{0}.ds'.format (name)
    print ('Adding {0} from {1} ...'.format (label, filename))
    v.add_data (name, filename, label=label)

add_data ('train_sig', 'Training signal sim')
add_data ('train_data', 'Training data')
add_data ('test_sig', 'Testing signal sim')
add_data ('test_data', 'Testing data')
add_data ('bg', 'Background sim')

print ('Configuring weighting...')

v.add_weighting ('weight', 'train_sig',
        color='cyan',
        )
v.add_weighting ('livetime', 'train_data',
        line=False, markers=True, marker='.', color='.5',
        errorbars=True,
        )

v.add_weighting ('weight', 'test_sig',
        color='blue',
        add_to_mc=True,
        )
v.add_weighting ('livetime', 'test_data',
        line=False, markers=True, marker='.', color='black',
        errorbars=True,
        use_as_data=True,
        )

v.setup_total_mc (
        color='green',
        )

v.add_weighting ('livetime', 'bg',
        color='purple',
        add_to_mc=True,
        )

save (v, 'output/fiveparam_final_hes.validator')
