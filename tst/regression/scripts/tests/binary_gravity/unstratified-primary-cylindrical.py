"""
unstratified-primary-cylindrical.py

A regression test for globAccDisk-src-primary_grav.hpp (source terms for gravity with a source at the coordinate center).

Current test:
 -- unstratified gravity (cylindrical-radial component only)
 -- cylindrical coordinate version
 -- a Keplerian disk is initialized, the test ensures that it remains stationary

"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data # noqa
import socket                                  # recognizing host machine


def prepare(**kwargs):

    # check which machine we're running on and configure accordingly
    if socket.gethostname() == 'ast1506-astro':
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-binary_gravity',
                     flux='roe',
                     eos='adiabatic',
                     coord='cylindrical',
                     hdf5_path='/usr/local/Cellar/hdf5/1.10.4',
                     **kwargs)
    else:
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-binary_gravity',
                     flux='roe',
                     eos='adiabatic',
                     coord='cylindrical',
                     **kwargs)


    athena.make()


def run(**kwargs):

    arguments = ['job/problem_id=unstratified-primary-cylindrical',
                 'output1/file_type=hdf5',
                 'output1/variable=prim',
                 'output1/dt=%.20f' % (2.5),
                 'time/cfl_number=0.3',
                 'time/tlim=%.20f' % (2.5),
                 'mesh/nx1=32',
                 'mesh/x1min=0.5',
                 'mesh/x1max=1.0',
                 'mesh/ix1_bc=reflecting',
                 'mesh/ox1_bc=reflecting',
                 'mesh/nx2=128',
                 'mesh/x2min=0.0',
                 'mesh/x2max=%.20f' % (2.*np.pi),
                 'mesh/ix2_bc=periodic',
                 'mesh/ox2_bc=periodic',
                 'mesh/nx3=32',
                 'mesh/x3min=-0.5',
                 'mesh/x3max=0.5',
                 'mesh/ix3_bc=reflecting',
                 'mesh/ox3_bc=reflecting',
                 'meshblock/nx1=16',
                 'meshblock/nx2=16',
                 'meshblock/nx3=16',
                 'hydro/gamma=1.1',
                 'hydro/dfloor=1.0e-3',
                 'hydro/pfloor=1.0e-3',
                 'problem/stratified=0',
                 'problem/binary_component=0',
                 'problem/GM1=1.0',
                 'problem/rho=1.0']

    athena.run('hydro/athinput.binary_gravity', arguments)

def analyze():

    initial_state = athena_read.athdf('bin/unstratified-primary-cylindrical.out1.00000.athdf')
    r, vphi_expected = [initial_state[key] for key in ['x1f', 'vel2']]

    # check if the frames are identical
    for i in range(1,2):
        test_results = athena_read.athdf('bin/unstratified-primary-cylindrical.out1.0000%i.athdf' % i)
        vphi = test_results['vel2']

        error_rel_vphi = np.sum(np.abs(vphi_expected - vphi)) / np.sum(vphi)
        print("[unstratified-primary-cylindrical]: error_rel_vphi = %.2e" % error_rel_vphi)

        if np.isnan(error_rel_vphi):
            return False
        if error_rel_vphi > 0.01:
            return False

    return True
