"""
barycenterRotation_centeredSphere-cylindrical.py

A regression test for globAccDisk-src-rotating_barycenter.hpp (source terms for a rotating frame of reference). In this example, GM < 1 is set, causing the rotation to occur around a point off the origin of the coordinate system.

Current test:
 -- a sphere centered around the rotation axis remains unchanged
 -- cylindrical coordinate version

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
                     prob='globAccDisk-test-rotating_barycenter',
                     flux='roe',
                     eos='adiabatic',
                     coord='cylindrical',
                     hdf5_path='/usr/local/Cellar/hdf5/1.10.4',
                     **kwargs)
    elif socket.gethostname()[:7] == 'perseus':
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-rotating_barycenter',
                     flux='roe',
                     eos='adiabatic',
                     coord='cylindrical',
                     cxx='icc',
                     **kwargs)
    else:
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-rotating_barycenter',
                     flux='roe',
                     eos='adiabatic',
                     coord='cylindrical',
                     **kwargs)


    athena.make()


def run(**kwargs):

    arguments = ['job/problem_id=barycenterRotation_centeredSphere-cylindrical',
                 'output1/file_type=hdf5',
                 'output1/variable=prim',
                 'output1/dt=%.20f' % (1.0),
                 'time/cfl_number=0.3',
                 'time/tlim=%.20f' % (1.0), # for longer times the pressure artifacts from rotating over boundary conditions distort the sphere
                 'mesh/nx1=64',
                 'mesh/x1min=0.02',
                 'mesh/x1max=0.62',
                 'mesh/ix1_bc=user',
                 'mesh/ox1_bc=user',
                 'mesh/nx2=256',
                 'mesh/x2min=0.0',
                 'mesh/x2max=%.20f' % (2.*np.pi),
                 'mesh/ix2_bc=periodic',
                 'mesh/ox2_bc=periodic',
                 'mesh/nx3=8',
                 'mesh/x3min=-0.15',
                 'mesh/x3max=0.15',
                 'mesh/ix3_bc=reflecting',
                 'mesh/ox3_bc=reflecting',
                 'meshblock/nx1=32',
                 'meshblock/nx2=32',
                 'meshblock/nx3=8',
                 'hydro/gamma=1.1',
                 'problem/GM1=0.7',
                 'problem/rho0=100.0',
                 'problem/rho1=110.0',
                 'problem/press=0.1',
                 'problem/obj_radius=0.3',
                 'problem/obj_phi0=0.0',
                 'problem/obj_size=0.1']

    athena.run('hydro/athinput.rotating_frame', arguments)

def analyze():

    initial_state = athena_read.athdf('bin/barycenterRotation_centeredSphere-cylindrical.out1.00000.athdf')
    r, rho_expected = [initial_state[key] for key in ['x1f', 'rho']]

    # check if the frames are identical
    for i in range(1,2):
        test_results = athena_read.athdf('bin/barycenterRotation_centeredSphere-cylindrical.out1.0000%i.athdf' % i)
        rho = test_results['rho']

        error_rel_rho = np.sum(np.abs(rho_expected - rho)) / np.sum(rho)
        print("[barycenterRotation_centeredSphere-cylindrical]: error_rel_rho = %.2e" % error_rel_rho)

        if np.isnan(error_rel_rho):
            return False
        if error_rel_rho > 0.01:
            return False

    return True
