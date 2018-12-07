"""
coordCenterRotation_offCenterSphere-cylindrical.py

A regression test for globAccDisk-src-rotating_barycenter.hpp (source terms for a rotating frame of reference). In this example, GM = 1 is set, causing the rotation to occur around the origin of the coordinate system.

Current test:
 -- a sphere off-center wrt the rotation axis
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

    arguments = ['job/problem_id=coordCenterRotation_offCenterSphere-cylindrical',
                 'output1/file_type=hdf5',
                 'output1/variable=prim',
                 'output1/dt=%.20f' % (0.5*np.pi),
                 'time/cfl_number=0.3',
                 'time/tlim=%.20f' % (2.*np.pi), # for longer times the pressure artifacts from rotating over boundary conditions distort the sphere
                 'mesh/nx1=32',
                 'mesh/x1min=0.02',
                 'mesh/x1max=0.62',
                 'mesh/ix1_bc=reflecting',
                 'mesh/ox1_bc=reflecting',
                 'mesh/nx2=128',
                 'mesh/x2min=0.0',
                 'mesh/x2max=%.20f' % (2.*np.pi),
                 'mesh/ix2_bc=periodic',
                 'mesh/ox2_bc=periodic',
                 'mesh/nx3=8',
                 'mesh/x3min=-0.3',
                 'mesh/x3max=0.3',
                 'mesh/ix3_bc=reflecting',
                 'mesh/ox3_bc=reflecting',
                 'meshblock/nx1=32',
                 'meshblock/nx2=32',
                 'meshblock/nx3=8',
                 'hydro/gamma=1.1',
                 'problem/GM1=1.0',
                 'problem/rho0=1.0',
                 'problem/rho1=1.1',
                 'problem/press=0.1',
                 'problem/obj_radius=0.3',
                 'problem/obj_phi0=0.0',
                 'problem/obj_size=0.2']

    athena.run('hydro/athinput.rotating_frame', arguments)

def in_sphere (r, phi, z, time):
    dist = np.sqrt(r*r + 0.3*0.3 - 2.*r*0.3*np.cos(phi-time) + z*z)
    return dist < 0.1

def analyze():

    initial_state = athena_read.athdf('bin/coordCenterRotation_offCenterSphere-cylindrical.out1.00000.athdf')
    r0, phi0, z0 = [initial_state[key] for key in ['x1v', 'x2v', 'x3v']]
    z, phi, r = np.meshgrid(z0, phi0, r0, indexing='ij')

    # check if the frames are identical
    for i in range(1,5):
        test_results = athena_read.athdf('bin/coordCenterRotation_offCenterSphere-cylindrical.out1.0000%i.athdf' % i)
        rho = test_results['rho']
        time = test_results['Time']

        rho_expected = np.ones(rho.shape) * 1.
        rho_expected[np.where(np.sqrt(r*r + 0.3*0.3 - 2.*r*0.3*np.cos(phi+time) + z*z) < 0.2)] = 1.1

        if False: # useful if the test fails
            import matplotlib.pyplot as plt
            plt.clf()
            plt.subplot(121)
            plt.contourf(r0, phi0, rho_expected[3])
            plt.subplot(122)
            plt.contourf(r0, phi0, rho[3])
            plt.show()

        error_rel_rho = np.sum(np.abs(rho_expected - rho)) / np.sum(rho_expected)
        print("[coordCenterRotation_offCenterSphere-cylindrical]: error_rel_rho = %.2e" % error_rel_rho)

        if np.isnan(error_rel_rho):
            return False
        if error_rel_rho > 0.01:
            return False

    return True
