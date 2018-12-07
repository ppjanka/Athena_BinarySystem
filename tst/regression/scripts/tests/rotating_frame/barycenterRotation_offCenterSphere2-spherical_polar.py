"""
barycenterRotation_offCenterSphere2-spherical_polar.py

A regression test for globAccDisk-src-rotating_barycenter.hpp (source terms for a rotating frame of reference). In this example, GM < 1 is set, causing the rotation to occur around a point off the origin of the coordinate system.

Current test:
 -- a sphere off-center wrt the rotation axis
 -- spherical_polar coordinate version

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
                     coord='spherical_polar',
                     hdf5_path='/usr/local/Cellar/hdf5/1.10.4',
                     **kwargs)
    elif socket.gethostname()[:7] == 'perseus':
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-rotating_barycenter',
                     flux='roe',
                     eos='adiabatic',
                     coord='spherical_polar',
                     cxx='icc',
                     **kwargs)
    else:
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-rotating_barycenter',
                     flux='roe',
                     eos='adiabatic',
                     coord='spherical_polar',
                     **kwargs)


    athena.make()


def run(**kwargs):

    arguments = ['job/problem_id=barycenterRotation_offCenterSphere2-spherical_polar',
                 'output1/file_type=hdf5',
                 'output1/variable=prim',
                 'output1/dt=%.20f' % (0.25),
                 'time/cfl_number=0.3',
                 'time/tlim=%.20f' % (1.0), # for longer times the pressure artifacts from rotating over boundary conditions distort the sphere
                 'mesh/nx1=32',
                 'mesh/x1min=0.02',
                 'mesh/x1max=0.62',
                 'mesh/ix1_bc=user',
                 'mesh/ox1_bc=user',
                 'mesh/nx2=8',
                 'mesh/x2min=%.20f' % (0.25*np.pi),
                 'mesh/x2max=%.20f' % (0.75*np.pi),
                 'mesh/ix2_bc=user',
                 'mesh/ox2_bc=user',
                 'mesh/nx3=128',
                 'mesh/x3min=0.0',
                 'mesh/x3max=%.20f' % (2.*np.pi),
                 'mesh/ix3_bc=periodic',
                 'mesh/ox3_bc=periodic',
                 'meshblock/nx1=32',
                 'meshblock/nx2=8',
                 'meshblock/nx3=32',
                 'hydro/gamma=1.1',
                 'problem/GM1=0.7',
                 'problem/rho0=1.0',
                 'problem/rho1=2.0',
                 'problem/press=0.01',
                 'problem/obj_radius=0.45',
                 'problem/obj_phi0=%.20f' % (1.3333333333 * np.pi),
                 'problem/obj_size=0.1']

    athena.run('hydro/athinput.rotating_frame', arguments)

def in_sphere (r, phi, z, time, r_init, phi_init, r_bary):
    # find the position of the sphere's center

    # transform initial sphere position to barycentered coordinates
    r_init_bary = np.sqrt(r_init**2 + r_bary**2 - 2.*r_init*r_bary*np.cos(phi_init))
    phi_init_bary = np.arccos(1. - (r_init*np.sin(phi_init) / r_init_bary)**2)
    if phi_init > np.pi:
        phi_init_bary = 2.*np.pi - phi_init_bary
    # advect the sphere center
    phi0_bary = phi_init_bary - time
    # transform back to the LAB coordinates
    r_exp = np.sqrt(r_init_bary**2 + r_bary**2 + 2.*r_init_bary*r_bary*np.cos(phi0_bary))
    phi_exp = np.arccos(np.sqrt(1. - (r_init_bary*np.sin(phi0_bary)/r_exp)**2))
    if phi0_bary > np.pi:
        phi0_bary = 2.*np.pi - phi0_bary
    # decide if we're inside the advected sphere
    return np.sqrt(r**2 + r_exp**2 - 2.*r*r_exp*np.cos(phi-phi_exp) + z*z) < 0.1


def analyze():

    initial_state = athena_read.athdf('bin/barycenterRotation_offCenterSphere2-spherical_polar.out1.00000.athdf')
    r0, theta0, phi0 = [initial_state[key] for key in ['x1v', 'x2v', 'x3v']]
    phi, theta, r = np.meshgrid(phi0, theta0, r0, indexing='ij')

    # transform initial sphere position to barycentered coordinates
    r_bary = 0.3
    r_init = 0.45
    phi_init = 1.3333333333 * np.pi
    r_init_bary = np.sqrt(r_init**2 + r_bary**2 - 2.*r_init*r_bary*np.cos(phi_init))
    phi_init_bary = np.arccos(0.5*(r_init**2 - r_bary**2 - r_init_bary**2) / (r_bary*r_init_bary))
    if phi_init > np.pi:
        phi_init_bary = 2.*np.pi - phi_init_bary

    # check if the sphere is correctly advected
    for i in range(1,5):
        test_results = athena_read.athdf('bin/barycenterRotation_offCenterSphere2-spherical_polar.out1.0000%i.athdf' % i)
        rho = test_results['rho']
        time = test_results['Time']

        # find the position of the sphere's center
        # advect the sphere center
        phi0_bary = phi_init_bary - time
        # transform back to the LAB coordinates
        r_exp = np.sqrt(r_init_bary**2 + r_bary**2 + 2.*r_init_bary*r_bary*np.cos(phi0_bary))
        phi_exp = np.arccos(0.5*(- r_init_bary**2 + r_exp**2 + r_bary**2) / (r_exp*r_bary))
        if phi0_bary > np.pi:
            phi_exp = 2.*np.pi - phi_exp
        # generate an image of the rotated sphere
        rho_expected = np.ones(rho.shape) * 1.
        rho_expected[np.where(np.sqrt(r**2 + r_exp**2 - 2.*r*np.sin(theta)*r_exp*np.cos(phi-phi_exp)) < 0.1)] = 2.5

        if False: # useful if the test fails
            import matplotlib.pyplot as plt
            plt.clf()
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            plt.subplot(121)
            plt.contourf(r0, phi0, rho_expected[:,3,:])
            plt.subplot(122)
            plt.contourf(r0, phi0, rho[:,3,:])
            plt.show()

        error_rel_rho = np.sum(np.abs(rho_expected - rho)) / np.sum(rho)
        print("[barycenterRotation_offCenterSphere2-spherical_polar]: error_rel_rho = %.2e" % error_rel_rho)

        if np.isnan(error_rel_rho):
            return False
        if error_rel_rho > 0.05:
            return False

    return True
