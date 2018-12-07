"""
unstratified-secondary-spherical_polar.py

A regression test for globAccDisk-src-primary_grav.hpp (source terms for gravity with a source away from the coordinate center).

Current test:
 -- stratified gravity
 -- spherical_polar coordinate version
 -- an immobile sphere is initialized, the test ensures its center of mass is properly accelerated towards the secondary body

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
                     coord='spherical_polar',
                     hdf5_path='/usr/local/Cellar/hdf5/1.10.4',
                     **kwargs)
    else:
        athena.configure('hdf5', 'mpi',
                     prob='globAccDisk-test-binary_gravity',
                     flux='roe',
                     eos='adiabatic',
                     coord='spherical_polar',
                     **kwargs)


    athena.make()


def run(**kwargs):

    arguments = ['job/problem_id=stratified-secondary-spherical_polar',
                 'output1/file_type=hdf5',
                 'output1/variable=prim',
                 'output1/ghost_zones=false',
                 'output1/dt=%.20f' % (0.01),
                 'time/cfl_number=0.3',
                 'time/tlim=%.20f' % (1.0),
                 'mesh/nx1=32',
                 'mesh/x1min=0.1',
                 'mesh/x1max=0.9',
                 'mesh/ix1_bc=outflow',
                 'mesh/ox1_bc=outflow',
                 'mesh/nx2=32',
                 'mesh/x2min=%.20f' % (0.25*np.pi),
                 'mesh/x2max=%.20f' % (0.5*np.pi),
                 'mesh/ix2_bc=outflow',
                 'mesh/ox2_bc=outflow',
                 'mesh/nx3=32',
                 'mesh/x3min=%.20f' % (0.0),
                 'mesh/x3max=%.20f' % (0.75*np.pi),
                 'mesh/ix3_bc=outflow',
                 'mesh/ox3_bc=outflow',
                 'meshblock/nx1=16',
                 'meshblock/nx2=16',
                 'meshblock/nx3=16',
                 'hydro/gamma=1.1',
                 'hydro/dfloor=1.0e-3',
                 'hydro/pfloor=1.0e-3',
                 'problem/stratified=1',
                 'problem/binary_component=1',
                 'problem/GM1=0.5',
                 'problem/rho=10.0',
                 'problem/rho0=1.0e-3',
                 'problem/press=1.0e-3',
                 'problem/obj_r0=0.6',
                 'problem/obj_phi0=1.5',
                 'problem/obj_z0=0.25',
                 'problem/obj_size=0.15']

    athena.run('hydro/athinput.binary_gravity', arguments)

visualize = False
movie_dir = ''
def analyze():

    initial_state = athena_read.athdf('bin/stratified-secondary-spherical_polar.out1.00000.athdf')
    r0, theta0, phi0, rho = [initial_state[key] for key in ['x1v', 'x2v', 'x3v', 'rho']]
    dr, dtheta, dphi = (r0[1]-r0[0]), (theta0[1]-theta0[0]), (phi0[1]-phi0[0])
    phi, theta, r = np.meshgrid(phi0, theta0, r0, indexing='ij')
    vol = dr * np.array(r) * dtheta * np.array(r) * np.sin(np.array(theta)) * dphi
    z = np.array(r) * np.cos(np.array(theta))
    x = np.array(r) * np.sin(np.array(theta)) * np.cos(np.array(phi))
    y = np.array(r) * np.sin(np.array(theta)) * np.sin(np.array(phi))

    #calculate the center of mass
    com_cartesian = np.array([np.sum(rho*vol*x), np.sum(rho*vol*y), np.sum(rho*vol*z)]) / np.sum(rho*vol)
    com_initial = 1.0 * com_cartesian

    # initialize integration of the CoM motion
    GM = 0.5
    time_curr = 0.
    com_curr = 1.*com_initial
    vtot = 0.

    if visualize:
        comm = []
        com_expecteds = []
        times = []
    # check if the sphere is correctly accelerated
    import matplotlib.pyplot as plt
    for i in range(1,101):
        test_results = athena_read.athdf('bin/stratified-secondary-spherical_polar.out1.%05i.athdf' % i)
        rho = test_results['rho']
        time = test_results['Time']

        # calculate the current center of mass
        com_cartesian = np.array([np.sum(rho*vol*x), np.sum(rho*vol*y), np.sum(rho*vol*z)]) / np.sum(rho*vol)

        # calculate the current acceleration - needs to take into account the extension of the object
        acc = np.sum(rho * vol * (1.-GM) / ((1.-x)**2 + y**2 + z**2)) / np.sum(rho*vol)

        dt = time - time_curr

        # calculate where the CoM should be
        rtot = np.sqrt((1.-com_curr[0])**2 + com_curr[1]**2 + com_curr[2]**2)
        dr = vtot * dt + 0.5 * acc * dt**2
        com_curr[0] += (1.-com_curr[0]) * dr / rtot
        com_curr[1] += - com_curr[1] * dr / rtot
        com_curr[2] += - com_curr[2] * dr / rtot

        vtot += acc * dt
        time_curr += dt

        error_tot_com = np.sqrt(np.sum((com_cartesian - com_curr)**2))
        print("[stratified-secondary-spherical_polar]: error_tot_com = %.2e (max %.2e)" % (error_tot_com, 0.01*time))

        if np.isnan(error_tot_com):
            return False
        if error_tot_com > 0.01 * time:
            return False

    return True
