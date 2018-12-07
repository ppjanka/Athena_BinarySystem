"""
unstratified-secondary-cylindrical.py

A regression test for globAccDisk-src-primary_grav.hpp (source terms for gravity with a source away from the coordinate center).

Current test:
 -- stratified gravity
 -- cylindrical coordinate version
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

    arguments = ['job/problem_id=stratified-secondary-cylindrical',
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
                 'mesh/x2min=%.20f' % (0.0),
                 'mesh/x2max=%.20f' % (0.75*np.pi),
                 'mesh/ix2_bc=outflow',
                 'mesh/ox2_bc=outflow',
                 'mesh/nx3=32',
                 'mesh/x3min=0.0',
                 'mesh/x3max=0.5',
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

    initial_state = athena_read.athdf('bin/stratified-secondary-cylindrical.out1.00000.athdf')
    r0, phi0, z0, rho = [initial_state[key] for key in ['x1v', 'x2v', 'x3v', 'rho']]
    dr, dphi, dz = (r0[1]-r0[0]), (phi0[1]-phi0[0]), (z0[1]-z0[0])
    z, phi, r = np.meshgrid(z0, phi0, r0, indexing='ij')
    vol = dr * np.array(r) * dphi * dz
    x = np.array(r) * np.cos(np.array(phi))
    y = np.array(r) * np.sin(np.array(phi))

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
        test_results = athena_read.athdf('bin/stratified-secondary-cylindrical.out1.%05i.athdf' % i)
        rho = test_results['rho']
        time = test_results['Time']

        # calculate the current center of mass
        com_cartesian = np.array([np.sum(rho*vol*x), np.sum(rho*vol*y), np.sum(rho*vol*z)]) / np.sum(rho*vol)

        if visualize or movie_dir != '':
            com_r = np.sqrt(com_cartesian[0]**2 + com_cartesian[1]**2)
            com_phi = np.arccos(com_cartesian[0] / com_r)
            if com_cartesian[1] < 0:
                com_phi = 2.*np.pi - com_phi
            com = np.array([com_r, com_phi, com_cartesian[2]])
            comm.append(com)
            times.append(time)

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

        if visualize or movie_dir != '':
            com_r = np.sqrt(com_curr[0]**2 + com_curr[1]**2)
            com_phi = np.arccos(com_curr[0] / com_r)
            if com_curr[1] < 0:
                com_phi = 2.*np.pi - com_phi
            com_expected = np.array([com_r, com_phi, com_curr[2]])
            com_expecteds.append(com_expected)
            print(com, com_expected)

        error_tot_com = np.sqrt(np.sum((com_cartesian - com_curr)**2))
        print("[stratified-secondary-cylindrical]: error_tot_com = %.2e (max %.2e)" % (error_tot_com, 0.05*time))

        if np.isnan(error_tot_com):
            return False
        if error_tot_com > 0.05 * time:
            return False

        if movie_dir != '':
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(121, projection='polar')
            plt.contour(phi0, r0, np.sum(rho, axis=0).transpose())
            ax.scatter([com[1],], [com[0],], color='r', label='actual')
            ax.scatter([com_expected[1],], [com_expected[0],], color='g', label='expected')
            plt.ylim(0.0, 1.0)
            plt.legend()
            ax = fig.add_subplot(122)
            ax.scatter([time,], [com[2],], color='r', label='actual')
            ax.scatter([time,], [com_expected[2],], color='g', label='expected')
            plt.savefig("%s/%05i.png" % (movie_dir, i))
            plt.close()

    if visualize:
        comm = np.array(comm).transpose()
        com_expecteds = np.array(com_expecteds).transpose()

        test_results = athena_read.athdf('bin/stratified-secondary-cylindrical.out1.00025.athdf')
        rho = test_results['rho']
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(121, projection='polar')
        plt.contour(phi0, r0, np.sum(rho, axis=0).transpose())
        ax.scatter(comm[1], comm[0], color='r', label='actual')
        ax.scatter(com_expecteds[1], com_expecteds[0], color='g', label='expected')
        plt.ylim(0.0, 1.0)
        plt.legend()
        ax = fig.add_subplot(122)
        ax.scatter(times, comm[2], color='r', label='actual')
        ax.scatter(times, com_expecteds[2], color='g', label='expected')
        plt.show()

    return True
