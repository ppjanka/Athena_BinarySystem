athena
======
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Athena++ radiation MHD code

======

Additions: rotating_frame, binary gravity.

------

Comments regarding the rotation src terms:

The source terms are handled in src/pgen/globAccDisk-src-rotating_barycenter.hpp. There's a simple testing problem generator in src/pgen/globAccDisk-test-rotating_barycenter.cpp (it creates a sphere and rotates the frame to see if it behaves as expected) and an example input file for it in inputs/hydro/athinput.rotating_frame.

A test suite to ensure the rotating frame is working correctly is located in tst/regression/scripts/tests/rotating_frame/. See Athena++ main page wiki for info on how to run them. Before running, the location of hdf5 installation must be changed in the regression scripts. Also, since my h5py has trouble with the new vis/python/athena_read.py, I am using the old version (included in the repo).

The implementation assumes rotation around the barycenter of a binary system with G(M1+M2) = 1, where GM1 is set in pmb->ruser_meshblock_data[0]\(0\) by the problem generator, and the frame is centered on the GM1 object. Note that G(M1+M2) = 1 results in the orbital angular frequency of 1. To avoid multiplication by 1 and SQR(1) in the code, the orbital frequency is hard-coded to 1 (i.e., omitted), to improve performance. I tried to comment it as I went (hence all the /*omega_frame*/), so it should be easy to put the frame rotation frequency back in.

------

Comments concerning binary gravity:

The gravity source terms for the primary (GM1) and the secondary (GM2) component of the binary are handled in src/pgen/globAccDisk-src-primary_grav.hpp and src/pgen/globAccDisk-src-secondary_grav.hpp, respectively. There is a simple testing problem generator in src/pgen/globAccDisk-test-binary_gravity.cpp with an input file in inputs/hydro/athinput.binary_gravity, and a testing suite tst/regression/scripts/tests/binary_gravity (see HDF5 caveat above).

Again, G(M1+M2) = 1, where GM1 is set in pmb->ruser_meshblock_data[0]\(0\).

------

Note:

 - One can only enroll a single function with EnrollUserExplicitSourceFunction(). If there are multiple functions handling source terms (e.g., two for primary and secondary gravity, one for rotating frame), they need to be called from a wrapper function, and then EnrollUserExplicitSourceFunction() should be given that wrapper function as argument.
