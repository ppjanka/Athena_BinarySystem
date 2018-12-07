athena
======
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Athena++ radiation MHD code

======

Additions: rotating_frame.
Comments for the rotation src terms:

The source terms are handled in src/pgen/globAccDisk-src-rotating_barycenter.hpp. There's a simple testing problem generator in src/pgen/globAccDisk-test-rotating_barycenter.cpp (it creates a sphere and rotates the frame to see if it behaves as expected) and an example input file for it in inputs/hydro/athinput.rotating_frame.

A test suite to ensure the rotating frame is working correctly is located in tst/regression/scripts/tests/rotating_frame/. See Athena++ main page wiki for info on how to run them. Before running, the location of hdf5 installation must be changed in the regression scripts. Also, since my h5py has trouble with the new vis/python/athena_read.py, I am using the old version (included in the repo).

The implementation assumes rotation around the barycenter of a binary system with G(M1+M2) = 1, where GM1 is set in pmb->ruser_meshblock_data[0](0) by the problem generator, and the frame is centered on the GM1 object. Note that G(M1+M2) = 1 results in the orbital angular frequency of 1. To avoid multiplication by 1 and SQR(1) in the code, the orbital frequency is hard-coded to 1 (i.e., omitted), to improve performance. I tried to comment it as I went (hence all the /*omega_frame*/), so it should be easy to put the frame rotation frequency back in.
