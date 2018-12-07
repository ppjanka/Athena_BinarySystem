//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globAccDisk-test-binary_gravity
//  \brief Unit test for the binary system's gravitational source terms. The test initializes a disk in Keplerian rotation; Patryk Pjanka

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"

// gravitation source terms
#include "globAccDisk-src-primary_grav.hpp"
#include "globAccDisk-src-secondary_grav.hpp"

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void __attribute__((weak)) Mesh::InitUserMeshData(ParameterInput *pin) {
  
  if (pin->GetInteger("problem", "binary_component") == 0) { // the primary
    if (pin->GetInteger("problem", "stratified") == 0)
      EnrollUserExplicitSourceFunction(PrimarySourceGrav_unstratified);
    else if (pin->GetInteger("problem", "stratified") == 1) {
      if (COORDINATE_SYSTEM == "cylindrical")
        EnrollUserExplicitSourceFunction(PrimarySourceGrav_stratified_cylindrical);
      else if (COORDINATE_SYSTEM == "spherical_polar")
        EnrollUserExplicitSourceFunction(PrimarySourceGrav_stratified_spherical);
    }
  } else if (pin->GetInteger("problem", "binary_component") == 1) { // the secondary
    if (pin->GetInteger("problem", "stratified") == 0)
      EnrollUserExplicitSourceFunction(SecondarySourceGrav_unstratified);
    else if (pin->GetInteger("problem", "stratified") == 1) {
      if (COORDINATE_SYSTEM == "cylindrical")
        EnrollUserExplicitSourceFunction(SecondarySourceGrav_stratified_cylindrical);
      else if (COORDINATE_SYSTEM == "spherical_polar")
        EnrollUserExplicitSourceFunction(SecondarySourceGrav_stratified_spherical);
    }
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//  used to initialize variables which are global to other functions in this file.
//  Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================

void __attribute__((weak)) MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(1);
  ruser_meshblock_data[0](0) = pin->GetReal("problem", "GM1");

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//========================================================================================

void __attribute__((weak)) MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real GM1 = pin->GetReal("problem", "GM1");
  Real GM2 = 1. - GM1;
  Real rho = pin->GetReal("problem", "rho"); // disk/ring density
  Real rho0 = pin->GetReal("problem", "rho0"); // ambient medium density
  #if NON_BAROTROPIC_EOS == 1
  Real press = pin->GetReal("problem", "press"); // pressure in the box
  Real adiab_idx = pin->GetReal("hydro", "gamma"); // adiabatic idx
  #endif
  // specify orientation of the inclined disk in hte stratified case
  Real position_angle, inclination;
  Real sin_incl, cos_incl, sin_pa, cos_pa;
  Real sin_ang_dist_max;
  if (pin->GetInteger("problem", "stratified") == 1
      && pin->GetInteger("problem", "binary_component") == 0) {
    inclination = pin->GetReal("problem", "inclination");
    sin_incl = sin(inclination); cos_incl = cos(inclination);
    position_angle =  pin->GetReal("problem", "position_angle");
    sin_pa = sin(position_angle); cos_pa = cos(position_angle);
    sin_ang_dist_max = sin(pin->GetReal("problem", "ang_dist_max"));
  }
  Real obj_r0, obj_phi0, obj_z0, obj_size, sqr_obj_size;
  if (pin->GetInteger("problem", "binary_component") == 1) {  
    obj_r0 = pin->GetReal("problem", "obj_r0"); // radius of the center of the object
    obj_phi0 = pin->GetOrAddReal("problem", "obj_phi0", 0.0); // initial azimuth of the center of the object
    obj_z0 = pin->GetOrAddReal("problem", "obj_z0", 0.0);
    obj_size = pin->GetReal("problem", "obj_size");
    sqr_obj_size = SQR(obj_size);
  }
  
  int i,j,k;
  Real x1, x2, x3, r_sph, r_cyl, z;
  Real phi, sin_phi, cos_phi;
  Real sin_theta, cos_theta;
  Real sin_ang_dist; // angular distance from the ring
  Real sin_xi, cos_xi; // angle between vec(v) and the direction to the Zenith
  Real mom_Kepl;
  Real cos_phi_obj, sqr_dist_obj;

  if (pin->GetInteger("problem", "binary_component") == 0) { // DISK AROUND THE PRIMARY (secondary gravity OFF)
    for (k=ks; k<=ke; k++) {   //z / phi
      x3 = pcoord->x3v(k);
      if (COORDINATE_SYSTEM == "spherical_polar") {
        phi = x3;
        cos_phi = cos(x3); sin_phi = sin(x3);
      }
      for (j=js; j<=je; j++) {  //phi / theta
        x2 = pcoord->x2v(j);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          sin_theta = sin(x2); cos_theta = cos(x2);
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          phi = x2;
          cos_phi = cos(x2); sin_phi = sin(x2);
        }
        #pragma omp simd
        for (i=is; i<=ie; i++) { //r_cyl / r
          x1 = pcoord->x1v(i);
          if (COORDINATE_SYSTEM == "spherical_polar") {
            r_sph = x1;
          }
          else if (COORDINATE_SYSTEM == "cylindrical") {
            r_sph = sqrt(SQR(x1) + SQR(x3));
            sin_theta = x1 / r_sph;
            cos_theta = x3 / r_sph;
          }
          mom_Kepl = rho * sqrt(GM1/pcoord->x1v(i));
          // fill out all parameters
          if (pin->GetInteger("problem", "stratified") == 0) {
            // unstratified Keplerian disk
            phydro->u(IDN, k,j,i) = rho;
            phydro->u(IM1, k,j,i) = 0.0;
            phydro->u(IM2, k,j,i) = mom_Kepl;
            phydro->u(IM3, k,j,i) = 0.0;
          } else if (pin->GetInteger("problem", "stratified") == 1) {
            // thin Keplerian disk inclined wrt the coordinate system
            sin_ang_dist = cos_incl*cos_theta + sin_incl*sin_theta*sin(position_angle - phi);
            if (abs(sin_ang_dist) < sin_ang_dist_max) { // inside the ring
              phydro->u(IDN, k,j,i) = rho;
              sin_xi = (cos_incl - cos_theta*sin_ang_dist) / (sin_theta*sqrt(1.-SQR(sin_ang_dist)));
              cos_xi = sqrt(1.-SQR(sin_xi));
              if (cos(position_angle - phi) < 0.)
                cos_xi = -cos_xi;
              if (COORDINATE_SYSTEM == "spherical_polar") {
                phydro->u(IM1, k,j,i) = 0.0;
                phydro->u(IM2, k,j,i) = -mom_Kepl * cos_xi;
                phydro->u(IM3, k,j,i) = mom_Kepl * sin_xi;
              } else if (COORDINATE_SYSTEM == "cylindrical") {
                phydro->u(IM1, k,j,i) = mom_Kepl * cos_xi * cos_theta;
                phydro->u(IM2, k,j,i) = mom_Kepl * sin_xi;
                phydro->u(IM3, k,j,i) = mom_Kepl * cos_xi * sin_theta;
              }
            } else { // outside the ring: void
                phydro->u(IDN, k,j,i) = rho0;
                phydro->u(IM1, k,j,i) = 0.0;
                phydro->u(IM2, k,j,i) = 0.0;
                phydro->u(IM3, k,j,i) = 0.0;
            }
          }
          #if NON_BAROTROPIC_EOS == 1
          phydro->u(IEN, k,j,i) = press/(adiab_idx-1.) + 0.5 * (SQR(phydro->u(IM1, k,j,i)) + SQR(phydro->u(IM2, k,j,i)) + SQR(phydro->u(IM3, k,j,i))) / rho;
          #endif
        }
      }
    }
  } else if (pin->GetInteger("problem", "binary_component") == 1) { // DROPPING SPHERES ONTO THE SECONDARY (primary gravity OFF)
    for (k=ks; k<=ke; k++) {   //z / phi
      x3 = pcoord->x3v(k);
      if (COORDINATE_SYSTEM == "cylindrical") {
        z = x3;
      }
      else if (COORDINATE_SYSTEM == "spherical_polar") {
        cos_phi = cos(x3); sin_phi = sin(x3);
        cos_phi_obj = cos(x3 - obj_phi0);// sin_phi_obj = sin(x3 - obj_phi0);
      }
      for (j=js; j<=je; j++) {  //phi / theta
        x2 = pcoord->x2v(j);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          sin_theta = sin(x2); cos_theta = cos(x2);
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          cos_phi = cos(x2); sin_phi = sin(x2);
          cos_phi_obj = cos(x2 - obj_phi0);// sin_phi_obj = sin(x2 - obj_phi0);
        }
        #pragma omp simd
        for (i=is; i<=ie; i++) { //r_cyl / r
          x1 = pcoord->x1v(i);
          if (COORDINATE_SYSTEM == "spherical_polar") {
            r_cyl = x1 * sin_theta;
            z = x1 * cos_theta;
          }
          else if (COORDINATE_SYSTEM == "cylindrical") {
            r_cyl = x1;
          }
          // distance to the object
          sqr_dist_obj = SQR(r_cyl) + SQR(obj_r0) - 2.*r_cyl*obj_r0*cos_phi_obj + SQR(z - obj_z0);
          if (sqr_dist_obj < sqr_obj_size) {
            phydro->u(IDN, k,j,i) = rho;
          } else {
            phydro->u(IDN, k,j,i) = rho0;
          }
          phydro->u(IM1, k,j,i) = 0.0;
          phydro->u(IM2, k,j,i) = 0.0;
          phydro->u(IM3, k,j,i) = 0.0;
          #if NON_BAROTROPIC_EOS == 1
          phydro->u(IEN, k,j,i) = press/(adiab_idx-1.) + 0.5 * (SQR(phydro->u(IM1, k,j,i)) + SQR(phydro->u(IM2, k,j,i)) + SQR(phydro->u(IM3, k,j,i))) / rho;
          #endif
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void __attribute__((weak)) MeshBlock::UserWorkInLoop(void) {
  // do nothing
  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief Function called before generating output files
//========================================================================================

void __attribute__((weak)) MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // do nothing
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Function called after main loop is finished for user-defined work.
//========================================================================================

void __attribute__((weak)) Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // do nothing
  return;
}
