//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globAccDisk-test-rotating_barycenter
//  \brief Unit test for the frame of reference in rotation around a binary system's barycenter: a static object as seen from a rotating frame of reference; Patryk Pjanka

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"

// rotating frame of reference
#include "globAccDisk-src-rotating_barycenter.hpp"

// declarations
void Boundary_ix1 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Boundary_ox1 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Boundary_ix2 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Boundary_ox2 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void __attribute__((weak)) Mesh::InitUserMeshData(ParameterInput *pin) {
  
  if (COORDINATE_SYSTEM == "cylindrical")
    EnrollUserExplicitSourceFunction(RotatingFrame_barycenter_cylindrical);
  else if (COORDINATE_SYSTEM == "spherical_polar")
    EnrollUserExplicitSourceFunction(RotatingFrame_barycenter_spherical);

  if (pin->GetString("mesh", "ix1_bc") == "user")
    EnrollUserBoundaryFunction(INNER_X1, Boundary_ix1);
  if (pin->GetString("mesh", "ox1_bc") == "user")
    EnrollUserBoundaryFunction(OUTER_X1, Boundary_ox1);
  if (pin->GetString("mesh", "ix2_bc") == "user")
    EnrollUserBoundaryFunction(INNER_X2, Boundary_ix2);
  if (pin->GetString("mesh", "ox2_bc") == "user")
    EnrollUserBoundaryFunction(OUTER_X2, Boundary_ox2);

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
  #if NON_BAROTROPIC_EOS == 1
  ruser_meshblock_data[0].NewAthenaArray(6);
  ruser_meshblock_data[0](5) = pin->GetReal("problem", "press");
  #else
  ruser_meshblock_data[0].NewAthenaArray(5);
  #endif
  ruser_meshblock_data[0](0) = pin->GetReal("problem", "GM1"); // mass of the central object (note: GM1+GM2 = 1 for the binary, GM < 1 required)
  ruser_meshblock_data[0](1) = pin->GetReal("problem", "obj_radius");
  ruser_meshblock_data[0](2) = pin->GetReal("problem", "obj_size");
  ruser_meshblock_data[0](3) = pin->GetReal("problem", "rho0");
  ruser_meshblock_data[0](4) = pin->GetReal("problem", "rho1");

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//========================================================================================

void __attribute__((weak)) MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rho0 = pin->GetReal("problem", "rho0"); // ambient medium density
  Real rho1 = pin->GetReal("problem", "rho1"); // density of the object
  #if NON_BAROTROPIC_EOS == 1
  Real press = pin->GetReal("problem", "press"); // pressure in the box
  Real adiab_idx = pin->GetReal("hydro", "gamma"); // adiabatic idx
  #endif
  Real obj_radius = pin->GetReal("problem", "obj_radius"); // radius of the center of the object
  Real obj_phi0 = pin->GetOrAddReal("problem", "obj_phi0", 0.0); // initial azimuth of the center of the object
  Real obj_size = pin->GetReal("problem", "obj_size");
  Real sqr_obj_size = SQR(obj_size);
  Real rho;

  Real GM = pin->GetReal("problem", "GM1");
  Real r_bary = 1. - GM;
  
  int i,j,k;
  Real x1, x2, x3, r_cyl, z;
  Real sin_theta, cos_theta;
  Real sin_phi, cos_phi;
  Real cos_phi_obj;
  Real sin_beta, cos_beta; // pi/2 - (angle coord.ctr.-curr.pos.-barycen.)
  Real sqr_dist_obj;
  Real sqr_dist_bary, dist_bary;
  Real mom_tot;
  Real x1max = pcoord->x1f(ie+1), x2max = pcoord->x2f(je+1), x3max = pcoord->x3f(ke+1);
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
        sqr_dist_obj = SQR(r_cyl) + SQR(obj_radius) - 2.*r_cyl*obj_radius*cos_phi_obj + SQR(z);
        if (sqr_dist_obj < sqr_obj_size) {
          rho = rho1;
        } else {
          rho = rho0;
        }
        // fill out other parameters
        phydro->u(IDN, k,j,i) = rho;
        // distance to the barycenter
        sqr_dist_bary = SQR(r_cyl) + SQR(r_bary) - 2.*r_cyl*r_bary*cos_phi;
        dist_bary = sqrt(sqr_dist_bary);
        sin_beta = 0.5 * (SQR(r_cyl)+sqr_dist_bary-SQR(r_bary)) / (r_cyl * dist_bary);
        cos_beta = r_bary * sin_phi / dist_bary;
        mom_tot = rho * /*omega_frame*/ dist_bary;
        if (COORDINATE_SYSTEM == "cylindrical") {
          phydro->u(IM1, k,j,i) = mom_tot * cos_beta;
          phydro->u(IM2, k,j,i) = -mom_tot * sin_beta;
          phydro->u(IM3, k,j,i) = 0.0;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          phydro->u(IM1, k,j,i) = mom_tot * cos_beta * sin_theta;
          phydro->u(IM2, k,j,i) = mom_tot * cos_beta * cos_theta;
          phydro->u(IM3, k,j,i) = -mom_tot * sin_beta;
        }
        #if NON_BAROTROPIC_EOS == 1
        phydro->u(IEN, k,j,i) = press/(adiab_idx-1.) + 0.5 * (SQR(phydro->u(IM1, k,j,i)) + SQR(phydro->u(IM2, k,j,i)) + SQR(phydro->u(IM3, k,j,i))) / rho;
        #endif
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

//========================================================================================
//! \fn void Boundary_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Inner radial boundary conditions for a reference frame rotating around a barycenter.
//========================================================================================

void Boundary_ix1 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real GM = pmb->ruser_meshblock_data[0](0);
  Real obj_radius = pmb->ruser_meshblock_data[0](1);
  Real obj_size = pmb->ruser_meshblock_data[0](2);
  Real sqr_obj_size = SQR(obj_size);
  Real rho0 = pmb->ruser_meshblock_data[0](3);
  Real rho1 = pmb->ruser_meshblock_data[0](4);
  #if NON_BAROTROPIC_EOS == 1
  Real press = pmb->ruser_meshblock_data[0](5);
  #endif

  Real r_bary = 1.-GM;
  Real rho;

  int i,j,k;
  Real x1, x2, x3, r_cyl;
  Real sin_theta, cos_theta;
  Real sin_phi, cos_phi;
  Real sin_beta, cos_beta; // pi/2 - (angle coord.ctr.-curr.pos.-barycen.)
  Real sqr_dist_obj;
  Real sqr_dist_bary, dist_bary;
  Real vel_tot;
  for (k=ks; k<=ke; ++k) {   //z / phi
    x3 = pco->x3v(k);
    if (COORDINATE_SYSTEM == "spherical_polar") {
      cos_phi = cos(x3); sin_phi = sin(x3);
    }
    for (j=js; j<=je; ++j) {  //phi / theta
      x2 = pco->x2v(j);
      if (COORDINATE_SYSTEM == "spherical_polar") {
        sin_theta = sin(x2); cos_theta = cos(x2);
      }
      else if (COORDINATE_SYSTEM == "cylindrical") {
        cos_phi = cos(x2); sin_phi = sin(x2);
      }
      #pragma omp simd
      for (i=1; i<=NGHOST; ++i) { //r_cyl / r
        x1 = pco->x1v(is-i);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          r_cyl = x1 * sin_theta;
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          r_cyl = x1;
        }
        // fill out other parameters
        rho = rho0;
        prim(IDN, k,j,is-i) = rho;
        // distance to the barycenter
        sqr_dist_bary = SQR(r_cyl) + SQR(r_bary) - 2.*r_cyl*r_bary*cos_phi;
        dist_bary = sqrt(sqr_dist_bary);
        sin_beta = 0.5 * (SQR(r_cyl)+sqr_dist_bary-SQR(r_bary)) / (r_cyl * dist_bary);
        cos_beta = r_bary * sin_phi / dist_bary;
        vel_tot = /*omega_frame*/ dist_bary;
        if (COORDINATE_SYSTEM == "cylindrical") {
          prim(IVX, k,j,is-i) = vel_tot * cos_beta;
          prim(IVY, k,j,is-i) = -vel_tot * sin_beta;
          prim(IVZ, k,j,is-i) = 0.0;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          prim(IVX, k,j,is-i) = vel_tot * cos_beta * sin_theta;
          prim(IVY, k,j,is-i) = vel_tot * cos_beta * cos_theta;
          prim(IVZ, k,j,is-i) = -vel_tot * sin_beta;
        }
        #if NON_BAROTROPIC_EOS == 1
        prim(IPR, k,j,is-i) = press;
        #endif
      }
    }
  }
}

// [OPT]: pre-calculate the boundary conditions and hold them in memory
//========================================================================================
//! \fn void Boundary_ox1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Outer radial boundary conditions for a reference frame rotating around a barycenter.
//========================================================================================

void Boundary_ox1 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real GM = pmb->ruser_meshblock_data[0](0);
  Real obj_radius = pmb->ruser_meshblock_data[0](1);
  Real obj_size = pmb->ruser_meshblock_data[0](2);
  Real sqr_obj_size = SQR(obj_size);
  Real rho0 = pmb->ruser_meshblock_data[0](3);
  Real rho1 = pmb->ruser_meshblock_data[0](4);
  #if NON_BAROTROPIC_EOS == 1
  Real press = pmb->ruser_meshblock_data[0](5);
  #endif

  Real r_bary = 1.-GM;
  Real rho;

  int i,j,k;
  Real x1, x2, x3, r_cyl;
  Real sin_theta, cos_theta;
  Real sin_phi, cos_phi;
  Real sin_beta, cos_beta; // pi/2 - (angle coord.ctr.-curr.pos.-barycen.)
  Real sqr_dist_obj;
  Real sqr_dist_bary, dist_bary;
  Real vel_tot;
  for (k=ks; k<=ke; ++k) {   //z / phi
    x3 = pco->x3v(k);
    if (COORDINATE_SYSTEM == "spherical_polar") {
      cos_phi = cos(x3); sin_phi = sin(x3);
    }
    for (j=js; j<=je; ++j) {  //phi / theta
      x2 = pco->x2v(j);
      if (COORDINATE_SYSTEM == "spherical_polar") {
        sin_theta = sin(x2); cos_theta = cos(x2);
      }
      else if (COORDINATE_SYSTEM == "cylindrical") {
        cos_phi = cos(x2); sin_phi = sin(x2);
      }
      #pragma omp simd
      for (i=1; i<=NGHOST; ++i) { //r_cyl / r
        x1 = pco->x1v(ie+i);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          r_cyl = x1 * sin_theta;
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          r_cyl = x1;
        }
        // fill out other parameters
        rho = rho0;
        prim(IDN, k,j,ie+i) = rho;
        sqr_dist_bary = SQR(r_cyl) + SQR(r_bary) - 2.*r_cyl*r_bary*cos_phi;
        dist_bary = sqrt(sqr_dist_bary);
        sin_beta = 0.5 * (SQR(r_cyl)+sqr_dist_bary-SQR(r_bary)) / (r_cyl * dist_bary);
        cos_beta = r_bary * sin_phi / dist_bary;
        vel_tot = /*omega_frame*/ dist_bary;
        if (COORDINATE_SYSTEM == "cylindrical") {
          prim(IVX, k,j,ie+i) = vel_tot * cos_beta;
          prim(IVY, k,j,ie+i) = -vel_tot * sin_beta;
          prim(IVZ, k,j,ie+i) = 0.0;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          prim(IVX, k,j,ie+i) = vel_tot * cos_beta * sin_theta;
          prim(IVY, k,j,ie+i) = vel_tot * cos_beta * cos_theta;
          prim(IVZ, k,j,ie+i) = -vel_tot * sin_beta;
        }
        #if NON_BAROTROPIC_EOS == 1
        prim(IPR, k,j,ie+i) = press;
        #endif
      }
    }
  }
}

//========================================================================================
//! \fn void Boundary_ix2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Inner polar boundary conditions for a reference frame rotating around a barycenter.
//========================================================================================

void Boundary_ix2 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real GM = pmb->ruser_meshblock_data[0](0);
  Real obj_radius = pmb->ruser_meshblock_data[0](1);
  Real obj_size = pmb->ruser_meshblock_data[0](2);
  Real sqr_obj_size = SQR(obj_size);
  Real rho0 = pmb->ruser_meshblock_data[0](3);
  Real rho1 = pmb->ruser_meshblock_data[0](4);
  #if NON_BAROTROPIC_EOS == 1
  Real press = pmb->ruser_meshblock_data[0](5);
  #endif

  Real r_bary = 1.-GM;
  Real rho;

  int i,j,k;
  Real x1, x2, x3, r_cyl;
  Real sin_theta, cos_theta;
  Real sin_phi, cos_phi;
  Real sin_beta, cos_beta; // pi/2 - (angle coord.ctr.-curr.pos.-barycen.)
  Real sqr_dist_obj;
  Real sqr_dist_bary, dist_bary;
  Real vel_tot;
  for (k=ks; k<=ke; ++k) {   //z / phi
    x3 = pco->x3v(k);
    if (COORDINATE_SYSTEM == "spherical_polar") {
      cos_phi = cos(x3); sin_phi = sin(x3);
    }
    for (j=1; j<=NGHOST; ++j) {  //phi / theta
      x2 = pco->x2v(js-j);
      if (COORDINATE_SYSTEM == "spherical_polar") {
        sin_theta = sin(x2); cos_theta = cos(x2);
      }
      else if (COORDINATE_SYSTEM == "cylindrical") {
        cos_phi = cos(x2); sin_phi = sin(x2);
      }
      #pragma omp simd
      for (i=is; i<=ie; ++i) { //r_cyl / r
        x1 = pco->x1v(i);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          r_cyl = x1 * sin_theta;
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          r_cyl = x1;
        }
        // fill out other parameters
        rho = rho0;
        prim(IDN, k,js-j,i) = rho;
        // distance to the barycenter
        sqr_dist_bary = SQR(r_cyl) + SQR(r_bary) - 2.*r_cyl*r_bary*cos_phi;
        dist_bary = sqrt(sqr_dist_bary);
        sin_beta = 0.5 * (SQR(r_cyl)+sqr_dist_bary-SQR(r_bary)) / (r_cyl * dist_bary);
        cos_beta = r_bary * sin_phi / dist_bary;
        vel_tot = /*omega_frame*/ dist_bary;
        if (COORDINATE_SYSTEM == "cylindrical") {
          prim(IVX, k,js-j,i) = vel_tot * cos_beta;
          prim(IVY, k,js-j,i) = -vel_tot * sin_beta;
          prim(IVZ, k,js-j,i) = 0.0;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          prim(IVX, k,js-j,i) = vel_tot * cos_beta * sin_theta;
          prim(IVY, k,js-j,i) = vel_tot * cos_beta * cos_theta;
          prim(IVZ, k,js-j,i) = -vel_tot * sin_beta;
        }
        #if NON_BAROTROPIC_EOS == 1
        prim(IPR, k,js-j,i) = press;
        #endif
      }
    }
  }
}

//========================================================================================
//! \fn void Boundary_ox2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Outer polar boundary conditions for a reference frame rotating around a barycenter.
//========================================================================================

void Boundary_ox2 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real GM = pmb->ruser_meshblock_data[0](0);
  Real obj_radius = pmb->ruser_meshblock_data[0](1);
  Real obj_size = pmb->ruser_meshblock_data[0](2);
  Real sqr_obj_size = SQR(obj_size);
  Real rho0 = pmb->ruser_meshblock_data[0](3);
  Real rho1 = pmb->ruser_meshblock_data[0](4);
  #if NON_BAROTROPIC_EOS == 1
  Real press = pmb->ruser_meshblock_data[0](5);
  #endif

  Real r_bary = 1.-GM;
  Real rho;

  int i,j,k;
  Real x1, x2, x3, r_cyl;
  Real sin_theta, cos_theta;
  Real sin_phi, cos_phi;
  Real sin_beta, cos_beta; // pi/2 - (angle coord.ctr.-curr.pos.-barycen.)
  Real sqr_dist_obj;
  Real sqr_dist_bary, dist_bary;
  Real vel_tot;
  for (k=ks; k<=ke; ++k) {   //z / phi
    x3 = pco->x3v(k);
    if (COORDINATE_SYSTEM == "spherical_polar") {
      cos_phi = cos(x3); sin_phi = sin(x3);
    }
    for (j=1; j<=NGHOST; ++j) {  //phi / theta
      x2 = pco->x2v(je+j);
      if (COORDINATE_SYSTEM == "spherical_polar") {
        sin_theta = sin(x2); cos_theta = cos(x2);
      }
      else if (COORDINATE_SYSTEM == "cylindrical") {
        cos_phi = cos(x2); sin_phi = sin(x2);
      }
      #pragma omp simd
      for (i=is; i<=ie; ++i) { //r_cyl / r
        x1 = pco->x1v(i);
        if (COORDINATE_SYSTEM == "spherical_polar") {
          r_cyl = x1 * sin_theta;
        }
        else if (COORDINATE_SYSTEM == "cylindrical") {
          r_cyl = x1;
        }
        // fill out other parameters
        rho = rho0;
        prim(IDN, k,je+j,i) = rho;
        // distance to the barycenter
        sqr_dist_bary = SQR(r_cyl) + SQR(r_bary) - 2.*r_cyl*r_bary*cos_phi;
        dist_bary = sqrt(sqr_dist_bary);
        sin_beta = 0.5 * (SQR(r_cyl)+sqr_dist_bary-SQR(r_bary)) / (r_cyl * dist_bary);
        cos_beta = r_bary * sin_phi / dist_bary;
        vel_tot = /*omega_frame*/ dist_bary;
        if (COORDINATE_SYSTEM == "cylindrical") {
          prim(IVX, k,je+j,i) = vel_tot * cos_beta;
          prim(IVY, k,je+j,i) = -vel_tot * sin_beta;
          prim(IVZ, k,je+j,i) = 0.0;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          prim(IVX, k,je+j,i) = vel_tot * cos_beta * sin_theta;
          prim(IVY, k,je+j,i) = vel_tot * cos_beta * cos_theta;
          prim(IVZ, k,je+j,i) = -vel_tot * sin_beta;
        }
        #if NON_BAROTROPIC_EOS == 1
        prim(IPR, k,je+j,i) = press;
        #endif
      }
    }
  }
}
