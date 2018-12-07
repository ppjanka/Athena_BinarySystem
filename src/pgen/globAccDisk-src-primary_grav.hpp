#ifndef GLOB_ACC_DISK_SRC_PRIMARY_GRAV
#define GLOB_ACC_DISK_SRC_PRIMARY_GRAV

//! \file globAccDisk-src-primary_grav.hpp
//  \brief Source terms for the gravitation of a binary system's primary at the coordinate center; Patryk Pjanka & James Stone
/*
 * Options for different coordinate systems enclosed in separate functions for optimization.
 * External variables:
    -- Real GM = pmb->ruser_meshblock_data[0](0): mass of the primary (note: GM_tot = 1)
*/
#include "globAccDisk.hpp" // includes and declarations

//========================================================================================
//! \fn void void PrimarySourceGrav_unstratified (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The unstratified version of the source term for the gravity of the primary body in a system. Note: cylindrical coordinates required.
//========================================================================================
void PrimarySourceGrav_unstratified (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, dM;
  for (k=pmb->ks; k<=pmb->ke; k++) {
    for (j=pmb->js; j<=pmb->je; j++) {
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) {

        r = pmb->pcoord->x1v(i);
        dM = - dt*cons(IDN,k,j,i)*GM / SQR(r); // total momentum change
        cons(IM1,k,j,i) += dM;

        #if NON_BAROTROPIC_EOS == 1
        cons(IEN,k,j,i) += (0.5*dM*dM - cons(IM1,k,j,i)*dM) / cons(IDN,k,j,i);
        #endif
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void void PrimarySourceGrav_stratified_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The stratified version of the source term for the gravity of the primary body in a system, cylindrical coordinates' version.
//========================================================================================
void PrimarySourceGrav_stratified_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, z;
  Real rtot, sqr_rtot; // total distance to the coordinate center
  Real dM, dMr, dMz;
  for (k=pmb->ks; k<=pmb->ke; k++) {
    z = pmb->pcoord->x3v(k);
    for (j=pmb->js; j<=pmb->je; j++) {
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) {

        r = pmb->pcoord->x1v(i);
        sqr_rtot = SQR(r) + SQR(z);
        rtot = sqrt(sqr_rtot);

        dM = - dt*cons(IDN,k,j,i)*GM / sqr_rtot; // total momentum change
        dMr = dM * r / rtot;
        dMz = dM * z / rtot;
        cons(IM1,k,j,i) += dMr;
        cons(IM3,k,j,i) += dMz;

        #if NON_BAROTROPIC_EOS == 1
        cons(IEN,k,j,i) += (0.5*dM*dM - cons(IM1,k,j,i)*dMr - cons(IM3,k,j,i)*dMz) / cons(IDN,k,j,i);
        #endif
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void void PrimarySourceGrav_stratified_spherical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The stratified version of the source term for the gravity of the primary body in a system, spherical coordinates' version.
//========================================================================================
void PrimarySourceGrav_stratified_spherical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, dMr;
  for (k=pmb->ks; k<=pmb->ke; k++) {
    for (j=pmb->js; j<=pmb->je; j++) {
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) {

        r = pmb->pcoord->x1v(i);

        dMr = - dt*cons(IDN,k,j,i)*GM / SQR(r); // total momentum change
        cons(IM1,k,j,i) += dMr;

        #if NON_BAROTROPIC_EOS == 1
        cons(IEN,k,j,i) += (0.5*dMr*dMr - cons(IM1,k,j,i)*dMr) / cons(IDN,k,j,i);
        #endif
      }
    }
  }

  return;
}

#endif //GLOB_ACC_DISK_SRC_PRIMARY_GRAV
