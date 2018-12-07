#ifndef GLOB_ACC_DISK_SRC_ROTATING_BARYCENTER
#define GLOB_ACC_DISK_SRC_ROTATING_BARYCENTER

//! \file globAccDisk-src-rotating_barycenter.hpp
//  \brief Source terms for the frame of reference rotating around a binary system's barycenter; Patryk Pjanka & James Stone
/*
 * Options for different coordinate systems enclosed in separate functions for optimization.
 * The code assumes omega_frame = 1., i.e., that angular velocity of the frame rotation equals orbital velocity of the binary (when GM_tot = 1.0)
 * External variables:
    -- Real GM = pmb->ruser_meshblock_data[0](0): mass of the primary (note: GM_tot = 1)
*/
#include "globAccDisk.hpp" // includes and declarations

//----------------------------------------------------------------------------------------
//! \fn void RotatingFrame_barycenter
//  \brief Calculates the source terms at single point due to rotation around the barycenter of the binary system, in cylindrical coordinates.

#if NON_BAROTROPIC_EOS == 1
inline void RotatingFrame_barycenter (Real dt, Real GM, Real r, Real phi, Real rho, Real Mr, Real Mphi, Real& dMr, Real& dMphi, Real& dE) {
#else
inline void RotatingFrame_barycenter (Real dt, Real GM, Real r, Real phi, Real rho, Real Mr, Real Mphi, Real& dMr, Real& dMphi) {
#endif

  Real r_bary, phi_bary, Mr_bary, Mphi_bary, Fr_bary, Fphi_bary; // variables in barycenter coord. sys.
  Real sin_alpha, cos_alpha;
  //Real a_prim = 1.-GM; // distance between LAB origin and the barycenter

  // move to the barycenter system of coordinates
  r_bary = sqrt(SQR(r) + SQR(1.-GM) - 2.*r*(1.-GM)*cos(phi));
  phi_bary = asin(r*sin(phi)/r_bary);
  if (r*cos(phi) < (1.-GM)) // [OPT]: any way we can avoid if??
    phi_bary = M_PI - phi_bary;
  // useful angle for coordinate transformations
  sin_alpha = sin(phi - phi_bary + 0.5*M_PI);
  cos_alpha = cos(phi - phi_bary + 0.5*M_PI);
  // momenta
  Mr_bary = Mr * sin_alpha + Mphi * cos_alpha;
  Mphi_bary = - Mr * cos_alpha + Mphi * sin_alpha;

  // calculate forces in the barycenter coord. sys.
  // hard-setting omega_frame = 1.0 (binary orbital speed)
  Fr_bary = rho * r_bary /*SQR(omega_frame)*/ //centripetal acceleration
      + 2.0 * Mphi_bary /*omega_frame*/; //Coriollis force
  Fphi_bary = -2.0 * Mr_bary /*omega_frame*/ ; //Coriollis force

  // move back to the LAB coordinates
  // multiply by dt to get momentum updates
  dMr = dt * (Fr_bary * sin_alpha - Fphi_bary * cos_alpha);
  dMphi = dt * (Fr_bary * cos_alpha + Fphi_bary * sin_alpha);

  #if NON_BAROTROPIC_EOS == 1
  // total energy update
  dE = - (dMr * (Mr + 0.5*dMr) + dMphi * (Mphi + 0.5*dMphi)) / rho;
  #endif

}

//----------------------------------------------------------------------------------------
//! \fn void RotatingFrame_barycenter_cylindrical
//  \brief Adds source terms due to rotation around the barycenter of the binary system, version for the cylindrical coordinate system.

void RotatingFrame_barycenter_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = pmb->ruser_meshblock_data[0](0);

  Real r, phi, rho, dMr, dMphi; // variables in LAB
  #if NON_BAROTROPIC_EOS == 1
  Real dE;
  #endif
  int i,j,k;
  for (k=pmb->ks; k<=pmb->ke; ++k) {    //z
    for (j=pmb->js; j<=pmb->je; ++j) {  //phi
      phi = pmb->pcoord->x2v(j);
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) { //r

          r = pmb->pcoord->x1v(i);
          rho = cons(IDN,k,j,i);

          //calculate the source terms
          #if NON_BAROTROPIC_EOS == 1
          RotatingFrame_barycenter (dt, GM, r, phi, rho, cons(IM1,k,j,i), cons(IM2,k,j,i), dMr, dMphi, dE);
          #else
          RotatingFrame_barycenter (dt, GM, r, phi, rho, cons(IM1,k,j,i), cons(IM2,k,j,i), dMr, dMphi);
          #endif

          //update momenta
          cons(IM1,k,j,i) += dMr;
          cons(IM2,k,j,i) += dMphi;

          #if NON_BAROTROPIC_EOS == 1
          //update total energy
          cons(IEN,k,j,i) += dE;
          #endif

      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RotatingFrame_barycenter_spherical
//  \brief Adds source terms due to rotation around the barycenter of the binary system, version for the spherical coordinate system.

void RotatingFrame_barycenter_spherical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = pmb->ruser_meshblock_data[0](0);

  // variables in LAB
  Real rho;
  Real theta, cos_theta, sin_theta; // polar angle in spherical coords
  Real r, phi; // cylindrical coordinates of a given point
  Real Mr, Mphi; // original momenta in CYLINDRICAL coords
  Real dMr, dMphi; // momentum updates in CYLINDRICAL coords
  #if NON_BAROTROPIC_EOS == 1
  Real dE;
  #endif

  int i,j,k;
  for (k=pmb->ks; k<=pmb->ke; ++k) {     // phi
    phi = pmb->pcoord->x3v(k);
    for (j=pmb->js; j<=pmb->je; ++j) {   // theta
      theta = pmb->pcoord->x2v(j);
      cos_theta = cos(theta);
      sin_theta = sin(theta);
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) { //r

          rho = cons(IDN,k,j,i);

          // transform to cylindrical coordinates
          r = pmb->pcoord->x1v(i) * sin_theta;
          Mr = cons(IM1,k,j,i) * sin_theta + cons(IM2,k,j,i) * cos_theta;
          Mphi = cons(IM3,k,j,i);

          //calculate the source terms
          #if NON_BAROTROPIC_EOS == 1
          RotatingFrame_barycenter (dt, GM, r, phi, rho, Mr, Mphi, dMr, dMphi, dE);
          #else
          RotatingFrame_barycenter (dt, GM, r, phi, rho, Mr, Mphi, dMr, dMphi);
          #endif

          //update momenta
          cons(IM1,k,j,i) += dMr * sin_theta;
          cons(IM2,k,j,i) += dMr * cos_theta;
          cons(IM3,k,j,i) += dMphi;

          #if NON_BAROTROPIC_EOS == 1
          //update total energy
          cons(IEN,k,j,i) += dE;
          #endif

      }
    }
  }

  return;
}

#endif //GLOB_ACC_DISK_SRC_ROTATING_BARYCENTER
