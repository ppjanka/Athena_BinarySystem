#ifndef GLOB_ACC_DISK_SRC_SECONDARY_GRAV
#define GLOB_ACC_DISK_SRC_SECONDARY_GRAV

//! \file globAccDisk-src-primary_grav.hpp
//  \brief Source terms for the gravitation of a binary system's primary at the coordinate center; Patryk Pjanka & James Stone
/*
 * Options for different coordinate systems enclosed in separate functions for optimization.
 * The code assumes omega_frame = 1., i.e., that angular velocity of the frame rotation equals orbital velocity of the binary (when GM_tot = 1.0)
 * External variables:
    -- Real GM = pmb->ruser_meshblock_data[0](0): mass of the primary (note: GM_tot = 1)
*/
#include "globAccDisk.hpp" // includes and declarations

//========================================================================================
//! \fn void void SecondarySourceGrav_unstratified_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The unstratified version of the source term for the gravity of the secondary body in a system. Note: cylindrical coordinates required.
//========================================================================================
void SecondarySourceGrav_unstratified (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = 1.0 - pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, sqr_rprim, rprim, z;
  Real sin_phi, cos_phi;
  Real dM, dMr, dMphi, dMz;
  Real sin_alpha, cos_alpha;

  for (k=pmb->ks; k<=pmb->ke; k++) {
    for (j=pmb->js; j<=pmb->je; j++) {
      sin_phi = sin(pmb->pcoord->x2v(j));
      cos_phi = cos(pmb->pcoord->x2v(j));
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) {

        r = pmb->pcoord->x1v(i);
        sqr_rprim = SQR(r) + 1.0 - 2.0 * r * cos_phi;
        rprim = sqrt(sqr_rprim);

        sin_alpha = sin_phi / rprim;
        cos_alpha = 0.5 * (SQR(r) + sqr_rprim - 1.0) / (r*rprim);

        dM = dt*cons(IDN,k,j,i)*GM / sqr_rprim; // total momentum change
        dMr = - cos_alpha * dM;
        dMphi = - sin_alpha * dM;

        cons(IM1,k,j,i) += dMr;
        cons(IM2,k,j,i) += dMphi;

        #if NON_BAROTROPIC_EOS == 1
        cons(IEN,k,j,i) += (0.5*dM*dM - cons(IM1,k,j,i)*dMr - cons(IM2,k,j,i)*dMphi) / cons(IDN,k,j,i);
        #endif

      }
    }
  }

  return;
}

//========================================================================================
//! \fn void void SecondarySourceGrav_stratified_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The stratified version of the source term for the gravity of the secondary body in a system, cylindrical coordinates' version.
//========================================================================================
void SecondarySourceGrav_stratified_cylindrical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = 1.0 - pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, rprim, /*rtot,*/ sqr_rprim, sqr_rtot, phi, z;
  Real dM, dM_bin, dMr, dMphi, dMz;
  Real sin_phi, cos_phi;
  Real sin_alpha, cos_alpha;

  for (k=pmb->ks; k<=pmb->ke; k++) {
    z = pmb->pcoord->x3v(k);
    for (j=pmb->js; j<=pmb->je; j++) {
      phi = pmb->pcoord->x2v(j);
      sin_phi = sin(phi); cos_phi = cos(phi);
      #pragma omp simd
      for (i=pmb->is; i<=pmb->ie; i++) {

        r = pmb->pcoord->x1v(i);
        sqr_rprim = SQR(r) + 1.0 - 2.0 * r * cos_phi;
        rprim = sqrt(sqr_rprim);
        sqr_rtot = sqr_rprim + SQR(z);
        //rtot = sqrt(sqr_rtot); // <- for some reason this crashes on Intel compilers (returns rtot = 0.), substituted explicitly in the code [Oct 3rd, 2018; Perseus]

        sin_alpha = sin_phi / rprim;
        cos_alpha = 0.5 * (SQR(r) + sqr_rprim - 1.0) / (r*rprim);

        dM = dt*cons(IDN,k,j,i)*GM / sqr_rtot; // total momentum change

        dMz = - dM * z / sqrt(sqr_rtot);
        dM_bin = dM * sqrt(sqr_rprim / sqr_rtot);   // momentum change projected onto the binary plane
        dMr = - cos_alpha * dM_bin;
        dMphi = - sin_alpha * dM_bin;

        cons(IM1,k,j,i) += dMr;
        cons(IM2,k,j,i) += dMphi;
        cons(IM3,k,j,i) += dMz;

        #if NON_BAROTROPIC_EOS == 1
        cons(IEN,k,j,i) += (0.5*SQR(dM) - cons(IM1,k,j,i)*dMr - cons(IM2,k,j,i)*dMphi - cons(IM3,k,j,i)*dMz) / cons(IDN,k,j,i);
        #endif
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void void SecondarySourceGrav_stratified_spherical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
// const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief The stratified version of the source term for the gravity of the secondary body in a system, spherical coordinates' version.
//========================================================================================
void SecondarySourceGrav_stratified_spherical (MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real GM = 1.0 - pmb->ruser_meshblock_data[0](0); //GM1 + GM2 = 1

  int i,j,k;
  Real r, sin_theta, cos_theta, sin_phi, cos_phi;
  Real dM, dM_tangent, dMr, dMtheta, dMphi;
  Real sin_xi, cos_xi, sin_eta, cos_eta, sin_nu, cos_nu;
  Real sqr_r, sqr_rtot, rtot;

  for (k=pmb->ks; k<=pmb->ke; k++) {     //phi
    sin_phi = sin(pmb->pcoord->x3v(k));
    cos_phi = cos(pmb->pcoord->x3v(k));
    for (j=pmb->js; j<=pmb->je; j++) {   //theta
      sin_theta = sin(pmb->pcoord->x2v(j));
      cos_theta = cos(pmb->pcoord->x2v(j));
      // angle xi: secondary star - coord. center - current point
      cos_xi = sin_theta*cos_phi;
      sin_xi = sin(acos(cos_xi)); // xi \in [0,pi], so this works
      if (sin_xi != 0.0) { // (r, theta, phi) off the binary axis
        #pragma omp simd
        for (i=pmb->is; i<=pmb->ie; i++) { //r

          r = pmb->pcoord->x1v(i);
          sqr_r = SQR(r);

          // sqr distance point - secondary star
          sqr_rtot = 1.0 + sqr_r - 2.0 * r * cos_xi;
          rtot = sqrt(sqr_rtot);

          // total momentum change
          dM = dt*cons(IDN,k,j,i)*GM / sqr_rtot;

          // angle eta: coord. center - curr. point - secondary star
          cos_eta = 0.5 * (sqr_r + sqr_rtot - 1.0) / (r * rtot);
          sin_eta = sin_xi / rtot;
          dMr = - dM * cos_eta;

          // momentum update component tangent to the sphere
          dM_tangent = dM * sin_eta;

          // angle nu: between the great circle arc from the current point towards the binary plane and from the current point to the point where the binary axis crosses the sphere (position angle of dM_tangent on the sphere)
          cos_nu = (cos_phi - cos_xi*sin_theta)
              / ( sin_xi * cos_theta );
          sin_nu = sin_phi / sin_xi;

          dMtheta = dM_tangent * cos_nu;
          dMphi = - dM_tangent * sin_nu;

          cons(IM1,k,j,i) += dMr;
          cons(IM2,k,j,i) += dMtheta;
          cons(IM3,k,j,i) += dMphi;

          #if NON_BAROTROPIC_EOS == 1
          cons(IEN,k,j,i) += (0.5*SQR(dM) - cons(IM1,k,j,i)*dMr - cons(IM2,k,j,i)*dMtheta - cons(IM3,k,j,i)*dMphi) / cons(IDN,k,j,i);
          #endif
        }
      } else { // (r, theta, phi) on the binary axis
        #pragma omp simd
        for (i=pmb->is; i<=pmb->ie; i++) { //r

          r = pmb->pcoord->x1v(i);
          sqr_r = SQR(r);

          // sqr distance point - secondary star
          sqr_rtot = 1.0 + sqr_r - 2.0 * r * cos_xi;

          // total momentum change
          dM = dt*cons(IDN,k,j,i)*GM / sqr_rtot;

          // angle eta: coord. center - curr. point - secondary star
          cos_eta = -cos_xi; // here: eta = pi - xi
          dMr = - dM * cos_eta;
          dMtheta = 0.;
          dMphi = 0.;

          cons(IM1,k,j,i) += dMr;
          cons(IM2,k,j,i) += dMtheta;
          cons(IM3,k,j,i) += dMphi;

          #if NON_BAROTROPIC_EOS == 1
          cons(IEN,k,j,i) += (0.5*SQR(dM) - cons(IM1,k,j,i)*dMr - cons(IM2,k,j,i)*dMtheta - cons(IM3,k,j,i)*dMphi) / cons(IDN,k,j,i);
          #endif
        }
      }
    }
  }

  return;
}

#endif //GLOB_ACC_DISK_SRC_SECONDARY_GRAV
