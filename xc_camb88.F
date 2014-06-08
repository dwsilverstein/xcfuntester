      Subroutine xc_camb88(cam_omega, n, tol_rho, fac, lfac, nlfac, 
     1                     rho, gamma, Amat, Amat2, Amat3, Cmat, 
     2                     Cmat2, Cmat11, Cmat3, Cmat21, Cmat12, nq, 
     3                     ipol, Ex, qwght,ldew,func, as)
c
C$Id: xc_camb88.F 23977 2013-04-05 23:11:55Z dwsilverstein $
c
c     Coulomb attenuated Becke88 functional
c
      implicit none
c      
c#include "dft2drv.fh"
c Daniel (4-4-13): Third derivative header file
c#include "dft3drv.fh"
c
      double precision as
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,3)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,2), Cmat(nq,3)
c
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,3), Cmat2(nq,6), Cmat11(nq,6)
      double precision A2tmp, C2tmp, C3tmp
c Daniel (4-4-13): Third derivatives
c
c     Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,4), Cmat3(nq,10), Cmat21(nq,9),
     1                 Cmat12(nq,12)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
c
      double precision BETA
      Parameter (BETA = 0.0042D0)
c
      double precision cam_omega
      double precision cam_alpha, cam_beta
c
c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************
c
      integer n
      double precision arcsinh, darcsinh, d2arcsinh
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg,
     &     dgdenom, t, Etmp, Atmp, Ctmp
c Daniel (4-4-13): Suggested variable to add
      double precision gdenom2
      double precision rhom23, d2g, d2gdenom
c Daniel (4-4-13): Suggested variable to add
      double precision gdenom3
c
      double precision rhom53, d3g, d3gdenom
      double precision gdenom4
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c Daniel (4-4-13): Second derivative of arcsinh
      d2arcsinh(x) = -x/dsqrt(1d0+x*x)**3
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
c         do 10 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
c            if (dsqrt(gamma).gt.tol_rho)then
c               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
c            else
c               x = 0d0
c            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
c Daniel (4-4-13): New variable for simplifying expressions
            gdenom2 = gdenom*gdenom
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = 2d0*rho43*C*fac
               Atmp = (4d0/3d0)*rho13*C*fac
            endif
c
c            if (nlfac) then
               Etmp = Etmp + 2d0*rho43*g*fac
               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
c            endif
c
c            if (x.gt.tol_rho) then
               Ctmp = 0.5d0 * dg / sqrt(gamma) * fac
c            endif
c
            A2tmp = 0d0
            C2tmp = 0d0
            C3tmp = 0d0
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
c Daniel (4-4-13): New variable for simplifying expressions
            gdenom3 = gdenom2*gdenom
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
     &           + BETA*x*x*d2gdenom/gdenom2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
c     &           + BETA*x*x*d2gdenom/gdenom**2
c     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
c rr
            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c rg
            C2tmp = - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
c gg
c            if (x.gt.tol_rho) then
               C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
c            endif
c Daniel (4-4-13): Third derivative stuff
            A3tmp = 0.0d0
            C4tmp = 0.0d0
            C5tmp = 0.0d0
            C6tmp = 0.0d0
c
            rhom53 = rhom23 / (0.5d0*rho(n,1))
c
            d3gdenom = 6.0d0*BETA*
     1          d2arcsinh(x)*( 3.0d0
     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
            gdenom4 = gdenom3*gdenom
c
            d3g = 6.0d0*BETA*dgdenom/gdenom2
     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
     5          + BETA*x*x*d3gdenom/gdenom2
c
c rrr
            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
     1                                    - 18.0d0*x*x*d2g
     2                                    - 8.0d0*x*x*x*d3g)*fac
c
c rrg
            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
     1                                           + 4.0d0*x*x*x*d3g)*fac
c
c rgg
            C5tmp = -(8.0d0/3.0d0)*(rhom23/rho(n,1)**3)/dsqrt(gamma)
     1               *d3g*fac
c
c ggg
c            if (x.gt.tol_rho) then
              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
     1                                              - 3.0d0*x*d2g
     2                                              + x*x*d3g)*fac
c            endif
c
c Daniel (4-4-13): Third derivatives
            cam_alpha = 0.19d0
            cam_beta = 0.46d0
            call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp,cam_omega,
     &           cam_alpha,cam_beta,as)
c
            Amat2(n,1) = Amat2(n,1) + A2tmp
            Cmat11(n,1) = Cmat2(n,1) + C2tmp
            Cmat2(n,1) = Cmat2(n,1) + C3tmp
c
            Amat3(n,1) = Amat3(n,1) + A3tmp
            Cmat21(n,1) = Cmat21(n,1) + C4tmp
            Cmat12(n,1) = Cmat12(n,1) + C5tmp
            Cmat3(n,1) = Cmat3(n,1) + C6tmp
c
            Ex = Ex + qwght(n)*Etmp
            if(ldew) func(n) = func(n) + Etmp
            Amat(n,1) = Amat(n,1) + Atmp
            Cmat(n,1) = Cmat(n,1) + Ctmp
c 10      continue
c
      else
cc
cc        ======> SPIN-UNRESTRICTED <======
cc
c         do 20 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 20
c            if (rho(n,2).lt.tol_rho) goto 25
cc
cc           Spin alpha:
cc
c            rho13 = rho(n,2)**(1.d0/3.d0)
c            rho43 = rho13*rho(n,2)
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
c            if (dsqrt(gamma).gt.tol_rho)then
c               x = dsqrt(gamma) / rho43
c            else
c               x = 0d0
c            endif
cc
c            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
c            g = -BETA*x*x / gdenom
c            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
cc Daniel (4-4-13): New variable for simplifying expressions
c            gdenom2 = gdenom*gdenom
c            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
cc            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
cc
c            Etmp = 0.d0
c            Atmp = 0.d0
c            Ctmp = 0.d0
c            if (lfac) then
c               Etmp = rho43*C*fac
c               Atmp = (4d0/3d0)*rho13*C*fac
c            endif
cc
c            if (nlfac) then
c               Etmp = Etmp + rho43*g*fac
c               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
c            endif
cc
c            if (x.gt.tol_rho) then
c               Ctmp = 0.5d0*dg / sqrt(gamma) * fac
c            endif
cc
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            if (lfac) g = g + C           ! Add local contribution back to g
c            rhom23 = rho13 / rho(n,2)
c            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
cc Daniel (4-4-13): New variable for simplifying expressions
c            gdenom3 = gdenom2*gdenom
c            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
c     &           + BETA*x*x*d2gdenom/gdenom2
c     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
cc            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
cc     &           + BETA*x*x*d2gdenom/gdenom**2
cc     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
cc
c            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c            C2tmp =  (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
c            if (x.gt.tol_rho) then
c               C3tmp = -0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
c            endif
c#endif
cc Daniel (11-14-12): Third derivative stuff
c#ifdef THIRD_DERIV
c            rhom53 = rhom23 / rho(n,2)
cc
c            d3gdenom = 6.0d0*BETA*
c     1          d2arcsinh(x)*( 3.0d0
c     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
cc
c            gdenom4 = gdenom3*gdenom
cc
c            d3g = 6.0d0*BETA*dgdenom/gdenom2
c     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
c     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
c     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
c     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
c     5          + BETA*x*x*d3gdenom/gdenom2
cc
c            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
c     1                                    - 18.0d0*x*x*d2g
c     2                                    - 8.0d0*x*x*x*d3g)*fac
cc
c            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
c     1                                           + 4.0d0*x*x*x*d3g)*fac
cc
c            C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,2)**3)/dsqrt(gamma)
c     1               *d3g*fac
cc
c            if (x.gt.tol_rho) then
c              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
c     1                                              - 3.0d0*x*d2g
c     2                                              + x*x*d3g)*fac
c            endif
c#endif
cc
cc Daniel (4-4-13): Third derivatives
c#ifdef THIRD_DERIV
c            call xc_att_xc_d3(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
cc
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
cc
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
c            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp
c            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp
c            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + C6tmp
c#elif defined(SECOND_DERIV)
c            call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,C2tmp,
c     &           C3tmp)
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
c#else
c           call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
c#endif
c            Ex = Ex + qwght(n)*Etmp
c            if(ldew) func(n) = func(n) + Etmp
c            Amat(n,1) = Amat(n,1) + Atmp
c            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
cc
c 25         continue
cc
cc           Spin beta:
cc
c            if (rho(n,3).lt.tol_rho) goto 20
cc
c            rho13 = rho(n,3)**(1.d0/3.d0)
c            rho43 = rho13*rho(n,3)
c            gamma = delrho(n,1,2)*delrho(n,1,2) +
c     &              delrho(n,2,2)*delrho(n,2,2) +
c     &              delrho(n,3,2)*delrho(n,3,2)
c            if (dsqrt(gamma).gt.tol_rho)then
c               x = dsqrt(gamma) / rho43
c            else
c               x = 0d0
c            endif
cc
c            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
c            g = -BETA*x*x / gdenom
c            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
cc Daniel (4-4-13): New variable for simplifying expressions
c            gdenom2 = gdenom*gdenom
c            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
cc            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
cc
c            Etmp = 0.d0
c            Atmp = 0.d0
c            Ctmp = 0.d0
c            if (lfac) then
c               Etmp = rho43*C*fac
c               Atmp = (4d0/3d0)*rho13*C*fac
c            endif
cc
c            if (nlfac) then
c               Etmp = Etmp + rho43*g*fac
c               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
c            endif
cc
c            if (x.gt.tol_rho) then
c               Ctmp = 0.5d0*dg / sqrt(gamma) * fac
c            endif
cc
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = 0d0
c            C2tmp = 0d0
c            C3tmp = 0d0
c            if(lfac) g = g + C           ! Add local contribution back to g
c            rhom23 = rho13 / rho(n,3)
c            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
cc Daniel (4-4-13): New variable for simplifying expressions
c            gdenom3 = gdenom2*gdenom
c            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
c     &           + BETA*x*x*d2gdenom/gdenom2
c     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
cc            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
cc     &           + BETA*x*x*d2gdenom/gdenom**2
cc     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
cc
c            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c            C2tmp = -(2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
c            if (x.gt.tol_rho) then
c               C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
c            endif
c#endif
cc Daniel (4-4-13): Third derivative stuff
c#ifdef THIRD_DERIV
c            rhom53 = rhom23 / rho(n,3)
cc
c            d3gdenom = 6.0d0*BETA*
c     1          d2arcsinh(x)*( 3.0d0
c     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
cc
c            gdenom4 = gdenom3*gdenom
cc
c            d3g = 6.0d0*BETA*dgdenom/gdenom2
c     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
c     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
c     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
c     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
c     5          + BETA*x*x*d3gdenom/gdenom2
cc
c            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
c     1                                    - 18.0d0*x*x*d2g
c     2                                    - 8.0d0*x*x*x*d3g)*fac
cc
c            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
c     1                                           + 4.0d0*x*x*x*d3g)*fac
cc
c            C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,3)**3)/dsqrt(gamma)
c     1               *d3g*fac
cc
c            if (x.gt.tol_rho) then
c              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
c     1                                              - 3.0d0*x*d2g
c     2                                              + x*x*d3g)*fac
c            endif
c#endif
cc
cc Daniel (4-4-13): Third derivatives
c#ifdef THIRD_DERIV
c            call xc_att_xc_d3(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
cc
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
c            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
cc
c            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
c            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp
c            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp
c            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) + C6tmp
c#elif defined(SECOND_DERIV)
c            call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,C2tmp,
c     &           C3tmp)
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
c            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
c#else
c            call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
c#endif
c            Ex = Ex + qwght(n)*Etmp
c            if(ldew) func(n) = func(n) + Etmp
c            Amat(n,2) = Amat(n,2) + Atmp
c            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp
c 20      continue
cc
      endif
c
      return
      end
