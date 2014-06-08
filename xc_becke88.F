      Subroutine xc_becke88(n, tol_rho, fac, lfac, nlfac, rho, gamma, 
     &                      Amat, Amat2, Amat3, Cmat, Cmat2, Cmat11,
     &                      Cmat3, Cmat21, Cmat12,
     &                      nq, ipol, Ex, qwght,ldew,func)
c
C$Id: xc_becke88.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      implicit none
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
c
c     Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,4), Cmat3(nq,10), Cmat21(nq,9),
     1                 Cmat12(nq,12)
c
      double precision BETA
      Parameter (BETA = 0.0042D0)
c
c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************
c
      integer n
      double precision arcsinh, darcsinh
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg,
     &     dgdenom, t
c Daniel (11-10-12): Variables to add
      double precision gdenom2
c Second derivative
      double precision rhom23, d2g, d2gdenom
c Daniel (11-10-12): Variables to add
      double precision gdenom3
c Third derivative
      double precision d2arcsinh
      double precision rhom53, d3g, d3gdenom
c Daniel (11-10-12): Variables to add
      double precision gdenom4
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c For third derivatives
      d2arcsinh(x) = -x/dsqrt(1.0d0+x*x)**3
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
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            gdenom2 = gdenom*gdenom
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            if (lfac) then
               Ex = Ex + 2d0*rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.0d0*rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
c            if (nlfac) then
               Ex = Ex + 2d0*rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
c            endif
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
c               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
               Cmat(n,1) = Cmat(n,1) + t
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
            endif
c
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            gdenom3 = gdenom2*gdenom
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
     &           + BETA*x*x*d2gdenom/gdenom2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c
            Amat2(n,1) = Amat2(n,1)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat11(n,1) = Cmat11(n,1)
     &           - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,1) = Cmat2(n,1)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
c Third derivative stuff
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
            Amat3(n,1) = Amat3(n,1)
     1           + (8.0d0/27.0d0)*rhom53*(-g + x*dg
     2                                   - 18.0d0*x*x*d2g
     3                                   - 8.0d0*x*x*x*d3g)*fac   
            Cmat21(n,1) = Cmat21(n,1)
     1           + (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
     2                                          + 4.0d0*x*x*x*d3g)*fac
            Cmat12(n,1) = Cmat12(n,1)
     1           - (8.0d0/3.0d0)*(rhom23/rho(n,1)**3)/dsqrt(gamma)
     2             *d3g*fac
            Cmat3(n,1) = Cmat3(n,1)
     1           + (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
     2                                           - 3.0d0*x*d2g
     3                                           + x*x*d3g)*fac
c
c 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c         do 20 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 20
c            if (rho(n,2).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13*rho(n,2)
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
c               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
               Cmat(n,1) = Cmat(n,1) + t * 0.5d0
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
            endif
c
            if (lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,2)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Amat2(n,1) = Amat2(n,1)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c     &           - (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            Cmat11(n,1) = Cmat11(n,1)
     &           - (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            if (x.gt.tol_rho) then
c               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
               Cmat2(n,1) = Cmat2(n,1)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
c
c 25         continue
c
c           Spin beta:
c
c            if (rho(n,3).lt.tol_rho) goto 20
c
            rho13 = rho(n,3)**(1.d0/3.d0)
            rho43 = rho13*rho(n,3)
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) +rho43*g*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
c               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
               Cmat(n,2) = Cmat(n,2) + t * 0.5d0
c               Cmat(n,1,2) = Cmat(n,1,2) + t*delrho(n,1,2)
c               Cmat(n,2,2) = Cmat(n,2,2) + t*delrho(n,2,2)
c               Cmat(n,3,2) = Cmat(n,3,2) + t*delrho(n,3,2)
            endif
c
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,3)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Amat2(n,3) = Amat2(n,3)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c     &           - (2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            Cmat11(n,6) = Cmat11(n,6)
     &           - (2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            if (x.gt.tol_rho) then
c               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
               Cmat2(n,6) = Cmat2(n,6)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
c
c 20      continue
c
      endif
c
      return
      end
