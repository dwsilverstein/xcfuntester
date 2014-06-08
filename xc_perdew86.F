      Subroutine xc_perdew86(n, tol_rho, fac, lfac, nlfac, rho, gamma, 
     1                       Amat, Amat2, Amat3, Cmat, Cmat2, Cmat11,
     2                       Cmat3, Cmat21, Cmat12,
     2                       nq, ipol, Ec, qwght, ldew, ffunc)
c
c$Id: xc_perdew86.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      implicit none
      double precision tol_rho, fac ! [input]
      integer nq, ipol              ! [input]
      double precision Ec           ! [input/output]
      logical lfac, nlfac, ldew
      double precision ffunc(nq)  ! value of the functional [output]
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(n,2), Cmat(n,3)
      double precision Amat2(n,3), Cmat2(n,6), Cmat11(n,6)
      double precision Amat3(n,4), Cmat21(n,9), Cmat12(n,12),
     1                 Cmat3(n,10)
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONE3, THREE, FOUR3, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
c Third derivative stuff
      double precision SEVEN3, EIGHT3
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOUR3 = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
c Third derivative stuff
      parameter (SEVEN3 = 7.0d0/3.0d0, EIGHT3 = 8.0d0/3.0d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision d1rs
      double precision d2rs
      double precision gamma, gam12
      double precision anum, aden, d1anum, d1aden, Cn, d1Cn,
     &     expfac, phi, d1phi(2), dlnphi, func, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision zeta, d1z(2), d, dm1, adp, d1d(2), t,
     &     dt12, d1dt12
      double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
     &     d2phi(3), d2dt12, d2Cn
c Daniel (4/1/13): New variables
      double precision d2lnphi
      double precision d2f(3)
      double precision d2lnfrho(3), d2lnfrg(2), d2lnfgam
      double precision arho2, arho3
      double precision aden2, aden3, aden4
c Third derivative stuff
      double precision d3lnphi
      double precision d3rs, d3anum, d3aden, d3Cn, d3phi(4)
      double precision d3lnfrho(4), d3lnfgam
      double precision d3f(3)
c
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
c         do 10 n = 1, nq
            rhoval = rho(n,1)
c            if (rhoval.lt.tol_rho) goto 10
            arho = 1.d0/rhoval
            rho13 = abs(rhoval)**ONE3
            rho43 = rhoval*rho13
            rho76 = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
            d2rs = -FOUR3*d1rs*arho
c Third derivatives
            d3rs = -SEVEN3*d2rs*arho
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = sqrt(abs(gamma))
c         
c           C(n)
c         
c Daniel (3-29-13): anum = P in my notes, aden = Q in my notes
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
c Daniel (3-29-13): pgamma is eta in my notes
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
c
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
c Third derivatives
            d3anum = 0.0d0
            d3aden = 6.0d0*beta10 
c     First compute rs derivative of Cn
c Daniel (4/1/13): New variable
            aden2 = aden*aden
c
c            d1Cn = d1anum/aden - anum*d1aden/aden**2
            d1Cn = d1anum/aden - anum*d1aden/aden2
c
c Daniel (4/1/13): New variable
            aden3 = aden2*aden
c
c            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
c     &           + 2d0*anum*d1aden**2/aden**3
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden2
     &           + 2d0*anum*d1aden**2/aden3
c Third derivatives
c Daniel (4/1/13): New variable
            aden4 = aden3*aden
c
            d3Cn = -( 3.0d0*d2anum*d1aden + 3.0d0*d1anum*d2aden 
     1              + anum*d3aden )/aden2
     2           + 6.0d0*( d1anum*d1aden**2
     3                   + anum*d2aden*d1aden )/aden3
     4           - 6.0d0*anum*d1aden**3/aden4
c     Convert to rho derivative
c
c Third derivatives
            d3Cn = d3Cn*d1rs*d1rs*d1rs
     1           + 3.0d0*d2Cn*d2rs*d1rs
     2           + d1Cn*d3rs
c
            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
c
            d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
            d1phi(1) = phi*dlnphi
c
c Daniel (4/1/13): Added a variable for the second derivative.  This
c makes it easier to write the third derivative also.
c Daniel (4/1/13): New variable
            arho2 = arho*arho
c
            d2lnphi = (d1Cn/Cn)**2 - d2Cn/Cn + SEV6*arho2
c            d2phi(1) = d1phi(1)*dlnphi
c     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
            d2phi(1) = d1phi(1)*dlnphi + phi*d2lnphi
c Third derivatives
c Daniel (4/1/13): New variable
            arho3 = arho2*arho
c
            d3lnphi = -2.0d0*(d1Cn/Cn)**3
     1              + 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     2              - d3Cn/Cn
     3              - SEVEN3*arho3
            d3phi(1) = d2phi(1)*dlnphi
     1               + 2.0d0*d1phi(1)*d2lnphi 
     2               + phi*d3lnphi
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43
            dlnfrho(1) = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
            d1f(1) = dlnfrho(1)*func
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            if (gam12.gt.TOLL)then
               d1phi(2) = phi / (2d0*gamma)
               dlnfgam = 1d0/gamma - d1phi(2)
               d1f(3) = func*dlnfgam
c               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
               Cmat(n,1) = Cmat(n,1) + d1f(3)*fac
c               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*fac
               Cmat(n,2) = Cmat(n,2) + d1f(3)*2D0*fac
c
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
c Daniel (11-19-12): d is a constant (it equals 1) for a restricted 
c calculation, since there is no spin-polarization.  Thus, the
c derivatives are zero.
c               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c     &              + (d1f(1)*dlnfrho(1)
c     &              + func*t)*fac
c               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c     &              + (d1f(1)*dlnfrho(1)
c     &              + func*t)*fac
               d2lnfrho(1) = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn 
     1                     + FOUR3*arho2
               d2f(1) = d1f(1)*dlnfrho(1)
     1                + func*d2lnfrho(1)
               Amat2(n,1) = Amat2(n,1)
     &              + d2f(1)*fac
               Amat2(n,2) = Amat2(n,2)
     &              + d2f(1)*fac
c
c               t = (d1f(1)*dlnfgam - func*d2phi(2))*fac
c               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
c               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
c               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
               d2lnfrg(1) = -d2phi(2)
               d2f(2) = (d1f(1)*dlnfgam + func*d2lnfrg(1))
               t = d2f(2)*fac
               Cmat11(n,1) = Cmat11(n,1) + t 
c               Cmat11(n,2) = Cmat11(n,2) + d2lnfrg(1)*2d0
               Cmat11(n,2) = Cmat11(n,2) + t
               Cmat11(n,3) = Cmat11(n,3) + t
c
c               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*fac
c               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
c               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
c               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
c               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
               d2lnfgam = -1.0d0/gamma**2 - d2phi(3) 
               d2f(3) = d1f(3)*dlnfgam + func*d2lnfgam
               t = d2f(3)*fac
               Cmat2(n,1) = Cmat2(n,1) + t
               Cmat2(n,2) = Cmat2(n,2) + t
               Cmat2(n,3) = Cmat2(n,3) + t*2d0
               Cmat2(n,4) = Cmat2(n,4) + t*4d0
c Third derivatives
c rrr terms
               d3lnfrho(1) = -d3phi(1) 
     1                     + 2.0d0*(d1Cn/Cn)**3
     2                     - 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     3                     + d3Cn/Cn
     4                     - EIGHT3*arho3
c
               d3f(1) = d2f(1)*dlnfrho(1)
     1                + 2.0d0*d1f(1)*d2lnfrho(1)
     2                + func*d3lnfrho(1)
c
               Amat3(n,1) = Amat3(n,1)
     1                    + d3f(1)*fac
c
               Amat3(n,2) = Amat3(n,2)
     1                    + d3f(1)*fac
c           
               Amat3(n,3) = Amat3(n,3)
     1                    + d3f(1)*fac
c rrg terms
               d3phi(2) = d2phi(2)*dlnphi + d1phi(2)*d2lnphi
c
               t = ( d2f(2)*dlnfrho(1)
     1             - d1f(1)*d2phi(2)
     2             + d1f(3)*d2lnfrho(1)
     3             - func*d3phi(2) )*fac
c
               Cmat21(n,1) = Cmat21(n,1) + t                  ! raragaa
c
               Cmat21(n,2) = Cmat21(n,2) + t*2.0d0            ! raragab
c
               Cmat21(n,3) = Cmat21(n,3) + t                  ! raragbb
c
               Cmat21(n,4) = Cmat21(n,4) + t                  ! rarbgaa
c
               Cmat21(n,5) = Cmat21(n,5) + t*2.0d0            ! rarbgab
c
               Cmat21(n,6) = Cmat21(n,6) + t                  ! rarbgbb
c
c               Cmat21(n,7) = Cmat21(n,7) + t                  ! rbrbgaa
cc
c               Cmat21(n,8) = Cmat21(n,8) + t*2.0d0            ! rbrbgab
cc
c               Cmat21(n,9) = Cmat21(n,9) + t                  ! rbrbgbb
c rgg terms
               d3phi(3) = -d2phi(3)*dlnphi 
c
               t = ( d2f(2)*dlnfgam 
     1             + d1f(1)*d2lnfgam
     2             + d1f(3)*d2lnfrg(1)
     3             + func*d3phi(3) )*fac
c
               Cmat12(n,1) = Cmat12(n,1) + t                   ! ragaagaa
c
               Cmat12(n,2) = Cmat12(n,2) + t*2.0d0             ! ragaagab
c
               Cmat12(n,3) = Cmat12(n,3) + t                   ! ragaagbb
c
               Cmat12(n,4) = Cmat12(n,4) + t*4.0d0             ! ragabgab
c
               Cmat12(n,5) = Cmat12(n,5) + t*2.0d0             ! ragabgbb
c
               Cmat12(n,6) = Cmat12(n,6) + t                   ! ragbbgbb
c
c ggg terms              
               d3phi(4) = -3.0d0*d2phi(3)/(2.0d0*gamma) 
               d3lnfgam = 2.0d0/gamma**3 - d3phi(4) 
c
               t = ( d2f(3)*dlnfgam
     1             + 2.0d0*d1f(3)*d2lnfgam
     2             + func*d3lnfgam )*fac
c
               Cmat3(n,1) = Cmat3(n,1) + t                      ! gaagaagaa
c
               Cmat3(n,2) = Cmat3(n,2) + t*2.0d0                ! gaagaagab
c
               Cmat3(n,3) = Cmat3(n,3) + t                      ! gaagaagbb
c
               Cmat3(n,4) = Cmat3(n,4) + t*4.0d0                ! gaagabgab
c
               Cmat3(n,5) = Cmat3(n,5) + t*2.0d0                ! gaagabgbb
c
               Cmat3(n,6) = Cmat3(n,6) + t                      ! gaagbbgbb
c
               Cmat3(n,7) = Cmat3(n,7) + t*8.0d0                ! gabgabgab
c
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c         do 20 n = 1, nq
c            rhoval = rho(n,1)
c            if (rhoval.lt.tol_rho) goto 20
c            arho=1.d0/rhoval
c            rho13  = abs(rhoval)**ONE3
c            rho43  = rhoval*rho13
c            rho76  = abs(rhoval)**SEV6
c            rs = rsfact/rho13
c            rs2 = rs*rs
c            rs3 = rs2*rs
c            d1rs = -ONE3*rs*arho
c#ifdef SECOND_DERIV
c            d2rs = -FOUR3*d1rs*arho
c#endif
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1) +
c     &              delrho(n,1,2)*delrho(n,1,2) +
c     &              delrho(n,2,2)*delrho(n,2,2) +
c     &              delrho(n,3,2)*delrho(n,3,2) +
c     &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c     &              delrho(n,2,1)*delrho(n,2,2) +
c     &              delrho(n,3,1)*delrho(n,3,2))
c            gam12 = sqrt(abs(gamma))
c            zeta = (rho(n,2) - rho(n,3))*arho
c            if(zeta.lt.-1d0) zeta=-1d0
c            if(zeta.gt.1d0) zeta=1d0
c            d1z(1) =  (1.d0 - zeta)*arho
c            d1z(2) = -(1.d0 + zeta)*arho
c#ifdef SECOND_DERIV
c            rrho2 = 2.d0*arho*arho
cc           1 = aa, 2 = ab, 3 = bb
c            d2z(1) =-rrho2*(1.d0-zeta)
c            d2z(2) = rrho2*zeta
c            d2z(3) = rrho2*(1.d0+zeta)
c#endif
cc
cc           d(zeta)
cc
c            dt12 = ((ONE+zeta)*.5d0)**FIVE3 + ((ONE-zeta)*.5d0)**FIVE3
c            d1dt12 = FIVE3*0.5d0*(
c     &           ((ONE+zeta)*.5d0)**TWO3 - ((ONE-zeta)*.5d0)**TWO3 )
c            d = 2.d0**ONE3*dsqrt(dt12)
c            dm1 = 1.d0/d
c            adp = 0.5d0*d/dt12*d1dt12
c            d1d(1) = adp*d1z(1)
c            d1d(2) = adp*d1z(2)
c#ifdef SECOND_DERIV
c            d2dt12 = FIVE3*TWO3*0.25d0*(
c     &         ((ONE+zeta)*.5d0)**(-ONE3) + ((ONE-zeta)*.5d0)**(-ONE3) )
c            dpp =-0.5d0*adp/dt12*d1dt12
c     &        + 2.d0**(-TWO3)*d2dt12/dsqrt(dt12)
c            d2d(1) = dpp*d1z(1)*d1z(1) + adp*d2z(1)
c            d2d(2) = dpp*d1z(1)*d1z(2) + adp*d2z(2)
c            d2d(3) = dpp*d1z(2)*d1z(2) + adp*d2z(3)
c#endif
cc         
cc           C(n)
cc         
c            anum = fff+alpha*rs+beta*rs2
c            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
c            Cn = zzz + anum/aden
c            d1anum = alpha + 2d0*beta*rs
c            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
c#ifdef SECOND_DERIV
c            d2anum = 2d0*beta
c            d2aden = 2d0*delta + 6d0*beta10*rs
c#endif
cc     First compute rs derivative
c            d1Cn = d1anum/aden - anum*d1aden/aden**2
c#ifdef SECOND_DERIV
c            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
c     &           + 2d0*anum*d1aden**2/aden**3
c#endif
cc     Convert to rho derivative
c#ifdef SECOND_DERIV
c            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
c#endif
c            d1Cn = d1Cn*d1rs
cc         
cc           phi(n,gradn)
cc         
c            expfac = 0.d0
c            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
c            if (phi.lt.EXPTOL) expfac = exp(-phi)
c            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
c            d1phi(1) = phi*dlnphi
c#ifdef SECOND_DERIV
c            d2phi(1) = d1phi(1)*dlnphi
c     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
c#endif
cc         
cc           functional
cc         
c            func = expfac*Cn*gamma/rho43*dm1
c            t = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
c            dlnfrho(1) = t - dm1*d1d(1)
c            dlnfrho(2) = t - dm1*d1d(2)
c            d1f(1) = dlnfrho(1)*func
c            d1f(2) = dlnfrho(2)*func
c            Amat(n,1) = Amat(n,1) + d1f(1)*fac
c            Amat(n,2) = Amat(n,2) + d1f(2)*fac
c            if (gam12.gt.TOLL)then
c               d1phi(2) = phi / (2d0*gamma)
c               dlnfgam = 1d0/gamma - d1phi(2)
c               d1f(3) = func*dlnfgam
c               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
c               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*fac
c               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(3)*fac
c#ifdef SECOND_DERIV
c               d2phi(2) = d1phi(2)*dlnphi
c               d2phi(3) =-d1phi(2)/(2d0*gamma)
cc
c               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c     &              + (d1f(1)*dlnfrho(1)
c     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*fac
c               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c     &              + (d1f(1)*dlnfrho(2)
c     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*fac
c               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c     &              + (d1f(2)*dlnfrho(2)
c     &              + func*(d1d(2)*d1d(2)*dm1**2-d2d(3)*dm1+t))*fac
cc
c               t = (d1f(1)*dlnfgam - func*d2phi(2))*fac
c               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
c               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
c               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c               t = (d1f(2)*dlnfgam - func*d2phi(2))*fac
c               Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + t
c               Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + t*2d0
c               Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + t
cc
c               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*fac
c               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
c               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
c               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + t
c               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
c               Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB) + t*2d0
c               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
c#endif
c            endif
c            Ec = Ec + func*qwght(n)*fac
c            if (ldew) ffunc(n)=ffunc(n)+func*fac
c   20    continue
      endif
      return
      end

c#ifndef SECOND_DERIV
c      Subroutine xc_p81(tol_rho, fac, lfac, nlfac, rho, Amat, nq, ipol, 
c     &                  Ec, qwght, ldew, func)
c#else
c#include "dft2drv.fh"
c      Subroutine xc_p81_d2(tol_rho, fac, lfac, nlfac, rho, Amat, Amat2,
c     &                     nq, ipol, Ec, qwght, ldew, func)
c#endif
cc
cc     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
cc
c      implicit none
cc
c      integer nq, ipol
c      logical lfac, nlfac, ldew
c      double precision func(*)  ! value of the functional [output]
c      double precision Ec, fac
cc
cc     Charge Density 
cc
c      double precision rho(nq,ipol*(ipol+1)/2)
cc
cc     Quadrature Weights
cc
c      double precision qwght(nq)
cc
cc     Sampling Matrices for the XC Potential & Energy
cc
c      double precision Amat(nq,ipol)
c#ifdef SECOND_DERIV
c      double precision Amat2(nq,*)
c#endif
c      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
c     &                 pi, tol_rho, ONE3, FOUR3, TWO3
c      save A, B, C, D, G, B1, B2
c      parameter (pi = 3.1415926535897932385d0)
c      Parameter (ONE3 = 1.d0/3.d0, FOUR3 = 4.D0/3.D0)
c      Parameter (TWO3 = 2.d0/3.d0)
c      integer n, i
c      double precision rhoval, rs, alnrs, d1rs, e(2), d1e(2), rden(2),
c     &                 d1den(2), d1zeta(2), d1ersz(2), d1edrho(2), eps,
c     &                 sqrtrs, fz, d1fz, zeta
c#ifdef SECOND_DERIV
c      double precision d2rs, d2e(2), d2den(2), d2zeta(3), d2ersz(3),
c     &                 d2edrho(3), d2fzeta, d2fz, rrho2
c#endif
c      double precision x, fzeta, d1fzeta, rsfact
c      fzeta(x) = ((1.d0+x)**FOUR3 + 
c     &            (1.d0-x)**FOUR3 - 2.d0) / (2.d0**FOUR3-2.d0)
c      d1fzeta(x) = FOUR3*((1.d0+x)**ONE3 - 
c     &                    (1.d0-x)**ONE3) / (2.d0**FOUR3-2.d0)
c#ifdef SECOND_DERIV
c      d2fzeta(x) = ONE3*FOUR3*((1.d0+x)**(-TWO3) + 
c     &                    (1.d0-x)**(-TWO3)) / (2.d0**FOUR3-2.d0)
c#endif
c      data A / 0.0311d0, 0.01555d0 /
c      data B / -0.048d0, -0.0269d0 /
c      data C / 0.0020d0, 0.0007d0 /
c      data D / -0.0116d0, -0.0048d0 /
c      data G / -.1423d0, -.0843d0 /
c      data B1 / 1.0529d0, 1.3981d0 /
c      data B2 / 0.3334d0, 0.2611d0 /
cc
c      rsfact = (0.75d0/pi)**ONE3
cc
cc     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
cc
c      do n = 1, nq
c         if (rho(n,1).gt.tol_rho)then
c            rhoval = rho(n,1)
c            if (ipol.eq.1) then
c               zeta = 0.0d0
c               d1zeta(1) = 1.d0/rhoval
c               d1zeta(2) =-1.d0/rhoval
c               fz = 0d0
c               d1fz = 0d0
c            else
c               zeta = (rho(n,2)-rho(n,3))/rhoval
c               if(zeta.lt.-1d0) zeta=-1d0
c               if(zeta.gt.1d0) zeta=1d0
c               fz = fzeta(zeta)
c               d1fz = d1fzeta(zeta)
c               d1zeta(1) = (1.d0-zeta)/rhoval
c               d1zeta(2) =-(1.d0+zeta)/rhoval
c            endif
c            rs = rsfact/abs(rhoval)**ONE3
c            d1rs = -ONE3*rs/rhoval
c#ifdef SECOND_DERIV
c            d2rs = -FOUR3*d1rs/rhoval
c            d2fz = d2fzeta(zeta)
c            rrho2 = 2.d0/(rhoval*rhoval)
cc           1 = aa, 2 = ab, 3 = bb
c            d2zeta(1) =-rrho2*(1.d0-zeta)
c            d2zeta(2) = rrho2*zeta
c            d2zeta(3) = rrho2*(1.d0+zeta)
c#endif
c            if (rs.lt.1.d0)then
c               alnrs = log(rs)
c               do i = 1, 2
c                  e(i) = A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs
c                  d1e(i) = A(i)/rs+C(i)*(alnrs+1d0)+D(i)
c#ifdef SECOND_DERIV
c                  d2e(i) = (C(i)-A(i)/rs)/rs
c#endif
c               enddo
c            else
c               sqrtrs = sqrt(rs)
c               do i = 1, 2
c                  rden(i) = 1.d0/(1.d0+B1(i)*sqrtrs+B2(i)*rs)
c                  d1den(i) = B1(i)/(2.d0*sqrtrs)+B2(i)
c                  e(i) = G(i)*rden(i)
c                  d1e(i) = -G(i)*d1den(i)*rden(i)**2
c#ifdef SECOND_DERIV
c                  d2den(i) = -B1(i)/(4.d0*rs*sqrtrs)
c                  d2e(i) = G(i)*rden(i)**2
c     &                 *(2.d0*d1den(i)**2*rden(i)-d2den(i))
c#endif
c               enddo
c            endif
c            eps = e(1) + fz*(e(2)-e(1))
c            d1ersz(1) = d1e(1) + fz*(d1e(2)-d1e(1))
c            d1ersz(2) = d1fz*(e(2)-e(1))
c            d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
c            d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
c            Ec = Ec + eps*qwght(n)*rhoval*fac
c            if (ldew) func(n) = func(n) + eps*rhoval*fac
c            Amat(n,1) = Amat(n,1) + (eps + rhoval*d1edrho(1))*fac
c            if (ipol.eq.2)
c     &      Amat(n,2) = Amat(n,2) + (eps + rhoval*d1edrho(2))*fac
c#ifdef SECOND_DERIV
cc           1 = rsrs, 2 = rsz, 3 = zz
c            d2ersz(1) = d2e(1) + fz*(d2e(2)-d2e(1))
c            d2ersz(2) = d1fz*(d1e(2)-d1e(1))
c            d2ersz(3) = d2fz*(e(2)-e(1))
cc           1 = aa, 2 = ab, 3 = bb
c            d2edrho(1) = d2ersz(1)*d1rs*d1rs
c     &                 + d2ersz(2)*d1rs*d1zeta(1)*2.d0
c     &                 + d2ersz(3)*d1zeta(1)*d1zeta(1)
c     &                 + d1ersz(1)*d2rs
c     &                 + d1ersz(2)*d2zeta(1)
c            d2edrho(2) = d2ersz(1)*d1rs*d1rs
c     &                 + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
c     &                 + d2ersz(3)*d1zeta(1)*d1zeta(2)
c     &                 + d1ersz(1)*d2rs
c     &                 + d1ersz(2)*d2zeta(2)
c            d2edrho(3) = d2ersz(1)*d1rs*d1rs
c     &                 + d2ersz(2)*d1rs*d1zeta(2)*2.d0
c     &                 + d2ersz(3)*d1zeta(2)*d1zeta(2)
c     &                 + d1ersz(1)*d2rs
c     &                 + d1ersz(2)*d2zeta(3)
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c     &           + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
c            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c     &           + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*fac
c            if (ipol.eq.2)
c     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c     &           + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*fac
c#endif
c         endif
c      enddo
c      return
c      end
cc
c#ifndef SECOND_DERIV
c#define SECOND_DERIV
cc
cc     Compile source again for the 2nd derivative case
cc
c#include "xc_perdew86.F"
c#endif
