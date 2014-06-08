c
c     Coulomb attenuated PBE exchange functional 
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
      Subroutine xc_camxpbe96(whichf, cam_omega, n,
     W     tol_rho, fac, lfac, nlfac, rho, gamma, 
     &     Amat, Amat2, Amat3, Cmat, Cmat2, Cmat11, Cmat3, Cmat21,
     &     Cmat12, nq, ipol, Ex, qwght,ldew,func,as)
c
c$Id: xc_camxpbe96.F 24038 2013-04-11 22:01:33Z dwsilverstein $
c
      implicit none
c
      double precision cam_omega
      double precision cam_alpha,cam_beta
      double precision as
c
      character*4 whichf
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
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
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,2), Cmat(nq,3)
      double precision Atmp, Ctmp, Etmp
c
      double precision Amat2(nq,3), Cmat2(nq,6), Cmat11(nq,6)
      double precision A2tmp, C2tmp, C3tmp
c Daniel (4-11-13): XC-3rd derivative stuff
      double precision Amat3(nq,4), Cmat3(nq,10)
      double precision Cmat21(nq,9), Cmat12(nq,12)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13
      double precision F73
      double precision F10d3, F23
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F73=7.d0/3.d0)
      parameter (F10d3=10.0d0/3.0d0, F23=2.0d0/3.0d0)
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
      double precision rhom23, d2s(3), gpp, d2g(3)
      double precision d3s(4), d3g(4), rhom53, gppp
      double precision gpbe0,gpbe1,gpbe2,gpbe3
c Original PBE
      gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
      gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
      gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
      gpbe3(s)= 24.0d0*umk*um*s*
     1  (2.0d0*umk*s*s/(1.0d0+umk*s*s)-1.0d0)/(1.0d0+umk*s*s)**3
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
c         do 10 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho
            rhom23 = F13*rho13*rrho
            rhom53 = F23*rhom23*rrho
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
c
            A2tmp = 0.0d0
c
            A3tmp = 0.0d0
            if (lfac) then
               Etmp = rho43
               Atmp = rho13
               A2tmp =  2d0*rhom23
               A3tmp = -4.0d0*rhom53
            endif
c
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 10
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
               g=gpbe0(s)
               gp=gpbe1(s)
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Etmp = Etmp + rho43*g
            Atmp = Atmp + (rho13*g+rho43*d1g(1))
            Ctmp = rho43*d1g(2)
c Daniel (4-11-13): I tried writing this in a smarter way (i.e. not
c copying and repasting code) using the preprocessor, but it lead to 
c unknown numerical issues that broke the second derivative code.  
c This seems to work correctly as written.
c Second derivatives
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
               gpp=gpbe2(s)
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            A2tmp = A2tmp
     &           +(rhom23*g + 2.d0*rho13*d1g(1) + rho43*d2g(1))
            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
            C3tmp = rho43*d2g(3)
c Third derivatives
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
               gppp = gpbe3(s)
c Third derivatives of the enhancement factor
            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
     1             + gppp*d1s(1)*d1s(1)*d1s(1)
            d3g(2) = gp*d3s(2)
     1             + gpp*d1s(2)*d2s(1)
     2             + 2.0d0*gpp*d1s(1)*d2s(2)
     3             + gppp*d1s(1)*d1s(1)*d1s(2)
            d3g(3) = gp*d3s(3)
     1             + gpp*d1s(1)*d2s(3)
     2             + 2.0d0*gpp*d1s(2)*d2s(2)
     3             + gppp*d1s(1)*d1s(2)*d1s(2)
            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
     1             + gppp*d1s(2)*d1s(2)*d1s(2)
c
            A3tmp = A3tmp
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))
            C4tmp = 
     1             (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))
            C5tmp =
     1             (rho13*d2g(3)
     2           +  rho43*d3g(3))
            C6tmp =
     1             (rho43*d3g(4))
c
            cam_alpha = 0.19
            cam_alpha = 0.46
            call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp,cam_omega,
     &           cam_alpha,cam_beta,as)
c
            Amat2(n,1) = Amat2(n,1) + A2tmp*fac
            Cmat11(n,1) = Cmat2(n,1) + C2tmp*fac
            Cmat2(n,1) = Cmat2(n,1) + C3tmp*fac
c
            Amat3(n,1) = Amat3(n,1) + A3tmp*fac
            Cmat21(n,1) = Cmat3(n,1) + C4tmp*fac
            Cmat12(n,1) = Cmat3(n,1) + C5tmp*fac
            Cmat3(n,1) = Cmat3(n,1) 
     1                              + C6tmp*fac
            Ex = Ex + qwght(n)*Etmp*fac
            if (ldew) func(n) = func(n) + Etmp*fac
            Amat(n,1) = Amat(n,1) + Atmp*fac
            Cmat(n,1) = Cmat(n,1) + Ctmp*fac
c 10      continue
c
      else
cc
cc        ======> SPIN-UNRESTRICTED <======
cc
c#ifdef IFCV81
cCDEC$ NOSWP
c#endif
c         do 20 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 20
cc
cc     Alpha
cc
c            if (rho(n,2).lt.tol_rho) goto 25
c            rho43 = C*(2d0*rho(n,2))**F43
c            rrho = 0.5d0/rho(n,2)
c            rho13 = F43*rho43*rrho
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            rhom23 = F13*rho13*rrho
c#endif
cc Daniel (4-11-13): Derivative of the electron density
c#ifdef THIRD_DERIV
c            rhom53 = F23*rhom23*rrho
c#endif
c            Etmp = 0.d0
c            Atmp = 0.d0
c            Ctmp = 0.d0
cc
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = 0.0d0
c#endif
cc
c#ifdef THIRD_DERIV
c            A3tmp = 0.0d0
c#endif
c            if (lfac) then
c               Etmp = rho43*0.5d0
c               Atmp = rho13
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c               A2tmp = 2d0*rhom23
c#endif
cc
c#ifdef THIRD_DERIV
c               A3tmp = -4.0d0*rhom53
c#endif
c            endif
cc
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 25
cc
c            s = Cs*gam12/rho43
c            d1s(1) = -F43*s*rrho
c            d1s(2) = 0.5d0*s/gamma
cc
cc     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
cc
c            if(whichf.eq.'revp') then
c               g=grevpbe0(s)
c               gp=grevpbe1(s)
c            elseif(whichf.eq.'rpbe') then
c               g=grpbe0(s)
c               gp=grpbe1(s)
c            else
c               g=gpbe0(s)
c               gp=gpbe1(s)
c            endif
cc
c            d1g(1) = gp*d1s(1)
c            d1g(2) = gp*d1s(2)
c            Etmp = Etmp + rho43*g*0.5d0
c            Atmp = Atmp + (rho13*g+rho43*d1g(1))
c            Ctmp = 0.5d0*rho43*d1g(2)
c#ifdef THIRD_DERIV
cc Second derivatives
c            d2s(1) = -F73*d1s(1)*rrho
c            d2s(2) = -F43*d1s(2)*rrho
c            d2s(3) = -0.5d0*d1s(2)/gamma
c            if(whichf.eq.'revp') then
c               gpp=grevpbe2(s)
c            elseif(whichf.eq.'rpbe') then
c               gpp=grpbe2(s)
c            else
c               gpp=gpbe2(s)
c            endif
c            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
c     &           + rho43*d2g(1))*2d0
c            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
c            C3tmp = rho43*d2g(3)*0.5d0
cc Third derivatives
cc 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c            d3s(1) = -F10d3*d2s(1)*rrho
c            d3s(2) = 0.5d0*d2s(1)/gamma
c            d3s(3) = -F43*d2s(3)*rrho
c            d3s(4) = -1.5d0*d2s(3)/gamma
c            if(whichf.eq.'revp') then
c               gppp = grevpbe3(s)
c            elseif(whichf.eq.'rpbe') then
c               gppp = grpbe3(s)
c            else
c               gppp = gpbe3(s)
c            endif
cc Third derivatives of the enhancement factor
c            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
c     1             + gppp*d1s(1)*d1s(1)*d1s(1)
c            d3g(2) = gp*d3s(2)
c     1             + gpp*d1s(2)*d2s(1)
c     2             + 2.0d0*gpp*d1s(1)*d2s(2)
c     3             + gppp*d1s(1)*d1s(1)*d1s(2)
c            d3g(3) = gp*d3s(3)
c     1             + gpp*d1s(1)*d2s(3)
c     2             + 2.0d0*gpp*d1s(2)*d2s(2)
c     3             + gppp*d1s(1)*d1s(2)*d1s(2)
c            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
c     1             + gppp*d1s(2)*d1s(2)*d1s(2)
cc
c            A3tmp = A3tmp 
c     1           + (-rhom53*g
c     2           +  3.0d0*rhom23*d1g(1)
c     3           +  3.0d0*rho13*d2g(1)
c     4           +  rho43*d3g(1))*4.0d0
c            C4tmp =
c     1             (rhom23*d1g(2)
c     2           +  2.0d0*rho13*d2g(2)
c     3           +  rho43*d3g(2))*2.0d0
c            C5tmp = 
c     1             (rho13*d2g(3)
c     2           +  rho43*d3g(3))
c            C6tmp = 
c     1             (rho43*d3g(4))*0.5d0
cc
c            call xc_att_xc_d3(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
cc
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
c            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
cc
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp*fac
c            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp*fac
c            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp*fac
c            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA)
c     1                              + C6tmp*fac
c#elif defined(SECOND_DERIV)
c            d2s(1) = -F73*d1s(1)*rrho
c            d2s(2) = -F43*d1s(2)*rrho
c            d2s(3) = -0.5d0*d1s(2)/gamma
c            if(whichf.eq.'revp') then
c               gpp=grevpbe2(s)
c            elseif(whichf.eq.'rpbe') then
c               gpp=grpbe2(s)
c            else
c               gpp=gpbe2(s)
c            endif
c            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
c     &           + rho43*d2g(1))*2d0
c            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
c            C3tmp = rho43*d2g(3)*0.5d0
cc
c            call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp)
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp*fac
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp*fac
c            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp*fac
c#else
c            call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
c#endif
c            Ex = Ex + qwght(n)*Etmp*fac
c            if (ldew) func(n) = func(n) + Etmp*fac
c            Amat(n,1) = Amat(n,1) + Atmp*fac
c            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp*fac
cc
cc     Beta
cc
c 25         continue
c            if (rho(n,3).lt.tol_rho) goto 20
c            rho43 = C*(2d0*rho(n,3))**F43
c            rrho = 0.5d0/rho(n,3)
c            rho13 = F43*rho43*rrho
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            rhom23 = F13*rho13*rrho
c#endif
cc Daniel (4-11-13): Derivative of the electron density
c#ifdef THIRD_DERIV
c            rhom53 = F23*rhom23*rrho
c#endif
c            Etmp = 0.d0
c            Atmp = 0.d0
c            Ctmp = 0.d0
cc
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            A2tmp = 0.0d0
c#endif
cc
c#ifdef THIRD_DERIV
c            A3tmp = 0.0d0
c#endif
c            if (lfac) then
c               Etmp = rho43*0.5d0
c               Atmp = rho13
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c               A2tmp= 2.d0*rhom23
c#endif
cc
c#ifdef THIRD_DERIV
c               A3tmp = -4.0d0*rhom53*fac
c#endif
c            endif
cc
c            gamma = delrho(n,1,2)*delrho(n,1,2) +
c     &              delrho(n,2,2)*delrho(n,2,2) +
c     &              delrho(n,3,2)*delrho(n,3,2)
c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 20
cc
c            s = Cs*gam12/rho43
c            d1s(1) = -F43*s*rrho
c            d1s(2) = 0.5d0*s/gamma
cc
cc     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
cc
c            if(whichf.eq.'revp') then
c               g=grevpbe0(s)
c               gp=grevpbe1(s)
c            elseif(whichf.eq.'rpbe') then
c               g=grpbe0(s)
c               gp=grpbe1(s)
c            else
c               g=gpbe0(s)
c               gp=gpbe1(s)
c            endif
cc
c            d1g(1) = gp*d1s(1)
c            d1g(2) = gp*d1s(2)
c            Etmp = Etmp + rho43*g*0.5d0
c            Atmp = Atmp + (rho13*g+rho43*d1g(1))
c            Ctmp = 0.5d0*rho43*d1g(2)
c#ifdef THIRD_DERIV
cc Second derivatives
c            d2s(1) = -F73*d1s(1)*rrho
c            d2s(2) = -F43*d1s(2)*rrho
c            d2s(3) = -0.5d0*d1s(2)/gamma
c            if(whichf.eq.'revp') then
c               gpp=grevpbe2(s)
c            elseif(whichf.eq.'rpbe') then
c               gpp=grpbe2(s)
c            else
c               gpp=gpbe2(s)
c            endif
c            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
cc
c            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
c     &           + rho43*d2g(1))*2d0
c            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
c            C3tmp =  rho43*d2g(3)*0.5d0
cc Third derivatives
cc 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
c            d3s(1) = -F10d3*d2s(1)*rrho
c            d3s(2) = 0.5d0*d2s(1)/gamma
c            d3s(3) = -F43*d2s(3)*rrho
c            d3s(4) = -1.5d0*d2s(3)/gamma
c            if(whichf.eq.'revp') then
c               gppp = grevpbe3(s)
c            elseif(whichf.eq.'rpbe') then
c               gppp = grpbe3(s)
c            else
c               gppp = gpbe3(s)
c            endif
cc Third derivatives of the enhancement factor
c            d3g(1) = gp*d3s(1) + 3.0d0*gpp*d1s(1)*d2s(1)
c     1             + gppp*d1s(1)*d1s(1)*d1s(1)
c            d3g(2) = gp*d3s(2)
c     1             + gpp*d1s(2)*d2s(1)
c     2             + 2.0d0*gpp*d1s(1)*d2s(2)
c     3             + gppp*d1s(1)*d1s(1)*d1s(2)
c            d3g(3) = gp*d3s(3)
c     1             + gpp*d1s(1)*d2s(3)
c     2             + 2.0d0*gpp*d1s(2)*d2s(2)
c     3             + gppp*d1s(1)*d1s(2)*d1s(2)
c            d3g(4) = gp*d3s(4) + 3.0d0*gpp*d1s(2)*d2s(3)
c     1             + gppp*d1s(2)*d1s(2)*d1s(2)
cc
c            A3tmp = A3tmp
c     1           + (-rhom53*g
c     2           +  3.0d0*rhom23*d1g(1)
c     3           +  3.0d0*rho13*d2g(1)
c     4           +  rho43*d3g(1))*4.0d0
c            C4tmp = 
c     1             (rhom23*d1g(2)
c     2           +  2.0d0*rho13*d2g(2)
c     3           +  rho43*d3g(2))*2.0d0
c            C5tmp =
c     1             (rho13*d2g(3)
c     2           +  rho43*d3g(3))
c            C6tmp = 
c     1             (rho43*d3g(4))*0.5d0
cc
c            call xc_att_xc_d3(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
cc
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*fac
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*fac
c            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*fac
cc
c            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp*fac
c            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp*fac
c            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp*fac
c            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB)
c     1                              + C6tmp*fac
c#elif defined(SECOND_DERIV)
c            d2s(1) = -F73*d1s(1)*rrho
c            d2s(2) = -F43*d1s(2)*rrho
c            d2s(3) = -0.5d0*d1s(2)/gamma
c            if(whichf.eq.'revp') then
c               gpp=grevpbe2(s)
c            elseif(whichf.eq.'rpbe') then
c               gpp=grpbe2(s)
c            else
c               gpp=gpbe2(s)
c            endif
c            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
c            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
c            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
cc
c            A2tmp = A2tmp + (rhom23*g + 2.d0*rho13*d1g(1)
c     &           + rho43*d2g(1))*2d0
c            C2tmp = (rho13*d1g(2) + rho43*d2g(2))
c            C3tmp =  rho43*d2g(3)*0.5d0
cc
c            call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp)
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp*fac
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp*fac
c            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp*fac
c#else
c            call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
c#endif
c            Ex = Ex + qwght(n)*Etmp*fac
c            if (ldew) func(n) = func(n) + Etmp*fac
c            Amat(n,2) = Amat(n,2) + Atmp*fac
c            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp*fac
cc
c 20      continue
      endif
c
      return
      end
