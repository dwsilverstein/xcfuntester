c     PBE exchange functional 
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
      Subroutine xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, gamma,
     &     Amat, Amat2, Amat3, Cmat, Cmat11, Cmat2, Cmat21, Cmat12, 
     &     Cmat3, nq, ipol, Ex, qwght, ldew, func)
c
c$Id: xc_xpbe96.F 19696 2012-07-27 15:13:42Z dws219 $
c
      implicit none
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
      double precision Amat(n,2), Cmat(n,3)
      double precision Amat2(n,3), Cmat2(n,6), Cmat11(n,6)
      double precision Amat3(n,4), Cmat21(n,9), Cmat12(n,12), 
     1                 Cmat3(n,10)
c
      double precision tol_rho, pi, um, uk, umk,ukrev,umkrev
      double precision C, Cs
      double precision F43, F13, F23
c Daniel (7-29-12): XC-3rd derivatives
      double precision F73
c Daniel (7-27-12): XC-3rd derivative information for d3s
      double precision F10d3
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0, F23=2.0d0/3.0d0)
c Daniel (7-29-12): XC-3rd derivatives
      parameter (F73=7.d0/3.d0)
c Daniel (7-27-12)
      parameter (F10d3=10.0d0/3.0d0)
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      d, g, gp, d1g(2)
c Daniel (7-29-12): XC-third derivatives
      double precision rhom23, d2s(3), gpp, d2g(3)
c Daniel (7-27-12): Variables needed for XC-third derivatives
      double precision d3s(4), d3g(4), rhom53, gppp
      double precision gpbe0, gpbe1, gpbe2, gpbe3
      double precision grpbe0, grpbe1, grpbe2, grpbe3
      double precision grevpbe0, grevpbe1, grevpbe2, grevpbe3
c Original PBE
      gpbe0(s)= uk*(1d0 - 1d0/(1d0+umk*s*s))
      gpbe1(s)= 2d0*um*s/(1d0+umk*s*s)**2
      gpbe2(s)= 2d0*um*(1d0-4d0*umk*s*s/(1d0+umk*s*s))/(1d0+umk*s*s)**2
      gpbe3(s)= 24.0d0*umk*um*s*
     1  (2.0d0*umk*s*s/(1.0d0+umk*s*s)-1.0d0)/(1.0d0+umk*s*s)**3
c revPBE by Zhang et al.
      grevpbe0(s)= ukrev*(1d0 - 1d0/(1d0+umkrev*s*s))
      grevpbe1(s)= 2d0*um*s/(1d0+umkrev*s*s)**2
      grevpbe2(s)= 2d0*um*(1d0-4d0*umkrev*s*s/(1d0+umkrev*s*s))/
     /     (1d0+umkrev*s*s)**2
      grevpbe3(s)= 24.0d0*umkrev*um*s*
     1  (2.0d0*umkrev*s*s/(1.0d0+umkrev*s*s)-1.0d0)/
     2  (1.0d0+umkrev*s*s)**3
c RPBE by Hammer et al.
      grpbe0(s)= uk*(1d0 - exp(-umk*s*s))
      grpbe1(s)= 2d0*um*s*exp(-umk*s*s)
      grpbe2(s)= 2d0*um*exp(-umk*s*s)*(1d0-2d0*umk*s*s)
      grpbe3(s)= -4.0d0*umk*um*s*exp(-umk*s*s)*(3d0-2d0*umk*s*s)
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
c Daniel (7-29-12): Needed for XC-third derivatives also 
            rhom23 = F13*rho13*rrho
c Daniel (7-27-12): Electron density for the XC-third derivatives
            rhom53 = F23*rhom23*rrho
c Daniel (7-27-12): This is used only if the local part coming from LDA
c isn't done already.
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
c               Amat(n,1) = Amat(n,1) + rho13*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
c Daniel (7-27-12): XC-third derivative, coming from LDA.
c               Amat3(n,D2_RA_RA_RA) = Amat3(n,D2_RA_RA_RA) 
c     1                              - 4.0d0*rhom53*fac
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
            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c
c Daniel (7-27-12): gp is the derivative of the rational function,
c or whatever the function in the revision is.
c First derivatives of the enhancement factor
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,1) = Cmat(n,1) + rho43*d1g(2)*fac
c            Cmat(n,1) = Cmat(n,1) + 2d0*rho43*d1g(2)*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
c Second derivatives of the enhancement factor
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,1) = Amat2(n,1)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac
            Cmat11(n,1) = Cmat11(n,1)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,1) = Cmat2(n,1)
     &           + rho43*d2g(3)*fac
c            Amat2(n,1) = Amat2(n,1)
c     &           +(rhom23*g
c     &           + 2.d0*rho13*d1g(1)
c     &           + rho43*d2g(1))*fac*2d0
c            Cmat11(n,1) = Cmat11(n,1)
c     &           +(rho13*d1g(2)
c     &           + rho43*d2g(2))*fac*4d0
c            Cmat2(n,1) = Cmat2(n,1)
c     &           + rho43*d2g(3)*fac*8d0
c Daniel (7-27-12): XC-third derivatives of different PBEs.
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
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
c Currently guessing the coefficients based on what I did for the
c Slater/Dirac functional
            Amat3(n,1) = Amat3(n,1)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac
            Cmat21(n,1) = Cmat21(n,1)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac
            Cmat12(n,1) = Cmat12(n,1)
     1           + (rho13*d2g(3) 
     2           +  rho43*d3g(3))*fac
            Cmat3(n,1) = Cmat3(n,1)
     1           + (rho43*d3g(4))*fac
c            Amat3(n,1) = Amat3(n,1)
c     1           + (rhom53*g
c     2           -  3.0d0*rhom23*d1g(1)
c     3           -  3.0d0*rho13*d2g(1)
c     4           +  rho43*d3g(1))*fac*4.0d0
c            Cmat21(n,1) = Cmat21(n,1)
c     1           + (rhom23*d1g(2)
c     2           -  2.0d0*rho13*d2g(2)
c     3           +  rho43*d3g(2))*fac*8.0d0
c            Cmat12(n,1) = Cmat12(n,1)
c     1           + (rho13*d2g(3)
c     2           +  rho43*d3g(3))*fac*8.0d0
c            Cmat3(n,1) = Cmat3(n,1)
c     1           + (rho43*d3g(4))*fac*16.0d0
c 10      continue
c 
c            write(6,*) 'Ex', Ex
c            write(6,*) 'Amat', Amat(n,1)
c            write(6,*) 'Cmat', Cmat(n,1)
c            write(6,*) 'Amat2', Amat2(n,1)
c            write(6,*) 'Cmat2', Cmat2(n,1)
c            write(6,*) 'Cmat11', Cmat11(n,1)
c            write(6,*) 'Amat3', Amat3(n,1)
c            write(6,*) 'Cmat3', Cmat3(n,1)
c            write(6,*) 'Cmat21', Cmat21(n,1)
c            write(6,*) 'Cmat12', Cmat12(n,1)
c            write(6,*) '=================='
c
      endif
c
      if (ipol.eq.0) then
c
c        ======> SPIN-UNRESTRICTED <======
c
c         do 20 n = 1, nq
c            if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c
c            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = C*(2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = F43*rho43*rrho
c Daniel (7-29-12): Needed for XC-third derivatives also
            rhom23 = F13*rho13*rrho
c Daniel (7-27-12): Derivative of the electron density
            rhom53 = F23*rhom23*rrho
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c               Amat(n,1) = Amat(n,1) + rho13*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*rhom23*fac
c Daniel (7-27-12): Third-derivative of the XC-functional (is this
c correct??  It seems like the factor of 4 here and 2 in the second
c derivative shouldn't be there for the unrestricted part).
c               Amat3(n,D2_RA_RA_RA) = Amat3(n,D2_RA_RA_RA) 
c     1                              - 4.0d0*rhom53*fac
            endif
c
c            gamma = delrho(n,1,1)*delrho(n,1,1) +
c     &              delrho(n,2,1)*delrho(n,2,1) +
c     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 25
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c

            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,1) = Cmat(n,1) + 0.5d0*rho43*d1g(2)*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
c Daniel (7-27-12): Factor of 2 is used to combat the factor of
c 0.5 in rrho
            Amat2(n,1) = Amat2(n,1)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat11(n,1) = Cmat11(n,1)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,1) = Cmat2(n,1)
     &           + rho43*d2g(3)*fac*0.5d0
c Daniel (7-27-12): XC-third derivatives of different PBEs.
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
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
c Currently guessing the coefficients based on what I did for the
c Slater/Dirac functional.  The coefficients seem related to the
c factor of 0.5 in rrho.  We need a factor of 0.5 left over.
            Amat3(n,1) = Amat3(n,1)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac*4.0d0
            Cmat21(n,1) = Cmat21(n,1)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac*2.0d0
            Cmat12(n,1) = Cmat12(n,1)
     1           + (rho13*d2g(3)
     2           +  rho43*d3g(3))*fac
            Cmat3(n,1) = Cmat3(n,1)
     1           + (rho43*d3g(4))*fac*0.5d0
c
c     Beta
c
c 25         continue
c            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
c Daniel (7-29-12): Needed for XC-third derivatives also
            rhom23 = F13*rho13*rrho
c Daniel (7-27-12): Derivative of the electron density
            rhom53 = F23*rhom23*rrho
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c               Amat(n,2) = Amat(n,2) + rho13*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
c               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*rhom23*fac
c Daniel (7-27-12): Third-derivative of the XC-functional (is this
c correct??  It seems like the factor of 4 here and 2 in the second
c derivative shouldn't be there for the unrestricted part).
c               Amat3(n,D2_RA_RA_RA) = Amat3(n,D2_RA_RA_RA)
c     1                              - 4.0d0*rhom53*fac
            endif
c
c            gamma = delrho(n,1,2)*delrho(n,1,2) +
c     &              delrho(n,2,2)*delrho(n,2,2) +
c     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 20
c
            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            if(whichf.eq.'revp') then
               g=grevpbe0(s)
               gp=grevpbe1(s)
            elseif(whichf.eq.'rpbe') then
               g=grpbe0(s)
               gp=grpbe1(s)
            else
               g=gpbe0(s)
               gp=gpbe1(s)
            endif
c
            d1g(1) = gp*d1s(1)
            d1g(2) = gp*d1s(2)
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,2) = Amat(n,2) + (rho13*g+rho43*d1g(1))*fac
            Cmat(n,3) = Cmat(n,3) + 0.5d0*rho43*d1g(2)*fac
c Daniel (7-29-12): Needed for XC-third derivatives also
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
            if(whichf.eq.'revp') then
               gpp=grevpbe2(s)
            elseif(whichf.eq.'rpbe') then
               gpp=grpbe2(s)
            else
               gpp=gpbe2(s)
            endif
            d2g(1) = gp*d2s(1) + gpp*d1s(1)*d1s(1)
            d2g(2) = gp*d2s(2) + gpp*d1s(1)*d1s(2)
            d2g(3) = gp*d2s(3) + gpp*d1s(2)*d1s(2)
            Amat2(n,3) = Amat2(n,3)
     &           +(rhom23*g
     &           + 2.d0*rho13*d1g(1)
     &           + rho43*d2g(1))*fac*2d0
            Cmat11(n,6) = Cmat11(n,6)
     &           +(rho13*d1g(2)
     &           + rho43*d2g(2))*fac
            Cmat2(n,6) = Cmat2(n,6)
     &           + rho43*d2g(3)*fac*0.5d0
c Daniel (7-27-12): XC-third derivatives of different PBEs.
c 1 = drdrdr, 2 = drdrdg, 3 = drdgdg, 4 = dgdgdg
            d3s(1) = -F10d3*d2s(1)*rrho
            d3s(2) = 0.5d0*d2s(1)/gamma
            d3s(3) = -F43*d2s(3)*rrho
            d3s(4) = -1.5d0*d2s(3)/gamma
            if(whichf.eq.'revp') then
               gppp = grevpbe3(s)
            elseif(whichf.eq.'rpbe') then
               gppp = grpbe3(s)
            else
               gppp = gpbe3(s)
            endif
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
c Currently guessing the coefficients based on what I did for the
c Slater/Dirac functional.  The coefficients seem related to the
c factor of 0.5 in rrho.  We need a factor of 0.5 left over.
            Amat3(n,4) = Amat3(n,4)
     1           + (-rhom53*g
     2           +  3.0d0*rhom23*d1g(1)
     3           +  3.0d0*rho13*d2g(1)
     4           +  rho43*d3g(1))*fac*4.0d0
            Cmat21(n,9) = Cmat21(n,9)
     1           + (rhom23*d1g(2)
     2           +  2.0d0*rho13*d2g(2)
     3           +  rho43*d3g(2))*fac*2.0d0
            Cmat12(n,12) = Cmat12(n,12)
     1           + (rho13*d2g(3)
     2           +  rho43*d3g(3))*fac
            Cmat3(n,9) = Cmat3(n,9)
     1           + (rho43*d3g(4))*fac*0.5d0
c
c 20      continue
      endif
c
c      return
      end
