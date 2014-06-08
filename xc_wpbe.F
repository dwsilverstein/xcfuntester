c     PBE exchange functional 
c
c     References:
c     [a] J.P. Perdew, K. Burke, and M. Ernzerhof, PRL 77, 3865 (1996).
c     [b] J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).;
c                                               40, 3399 (1989) (E).
c     Hammer, Hansen and Norskov, PRB 59, 7413 (1999) [RPBE]
c     Zhang and Yang, PRL 80, 890 (1998) [RevPBE]
c
      Subroutine xc_wpbe(cam_omega,n,
     &     tol_rho, fac, lfac, nlfac, rho, gamma, 
     &     Amat, Amat2, Cmat, Cmat11, Cmat2, Amat3, Cmat21, Cmat12,
     &     Cmat3, nq, ipol, Ex,
     &     qwght,ldew,func)
c
c$Id: xc_wpbe96.F,v 1.10 2006/01/11 02:05:01 edo Exp $
c
      implicit none
c
c      
c      character*4 whichf
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
      double precision delrho(1,3,ipol)
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
      double precision F43, F13
      double precision F73
      double precision F23
      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter(ukrev=1.245d0, umkrev=um/ukrev)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F73=7.d0/3.d0)
      parameter (F23=2.d0/3.d0)
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, d1s(2),
     &      g, gp1, d1gs(2), gp2, d1gnu
      double precision nu, kf, d1nu
      double precision Fxhse,Fxhse10,Fxhse01
      double precision rhom23, gpp1, d2s(3),d2gs(3), gpp2, d2gnu
      double precision Fxhse20,Fxhse02,d2nu,gpp11,d2gsnu,Fxhse11
      double precision Fxhse30, Fxhse21, Fxhse12, Fxhse03
      double precision d3s(4),d3nu,d3gs(4),rhom53
      double precision d1nud1rhonu,d1nud1rhos,d1sd1rhos,d1sd1rhonu
c
      double precision cam_omega 
c
      Ex     = 0d0
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cs = Cs * C               ! account for including C in rho43
c
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
            kf = -1d0*pi*rho13 
            rhom23 = F13*rho13*rrho
            rhom53 = -F23*rhom23*rrho
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac
c               Amat(n,1) = Amat(n,1) + rho13*fac
c               Amat2(n,1) = Amat2(n,1) + rhom23*fac
c               Amat2(n,1) = Amat2(n,1) + 2d0*rhom23*fac
            endif
c
c            gamma = delrho(1,1,1)*delrho(1,1,1) +
c     &              delrho(1,2,1)*delrho(1,2,1) +
c     &              delrho(1,3,1)*delrho(1,3,1)
            gam12 = dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho**2)) goto 10A
c

            s = Cs*gam12/rho43
            d1s(1) = -F43*s*rrho
            d1s(2) = 0.5d0*s/gamma
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
c     Evaluate the GC part of F(s), i.e. g(s) = F(s) - 1
c
            call HSE08Fx(rho(n,1),s,nu,fxhse,fxhse10,fxhse01,
     &                   fxhse20,fxhse02,fxhse11,fxhse30,
     &                   fxhse21,fxhse12,fxhse03)

            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*(d1gs(1)+d1gnu))*fac
c            Cmat(n,1) = Cmat(n,1) + 2d0*rho43*d1gs(2)*fac
            Cmat(n,1) = Cmat(n,1) + rho43*d1gs(2)*fac
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*rrho*d1nu
c
            d1sd1rhos = -F43*rrho
            d1sd1rhonu = -F13*2d0*cam_omega/gam12
            d1nud1rhonu = -F13*rrho
            d1nud1rhos = -F43*0.5d0*gam12*rrho**2d0/cam_omega
c
            gpp1=fxhse20
            gpp2=fxhse02
            gpp11=fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c
            Amat2(n,1) = Amat2(n,1)
     &           + (rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac
            Cmat11(n,1) = Cmat11(n,1)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac
            Cmat2(n,1) = Cmat2(n,1)
     &           + rho43*d2gs(3)*fac
c
c d3s(1)=d3s_rho, d3s(2)=d3s_rho2_gamma, d3s(3)=d3s_rho_gamma2
c d3s(4)=d3s_gamma, d3nu=d3nu_rho

      d3s(1) = -280d0*s/(27d0*rho(n,1)**3)
      d3s(2) = 14d0*s/(9d0*rho(n,1)**2*gamma)
      d3s(3) = s/(3d0*rho(n,1)*gamma**2)
      d3s(4) = 3d0*s/(8d0*gamma**3)
      d3nu = -28d0*nu/(27d0*rho(n,1)**3)

      d3gs(1) = d3nu*Fxhse01 + 2d0*d2nu*(d1nu*Fxhse02 + 
     &  d1s(1)*Fxhse11) + 2d0*(d1nu*Fxhse11 + d1s(1)*
     &  Fxhse20)*d2s(1) + d1nu*(d2nu*Fxhse02 + d1nu*(d1nu*
     &  Fxhse03 + d1s(1)*Fxhse12) + Fxhse11*d2s(1) + d1s(1)*
     &  (d1nu*Fxhse12 + d1s(1)*Fxhse21)) + d1s(1)*(d2nu*
     &  Fxhse11 + Fxhse20*d2s(1) + d1nu*(d1nu*Fxhse12 + 
     &  d1s(1)*Fxhse21) + d1s(1)*(d1nu*Fxhse21 + d1s(1)*
     &  Fxhse30)) + Fxhse10*d3s(1)

      d3gs(2) =d2nu*d1s(2)*Fxhse11 + d2s(2)*(d1nu*Fxhse11 + 
     &  d1s(1)*Fxhse20) + d1s(2)*Fxhse20*d2s(1) + d1nu*
     &  (Fxhse11*d2s(2) + d1nu*d1s(2)*Fxhse12 + d1s(2)*d1s(1)*
     &  Fxhse21) + Fxhse10*d3s(2) + d1s(1)*(d2s(2)*Fxhse20 + 
     &  d1nu*d1s(2)*Fxhse21 + d1s(2)*d1s(1)*Fxhse30)

      d3gs(3) = Fxhse10*d3s(3) + 2d0*d1s(2)*d2s(2)*Fxhse20 + 
     &  d1nu*(d2s(3)*Fxhse11 + d1s(2)**2d0*Fxhse21) + 
     &  d1s(1)*(d2s(3)*Fxhse20 + d1s(2)**2d0*Fxhse30)

      d3gs(4) = Fxhse30*d1s(2)**3d0 + 3d0*Fxhse20*d1s(2)*d2s(3)
     &        + Fxhse10*d3s(4)

c
      Amat3(n,1) = Amat3(n,1) + ((rhom53*g
     &      + 3d0*(rhom23*(d1gs(1)+d1gnu))
     &      + 3d0*(rho13*(d2gs(1)+d2gnu+d2gsnu))
     &      + rho43*d3gs(1))*fac)
c      Amat3 = Amat3 + K*((-(8d0/27d0)*rhom53*g
c     &      + 3d0*((4d0/9d0)*rhom23*(d1gs(1)+d1gnu))
c     &      + 3d0*((4d0/3d0)*rho13*(d2gs(1)+d2gnu+d2gsnu)) 
c     &      + rho43*d3gs(1))*fac)

      Cmat21(n,1) = Cmat21(n,1) + ((rhom23*d1gs(2)
     &      + 2d0*rho13*(d2gs(2))
     &      + rho43*d3gs(2))*fac)
c      Cmat21 = Cmat21+ K*(((4d0/9d0)*rhom23*d1gs(2)
c     &      + 2d0*(4d0/3d0)*rho13*(d2gs(2))
c     &      + rho43*d3gs(2))*fac)

      Cmat12(n,1) = Cmat12(n,1) + ((rho43*d3gs(3)
     &      + rho13*d2gs(3))*fac)
c      Cmat12 = Cmat12 + K*((rho43*d3gs(3)
c     &      + (4d0/3d0)*rho13*d2gs(3))*fac)

      Cmat3(n,1) = Cmat3(n,1) + ((rho43*d3gs(4))*fac)
c      Cmat3 = Cmat3 + K*((rho43*d3gs(4))*fac)

c 10      continue
c
      endif
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
            kf = -1d0*pi*rho13
            rhom23 = F13*rho13*rrho
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c               Amat = Amat + rho13*fac
c               Amat2 = Amat2 + 2d0*rhom23*fac
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
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c            call HSE08Fx(rho(n,2),s,nu,fxhse,fxhse10,fxhse01,
c     &                     fxhse20,fxhse02,fxhse11)

            g = fxhse - 1.0d0
            gp1 = fxhse10
            gp2 = fxhse01 
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*(d1gs(1)+d1gnu))*fac
            Cmat(n,1) = Cmat(n,1) + 0.5d0*rho43*d1gs(2)*fac
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c     
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2) +
     &                gpp11*d1nu*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c            Amat2(n,1) = Amat2(n,1)
c     &           +(rhom23*g
c     &           + 2.d0*rho13*(d1gs(1)+d1nu)
c     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac*2d0
c            Cmat11(n,1) = Cmat11(n,1)
c     &           +(rho13*d1gs(2)
c     &           + rho43*d2gs(2))*fac
c            Cmat2(n,1) = Cmat2(n,1)
c     &           + rho43*d2gs(3)*fac*0.5d0
            Amat2(n,1) = Amat2(n,1)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac
            Cmat11(n,1) = Cmat11(n,1)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac
            Cmat2(n,1) = Cmat2(n,1)
     &           + rho43*d2gs(3)*fac*0.5d0

      endif
      if (ipol.eq.0)then

c
c     Beta
c
c 25         continue
c            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho
            kf = -1d0*pi*rho13
            rhom23 = F13*rho13*rrho
            if (lfac) then
c               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
c               Amat(n,1) = Amat(n,1) + rho13*fac
c               Amat2(n,1) = Amat2(n,1) + 2d0*rhom23*fac
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
c           Note that there is an extra factor of 0.5 in rrho, so a
c           factor of 2 needs to be multiplied in to have the correct
c           Fermi wavevector.
c
            nu = cam_omega/kf
            d1nu = -F13*nu*rrho
c
c            call HSE08Fx(rho(n,3),s,nu,fxhse,fxhse10,fxhse01,
c     &                      fxhse20,fxhse02,fxhse11)

            g = fxhse - 1.0d0
            gp1= fxhse10
            gp2 = fxhse01
c
            d1gs(1) = gp1*d1s(1)
            d1gs(2) = gp1*d1s(2)
            d1gnu = gp2*d1nu
            Ex = Ex + rho43*g*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*g*fac*0.5d0
            Amat(n,1) = Amat(n,1) + (rho13*g+rho43*(d1gs(1)+d1nu))*fac
            Cmat(n,1) = Cmat(n,1) + 0.5d0*rho43*d1gs(2)*fac
            d2s(1) = -F73*d1s(1)*rrho
            d2s(2) = -F43*d1s(2)*rrho
            d2s(3) = -0.5d0*d1s(2)/gamma
c
            d2nu = -F43*d1nu*rrho
c
            gpp1 = fxhse20
            gpp2 = fxhse02
            gpp11 = fxhse11
c
            d2gs(1) = gp1*d2s(1) + gpp1*d1s(1)*d1s(1)
            d2gs(2) = gp1*d2s(2) + gpp1*d1s(1)*d1s(2)
            d2gs(3) = gp1*d2s(3) + gpp1*d1s(2)*d1s(2)
            d2gnu = gp2*d2nu + gpp2*d1nu*d1nu
            d2gsnu = 2d0*gpp11*d1nu*d1s(1)
c            Amat2(n,1) = Amat2(n,1)
c     &           +(rhom23*g
c     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
c     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac*2d0
c            Cmat11(n,1) = Cmat11(n,1)
c     &           +(rho13*d1gs(2)
c     &           + rho43*d2gs(2))*fac
c            Cmat2(n,1) = Cmat2(n,1)
c     &           + rho43*d2gs(3)*fac*0.5d0
            Amat2(n,3) = Amat2(n,3)
     &           +(rhom23*g
     &           + 2.d0*rho13*(d1gs(1)+d1gnu)
     &           + rho43*(d2gs(1)+d2gnu+d2gsnu))*fac
            Cmat11(n,6) = Cmat11(n,6)
     &           +(rho13*d1gs(2)
     &           + rho43*d2gs(2))*fac
            Cmat2(n,6) = Cmat2(n,6)
     &           + rho43*d2gs(3)*fac*0.5d0

c
c 20      continue
      endif
c
c      return
      end
