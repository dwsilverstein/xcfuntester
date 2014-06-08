      Subroutine xc_lyp88(n, tol_rho, fac,  rho, gaa, 
     &                    Amat, Amat2, Amat3, Cmat, 
     &                    Cmat2, Cmat11, Cmat3, Cmat21, Cmat12,
     &                    nq, ipol, Ec, qwght, ldew, func)
c
C$Id: xc_lyp88.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      implicit none 
c
      double precision fac ! [input]
      integer nq
      integer ipol
      double precision Ec
      logical ldew
      double precision func(nq)  ! value of the functional [output]
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
      double precision F13, F43, F113, F83, F53, F19, F79, P1, tol_rho,
     &                 A, B, C, D
c
c      Parameter (F13 = 1.0D0/3.0D0, F43 = 4.0D0*F13, 
c     &           F113 = 11.0D0*F13, F83 = 8.0D0*F13, 
c     &           F53 = 5.0D0*F13, F19 = 1.0D0/9.0D0,
c     &           F79 = 7.0D0*F19)
      Parameter (F13 = 1.D0/3.D0, F43 = 4.D0*F13, F113 = 11.D0*F13,
     &           F83 = 8.D0*F13, F53 = 5.D0*F13, F19 = 1.D0/9.D0,
     &           F79 = 7.D0*F19)
c Third derivatives
      double precision F23, F73, F223
      Parameter (F23 = 2.0D0*F13, F73 = 7.0D0*F13, F223 = 22.0*F13)
c
c     P1 = 2**(11/3)*(3/10)*(3*PI**2)**(2/3)
c
      Parameter (P1 = 0.3646239897876487D+02)
c
c     Colle-Salvetti Empirical Parameters
c
      Parameter (A = 0.04918D0)
      Parameter (B = 0.13200D0)
      Parameter (C = 0.25330D0)
      Parameter (D = 0.34900D0)
c
c     Compute the partial derivatives of the correlation functional of 
c     Lee, Yang and Parr.
c
c     References:
c
c     Colle & Salvetti, Theor. Chim. Acta 37, 329 (1975)
c     Lee, Yang & Parr, Phys. Rev. B 37, 785 (1988)
c     Miehlich, Savin, Stoll & Preuss, Chem. Phys. Lett. 157, 200 (1989)
c     Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
      integer n
      double precision c1, c2, ab, ratcd
      double precision rrho, rhoa, rhob, rrho2, rhoa2, rhob2, rhoab,
     &     rhoa53, rhob53, rhoa83, rhob83, rho2,
     &     h1, h2, h3, d1h1, d1h2, d1h3, om, d1om, de, d1de, de11, de47,
     &     t, d1t(2)
c Suggested variable to add
      double precision xrarb, d1xrarb(2)
      double precision tm_in, d1tm_in(2)
c
      double precision gaa, gab, gbb
      double precision f1, f2, d1f1(2), d1f2(2), f, d1f(5),
     &     d2fgaa(2), d2fgab(2), d2fgbb(2)
c Second derivatives
      double precision d2h1, d2h2, d2h3, d2om, d2de, rrho3, d2f1(3),
     &     d2f2(3), rhoa113, rhob113, d3fgaa(3), d3fgab(3), d3fgbb(3),
     &     d2t(3), d2f(3)
      double precision dt
      double precision d2xrarb(3)
      double precision d2tm_in(3)
c Third derivatives (the 4's in the array dimensions correspond to
c permutations of spin-densities)
      double precision rrho4, rhoa23, d3h1, d3h2, d3h3, d3om, d3de, 
     1     d3f1(4), d3f2(4), d4fgaa(4), d4fgab(4), d4fgbb(4), d3t(4), 
     2     d3f(4)
      double precision ddt
      double precision d3xrarb(4)
      double precision d3tm_in(4)
c
c     Coefficients of first two terms in LYP functional and other
c     commonly occurring factors
c
      c1 = -4.0d0*a
c      c1 = -4d0*a
      c2 = -P1*a*b
      ab = a*b
      ratcd = c/d
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
c         do 10 n = 1, nq
c            if (rho(n,1).lt.tol_rho)goto 10
            rrho = 1.0d0/rho(n,1)
            rhoa = 0.5d0*rho(n,1)
            rrho2 = rrho*rrho
            rho2 = 1.0d0/rrho2
            rhoa2 = rhoa*rhoa
            rhoab = rhoa2
            rhoa53 = abs(rhoa)**F53*sign(1.0d0,rhoa)
            rhoa83 = abs(rhoa)**F83*sign(1.0d0,rhoa)
c Second derivatives
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
c Third derivatives
            rrho4 = rrho*rrho3
            rhoa23 = abs(rhoa)**F23*sign(1.0d0,rhoa)
c
            h2 = d*abs(rrho)**F13
            d1h2 = -F13*h2*rrho
c Second derivatives
            d2h2 = -F43*d1h2*rrho
c Third derivatives
            d3h2 = -F73*d2h2*rrho 
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
c Second derivatives
            d2h3 = ratcd*d2h2
c Third derivatives
            d3h3 = ratcd*d3h2
c
            h1 = 1.0d0/(1.0d0+h2)
            d1h1 = -h1*h1*d1h2
c Second derivatives
            d2h1 = -(2.0d0*h1*d1h1*d1h2 + h1*h1*d2h2)
c Third derivatives
            d3h1 = -6.0d0*d1h1*d1h1*d1h2
     1           - 6.0d0*h1*d2h2*d1h1 
     2           - h1*h1*d3h2
c
!            om = exp(-h3)*h1*rrho**F113
            om = exp(-h3)*h1*abs(rrho)**F113
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
c Second derivative
            dt = d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2
c            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
            d2om = -(d1om*t+om*dt)
c Third derivatives
            ddt = d3h3 + d2h1*d1h2 + 2.0d0*d1h1*d2h2
     1          + h1*d3h2 + F223*rrho3
            d3om = -(ddt*om + 2.0d0*d1om*dt + d2om*t) 
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
c Second derivative
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2.0d0*d1h1*d1h2
c Third derivative
            d3de = d3h3 + d3h1*h2 + 3.0d0*d2h1*d1h2 
     1           + 3.0d0*d1h1*d2h2 + h1*d3h2
c
c Daniel (10-24-12): f1 is what I call J in my notes (sans constants).
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(1) = d1f1(1) + h1*rhoa*rrho
            d1f1(2) = d1f1(1)
c Second derivatives (1 and 2 refer to rara and rarb, respectively)
c I'd suggest writing the expressions in a simpler form than what
c was originally present.  The forms given should be easier to 
c evaluate.
c            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
c     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
            d2f1(1) = d2h1*rhoab*rrho 
     1              + d1h1*(1.0d0 - rhoa*rrho)
     2              + h1*(-rrho + rhoa*rrho2)
c            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
c     &           + 2d0*h1*rhoab*rrho3
            d2f1(2) = d2h1*rhoab*rrho 
     1              + d1h1*(1.0d0-rhoa*rrho)
     2              + h1*rhoa*rrho2
            d2f1(3) = d2f1(1)
c Third derivatives (keep in mind rhoa = rhob)
c 1 = rarara, 2 = rararb, 3 = rarbrb
            d3f1(1) = d3h1*rhoab*rrho
     1              + d2h1*( 1.0d0 - rhoab*rrho2)
     2              + d1h1*(-1.5d0*rrho)
     3              + h1*(1.5d0*rrho2)
            d3f1(2) = d3h1*rhoab*rrho
     1              + d2h1*(1.0d0 - rhoab*rrho2)
     2              + d1h1*(0.50d0*rrho)
     3              + h1*(-0.50d0*rrho2)
            d3f1(3) = d3f1(2)
            d3f1(4) = d3f1(1)
c
c Daniel (10-25-12): Define xrarb here
            xrarb = rhoab*(rhoa83 + rhoa83)
c
c            f2 = om*rhoab*(rhoa83+rhoa83)
            f2 = om*xrarb
c
c Daniel (10-25-12): Define d1xrarb here
            d1xrarb(1) = rhoa*(F113*rhoa83+rhoa83)
c
c            d1f2(1) = d1om*rhoab*(rhoa83+rhoa83)
c            d1f2(1) = d1f2(1) + om*rhoa*(F113*rhoa83+rhoa83)
            d1f2(1) = d1om*xrarb + om*d1xrarb(1)
            d1f2(2) = d1f2(1)
c Second derivatives
c Daniel (10-25-12): Define d2xrarb here
            d2xrarb(1) = rhoa*F113*F83*rhoa53
            d2xrarb(2) = F113*(rhoa83+rhoa83)
c
c            d2f2(1) = d2om*rhoab*(rhoa83+rhoa83)
c     &          + 2d0*d1om*rhoa*(F113*rhoa83+rhoa83)
c     &          +       om*rhoa*F113*F83*rhoa53
c            d2f2(2) = d2om*rhoab*(rhoa83+rhoa83)
c     &           + d1om*(rhoa113+rhoa113+F113*(rhoa*rhoa83+rhoa*rhoa83))
c     &           +   om*F113*(rhoa83+rhoa83)
            d2f2(1) = d2om*xrarb
     &          + 2.0d0*d1om*d1xrarb(1)
     &          +       om*d2xrarb(1)
            d2f2(2) = d2om*xrarb
     &           + 2.0d0*d1om*d1xrarb(1)
     &           +   om*d2xrarb(2)
            d2f2(3) = d2f2(1)
c Third derivatives
c Daniel (10-25-12): Define d2xrarb here
            d3xrarb(1) = rhoa*F113*F83*F53*rhoa23
            d3xrarb(2) = F113*F83*rhoa53
            d3xrarb(3) = F113*F83*rhoa53
c
            d3f2(1) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(1)
     2              + 3.0d0*d1om*d2xrarb(1)
     3              + om*d3xrarb(1)
            d3f2(2) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(1)
     2              + d1om*(d2xrarb(1) + 2.0d0*d2xrarb(2))
     3              + om*d3xrarb(2)
c            d3f2(3) = d3om*xrarb
c     1              + 3.0d0*d2om*d1xrarb(1)
c     2              + d1om*(d2xrarb(1) + 2.0d0*d2xrarb(2))
c     3              + om*d3xrarb(3)
            d3f2(3) = d3f2(2)
            d3f2(4) = d3f2(1)
c
c This is read in for xcfuntester, so I comment it out.  They use the
c total density for a closed shell calculation, so gaa is multiplied by
c one-half twice.
c            gaa =(delrho(n,1,1)*delrho(n,1,1) +
c     &            delrho(n,2,1)*delrho(n,2,1) +
c     &            delrho(n,3,1)*delrho(n,3,1))*0.25d0
c
            de11 = de - 11.0d0
            de47 = 47.0d0 - 7.0d0*de
c            de11 = de - 11d0
c            de47 = 47d0 - 7d0*de
c
c Daniel (10-23-12): "t" is what I call Q or S.  I'm adding what I call
c Qi (tm_in) which allows me to write the expressions in a simpler way.
            tm_in = 1.0d0-3.0d0*de-de11*rhoa*rrho
            t = F19*rhoab*tm_in-rhoa2
c            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhoa2
c Daniel (10-23-12): d1t(1) is the derivative with respect to rhoa,
c and d1t(2) is the derivative with respect to rhob.
c
c            d1tm_in(1) = (3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2
            d1tm_in(1) = -(3.0d0+rhoa*rrho)*d1de-de11*rhoa*rrho2
            d1tm_in(2) = -(3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2 
c These are numerically unstable, although they are equivalent to
c the above two expressions.
c            d1tm_in(1) = rhoa*de11*rrho2 
c     1                 - de11*rrho 
c     2                 - 3.0d0*d1de - rhoa*d1de*rrho
c            d1tm_in(2) = rhoa*de11*rrho2 - 3.0d0*d1de 
c     1                 - rhoa*d1de*rrho
c
c            d1t(1) = F19*( rhoa*tm_in - rhoab*d1tm_in(1) )
            d1t(1) = F19*( rhoa*tm_in + rhoab*d1tm_in(1) )
            d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )    
     &             - 2.0d0*rhoa
c            d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )    
c            d1t(1) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
c     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
c            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
c     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
c     &             - 2d0*rhoa
c Daniel (10-23-12): "d1f(3)" represents the third term 
c (what I call L in my notes).  It is a derivative with respect to
c gamma, not rho.
            d1f(3) = -ab*om*t
            d2fgaa(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgaa(2) = -ab*( d1om*t + om*d1t(2) )
c Second derivatives
c
c            d2tm_in(1) = (3.0d0+rhoa*rrho)*d2de
c     1                 + 2.0d0*d1de*rhoa*rrho2
c     2                 - 2.0d0*de11*rhoa*rrho3
            d2tm_in(1) = -(3.0d0+rhoa*rrho)*d2de
     1                 - 2.0d0*d1de*rhoa*rrho2
     2                 + 2.0d0*de11*rhoa*rrho3
c            d2tm_in(2) = (3.0d0+rhoa*rrho) ! Written without d2de
            d2tm_in(2) = -(3.0d0+rhoa*rrho) ! Written without d2de
c            d2tm_in(3) = (3.0d0+rhoa*rrho)*d2de
c     1                 - 2.0d0*d1de*rhoa*rrho2
c     2                 + 2.0d0*de11*rhoa*rrho3
            d2tm_in(3) = -(3.0d0+rhoa*rrho)*d2de
     1                 + 2.0d0*d1de*rhoa*rrho2
     2                 - 2.0d0*de11*rhoa*rrho3
c
c            d2t(1) = -F19*( 2.0d0*rhoa*d1tm_in(1)
c     1                    + rhoab*d2tm_in(1) )
            d2t(1) = F19*( 2.0d0*rhoa*d1tm_in(1)
     1                   + rhoab*d2tm_in(1) )
c            d2t(2) = F19*( tm_in
c     1                   - rho(n,1)*d2tm_in(2)*d1de
c     2                   - rhoab*d2tm_in(2)*d2de )
            d2t(2) = F19*( tm_in
     1                   + rho(n,1)*d2tm_in(2)*d1de
     2                   + rhoab*d2tm_in(2)*d2de )
c            d2t(3) = -F19*( -2.0d0*rhoa*d1tm_in(2) 
c     1                    + rhoab*d2tm_in(3) )
c     2             - 2.0d0
            d2t(3) = F19*( 2.0d0*rhoa*d1tm_in(2) 
     1                   + rhoab*d2tm_in(3) )
     2             - 2.0d0
c            d2t(1) = -F19*(
c     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhoa*rrho2
c     &                                         -2d0*de11*rhoa*rrho3))
c            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
c     &           - rho(n,1)*(3d0+rhoa*rrho)*d1de
c     &           - rhoab*((3d0+rhoa*rrho)*d2de))
c            d2t(3) = -F19*(
c     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
c     &                                         +2d0*de11*rhoa*rrho3))
c     &           - 2d0
c
            d3fgaa(1) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(1) )
            d3fgaa(2) = -ab*( d2om*t + d1om*( d1t(1) + d1t(2) ) 
     1                      + om*d2t(2) )
            d3fgaa(3) = -ab*( d2om*t + 2.0d0*d1om*d1t(2) + om*d2t(3) )
c Third derivatives 
c 
            d3tm_in(1) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - 3.0d0*d2de*rhoa*rrho2
     2                 + 6.0d0*d1de*rhoa*rrho3
     3                 - 6.0d0*de11*rhoa*rrho4
c            d3tm_in(1) = ( 3.0d0 + rhoa*rrho )*d3de
c     1                 + 3.0d0*d2de*rhoa*rrho2
c     2                 - 6.0d0*d1de*rhoa*rrho3
c     3                 + 6.0d0*de11*rhoa*rrho4
            d3tm_in(2) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - d2de*rhoa*rrho2
     2                 + 2.0d0*d1de*rhoa*rrho3
     3                 - 2.0d0*de11*rhoa*rrho4
c            d3tm_in(2) = ( 3.0d0 + rhoa*rrho )*d3de
c     1                 + d2de*rhoa*rrho2
c     2                 - 2.0d0*d1de*rhoa*rrho3
c     3                 + 2.0d0*de11*rhoa*rrho4
            d3tm_in(3) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + d2de*rhoa*rrho2
     2                 - 2.0d0*d1de*rhoa*rrho3
     3                 + 2.0d0*de11*rhoa*rrho4
c            d3tm_in(3) = ( 3.0d0 + rhoa*rrho )*d3de
c     1                 - d2de*rhoa*rrho2
c     2                 + 2.0d0*d1de*rhoa*rrho3
c     3                 - 2.0d0*de11*rhoa*rrho4
            d3tm_in(4) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + 3.0d0*d2de*rhoa*rrho2
     2                 - 6.0d0*d1de*rhoa*rrho3
     3                 + 6.0d0*de11*rhoa*rrho4
c            d3tm_in(4) = ( 3.0d0 + rhoa*rrho )*d3de
c     1                 - 3.0d0*d2de*rhoa*rrho2
c     2                 + 6.0d0*d1de*rhoa*rrho3
c     3                 - 6.0d0*de11*rhoa*rrho4
c
            d3t(1) = F19*( 3.0d0*rhoa*d2tm_in(1)
     1                   + rhoab*d3tm_in(1) )
c
            d3t(2) = F19*( 2.0d0*d1tm_in(1)
     1                   + rhoa*d2tm_in(1)
     2                   + 2.0d0*rhoa*d2tm_in(2)*d2de
     3                   + rhoab*d3tm_in(2) )
c
            d3t(3) = F19*( 2.0d0*d1tm_in(2)
     2                   + 2.0d0*rhoa*d2tm_in(2)*d2de
     3                   + rhoa*d2tm_in(3)
     4                   + rhoab*d3tm_in(3) )
c
            d3t(4) = F19*( 3.0d0*rhoa*d2tm_in(3)
     1                   + rhoab*d3tm_in(4) )
c Here we collect the derivatives of L
            d4fgaa(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) ) 
            d4fgaa(2) = -ab*( d3om*t + d2om*(2.0d0*d1t(1) + d1t(2))
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2)) + om*d3t(2) )
            d4fgaa(3) = -ab*( d3om*t + d2om*(d1t(1) + 2.0d0*d1t(2))
     1                      + d1om*(2.0d0*d2t(2) + d2t(3)) + om*d3t(3) )
            d4fgaa(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(2)
     1                      + 3.0d0*d1om*d2t(3) + om*d3t(4) ) 
c
c Daniel (10-23-12): "t" is what I call R.
            t = F19*rhoab*de47-F43*rho2
            d1t(1) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,1)
            d1t(2) = d1t(1)
c Daniel (10-25-12): "d1f(4)" represents the fourth term (what I call M
c in my notes).  It is a derivative with respect to gamma, not rho.
            d1f(4) = -ab*om*t
            d2fgab(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgab(2) = -ab*( d1om*t + om*d1t(2) )
c Second derivatives
            d2t(1) = -F79*( 2.0d0*rhoa*d1de + rhoab*d2de ) - F83
            d2t(2) = F19*de47 - F79*( rho(n,1)*d1de + rhoab*d2de ) 
     1             - F83
            d2t(3) = d2t(1)
c            d2t(1) = -F79*(2d0*rhoa*d1de+rhoab*d2de) - F83
c            d2t(2) = F19*de47 - F79*(rho(n,1)*d1de+rhoab*d2de) - F83
            d3fgab(1) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(1) )
            d3fgab(2) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(2) )
            d3fgab(3) = -ab*( d2om*t + 2.0d0*d1om*d1t(2) + om*d2t(3) )
c            d3fgab(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
c            d3fgab(2) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(2))
c Third derivatives (1 = rarara, 2 = rararb, 3 = rarbrb)
            d3t(1) = -F79*( 3.0d0*rhoa*d2de + rhoab*d3de)
            d3t(2) = -F79*( 2.0d0*d1de + 3.0d0*rhoa*d2de 
     1                    + rhoab*d3de)
            d3t(3) = d3t(2)
            d3t(4) = d3t(1)
            d4fgab(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1))
            d4fgab(2) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(2))
            d4fgab(3) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(3))
            d4fgab(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(4))
c
c Daniel (10-26-12): This is where the variable N is evaluated 
            d2fgbb(1) = d2fgaa(2)
            d2fgbb(2) = d2fgaa(1)
c Second derivatives
            d3fgbb(1) = d3fgaa(3)
            d3fgbb(2) = d3fgaa(2)
            d3fgbb(3) = d3fgaa(1)
c Third derivatives
            d4fgbb(1) = d4fgaa(4)
            d4fgbb(2) = d4fgaa(3)
            d4fgbb(3) = d4fgaa(2)
            d4fgbb(4) = d4fgaa(1)
c
c Daniel (10-25-12): It is tricky that this is written with derivatives
c with respect to gamma.
            f = c1*f1 + c2*f2 + gaa*(2d0*d1f(3) + d1f(4))
c First derivatives
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*(d2fgaa(1) + d2fgab(1) + d2fgbb(1))
            d1f(2) = c1*d1f1(2) + c2*d1f2(2)
     &             + gaa*(d2fgaa(2) + d2fgab(2) + d2fgbb(2))
c Second derivatives
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*(d3fgaa(1) + d3fgab(1) + d3fgbb(1))
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*(d3fgaa(2) + d3fgab(2) + d3fgbb(2))
            d2f(3) = c1*d2f1(3) + c2*d2f2(3)
     &             + gaa*(d3fgaa(3) + d3fgab(3) + d3fgbb(3))
c Third derivatives (1 = rarara, 2 = rararb, 3 = rarbrb)
            d3f(1) = c1*d3f1(1) + c2*d3f2(1)
     1             + gaa*(d4fgaa(1) + d4fgab(1) + d4fgbb(1)) 
            d3f(2) = c1*d3f1(2) + c2*d3f2(2)
     1             + gaa*(d4fgaa(2) + d4fgab(2) + d4fgbb(2)) 
            d3f(3) = c1*d3f1(3) + c2*d3f2(3)
     1             + gaa*(d4fgaa(3) + d4fgab(3) + d4fgbb(3)) 
            d3f(4) = c1*d3f1(4) + c2*d3f2(4)
     1             + gaa*(d4fgaa(4) + d4fgab(4) + d4fgbb(4)) 
c
            Ec = Ec + f*fac*qwght(n)
            if (ldew) func(n) = func(n) + f*fac
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            Amat(n,2) = Amat(n,2) + d1f(2)*fac
c            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
c            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*fac
c            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(3)*fac
            Cmat(n,1) = Cmat(n,1) + d1f(3)*fac
            Cmat(n,2) = Cmat(n,2) + d1f(4)*fac
            Cmat(n,3) = Cmat(n,3) + d1f(3)*fac
c Second derivatives
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*fac
c            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*fac
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + d2f(2)*fac
            Amat2(n,1) = Amat2(n,1) + d2f(1)*fac
            Amat2(n,2) = Amat2(n,2) + d2f(2)*fac
            Amat2(n,3) = Amat2(n,3) + d2f(3)*fac
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*fac
c            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*fac
c            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*fac
            Cmat11(n,1) = Cmat11(n,1) + d2fgaa(1)*fac
            Cmat11(n,2) = Cmat11(n,2) + d2fgab(1)*fac
            Cmat11(n,3) = Cmat11(n,3) + d2fgbb(1)*fac
c Third derivatives
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + d3f(1)*fac
c            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + d3f(2)*fac
c            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + d3f(3)*fac
c            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + d3f(4)*fac
            Amat3(n,1) = Amat3(n,1) + d3f(1)*fac
            Amat3(n,2) = Amat3(n,2) + d3f(2)*fac
            Amat3(n,3) = Amat3(n,3) + d3f(3)*fac
            Amat3(n,4) = Amat3(n,4) + d3f(4)*fac
c            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + d3fgaa(1)*fac
c            Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB) + d3fgab(1)*fac
c            Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB) + d3fgbb(1)*fac
            Cmat21(n,1) = Cmat21(n,1) + d3fgaa(1)*fac
            Cmat21(n,2) = Cmat21(n,2) + d3fgab(1)*fac
            Cmat21(n,3) = Cmat21(n,3) + d3fgbb(1)*fac
            Cmat21(n,4) = Cmat21(n,4) + d3fgaa(2)*fac
            Cmat21(n,5) = Cmat21(n,5) + d3fgab(2)*fac
            Cmat21(n,6) = Cmat21(n,6) + d3fgbb(2)*fac
            Cmat21(n,7) = Cmat21(n,7) + d3fgaa(3)*fac
            Cmat21(n,8) = Cmat21(n,8) + d3fgab(3)*fac
            Cmat21(n,9) = Cmat21(n,9) + d3fgbb(3)*fac
c
c   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c         do 20 n = 1,nq
c            if (rho(n,1).lt.tol_rho)goto 20
            rrho = 1d0/rho(n,1)
            rhoa = max(0.0d0,rho(n,2))
            rhob = max(0.0d0,rho(n,3))
            rrho2 = rrho*rrho
            rho2 = 1d0/rrho2
            rhoa2 = rhoa*rhoa
            rhob2 = rhob*rhob
            rhoab = rhoa*rhob
            rhoa53 = rhoa**F53
            rhob53 = rhob**F53
csedo            rhoa53 = abs(rhoa)**F53
csedo            rhob53 = abs(rhob)**F53
            rhoa83 = rhoa*rhoa53
            rhob83 = rhob*rhob53
c
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
            rhob113 = rhob*rhob83
c
cedo            h2 = d*abs(rrho)**F13*sign(1d0,rrho)
            h2 = d*rrho**F13
            d1h2 = -F13*h2*rrho
c
            d2h2 = -F43*d1h2*rrho
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
c
            d2h3 = ratcd*d2h2
c
            h1 = 1d0/(1d0+h2)
            d1h1 = -h1*h1*d1h2
c
            d2h1 = -(2d0*h1*d1h1*d1h2 + h1*h1*d2h2)
c
            om = exp(-h3)*h1*rrho**F113
cedo            om = exp(-h3)*h1*abs(rrho)**F113*sign(1d0,rrho)
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
c
            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
c
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2d0*d1h1*d1h2
c
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(2) = d1f1(1)
            d1f1(1) = d1f1(1) + h1*rhob*rrho
            d1f1(2) = d1f1(2) + h1*rhoa*rrho
c
            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhob*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhob*rrho2+rhoab*rrho3)
            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
     &           + 2d0*h1*rhoab*rrho3
            d2f1(3) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
c
            f2 = om*rhoab*(rhoa83+rhob83)
            d1f2(1) = d1om*rhoab*(rhoa83+rhob83)
            d1f2(2) = d1f2(1)
            d1f2(1) = d1f2(1) + om*rhob*(F113*rhoa83+rhob83)
            d1f2(2) = d1f2(2) + om*rhoa*(F113*rhob83+rhoa83)
c
            d2f2(1) = d2om*rhoab*(rhoa83+rhob83)
     &          + 2d0*d1om*rhob*(F113*rhoa83+rhob83)
     &          +       om*rhob*F113*F83*rhoa53
            d2f2(2) = d2om*rhoab*(rhoa83+rhob83)
     &           + d1om*(rhoa113+rhob113+F113*(rhob*rhoa83+rhoa*rhob83))
     &           +   om*F113*(rhoa83+rhob83)
            d2f2(3) = d2om*rhoab*(rhoa83+rhob83)
     &          + 2d0*d1om*rhoa*(F113*rhob83+rhoa83)
     &          +       om*rhoa*F113*F83*rhob53
c
            gaa = delrho(n,1,1)*delrho(n,1,1) +
     &            delrho(n,2,1)*delrho(n,2,1) +
     &            delrho(n,3,1)*delrho(n,3,1)
            gab = delrho(n,1,1)*delrho(n,1,2) +
     &            delrho(n,2,1)*delrho(n,2,2) +
     &            delrho(n,3,1)*delrho(n,3,2)
            gbb = delrho(n,1,2)*delrho(n,1,2) +
     &            delrho(n,2,2)*delrho(n,2,2) +
     &            delrho(n,3,2)*delrho(n,3,2)
c
            de11 = de - 11d0
            de47 = 47d0 - 7d0*de
c
            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhob2
            d1t(1) = F19*(rhob*(1d0-3d0*de-de11*rhoa*rrho)
     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2))
            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
     &             - 2d0*rhob
            d1f(3) = -ab*om*t
            d2fgaa(1) = -ab*(d1om*t+om*d1t(1))
            d2fgaa(2) = -ab*(d1om*t+om*d1t(2))
c
            d2t(1) = -F19*(
     &           2d0*rhob*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhob*rrho2
     &                                         -2d0*de11*rhob*rrho3))
            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
     &           - rho(n,1)*(3d0+rhoa*rrho)*d1de
     &           - rhoab*((3d0+rhoa*rrho)*d2de-d1de*(rhoa-rhob)*rrho2
     &                                        +de11*(rhoa-rhob)*rrho3))
            d2t(3) = -F19*(
     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
     &                                         +2d0*de11*rhoa*rrho3))
     &           - 2d0
            d3fgaa(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgaa(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgaa(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
c
            t = F19*rhoab*de47-F43*rho2
            d1t(1) = F19*rhob*de47 - F79*rhoab*d1de - F83*rho(n,1)
            d1t(2) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,1)
            d1f(4) = -ab*om*t
            d2fgab(1) = -ab*(d1om*t+om*d1t(1))
            d2fgab(2) = -ab*(d1om*t+om*d1t(2))
c
            d2t(1) = -F79*(2d0*rhob*d1de+rhoab*d2de) - F83
            d2t(2) = F19*de47 - F79*(rho(n,1)*d1de+rhoab*d2de) - F83
            d2t(3) = -F79*(2d0*rhoa*d1de+rhoab*d2de) - F83
            d3fgab(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgab(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgab(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
c
            t = F19*rhoab*(1d0-3d0*de-de11*rhob*rrho)-rhoa2
            d1t(1) = F19*(rhob*(1d0-3d0*de-de11*rhob*rrho)
     &             + rhoab*(-(3d0+rhob*rrho)*d1de+de11*rhob*rrho2))
     &             - 2d0*rhoa
            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhob*rrho)
     &             - rhoab*((3d0+rhob*rrho)*d1de+de11*rhoa*rrho2))
            d1f(5) = -ab*om*t
            d2fgbb(1) = -ab*(d1om*t+om*d1t(1))
            d2fgbb(2) = -ab*(d1om*t+om*d1t(2))
c
            d2t(1) = -F19*(
     &           2d0*rhob*((3d0+rhob*rrho)*d1de-de11*rhob*rrho2)
     &           +  rhoab*((3d0+rhob*rrho)*d2de-2d0*d1de*rhob*rrho2
     &                                         +2d0*de11*rhob*rrho3))
     &           - 2d0
            d2t(2) = F19*(1d0-3d0*de-de11*rhob*rrho
     &           - rho(n,1)*(3d0+rhob*rrho)*d1de
     &           - rhoab*((3d0+rhob*rrho)*d2de+d1de*(rhoa-rhob)*rrho2
     &                                        -de11*(rhoa-rhob)*rrho3))
            d2t(3) = -F19*(
     &           2d0*rhoa*((3d0+rhob*rrho)*d1de+de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhob*rrho)*d2de+2d0*d1de*rhoa*rrho2
     &                                         -2d0*de11*rhoa*rrho3))
            d3fgbb(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgbb(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgbb(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
c
            f = c1*f1 + c2*f2 + gaa*d1f(3) + gab*d1f(4) + gbb*d1f(5)
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*d2fgaa(1) + gab*d2fgab(1) + gbb*d2fgbb(1)
            d1f(2) = c1*d1f1(2) + c2*d1f2(2)
     &             + gaa*d2fgaa(2) + gab*d2fgab(2) + gbb*d2fgbb(2)
c
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*d3fgaa(1) + gab*d3fgab(1) + gbb*d3fgbb(1)
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*d3fgaa(2) + gab*d3fgab(2) + gbb*d3fgbb(2)
            d2f(3) = c1*d2f1(3) + c2*d2f2(3)
     &             + gaa*d3fgaa(3) + gab*d3fgab(3) + gbb*d3fgbb(3)
c
            Ec = Ec + f*fac*qwght(n)
            if (ldew) func(n) = func(n) + f*fac
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            Amat(n,2) = Amat(n,2) + d1f(2)*fac
c            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
c            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*fac
c            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(5)*fac
cc
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*fac
c            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*fac
c            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + d2f(3)*fac
c            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*fac
c            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*fac
c            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*fac
c            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + d2fgaa(2)*fac
c            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + d2fgab(2)*fac
c            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + d2fgbb(2)*fac
c
c   20    continue
      endif
      return
      end
