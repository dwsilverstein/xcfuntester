      Subroutine xc_dirac(n, tol_rho, fac, lfac, nlfac, rho, Amat,
     &                       Amat2, Amat3, nq, ipol, Ex, qwght, ldew,
     &                       func)
c
C$Id: xc_dirac.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      Implicit none
c
      integer nq, ipol
      double precision fac, Ex
      logical ldew, lfac, nlfac
      double precision func(nq)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,3)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(n,2)
      double precision Amat2(n,3)
      double precision Amat3(n,4)
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
c
c Daniel (5-23-12): P1 and P2 are for restricted calculations, P3 
c and P4 are for unrestricted calculations.
c
c     P1 =       -(3/PI)**(1/3)
c     P2 = -(3/4)*(3/PI)**(1/3)
c     P3 =       -(6/PI)**(1/3)
c     P4 = -(3/4)*(6/PI)**(1/3)
c
      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      double precision rho13, rho32, rho33, one_third, term
      double precision two_ninth, rhom23
      Parameter (one_third = 1.d0/3.d0)
      Parameter (two_ninth = 2.d0/9.d0)
      integer n
c
      term = 0.d0
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
        rho13=rho(n,1)**one_third
        Amat(n,1) = Amat(n,1) + rho13
        Amat2(n,1) = Amat2(n,1) + rho13/rho(n,1)
        rhom23 = rho13/rho(n,1) 
        Amat3(n,1) = Amat3(n,1) + rhom23/rho(n,1)
        term = term + rho(n,1)*rho13*qwght(n)
        if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
        Amat(n,1) = Amat(n,1)*P1*fac
        Amat2(n,1) = Amat2(n,1)*P1*fac*one_third
        Amat3(n,1) = Amat3(n,1)*P1*fac*-two_ninth
        Ex = Ex + P2*term*fac
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
        rho32=0.0d0
        rho33=0.0d0
        if (rho(n,2).gt.tol_rho) rho32=rho(n,2)**one_third
        if (rho(n,3).gt.tol_rho) rho33=rho(n,3)**one_third
        Amat(n,1) = Amat(n,1) + P3*rho32*fac
        Amat(n,2) = Amat(n,2) + P3*rho33*fac
        Amat2(n,1) = Amat2(n,1) + one_third*P3*rho32/rho(n,2)*fac
        Amat2(n,3) = Amat2(n,3) + one_third*P3*rho33/rho(n,3)*fac
        Amat3(n,1) = Amat3(n,1) - two_ninth*P3*rho32/(rho(n,2)**2)*fac
        Amat3(n,4) = Amat3(n,4) - two_ninth*P3*rho33/(rho(n,3)**2)*fac
        if (ldew)func(n) = func(n) + ( rho32*rho(n,2) + 
     1                                 rho33*rho(n,3)   )*P4*fac
        term = term + (rho32*rho(n,2)+rho33*rho(n,3))*qwght(n)
        Ex = Ex + P4*term*fac
c
      endif
      return
      end
