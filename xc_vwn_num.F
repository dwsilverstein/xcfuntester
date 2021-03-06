c
c     General implementation of VWN functionals
c
c     BGJ (4/00)
*
* $Id: xc_vwn_a.fh 19696 2010-10-29 16:53:42Z d3y133 $
*
c     Meanings of macros used in constructing different VWN parameterizations
c     The determining macro FNL_VWN_xxx is set in the .F file
c
c     VWN_NAME - name of subroutine
c     USE_VWN_DEC - defined if VWN form for deltaEc(rs,zeta), Eq. 3.2, is used
c     NEED_RPA_RS1 - defined if the RPA epsilon difference deltaERPA(rs,1) is
c                    needed in addition to the Monte Carlo value deltaEc(rs,1)
c     The following are determined based on how the above are set
c     N_EPS - number of different parameterizations of epsilon involved
c     BETA_SET_RS1 - governs which set of parameters to use for the epsilon
c                    difference in beta in Eq. 3.2
      Subroutine xc_vwn_num(n, tol_rho, fac, rho,
     &     Amat, Amat2, Amat3, nq, ipol, Ec, qwght, ldew, func,
     &     functional)

c Daniel (5-23-12): Which VWN is being used?
c For LDA, we use VWN 5
c#if defined(FNL_VWN_1)
c#define VWN_NAME xc_vwn_1

c#elif defined(FNL_VWN_1_RPA)
c#define VWN_NAME xc_vwn_1_rpa

c#elif defined(FNL_VWN_2)
c#define VWN_NAME xc_vwn_2
c#define USE_VWN_DEC
c#define NEED_RPA_RS1

c#elif defined(FNL_VWN_3)
c#define VWN_NAME xc_vwn_3
c#define USE_VWN_DEC
c#define NEED_RPA_RS1

c#elif defined(FNL_VWN_4)
c#define VWN_NAME xc_vwn_4
c#define USE_VWN_DEC

c#elif defined(FNL_VWN_5)
c#define VWN_NAME xc_vwn_5
c#define USE_VWN_DEC
c#endif

c     Set how many different parameterizations of epsilon are involved
c#ifndef USE_VWN_DEC
c#define N_EPS 2
c#else
c#ifndef NEED_RPA_RS1
c#define N_EPS 3
c#else
c#define N_EPS 5
c#endif
c#endif

c#ifdef USE_VWN_DEC
cc     Set which parameters to use for epsilon difference in beta(rs)
c#ifdef NEED_RPA_RS1
c#define BETA_SET_RS1 4
c#else
c#define BETA_SET_RS1 1
c#endif
c#endif

c      Subroutine VWN_NAME(n, tol_rho, fac, rho,
c     &     Amat, Amat2, Amat3, nq, ipol, Ec, qwght, ldew, func)
      implicit none
c
      integer nq, ipol
      double precision fac, Ec, tol_rho
      logical ldew
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
c     Partial Derivatives of the Correlation Energy Functional
c
      double precision Amat(n,2)
      double precision Amat2(n,3)
      double precision Amat3(n,4)
c
      double precision onethird, fourthirds, twothirds, pi
      double precision seventhirds, threehalf
      double precision fivethirds
      Parameter (onethird = 1.D0/3.D0, fourthirds = 4.D0/3.D0)
      Parameter (twothirds = 2.D0/3.D0)
      Parameter (fivethirds = 5.D0/3.D0)
      Parameter (seventhirds = 7.D0/3.D0)
      Parameter (threehalf = 3.D0/2.D0)
      Parameter (pi = 3.1415926535898D0)
c
c     Functional Parameters (different flavors of VWN)
c     1 = EcP - Monte Carlo
c     2 = EcF - Monte Carlo
c     3 = alpha - Monte Carlo
c     4 = EcP - RPA
c     5 = EcF - RPA
c     6 = alpha - RPA
c

c See text below 
      double precision A(6), b(6), c(6), x0(6), Q(6),
     &     ct1(6), ct2(6), ct3(6)
      save A, b, c, x0, Q, ct1, ct2, ct3
      data A / 0.0621814D0, 0.D0, 0.D0,
     &         0.0621814D0, 0.D0, 0.D0 / ! see initialization below
      data b / 3.72744D0, 7.06042D0, 1.13107D0,
     &        13.0720D0, 20.1231D0, 1.06835D0 /
      data c / 12.9352D0, 18.0578D0, 13.0045D0,
     &         42.7198D0, 101.578D0, 11.4813D0 /
      data x0 / -0.10498D0, -0.32500D0, -0.00475840D0,
     &          -0.409286D0, -0.743294D0, -0.228344D0 /
c
      double precision e(5), d1e(5), rhoval, rs, d1rs, x, d1x,
     &     rchi, d1chi, t1, t2, t3, d1t1, d1t2, d1t3, s, d1s,
     &     d1zeta(2), d1ersz(2), d1edrho(2), zeta, fz, d1fz, eps,
     &     dec_rs1, dec_rsz, d1dec_rs1, d1dec_rsz(2)
      double precision d2fz0
c#ifdef NEED_RPA_RS1
      double precision derpa_rs1, d1derpa_rs1
      double precision d2derpa_rs1
c Daniel (5-29-12): Third order derivative stuff for Delta[epsilon_RPA]
      double precision d3derpa_rs1
c#endif
c#ifdef USE_VWN_DEC
      double precision devwn_rsz, d1devwn_rsz(2), zeta2, zeta3, zeta4,
     &     beta_rs1, d1beta_rs1, t_vwn, d1t_vwn
      double precision d2beta_rs1, d2t_vwn, d2devwn_rsz(3)
c Daniel (5-30-12): Third derivative stuff
      double precision d3beta_rs1, d3t_vwn, d3devwn_rsz(4)
c#endif
      double precision d2rs, d2x, d2t1, d2t2, d2t3, d2s,
     &     d2e(5), d2zeta(3), d2dec_rs1, d2dec_rsz(3),
     &     d2ersz(3), d2edrho(3), d2fz, rrho2
c Daniel (5-23-12): Third derivative stuff
      double precision d3rs, d3x, d3t1, d3t2, d3t3, d3s,
     &     d3e(5), d3fz, d3zeta(4), rrho3, d3dec_rs1,
     &     d3dec_rsz(4), d3ersz(4), d3edrho(4)
      double precision p0, p1, p2, p3
c Daniel (5-29-12): Third derivative coefficient
      double precision p4
c
      integer i, n, initial
      save initial
      data initial /1/
c Daniel (6-10-12): Annoying C-preprocessor variables
      logical USE_VWN_DEC, NEED_RPA_RS1
      integer N_EPS, BETA_SET_RS1
c
      character*10 functional
c Define the preprocessor variables.
      if (functional.eq.'cvwn1') then
        USE_VWN_DEC = .false.
        NEED_RPA_RS1 = .false.
        N_EPS = 2
      else if (functional.eq.'cvwn2') then
        USE_VWN_DEC = .true.
        NEED_RPA_RS1 = .true.
        N_EPS = 5
        BETA_SET_RS1 = 4
      else if (functional.eq.'cvwn3') then
        USE_VWN_DEC = .true.
        NEED_RPA_RS1 = .true.
        N_EPS = 5
        BETA_SET_RS1 = 4
      else if (functional.eq.'cvwn4') then
        USE_VWN_DEC = .true.
        NEED_RPA_RS1 = .false.
        N_EPS = 3
        BETA_SET_RS1 = 1
      else if (functional.eq.'cvwn5') then
        USE_VWN_DEC = .true.
        NEED_RPA_RS1 = .false.
        N_EPS = 3
        BETA_SET_RS1 = 1
      else if (functional.eq.'cvwn1rpa') then
        USE_VWN_DEC = .false.
        NEED_RPA_RS1 = .false.
        N_EPS = 2
      end if
        

c
c     Define miscellaneous parameters.
c
c Daniel (5-23-12): 
c p0 is for defining the Wigner-Seitz radius
c p1 is for the function, f[zeta]
c p2 is for the first derivative of f[zeta]
c p3 is for the second derivative of f[zeta] 
c d2fz0 is the second derivative of f[zeta], zeta = 0
      p0 = (1.0d0/(fourthirds*pi))**onethird
      p1 = 0.5D0/(2.d0**onethird - 1.d0)
      p2 = fourthirds*p1
      p3 = onethird*p2
      p4 = -twothirds*p3
      d2fz0 = 2.d0*p3
      if (initial.eq.1)then
         initial = 0
c
c     Initialize data values that are expressions here
c     instead of in data statement
c
         A(2) = 0.5D0*0.0621814D0
         A(3) = -1.D0/(3.D0*pi*pi)
         A(5) = A(2)
         A(6) = A(3)
c
c Daniel (5-23-12): Analytical formulas for the correlation functional
c from Eq. [4.4] in VWN.  This is done for all VWN functionals, then
c the parameters of the functional requested by the user are employed
c in the calculation.
c ct1 = constant for term 1
c ct2 = constant for term 2
c ct3 = constant for term 3
         do i = 1, 6
            Q(i) = sqrt(4.d0*c(i)-b(i)*b(i))
            ct1(i) = A(i)
            ct2(i) =-A(i)*b(i)*x0(i)/(x0(i)*x0(i)+b(i)*x0(i)+c(i))
            ct3(i) = A(i)*(2.d0*b(i)/Q(i))
     &           * ((c(i)-x0(i)*x0(i))/(x0(i)*x0(i)+b(i)*x0(i)+c(i)))
c           Convert interpolation formulas from Rydberg to Hartree units
            ct1(i) = 0.5D0*ct1(i)
            ct2(i) = 0.5D0*ct2(i)
            ct3(i) = 0.5D0*ct3(i)
         enddo
         if (functional.eq.'cvwn1rpa') then
c        Copy RPA over Monte Carlo for all parameters
           do i = 1, 3       
             A(i) = A(4)
             b(i) = b(4)
             x0(i) = x0(4)
             Q(i) = Q(4)
             ct1(i) = ct1(4)
             ct2(i) = ct2(4)
             ct3(i) = ct3(4)
           enddo
c         call dcopy(3, A(4), 1, A(1), 1)
c         call dcopy(3, b(4), 1, b(1), 1)
c         call dcopy(3, c(4), 1, c(1), 1)
c         call dcopy(3, x0(4), 1, x0(1), 1)
c         call dcopy(3, Q(4), 1, Q(1), 1)
c         call dcopy(3, ct1(4), 1, ct1(1), 1)
c         call dcopy(3, ct2(4), 1, ct2(1), 1)
c         call dcopy(3, ct3(4), 1, ct3(1), 1)
         endif
c        The appropriate spin stiffness parameters are always assumed
c        to be in slot 3 - these start out as the Monte Carlo values,
c        so copy RPA values here when RPA is actually used for this
         if ((functional.eq.'cvwn2').or.
     &       (functional.eq.'cvwn3').or.
     &       (functional.eq.'cvwn4')) then
           A(3) = A(6)
           b(3) = b(6)
           c(3) = c(6)
           x0(3) = x0(6)
           Q(3) = Q(6)
           ct1(3) = ct1(6)
           ct2(3) = ct2(6)
           ct3(3) = ct3(6)
         endif
      endif
c
c     ======> BOTH SPIN-RESTRICTED AND SPIN-UNRESTRICTED <======
c
c      do 200 n = 1, nq
c         if (rho(n,1).lt.tol_rho)goto 200
c        if (abs(rho(n,1)).lt.tol_rho)goto 200
c
c Daniel (5-23-12): Check the density at the quadrature point against
c the tolerance.  rs => Wigner-Seitz radius 
         rhoval = rho(n,1)
         if(rhoval.gt.tol_rho) then
            rs = p0*rhoval**(-onethird)
         else
            rs=0d0
         endif
c Daniel (5-23-12): x = rs**{1/2} is defined under Eq. [4.3] in VWN.
c d1rs => 1st derivative of rs with respect to rho
c d1x => 1st derivative of x with respect to rs
         d1rs = -onethird*rs/rhoval
         x = sqrt(rs)
         d1x = 0.5d0/x
c Daniel (5-23-12)
c d2rs => 2nd derivative of rs with respect to rho
c d2x => 2nd derivative of x with respect to rs
         d2rs = -fourthirds*d1rs/rhoval
         d2x = -0.5d0*d1x/rs
c Daniel (5-23-12): Third order derivative stuff
         d3rs = -seventhirds*d2rs/rhoval 
         d3x = -threehalf*d2x/rs 
c
c        Evaluate the individual correlation energy formulas
c
c Daniel (5-28-12): Road map to the variables
c rchi = reciprocal of chi
c d1chi = derivative of chi with respect to x
c s = stuff inside of inverse tangent
c t1 = first logarithm
c t2 = second logarithm
c t3 = both inverse tangents
         do i = 1, N_EPS
            rchi = 1.0d0/(rs + b(i)*x + c(i))
            d1chi = 2.0d0*x + b(i)
            s = Q(i)/d1chi
            d1s = -2.0d0*s/d1chi
            t1 = log(abs(rs*rchi))
            d1t1=(x - x0(i))*(x - x0(i))*rchi
            t2 = log(abs(d1t1))
            t3 = atan(s)
c d1t1 = derivative of the first term
c d1t2 = derivative of the second term
c d1t3 = derivative of the third term
            d1t1 = 2.0d0/x - d1chi*rchi
            d1t2 = 2.0d0/(x-x0(i)) - d1chi*rchi
            d1t3 = d1s/(1.0d0 + s*s)
c e(i) is the correlation energy per particle
c d1e(i) is the derivative of the correlation energy per particle
            e(i) = ct1(i)*t1 + ct2(i)*t2 + ct3(i)*t3
            d1e(i) = ct1(i)*d1t1 + ct2(i)*d1t2 + ct3(i)*d1t3
            d2s = -4.d0*d1s/d1chi
            d2t1 =-2.0d0/rs + (d1chi*rchi)**2 - 2.d0*rchi
            d2t2 =-2.0d0/(x-x0(i))**2 + (d1chi*rchi)**2 - 2.d0*rchi
            d2t3 = (d2s - 2.d0*s*d1s*d1t3)/(1.0d0 + s*s)
c d2e(i) is the second derivative of the correlation energy per particle
            d2e(i) = ct1(i)*d2t1 + ct2(i)*d2t2 + ct3(i)*d2t3
c Daniel (5-29-12): Third order derivative stuff
c            d3s = -6.0d0*d2s/d1chi
c Might be better to implement this using the following, so that we 
c don't rely on the second derivative being defined:
            d3s = 24.0d0*d1s/(d1chi*d1chi)
            d3t1 = 4.0d0/x**3 + 6.0d0*d1chi*rchi*rchi
     &        - 2.0d0*(d1chi*rchi)**3
            d3t2 = 4.0d0/(x-x0(i))**3 + 6.0d0*d1chi*rchi*rchi
     &        - 2.0d0*(d1chi*rchi)**3
            d3t3 = (d3s + (6.0d0*s*s-2.0d0)*d1s*d1t3*d1t3 
     &        - 6.0d0*s*d2s*d1t3)/(1.0d0 + s*s)
            d3e(i) = ct1(i)*d3t1 + ct2(i)*d3t2 + ct3(i)*d3t3    
c
c           Transform derivatives wrt x to derivatives wrt rs
c
c Daniel (5-29-12): Third order derivative stuff
c Need to do this first because we need the first and second order
c derivatives with respect to x aren't replaced.
            d3e(i) = d3e(i)*d1x*d1x*d1x + 3.0d0*d2e(i)*d1x*d2x
     &        + d1e(i)*d3x            
c           Do 2nd derivative first so the x first derivative in d1e
c           is not lost
            d2e(i) = d2e(i)*d1x*d1x + d1e(i)*d2x
            d1e(i) = d1e(i)*d1x
         enddo
c
c        Compute the polarization function and its derivatives
c
c Daniel (5-29-12): The spin polarization function is defined as
c the alpha-density minus the beta-density divided by the total density
c This is nonzero only if an open-shell calculation is performed.
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,2) - rho(n,3))/rhoval
         endif
c By definition, zeta should be between -1 and 1.
         if (zeta.gt.1.d0)then
            zeta = 1.d0
         elseif (zeta.lt.-1.d0)then
            zeta =-1.d0
         endif
c Daniel (5-29-12): This is [2.3] in VWN.
         fz = ((1.d0+zeta)**fourthirds + 
     &         (1.d0-zeta)**fourthirds - 2.d0)*p1
c Daniel (5-29-12): Derivative of f[zeta] with respect
c to zeta.
         d1fz = ((1.d0+zeta)**onethird - 
     &           (1.d0-zeta)**onethird)*p2
c Daniel (5-29-12): Derivatives of zeta with respect to the 
c spin-densities (1 = alpha, 2 = beta)
         d1zeta(1) = (1.d0-zeta)/rhoval
         d1zeta(2) =-(1.d0+zeta)/rhoval
c Daniel (5-29-12): If zeta is less than the tolerance on the
c electron density, we insert zeta=0 into the equation.
         if(dabs(zeta).lt.tol_rho) then
            d2fz = d2fz0 
         else
c Daniel (5-29-12): If zeta is ~ -1 or ~ +1, we can have a problem with
c numerical inaccuracy.  This part screens for that situation.
            if (dabs(zeta+1.d0).le.tol_rho) then
              d2fz = ((1.d0-zeta)**(-twothirds))*p3
            else if (dabs(zeta-1.d0).le.tol_rho) then
              d2fz = ((1.d0+zeta)**(-twothirds))*p3
            else
              d2fz = ((1.d0+zeta)**(-twothirds) +
     &                (1.d0-zeta)**(-twothirds))*p3
            endif
         endif
         rrho2 = 2.d0/(rhoval*rhoval)
c        1 = aa, 2 = ab, 3 = bb
c Daniel (5-29-12): Second derivatives of zeta with respect to the 
c spin density.
         d2zeta(1) =-rrho2*(1.d0-zeta)
         d2zeta(2) = rrho2*zeta
         d2zeta(3) = rrho2*(1.d0+zeta)
c
c Daniel (5-29-12): Third order derivatives of spin polarization stuff
         if (dabs(zeta+1.0d0).le.tol_rho) then
           d3fz = (-(1.0d0-zeta)**(-fivethirds))*p4 
         else if (dabs(zeta-1.0d0).le.tol_rho) then
           d3fz = ((1.0d0+zeta)**(-fivethirds))*p4
         else
           d3fz = ((1.0d0+zeta)**(-fivethirds) -
     &             (1.0d0-zeta)**(-fivethirds))*p4
         end if          
         rrho3 = 1.0d0/(rhoval*rhoval*rhoval) 
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
         d3zeta(1) = 6.0d0*(1.0d0-zeta)*rrho3
         d3zeta(2) = 2.0d0*(1.0d0-3.0d0*zeta)*rrho3
         d3zeta(3) = -2.0d0*(1.0d0+3.0d0*zeta)*rrho3
         d3zeta(4) = -6.0d0*(1.0d0+3.0d0*zeta)*rrho3
c
         dec_rs1 = e(2)-e(1)
         d1dec_rs1 = d1e(2)-d1e(1)
         d2dec_rs1 = d2e(2)-d2e(1)
c Daniel (5-29-12): Third order derivative of Delta[epsilon_c]
         d3dec_rs1 = d3e(2)-d3e(1)         
         if (NEED_RPA_RS1) then
           derpa_rs1 = e(5)-e(4)
           d1derpa_rs1 = d1e(5)-d1e(4)
           d2derpa_rs1 = d2e(5)-d2e(4)
c Daniel (5-29-12): Third order derivative of Delta[epsilon_RPA]
           d3derpa_rs1 = d3e(5)-d3e(4)
         endif
c
c Daniel (5-29-12): This part deals with the functional form of
c Delta[epsilon] 
         if (USE_VWN_DEC) then
c Delta[espilon](rs,1) in Eq. [4.8] in VWN and its derivative
           beta_rs1 = e(BETA_SET_RS1+1)-e(BETA_SET_RS1)
           d1beta_rs1 = d1e(BETA_SET_RS1+1)-d1e(BETA_SET_RS1)
c Various zetas
           zeta2 = zeta*zeta
           zeta3 = zeta2*zeta
           zeta4 = zeta3*zeta
c Eq. [4.8] in VWN and its derivative with respect to rs 
           t_vwn = d2fz0*beta_rs1-e(3)
           d1t_vwn = d2fz0*d1beta_rs1-d1e(3)
c Delta[epsilon_{VWN}](rs,zeta)
           devwn_rsz = fz/d2fz0*(e(3)+t_vwn*zeta4)
c Derivatives: 1 = drs, 2 = dzeta
           d1devwn_rsz(1) = fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
           d1devwn_rsz(2) = d1fz/d2fz0*(e(3)+t_vwn*zeta4)
     &          + fz/d2fz0*t_vwn*4.d0*zeta3
           d2beta_rs1 = d2e(BETA_SET_RS1+1)-d2e(BETA_SET_RS1)
           d2t_vwn = d2fz0*d2beta_rs1-d2e(3)
c Derivatives: 1 = drsdrs, 2 = drsdzeta, 3 = dzetadzeta
           d2devwn_rsz(1) = fz/d2fz0*(d2e(3)+d2t_vwn*zeta4)
           d2devwn_rsz(2) = d1fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
     &          + fz/d2fz0*d1t_vwn*4.d0*zeta3
           d2devwn_rsz(3) = d2fz/d2fz0*(e(3)+t_vwn*zeta4)
     &          + d1fz/d2fz0*t_vwn*8.d0*zeta3
     &          + fz/d2fz0*t_vwn*12.d0*zeta2
c Daniel (5-30-12): Third order derivative stuff
           d3beta_rs1 = d3e(BETA_SET_RS1+1)-d3e(BETA_SET_RS1) 
           d3t_vwn = d2fz0*d3beta_rs1-d3e(3)
c Derivatives: 1 = drsdrsdrs, 2 = drsdrsdzeta, 3 = drsdzetadzeta,
c              4 = dzetadzetadzeta
           d3devwn_rsz(1) = fz/d2fz0*(d3e(3)+d3t_vwn*zeta4)
           d3devwn_rsz(2) = d1fz/d2fz0*(d2e(3)+d2t_vwn*zeta4)
     &          + fz/d2fz0*d2t_vwn*4.0d0*zeta3
           d3devwn_rsz(3) = d2fz/d2fz0*(d1e(3)+d1t_vwn*zeta4)
     &          + d1fz/d2fz0*d1t_vwn*8.0d0*zeta3
     &          + fz/d2fz0*d1t_vwn*12.0d0*zeta2
           d3devwn_rsz(4) = d3fz/d2fz0*(e(3)+t_vwn*zeta4)
     &          + d2fz/d2fz0*t_vwn*12.0d0*zeta3
     &          + d1fz/d2fz0*t_vwn*36.0d0*zeta2
     &          + fz/d2fz0*t_vwn*24.0d0*zeta
         endif
c
c --- Functional-dependent unrestricted part, deltaEc(rs,zeta),
c --- goes in the .F file after inclusion of this header
      if (functional.eq.'cvwn1') then
        dec_rsz = dec_rs1*fz
        d1dec_rsz(1) = d1dec_rs1*fz
        d1dec_rsz(2) = dec_rs1*d1fz
        d2dec_rsz(1) = d2dec_rs1*fz
        d2dec_rsz(2) = d1dec_rs1*d1fz
        d2dec_rsz(3) = dec_rs1*d2fz
        d3dec_rsz(1) = d3dec_rs1*fz
        d3dec_rsz(2) = d2dec_rs1*d1fz
        d3dec_rsz(3) = d1dec_rs1*d2fz
        d3dec_rsz(4) = dec_rs1*d3fz
      else if (functional.eq.'cvwn2') then
        dec_rsz = devwn_rsz + (dec_rs1-derpa_rs1)*fz
        d1dec_rsz(1) = d1devwn_rsz(1) + (d1dec_rs1-d1derpa_rs1)*fz
        d1dec_rsz(2) = d1devwn_rsz(2) + (dec_rs1-derpa_rs1)*d1fz
        d2dec_rsz(1) = d2devwn_rsz(1) + (d2dec_rs1-d2derpa_rs1)*fz
        d2dec_rsz(2) = d2devwn_rsz(2) + (d1dec_rs1-d1derpa_rs1)*d1fz
        d2dec_rsz(3) = d2devwn_rsz(3) + (dec_rs1-derpa_rs1)*d2fz
        d3dec_rsz(1) = d3devwn_rsz(1) + (d3dec_rs1-d3derpa_rs1)*fz
        d3dec_rsz(2) = d3devwn_rsz(2)
     &               + (d2dec_rs1-d2derpa_rs1)*d1fz
        d3dec_rsz(3) = d3devwn_rsz(3)
     &               + (d1dec_rs1-d1derpa_rs1)*d2fz
        d3dec_rsz(4) = d3devwn_rsz(4)
     &               + (dec_rs1-derpa_rs1)*d3fz
      else if (functional.eq.'cvwn3') then
        dec_rsz = (dec_rs1/derpa_rs1)*devwn_rsz
        d1dec_rsz(1) = (d1dec_rs1/derpa_rs1)*devwn_rsz
     &       - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
     &       + (dec_rs1/derpa_rs1)*d1devwn_rsz(1)
        d1dec_rsz(2) = (dec_rs1/derpa_rs1)*d1devwn_rsz(2)
        d2dec_rsz(1) = (dec_rs1/derpa_rs1)*d2devwn_rsz(1)
     &       + (d2dec_rs1/derpa_rs1)*devwn_rsz
     &       - (dec_rs1*d2derpa_rs1/derpa_rs1**2)*devwn_rsz
     &       + 2.d0*((d1dec_rs1/derpa_rs1)*d1devwn_rsz(1)
     &       - (d1dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
     &       - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*d1devwn_rsz(1)
     &       + (dec_rs1*d1derpa_rs1**2/derpa_rs1**3)*devwn_rsz)
        d2dec_rsz(2) = (dec_rs1/derpa_rs1)*d2devwn_rsz(2)
     &       + ((d1dec_rs1/derpa_rs1)
     &       - (dec_rs1*d1derpa_rs1/derpa_rs1**2))*d1devwn_rsz(2)
        d2dec_rsz(3) = (dec_rs1/derpa_rs1)*d2devwn_rsz(3)
        d3dec_rsz(1) = (dec_rs1/derpa_rs1)*d3devwn_rsz(1)
     &       + 3.0d0*(derpa_rs1*d1dec_rs1
     &          - dec_rs1*d1derpa_rs1)*d2devwn_rsz(1)/derpa_rs1**2
     &       + 3.0d0*(derpa_rs1**2*d2dec_rs1
     &         - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
     &          +dec_rs1*d2derpa_rs1)
     &          +2.0d0*dec_rs1*d1derpa_rs1**2)
     &          *d1devwn_rsz(1)/derpa_rs1**3
     &       + (d3dec_rs1*derpa_rs1**3
     &          + 6.0d0*derpa_rs1*d1derpa_rs1
     &           *(d1dec_rs1*d1derpa_rs1+dec_rs1*d2derpa_rs1)
     &          - derpa_rs1**2*(3.0d0*d2dec_rs1*d1derpa_rs1
     &           + 3.0d0*d1dec_rs1*d2derpa_rs1
     &           + dec_rs1*d3derpa_rs1)
     &          - 6.0d0*dec_rs1*d1derpa_rs1**3)
     &          *devwn_rsz/derpa_rs1**4
        d3dec_rsz(2) = (dec_rs1/derpa_rs1)*d3devwn_rsz(2)
     &       + 2.0d0*(derpa_rs1*d1dec_rs1
     &          - dec_rs1*d1derpa_rs1)*d2devwn_rsz(2)/derpa_rs1**2
     &       + (derpa_rs1**2*d2dec_rs1
     &         - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
     &          +dec_rs1*d2derpa_rs1)
     &          +2.0d0*dec_rs1*d1derpa_rs1**2)
     &          *d1devwn_rsz(2)/derpa_rs1**3
        d3dec_rsz(3) = (dec_rs1/derpa_rs1)*d3devwn_rsz(3)
     &       + (derpa_rs1*d1dec_rs1
     &          - dec_rs1*d1derpa_rs1)*d2devwn_rsz(3)/derpa_rs1**2
        d3dec_rsz(4) = (dec_rs1/derpa_rs1)*d3devwn_rsz(4)
      else if (functional.eq.'cvwn4') then
        dec_rsz = devwn_rsz
        d1dec_rsz(1) = d1devwn_rsz(1)
        d1dec_rsz(2) = d1devwn_rsz(2)
        d2dec_rsz(1) = d2devwn_rsz(1)
        d2dec_rsz(2) = d2devwn_rsz(2)
        d2dec_rsz(3) = d2devwn_rsz(3)
        d3dec_rsz(1) = d3devwn_rsz(1)
        d3dec_rsz(2) = d3devwn_rsz(2)
        d3dec_rsz(3) = d3devwn_rsz(3)
        d3dec_rsz(4) = d3devwn_rsz(4)
      else if (functional.eq.'cvwn5') then
        dec_rsz = devwn_rsz
        d1dec_rsz(1) = d1devwn_rsz(1)
        d1dec_rsz(2) = d1devwn_rsz(2)
        d2dec_rsz(1) = d2devwn_rsz(1)
        d2dec_rsz(2) = d2devwn_rsz(2)
        d2dec_rsz(3) = d2devwn_rsz(3)
        d3dec_rsz(1) = d3devwn_rsz(1)
        d3dec_rsz(2) = d3devwn_rsz(2)
        d3dec_rsz(3) = d3devwn_rsz(3)
        d3dec_rsz(4) = d3devwn_rsz(4)
      else if (functional.eq.'cvwn1rpa') then
        dec_rsz = dec_rs1*fz
        d1dec_rsz(1) = d1dec_rs1*fz
        d1dec_rsz(2) = dec_rs1*d1fz
        d2dec_rsz(1) = d2dec_rs1*fz
        d2dec_rsz(2) = d1dec_rs1*d1fz
        d2dec_rsz(3) = dec_rs1*d2fz
        d3dec_rsz(1) = d3dec_rs1*fz
        d3dec_rsz(2) = d2dec_rs1*d1fz
        d3dec_rsz(3) = d1dec_rs1*d2fz
        d3dec_rsz(4) = dec_rs1*d3fz
      end if

c
c     General implementation of VWN functionals (continued)
c
c     BGJ (4/00)
*
* $Id: xc_vwn_b.fh 19696 2010-10-29 16:53:42Z d3y133 $
*
c     Finish off the unrestricted case:
c     Assemble the entire functional and its derivatives given the
c     parameterization-dependent part deltaEc(rs,zeta) and its derivatives
c
c Daniel (5-30-12): Actually calculate the correlation energy.
         eps = e(1) + dec_rsz
c Derivative with respect to rs
         d1ersz(1) = d1e(1) + d1dec_rsz(1)
c Derivative with respect to zeta
         d1ersz(2) = d1dec_rsz(2)
c Convert the derivatives to derivatives with respect to the density.
c 1 = drhoa, 2 = drhob
         d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
         d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
         Ec = Ec + eps*qwght(n)*rhoval*fac
         if (ldew) func(n) = func(n) + eps*rhoval*fac
         Amat(n,1) = Amat(n,1) + (eps + rhoval*d1edrho(1))*fac
         if (ipol.eq.2)
     &      Amat(n,2) = Amat(n,2) + (eps + rhoval*d1edrho(2))*fac
c        1 = rsrs, 2 = rsz, 3 = zz
         d2ersz(1) = d2e(1) + d2dec_rsz(1)
         d2ersz(2) = d2dec_rsz(2)
         d2ersz(3) = d2dec_rsz(3)
c        1 = aa, 2 = ab, 3 = bb
c Daniel (5-30-12): rs is independent of spin density, but zeta is by
c definition dependent on the spin density.
         d2edrho(1) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*d1zeta(1)*2.d0
     &              + d2ersz(3)*d1zeta(1)*d1zeta(1)
     &              + d1ersz(1)*d2rs 
     &              + d1ersz(2)*d2zeta(1)
         d2edrho(2) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
     &              + d2ersz(3)*d1zeta(1)*d1zeta(2)
     &              + d1ersz(1)*d2rs
     &              + d1ersz(2)*d2zeta(2)
         d2edrho(3) = d2ersz(1)*d1rs*d1rs
     &              + d2ersz(2)*d1rs*d1zeta(2)*2.d0
     &              + d2ersz(3)*d1zeta(2)*d1zeta(2)
     &              + d1ersz(1)*d2rs
     &              + d1ersz(2)*d2zeta(3)
         Amat2(n,1) = Amat2(n,1)
     &        + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
         Amat2(n,2) = Amat2(n,2)
     &        + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*fac
         Amat2(n,3) = Amat2(n,3)
     &        + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*fac
c Daniel (5-31-12): Third derivative stuff
c Daniel (5-31-12): Note that d3dec_rsz(#) is defined in the routine
c calling the header file.  It comes from d3devwn_rsz.
c 1 = rsrsrs, 2 = rsrsz, 3 = rszz, 4 = zzz
         d3ersz(1) = d3e(1) + d3dec_rsz(1) 
         d3ersz(2) = d3dec_rsz(2) 
         d3ersz(3) = d3dec_rsz(3)
         d3ersz(4) = d3dec_rsz(4)
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
         d3edrho(1) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*d1zeta(1)*d1zeta(1)*3.0d0 ! Bug?
     &              + d2ersz(2)*d1rs*d2zeta(1)*3.0d0
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d1zeta(1)*d2rs*3.0d0
     &              + d3ersz(2)*d1zeta(1)*d1rs*d1rs*3.0d0
     &              + d3ersz(4)*d1zeta(1)*d1zeta(1)*d1zeta(1)
     &              + d2ersz(3)*d1zeta(1)*d2zeta(1)*3.0d0
     &              + d1ersz(2)*d3zeta(1)
         d3edrho(2) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*(d1zeta(1)*d1zeta(1)
     &                              + d1zeta(1)*d1zeta(2)*2.0d0)
     &              + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     &                              + d2zeta(1))
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d2rs*(d1zeta(1)*2.0d0
     &                              + d1zeta(2))
     &              + d3ersz(2)*d1rs*d1rs*(d1zeta(2)
     &                                   + d1zeta(1)*2.0d0)
     &              + d3ersz(4)*d1zeta(2)*d1zeta(1)*d1zeta(1)
     &              + d2ersz(3)*(d1zeta(1)*d2zeta(2)*2.0d0
     &                         + d1zeta(2)*d2zeta(1))
     &              + d1ersz(2)*d3zeta(2) 
         d3edrho(3) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*(d1zeta(2)*d1zeta(2) 
     &                              + d1zeta(2)*d1zeta(1)*2.0d0)
     &              + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     &                              + d2zeta(3))
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d2rs*(d1zeta(2)*2.0d0
     &                              + d1zeta(1))
     &              + d3ersz(2)*d1rs*d1rs*(d1zeta(1)
     &                                   + d1zeta(2)*2.0d0)
     &              + d3ersz(4)*d1zeta(1)*d1zeta(2)*d1zeta(2)
     &              + d2ersz(3)*(d1zeta(2)*d2zeta(2)*2.0d0
     &                         + d1zeta(1)*d2zeta(3))
     &              + d1ersz(2)*d3zeta(3) 
         d3edrho(4) = d3ersz(1)*d1rs*d1rs*d1rs
     &              + d2ersz(1)*d1rs*d2rs*3.0d0
     &              + d3ersz(3)*d1rs*d1zeta(2)*d1zeta(2)*3.0d0
     &              + d2ersz(2)*d1rs*d2zeta(3)*3.0d0
     &              + d1ersz(1)*d3rs
     &              + d2ersz(2)*d1zeta(2)*d2rs*3.0d0
     &              + d3ersz(2)*d1zeta(2)*d1rs*d1rs*3.0d0
     &              + d3ersz(4)*d1zeta(2)*d1zeta(2)*d1zeta(2)
     &              + d2ersz(3)*d1zeta(2)*d2zeta(3)*3.0d0
     &              + d1ersz(2)*d3zeta(4)
         Amat3(n,1) = Amat3(n,1)
     &        + (3.0d0*d2edrho(1) + rhoval*d3edrho(1))*fac
         Amat3(n,2) = Amat3(n,2)
     &        + (d2edrho(1) + 2.0d0*d2edrho(2) + rhoval*d3edrho(2))*fac
         Amat3(n,3) = Amat3(n,3)
     &        + (2.0d0*d2edrho(2) + d2edrho(3) + rhoval*d3edrho(3))*fac
         Amat3(n,4) = Amat3(n,4)
     &        + (3.0d0*d2edrho(3) + rhoval*d3edrho(4))*fac
c  200 continue
c
      return
      end
