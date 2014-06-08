      subroutine vwn_wrapper(Amat, Amat2, Amat3,
     &  d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)

      implicit none

      integer n,ipol,nq
      logical lfac,nlfac,ldew
      double precision tol_rho,fac
      double precision rho(1,3),delrho(1,3,2)
      double precision Amat(1,2)
      double precision Amat2(1,3)
      double precision Ec
      double precision qwght(1),func(1)
      double precision Amat3(1,4)
      double precision tAmat3(1,4)
c
      double precision gamma
      double precision cam_omega
c
      double precision prho(1,3),mrho(1,3),p2rho(1,3),m2rho(1,3),
     &  drho(1,3)
c
      double precision p10Amat(1,2),p10Amat2(1,3)
      double precision m10Amat(1,2),m10Amat2(1,3)
      double precision p20Amat(1,2),p20Amat2(1,3)
      double precision m20Amat(1,2),m20Amat2(1,3)
      double precision p10Ec,m10Ec,p20Ec,m20Ec,d1E_rho(1,2),
     &                 d2E_rho(1,3)
c
      double precision dgamma,pgamma,mgamma,p2gamma,m2gamma
c
      double precision d3E_rho(1,4)
c
      character*10 functional

c --------------------------------------------------------------------
c Subroutine vwn_wrapper
c
c Purpose: Evaluate the derivatives of the VWN correlation functionals
c --------------------------------------------------------------------

c Initialize variables

      call var_init(n, ipol, nq, lfac, nlfac, ldew, tol_rho,
     &  fac, rho, delrho, qwght, gamma, cam_omega, prho, mrho,
     &  p2rho, m2rho, drho, dgamma, pgamma, mgamma, p2gamma, m2gamma) 

c ----------------------
c Analytical Derivatives
c ----------------------
c     Calculate the analytical derivatives of the correlation energy
c     functional of VWN with respect to the density and gradient of
c     the density.
c
      call xc_vwn_num(n, tol_rho, fac,  rho,
     &     Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &     ldew, func, functional)

c To correctly build the analytical derivative, we need to account
c for spin polarization.
c The factor of one half is needed for the second derivative to make
c the numerical result consistent with the analytical result.  This
c is a consequence of only calculating the first derivative for alpha
c spin electrons in xc_vwn_num.
      Amat2(1,1) = Amat2(1,1) 
     &           + 2.0d0*Amat2(1,2) 
     &           + Amat2(1,3)
      Amat2(1,1) = 0.50d0*Amat2(1,1)
 
      Amat3(1,1) = Amat3(1,1) 
     &           + 3.0d0*Amat3(1,2) 
     &           + 3.0d0*Amat3(1,3)
     &           + Amat3(1,4)
      Amat3(1,1) = Amat3(1,1)

c ---------------------
c Numerical Derivatives
c ---------------------
c     Evaluate the numerical derivatives of the correlation energy 
c     functional of VWN with respect to the density.
c    
c     Find the first, second, and third derivatives of E with respect to rho
c
c Step +1 rho
      call numder_init_lda(func, p10Ec, p10Amat, p10Amat2, tAmat3)
c
      call xc_vwn_num(n, tol_rho, fac, prho, p10Amat,
     &     p10Amat2, tAmat3, nq, ipol, p10Ec, qwght, ldew,
     &     func, functional)
c Step -1 rho
      call numder_init_lda(func, m10Ec, m10Amat, m10Amat2, tAmat3)
c
      call xc_vwn_num(n, tol_rho, fac, mrho, m10Amat,
     &     m10Amat2, tAmat3, nq, ipol, m10Ec, qwght, ldew,
     &     func, functional)
c Step +2 rho
      call numder_init_lda(func, p20Ec, p20Amat, p20Amat2, tAmat3)
c
      call xc_vwn_num(n, tol_rho, fac, p2rho, p20Amat,
     &     p20Amat2, tAmat3, nq, ipol, p20Ec, qwght, ldew,
     &     func, functional)
c Step -2 rho
      call numder_init_lda(func, m20Ec, m20Amat, m20Amat2, tAmat3)
c
      call xc_vwn_num(n, tol_rho, fac, m2rho, m20Amat,
     &     m20Amat2, tAmat3, nq, ipol, m20Ec, qwght, ldew,
     &     func, functional)
c
c First derivative is simple because it only uses the correlation 
c energy.  Here the derivative is with respect to the total density.
      d1E_rho(1,1) = -p20Ec + 8d0*p10Ec - 8d0*m10Ec + m20Ec
      d1E_rho(1,1) = d1E_rho(1,1)/(12d0*drho(1,1))
c
c To be consistent with the analytical expression, we must calculate
c the second derivative with respect to the total electron density
c by summing the contribution from the alpha- and beta-spin electrons.
c Since we are using lower order derivatives to build the numerical 
c derivatives, the coefficients are not expected to match those of
c the analytical derivatives.
c dradra
      d2E_rho(1,2) = -p20Amat(1,1) + 8d0*p10Amat(1,1) 
     &               - 8d0*m10Amat(1,1) + m20Amat(1,1)
      d2E_rho(1,2) = d2E_rho(1,2)/(12d0*drho(1,2))
c drbdrb (this is zero for restricted calculations)
      d2E_rho(1,3) = -p20Amat(1,2) + 8d0*p10Amat(1,2) 
     &               - 8d0*m10Amat(1,2) + m20Amat(1,2)
      d2E_rho(1,3) = d2E_rho(1,3)/(12d0*drho(1,3))
c dradra + drbdrb
      d2E_rho(1,1) = d2E_rho(1,2) 
     &             + d2E_rho(1,3)
c
c Third derivatives require the mixed partial derivative with respect
c to alpha- and beta-spin densities.
c dradradra
      d3E_rho(1,2) = -p20Amat2(1,1) + 8d0*p10Amat2(1,1) 
     &               - 8d0*m10Amat2(1,1) + m20Amat2(1,1)
      d3E_rho(1,2) = d3E_rho(1,2)/(12d0*drho(1,2))
c dradrbdra 
      d3E_rho(1,3) = -p20Amat2(1,2) + 8d0*p10Amat2(1,2)
     &               - 8d0*m10Amat2(1,2) + m20Amat2(1,2)
      d3E_rho(1,3) = d3E_rho(1,3)/(12d0*drho(1,2))
c drbdrbdrb
      d3E_rho(1,4) = -p20Amat2(1,3) + 8d0*p10Amat2(1,3)
     &               - 8d0*m10Amat2(1,3) + m20Amat2(1,3)
      d3E_rho(1,4) = d3E_rho(1,4)/(12d0*drho(1,3))   
c Since the calculation is submitted like a closed shell system, ra=rb,
c we can double the contribution from the mixed partial derivative.
      d3E_rho(1,1) = d3E_rho(1,2) 
     &             + 2.0d0*d3E_rho(1,3)
     &             + d3E_rho(1,4)

      end
