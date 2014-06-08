      subroutine bnl_wrapper(Amat, Amat2, Amat3, 
     1  d1E_rho, d2E_rho, d3E_rho, rho, gamma)

      implicit none

      integer n,ipol,nq
      logical lfac,nlfac,ldew
      double precision tol_rho,fac
      double precision rho(1,3),delrho(1,3,2)
      double precision Amat(1,2),Amat2(1,3),Ex
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
      double precision p10Ex,m10Ex,p20Ex,m20Ex,d1E_rho(1,2),
     &                 d2E_rho(1,3)
c
      double precision dgamma,pgamma,mgamma,p2gamma,m2gamma
c
      double precision d3E_rho(1,4)

c --------------------------------------------------------------------
c Subroutine bnl_wrapper
c
c Purpose: Evaluate the derivatives of the BNL exchange functional
c --------------------------------------------------------------------

c Initialize variables

      call var_init(n, ipol, nq, lfac, nlfac, ldew, tol_rho,
     &  fac, rho, delrho, qwght, gamma, cam_omega, prho, mrho,
     &  p2rho, m2rho, drho, dgamma, pgamma, mgamma, p2gamma, m2gamma) 

c ----------------------
c Analytical Derivatives
c ----------------------
c     Calculate the analytical derivatives of the exchange energy
c     functional of BNL with respect to the density and gradient of
c     the density.
c
      call xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, rho, 
     1     Amat, Amat2, Amat3, nq, ipol, Ex, qwght, ldew, func)

c ---------------------
c Numerical Derivatives
c ---------------------
c     Evaluate the numerical derivatives of the exchange energy 
c     functional of BNL with respect to the density. 
c    
c     Find the first, second, and third derivatives of E with respect to rho
c
c Step +1 rho
      call numder_init_lda(func, p10Ex, p10Amat, p10Amat2, tAmat3)
c
      call xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, prho, 
     1     p10Amat, p10Amat2, tAmat3, nq, ipol, p10Ex, qwght, ldew,
     2     func)
c Step -1 rho
      call numder_init_lda(func, m10Ex, m10Amat, m10Amat2, tAmat3)
c
      call xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, mrho, 
     1     m10Amat, m10Amat2, tAmat3, nq, ipol, m10Ex, qwght, ldew,
     2     func)
c Step +2 rho
      call numder_init_lda(func, p20Ex, p20Amat, p20Amat2, tAmat3)
c
      call xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, p2rho, 
     1     p20Amat, p20Amat2, tAmat3, nq, ipol, p20Ex, qwght, ldew,
     2     func)
c Step -2 rho
      call numder_init_lda(func, m20Ex, m20Amat, m20Amat2, tAmat3)
c
      call xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, m2rho, 
     1     m20Amat, m20Amat2, tAmat3, nq, ipol, m20Ex, qwght, ldew,
     2     func)
c
      d1E_rho(1,1) = -p20Ex + 8d0*p10Ex - 8d0*m10Ex + m20Ex
      d1E_rho(1,1) = d1E_rho(1,1)/(12d0*drho(1,1))
c
      d2E_rho(1,1) = -p20Amat(1,1) + 8d0*p10Amat(1,1)
     &               - 8d0*m10Amat(1,1) + m20Amat(1,1)
      d2E_rho(1,1) = d2E_rho(1,1)/(12d0*drho(1,1))
c
      d3E_rho(1,1) = -p20Amat2(1,1) + 8d0*p10Amat2(1,1)
     &               - 8d0*m10Amat2(1,1) + m20Amat2(1,1)
      d3E_rho(1,1) = d3E_rho(1,1)/(12d0*drho(1,1))

      end
