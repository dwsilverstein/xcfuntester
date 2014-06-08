      subroutine xpbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11, 
     &  Amat3, Cmat3, Cmat21, Cmat12,
     &  d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &  d3E_rho, d3E_gamma, d21E_rho_rho_gamma, 
     &  d12E_rho_gamma_gamma, rho, gamma, whichf)

      implicit none

      integer n,ipol,nq
      logical lfac,nlfac,ldew
      character*4 whichf
      double precision tol_rho,fac
      double precision rho(1,3),delrho(1,3,2)
      double precision Amat(1,2),Amat2(1,3),Cmat(1,3),Cmat2(1,6),Ex
      double precision qwght(1),func(1)
      double precision Cmat11(1,6)
      double precision Amat3(1,4),Cmat21(1,9),Cmat12(1,12),Cmat3(1,10)
      double precision tAmat3(1,4),tCmat21(1,9),tCmat12(1,12),
     &                 tCmat3(1,10)
c
      double precision gamma
      double precision cam_omega
c
      double precision prho(1,3),mrho(1,3),p2rho(1,3),m2rho(1,3),
     &  drho(1,3)
c
      double precision p10Amat(1,2),p10Amat2(1,3),p10Cmat(1,3),
     &                 p10Cmat2(1,6),p10Cmat11(1,6)
      double precision m10Amat(1,2),m10Amat2(1,3),m10Cmat(1,3),
     &                 m10Cmat2(1,6),m10Cmat11(1,6)
      double precision p20Amat(1,2),p20Amat2(1,3),p20Cmat(1,3),
     &                 p20Cmat2(1,6),p20Cmat11(1,6)
      double precision m20Amat(1,2),m20Amat2(1,3),m20Cmat(1,3),
     &                 m20Cmat2(1,6),m20Cmat11(1,6)
      double precision p10Ex,m10Ex,p20Ex,m20Ex,d1E_rho(1,2),
     &                 d2E_rho(1,3)
c
      double precision p01Amat(1,2),p01Amat2(1,3),p01Cmat(1,3),
     &                 p01Cmat2(1,6),p01Cmat11(1,6)
      double precision m01Amat(1,2),m01Amat2(1,3),m01Cmat(1,3),
     &                 m01Cmat2(1,6),m01Cmat11(1,6)
      double precision p02Amat(1,2),p02Amat2(1,3),p02Cmat(1,3),
     &                 p02Cmat2(1,6),p02Cmat11(1,6)
      double precision m02Amat(1,2),m02Amat2(1,3),m02Cmat(1,3),
     &                 m02Cmat2(1,6),m02Cmat11(1,6)
      double precision p01Ex,m01Ex,p02Ex,m02Ex,d1E_gamma(1,3),
     &                 d2E_gamma(1,6)
c
      double precision dgamma,pgamma,mgamma,p2gamma,m2gamma
c
      double precision d11E_rho_gamma(1,6)
c
      double precision d3E_rho(1,4),d3E_gamma(1,10),
     &                 d21E_rho_rho_gamma(1,9)
      double precision d12E_rho_gamma_gamma(1,12)

c --------------------------------------------------------------------
c Subroutine xpbe96_wrapper
c
c Purpose: Evaluate the derivatives of the PBE 96 exchange functional
c --------------------------------------------------------------------

c Initialize variables

      call var_init(n, ipol, nq, lfac, nlfac, ldew, tol_rho,
     &  fac, rho, delrho, qwght, gamma, cam_omega, prho, mrho,
     &  p2rho, m2rho, drho, dgamma, pgamma, mgamma, p2gamma, m2gamma) 

c ----------------------
c Analytical Derivatives
c ----------------------
c     Calculate the analytical derivatives of the exchange energy
c     functional of wPBE with respect to the density and gradient of
c     the density.
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, gamma,
     &     Amat, Amat2, Amat3, Cmat, Cmat11, Cmat2, Cmat21, Cmat12,
     &     Cmat3, nq, ipol, Ex, qwght, ldew, func)

c ---------------------
c Numerical Derivatives
c ---------------------
c     Evaluate the numerical derivatives of the exchange energy 
c     functional of PBE 96 with respect to the density and gradient of 
c     the density.
c    
c     Find the first, second, and third derivatives of E with respect to rho
c
c Step +1 rho
      call numder_initialize(func, p10Ex, p10Amat, p10Cmat,
     1  p10Amat2, p10Cmat2, p10Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, prho, gamma,
     &     p10Amat, p10Amat2, tAmat3, p10Cmat, p10Cmat11, p10Cmat2, 
     &     tCmat21, tCmat12, tCmat3, nq, ipol, p10Ex, qwght, ldew, func)
c Step -1 rho
      call numder_initialize(func, m10Ex, m10Amat, m10Cmat,
     1  m10Amat2, m10Cmat2, m10Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, mrho, gamma,
     &     m10Amat, m10Amat2, tAmat3, m10Cmat, m10Cmat11, m10Cmat2, 
     &     tCmat21, tCmat12, tCmat3, nq, ipol, m10Ex, qwght, ldew, func)
c Step +2 rho
      call numder_initialize(func, p20Ex, p20Amat, p20Cmat,
     1  p20Amat2, p20Cmat2, p20Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, p2rho, gamma,
     &     p20Amat, p20Amat2, tAmat3, p20Cmat, p20Cmat11, p20Cmat2, 
     &     tCmat21, tCmat12, tCmat3, nq, ipol, p20Ex, qwght, ldew, func)
c Step -2 rho
      call numder_initialize(func, m20Ex, m20Amat, m20Cmat,
     1  m20Amat2, m20Cmat2, m20Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, m2rho, gamma,
     &     m20Amat, m20Amat2, tAmat3, m20Cmat, m20Cmat11, m20Cmat2, 
     &     tCmat21, tCmat12, tCmat3, nq, ipol, m20Ex, qwght, ldew, func)
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
c**********************************************************************
c This part is for matching the weird PBE 96 factors.
c      d2E_rho(1,1) = d2E_rho(1,1)*2.0d0
c      d3E_rho(1,1) = d3E_rho(1,1)*4.0d0*0.50d0
c**********************************************************************
c
c     Find the first, second, and third derivatives of E with respect to 
c     gamma (the density gradient invariant)
c
c Step +1 gamma
      call numder_initialize(func, p01Ex, p01Amat, p01Cmat,
     1  p01Amat2, p01Cmat2, p01Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, pgamma,
     &     p01Amat, p01Amat2, tAmat3, p01Cmat, p01Cmat11, p01Cmat2,
     &     tCmat21, tCmat12, tCmat3, nq, ipol, p01Ex, qwght, ldew, func)
c Step -1 gamma
      call numder_initialize(func, m01Ex, m01Amat, m01Cmat,
     1  m01Amat2, m01Cmat2, m01Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, mgamma,
     &     m01Amat, m01Amat2, tAmat3, m01Cmat, m01Cmat11, m01Cmat2,
     &     tCmat21, tCmat12, tCmat3, nq, ipol, m01Ex, qwght, ldew, func)
c Step +2 gamma
      call numder_initialize(func, p02Ex, p02Amat, p02Cmat,
     1  p02Amat2, p02Cmat2, p02Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, p2gamma,
     &     p02Amat, p02Amat2, tAmat3, p02Cmat, p02Cmat11, p02Cmat2,
     &     tCmat21, tCmat12, tCmat3, nq, ipol, p02Ex, qwght, ldew, func)
c Step -2 gamma
      call numder_initialize(func, m02Ex, m02Amat, m02Cmat,
     1  m02Amat2, m02Cmat2, m02Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_xpbe96(whichf, n,
     &     tol_rho, fac, lfac, nlfac, rho, m2gamma,
     &     m02Amat, m02Amat2, tAmat3, m02Cmat, m02Cmat11, m02Cmat2,
     &     tCmat21, tCmat12, tCmat3, nq, ipol, m02Ex, qwght, ldew, func)
c
      d1E_gamma(1,1) = -p02Ex + 8d0*p01Ex - 8d0*m01Ex + m02Ex
      d1E_gamma(1,1) = d1E_gamma(1,1)/(12d0*dgamma)
c
      d2E_gamma(1,1) = -p02Cmat(1,1) + 8d0*p01Cmat(1,1) 
     &                 - 8d0*m01Cmat(1,1) + m02Cmat(1,1)
      d2E_gamma(1,1) = d2E_gamma(1,1)/(12d0*dgamma)
c
      d3E_gamma(1,1) = -p02Cmat2(1,1) + 8d0*p01Cmat2(1,1) 
     &                 - 8d0*m01Cmat2(1,1) + m02Cmat2(1,1)
      d3E_gamma(1,1) = d3E_gamma(1,1)/(12d0*dgamma)
c
c     Find the mixed second functional derivative of E with respect to
c     gamma and rho (take the derivative with respect to gamma first,
c     even though order is not important).
c
      d11E_rho_gamma(1,1) = -p20Cmat(1,1) + 8d0*p10Cmat(1,1) 
     &                      - 8d0*m10Cmat(1,1) + m20Cmat(1,1)
      d11E_rho_gamma(1,1) = d11E_rho_gamma(1,1)/(12d0*drho(1,1))

c
c     Find the mixed third order functional derivative of E: 
c     dgamma drho^2. 
c
      d21E_rho_rho_gamma(1,1) = -p20Cmat11(1,1) + 8d0*p10Cmat11(1,1) 
     &                          - 8d0*m10Cmat11(1,1) + m20Cmat11(1,1)
      d21E_rho_rho_gamma(1,1) = d21E_rho_rho_gamma(1,1)/(12d0*drho(1,1))

c
c     Find the mixed third order functional derivative of E: 
c     dgamma^2 drho.  There are 2 ways to calculate this.
c
      d12E_rho_gamma_gamma(1,1) = -p20Cmat2(1,1) + 8d0*p10Cmat2(1,1) 
     &                            - 8d0*m10Cmat2(1,1) + m20Cmat2(1,1)
      d12E_rho_gamma_gamma(1,1) = d12E_rho_gamma_gamma(1,1)/
     &  (12d0*drho(1,1))

      end
