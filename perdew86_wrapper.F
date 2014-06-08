      subroutine perdew86_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11, 
     &  Amat3, Cmat3, Cmat21, Cmat12,
     &  d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &  d3E_rho, d3E_gamma, d21E_rho_rho_gamma, 
     &  d12E_rho_gamma_gamma, rho, gamma)

      implicit none

      integer n,ipol,nq,i
      logical lfac,nlfac,ldew
      double precision tol_rho,fac
      double precision rho(1,3),delrho(1,3,2)
      double precision Amat(1,2),Amat2(1,3),Cmat(1,3),Cmat2(1,6),Ec
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
      double precision p10Ec,m10Ec,p20Ec,m20Ec
      double precision p01Ec,m01Ec,p02Ec,m02Ec
      double precision p10Amat(1,2),p10Amat2(1,3),p10Cmat(1,3),
     &                 p10Cmat2(1,6),p10Cmat11(1,6)
      double precision m10Amat(1,2),m10Amat2(1,3),m10Cmat(1,3),
     &                 m10Cmat2(1,6),m10Cmat11(1,6)
      double precision p20Amat(1,2),p20Amat2(1,3),p20Cmat(1,3),
     &                 p20Cmat2(1,6),p20Cmat11(1,6)
      double precision m20Amat(1,2),m20Amat2(1,3),m20Cmat(1,3),
     &                 m20Cmat2(1,6),m20Cmat11(1,6)
      double precision d1E_rho(1,2),d2E_rho(1,3)
c
      double precision p01Amat(1,2),p01Amat2(1,3),p01Cmat(1,3),
     &                 p01Cmat2(1,6),p01Cmat11(1,6)
      double precision m01Amat(1,2),m01Amat2(1,3),m01Cmat(1,3),
     &                 m01Cmat2(1,6),m01Cmat11(1,6)
      double precision p02Amat(1,2),p02Amat2(1,3),p02Cmat(1,3),
     &                 p02Cmat2(1,6),p02Cmat11(1,6)
      double precision m02Amat(1,2),m02Amat2(1,3),m02Cmat(1,3),
     &                 m02Cmat2(1,6),m02Cmat11(1,6)
      double precision d1E_gamma(1,3), d2E_gamma(1,6)
c
      double precision dgamma,pgamma,mgamma,p2gamma,m2gamma
c
      double precision d11E_rho_gamma(1,6)
c
      double precision d3E_rho(1,4),d3E_gamma(1,10),
     &                 d21E_rho_rho_gamma(1,9)
      double precision d12E_rho_gamma_gamma(1,12)
c Testing 7-point differentiation
c      double precision p3rho(1,3), m3rho(1,3)
c      double precision p30Ec,m30Ec
c      double precision p30Amat(1,2),p30Amat2(1,3),p30Cmat(1,3),
c     &                 p30Cmat2(1,6),p30Cmat11(1,6)
c      double precision m30Amat(1,2),m30Amat2(1,3),m30Cmat(1,3),
c     &                 m30Cmat2(1,6),m30Cmat11(1,6)
cc
c      double precision p3gamma, m3gamma
c      double precision p03Ec,m03Ec
c      double precision p03Amat(1,2),p03Amat2(1,3),p03Cmat(1,3),
c     &                 p03Cmat2(1,6),p03Cmat11(1,6)
c      double precision m03Amat(1,2),m03Amat2(1,3),m03Cmat(1,3),
c     &                 m03Cmat2(1,6),m03Cmat11(1,6)

c --------------------------------------------------------------------
c Subroutine perdew86_wrapper
c
c Purpose: Evaluate the derivatives of the Perdew 86 correlation functional
c --------------------------------------------------------------------

c Initialize variables

      call var_init(n, ipol, nq, lfac, nlfac, ldew, tol_rho,
     &  fac, rho, delrho, qwght, gamma, cam_omega, prho, mrho,
     &  p2rho, m2rho, drho, dgamma, pgamma, mgamma, p2gamma, m2gamma)

c ----------------------
c Analytical Derivatives
c ----------------------
c     Calculate the analytical derivatives of the exchange energy
c     functional of Perdew 86 with respect to the density and gradient of
c     the density.
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  rho, gamma, Amat, Amat2, Amat3,
     &  Cmat, Cmat2, Cmat11, Cmat3, Cmat21, Cmat12,
     &  nq, ipol, Ec, qwght, ldew, func)

c To correctly build the analytical derivative, we need to account
c for the total spin density.

      Amat2(1,3) = Amat2(1,1)

      Amat2(1,1) = Amat2(1,1)
     &           + 2.0d0*Amat2(1,2)
     &           + Amat2(1,3)
c      Amat2(1,1) = Amat2(1,1)
c     &           + Amat2(1,2)

      Amat3(1,4) = Amat3(1,1)

      Amat3(1,1) = Amat3(1,1)
     &           + 3.0d0*Amat3(1,2)
     &           + 3.0d0*Amat3(1,3)
     &           + Amat3(1,4)
c      Amat3(1,1) = Amat3(1,1)
c     &           + 2.0d0*Amat3(1,2)
c     &           + Amat3(1,3)

      Cmat21(1,9) = Cmat21(1,1)

      Cmat21(1,1) = Cmat21(1,1)
     &            + 2.0d0*Cmat21(1,4)
     &            + Cmat21(1,9)

c ---------------------
c Numerical Derivatives
c ---------------------
c     Evaluate the numerical derivatives of the correlation energy 
c     functional of Perdew 86 with respect to the density and gradient of 
c     the density.
c    
c     Find the first, second, and third derivatives of E with respect to rho
c
c Step +1 rho
      call numder_initialize(func, p10Ec, p10Amat, p10Cmat,
     1  p10Amat2, p10Cmat2, p10Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  prho, gamma, p10Amat, p10Amat2, tAmat3,
     &  p10Cmat, p10Cmat2, p10Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, p10Ec, qwght, ldew, func)
c Step -1 rho
      call numder_initialize(func, m10Ec, m10Amat, m10Cmat,
     1  m10Amat2, m10Cmat2, m10Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  mrho, gamma, m10Amat, m10Amat2, tAmat3,
     &  m10Cmat, m10Cmat2, m10Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, m10Ec, qwght, ldew, func)
c Step +2 rho
      call numder_initialize(func, p20Ec, p20Amat, p20Cmat,
     1  p20Amat2, p20Cmat2, p20Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  p2rho, gamma, p20Amat, p20Amat2, tAmat3,
     &  p20Cmat, p20Cmat2, p20Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, p20Ec, qwght, ldew, func)
c Step -2 rho
      call numder_initialize(func, m20Ec, m20Amat, m20Cmat,
     1  m20Amat2, m20Cmat2, m20Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  m2rho, gamma, m20Amat, m20Amat2, tAmat3,
     &  m20Cmat, m20Cmat2, m20Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, m20Ec, qwght, ldew, func)
c
c 5-point derivative
      d1E_rho(1,1) = -p20Ec + 8d0*p10Ec - 8d0*m10Ec + m20Ec
      d1E_rho(1,1) = d1E_rho(1,1)/(12d0*drho(1,1))
c
c To be consistent with the analytical expression, we must calculate
c the second derivative with respect to the total electron density
c by summing the contribution from the alpha- and beta-spin electrons.
c dradra, 5-point
      d2E_rho(1,2) = -p20Amat(1,1) + 8d0*p10Amat(1,1)
     &               - 8d0*m10Amat(1,1) + m20Amat(1,1)
      d2E_rho(1,2) = d2E_rho(1,2)/(12d0*drho(1,2))
c drbdrb, 5-point (this is zero for restricted calculations)
c For restricted
      d2E_rho(1,3) = d2E_rho(1,2)
c      d2E_rho(1,3) = -p20Amat(1,2) + 8d0*p10Amat(1,2)
c     &               - 8d0*m10Amat(1,2) + m20Amat(1,2)
c      d2E_rho(1,3) = d2E_rho(1,3)/(12d0*drho(1,3))
c dradra + drbdrb
      d2E_rho(1,1) = d2E_rho(1,2)
     &             + d2E_rho(1,3)
c
c Third derivatives require the mixed partial derivative with respect
c to alpha- and beta-spin densities.
c dradradra, 5-point
      d3E_rho(1,2) = -p20Amat2(1,1) + 8d0*p10Amat2(1,1)
     &               - 8d0*m10Amat2(1,1) + m20Amat2(1,1)
      d3E_rho(1,2) = d3E_rho(1,2)/(12d0*drho(1,2))
c dradrbdra, 5-point 
      d3E_rho(1,3) = -p20Amat2(1,2) + 8d0*p10Amat2(1,2)
     &               - 8d0*m10Amat2(1,2) + m20Amat2(1,2)
      d3E_rho(1,3) = d3E_rho(1,3)/(12d0*drho(1,2))
c drbdrbdrb, 5-point
c For a restricted calculation
      d3E_rho(1,4) = d3E_rho(1,2)
c      d3E_rho(1,4) = -p20Amat2(1,3) + 8d0*p10Amat2(1,3)
c     &               - 8d0*m10Amat2(1,3) + m20Amat2(1,3)
c      d3E_rho(1,4) = d3E_rho(1,4)/(12d0*drho(1,3))
c Since the calculation is submitted like a closed shell system, ra=rb,
c we can double the contribution from the mixed partial derivative.
      d3E_rho(1,1) = d3E_rho(1,2)
     &             + 2.0d0*d3E_rho(1,3)
     &             + d3E_rho(1,4)

c
c     Find the first, second, and third derivatives of E with respect to 
c     gamma (the density gradient invariant)
c
c Step +1 gamma
      call numder_initialize(func, p01Ec, p01Amat, p01Cmat,
     1  p01Amat2, p01Cmat2, p01Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  rho, pgamma, p01Amat, p01Amat2, tAmat3,
     &  p01Cmat, p01Cmat2, p01Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, p01Ec, qwght, ldew, func)
c Step -1 gamma
      call numder_initialize(func, m01Ec, m01Amat, m01Cmat,
     1  m01Amat2, m01Cmat2, m01Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  rho, mgamma, m01Amat, m01Amat2, tAmat3,
     &  m01Cmat, m01Cmat2, m01Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, m01Ec, qwght, ldew, func)
c Step +2 gamma
      call numder_initialize(func, p02Ec, p02Amat, p02Cmat,
     1  p02Amat2, p02Cmat2, p02Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  rho, p2gamma, p02Amat, p02Amat2, tAmat3,
     &  p02Cmat, p02Cmat2, p02Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, p02Ec, qwght, ldew, func)
c Step -2 gamma
      call numder_initialize(func, m02Ec, m02Amat, m02Cmat,
     1  m02Amat2, m02Cmat2, m02Cmat11, tAmat3, tCmat3, tCmat21,
     2  tCmat12)
c
      call xc_perdew86(n, tol_rho, fac, lfac, nlfac,
     &  rho, m2gamma, m02Amat, m02Amat2, tAmat3,
     &  m02Cmat, m02Cmat2, m02Cmat11, tCmat3, tCmat21, tCmat12,
     &  nq, ipol, m02Ec, qwght, ldew, func)
c
c 5-point derivative 
      d1E_gamma(1,1) = -p02Ec + 8d0*p01Ec - 8d0*m01Ec + m02Ec
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
c dradradg
      d21E_rho_rho_gamma(1,2) = -p20Cmat11(1,1) + 8d0*p10Cmat11(1,1)
     &                        - 8d0*m10Cmat11(1,1) + m20Cmat11(1,1)
      d21E_rho_rho_gamma(1,2) = d21E_rho_rho_gamma(1,2)/
     &                          (12d0*drho(1,2))
c drbdrbdg 
      d21E_rho_rho_gamma(1,4) = -p20Cmat11(1,2) + 8d0*p10Cmat11(1,2)
     &                        - 8d0*m10Cmat11(1,2) + m20Cmat11(1,2)
      d21E_rho_rho_gamma(1,4) = d21E_rho_rho_gamma(1,4)/
     &                          (12d0*drho(1,3))
c dradradg + drbdrbdg
      d21E_rho_rho_gamma(1,1) = d21E_rho_rho_gamma(1,2)
     &                        + d21E_rho_rho_gamma(1,4)
c
c     Find the mixed third order functional derivative of E: 
c     dgamma^2 drho.  There are 2 ways to calculate this.
c
      d12E_rho_gamma_gamma(1,1) = -p20Cmat2(1,1) + 8d0*p10Cmat2(1,1) 
     &                            - 8d0*m10Cmat2(1,1) + m20Cmat2(1,1)
      d12E_rho_gamma_gamma(1,1) = d12E_rho_gamma_gamma(1,1)/
     &  (12d0*drho(1,1))

      end
