      subroutine var_init(n, ipol, nq, lfac, nlfac, ldew, tol_rho,
     &  fac, rho, delrho, qwght, gamma, cam_omega, prho, mrho, 
     &  p2rho, m2rho, drho, dgamma, pgamma, mgamma, p2gamma, m2gamma)
c     &  p3rho, m3rho, p3gamma, m3gamma)

      implicit none

      integer n,ipol,nq,i
      logical lfac,nlfac,ldew
      double precision tol_rho,fac
      double precision rho(1,3),delrho(1,3,2)
      double precision qwght(1)
c
      double precision gamma
      double precision cam_omega
c
c Variables needed for 5-point derivatives
      double precision prho(1,3),mrho(1,3),p2rho(1,3),m2rho(1,3),
     &  drho(1,3)
c Variables needed for 7-point derivatives
c      double precision p3rho(1,3), m3rho(1,3)
c
c Variables needed for 5-point derivatives
      double precision dgamma,pgamma,mgamma,p2gamma,m2gamma
c Variables needed for 7-point derivatives
c      double precision p3gamma, m3gamma

c --------------------------------------------------------------------
c Subroutine var_init
c
c Purpose: Initialize common variables to the exchange and correlation
c functionals.
c --------------------------------------------------------------------

c -----------------
c Initialize Values
c -----------------
c
c     Define the iteration variables
c
      n = 1
      nq = 1
      ipol = 1
      lfac = .false.
      nlfac = .false.
      ldew = .false.
c
c      tol_rho = 0.01d0
      tol_rho = 1.0d-10
      qwght = 1.0d0
      fac = 1.0d0
      cam_omega = 0.5d0
c
c     Define values for the electron density and electron density
c     gradient.  
c
c For rho, the 3 dimensions are: 
c 1 = total charge density 
c 2 = alpha charge density
c 3 = beta spin charge density 
c --------------------------------------------------------------
c Values are taken from an N2 calculation with PBEx + LYP/STO-3G
c --------------------------------------------------------------
c These are numbered based on size, where 1 is for the largest density.
c **Note: I have found the numerical derivatives to be inaccurate for 
c values of the density smaller than 8e-7 (number 12) 
c 1. Extremely large values (1e2)
c      rho(n,1) = 123.32015832293358d0
c      rho(n,2) =  61.66007916146679d0
c      rho(n,3) =  61.66007916146679d0
c 2. Very large values (1e1) (Small a_sigma)
      rho(n,1) = 26.426613856209315d0
      rho(n,2) = 13.213306928104657d0
      rho(n,3) = 13.213306928104657d0
c 3. Large values (1e0)
c      rho(n,1) = 1.5530607745547451d0
c      rho(n,2) = 0.77653038727737256d0
c      rho(n,3) = 0.77653038727737256d0
c 4. Medium-large values (8e-1)
c      rho(n,1) = 0.78189665830502308d0
c      rho(n,2) = 0.39094832915251154d0
c      rho(n,3) = 0.39094832915251154d0
c 5. Medium values (3e-1)
c      rho(n,1) = 0.26112058424563056d0
c      rho(n,2) = 0.13056029212281528d0
c      rho(n,3) = 0.13056029212281528d0
c 6. Medium-small values (5e-2)
c      rho(n,1) = 0.0534721665082183856d0
c      rho(n,2) = 0.026736083254109193d0
c      rho(n,3) = 0.026736083254109193d0
c 7. Small values (3e-2)
c      rho(n,1) = 0.0260143192338825671d0
c      rho(n,2) = 0.013007159616941284d0
c      rho(n,3) = 0.013007159616941284d0
c 8. Small values (2e-3) (Intermediate a_sigma)
c      rho(n,1) = 0.00245485286886596311d0
c      rho(n,2) = 0.0012274264344329816d0
c      rho(n,3) = 0.0012274264344329816d0
c 9. Small values (1e-4)
c      rho(n,1) = 1.41619012840620169d-004 
c      rho(n,2) = 7.0809506420310085d-05
c      rho(n,3) = 7.0809506420310085d-05
c 10. Small values (3e-5)
c      rho(n,1) = 3.18273891584200803d-005 
c      rho(n,2) = 1.591369457921004d-05
c      rho(n,3) = 1.591369457921004d-05 
c 11. Small values (6e-6)
c      rho(n,1) = 6.42716971745045573d-006 
c      rho(n,2) = 3.2135848587252279d-06
c      rho(n,3) = 3.2135848587252279d-06
c 12. Small values (8e-7) (Large a_sigma)
c      rho(n,1) = 8.60549076595114122d-007 
c      rho(n,2) = 4.3027453829755706d-07
c      rho(n,3) = 4.3027453829755706d-07
c 13. Small values (8e-8)
c      rho(n,1) = 8.02056723730010681d-008 
c      rho(n,2) = 4.0102836186500534d-08 
c      rho(n,3) = 4.0102836186500534d-08
c 14. Small values (4e-9)
c      rho(n,1) = 4.31552350124199254d-009 
c      rho(n,2) = 2.1577617506209963d-09
c      rho(n,3) = 2.1577617506209963d-09
c 15. Smallest values (1e-10)
c      rho(n,1) = 1.14162379685153242d-010 
c      rho(n,2) = 5.7081189842576621d-11 
c      rho(n,3) = 5.7081189842576621d-11 
c For delrho, the 3 dimensions come from having 3 Cartesian 
c coordinates.  The array dimension 2 comes from 2 spins (alpha and
c beta).
c Large values
c      delrho(1,1,1) = 1.5d0
c      delrho(1,2,1) = 1.5d0
c      delrho(1,3,1) = 1.5d0
c      delrho(1,1,2) = 1.5d0
c      delrho(1,2,2) = 1.5d0
c      delrho(1,3,2) = 1.5d0
c Medium-large values
c      delrho(1,1,1) = 1.0d0
c      delrho(1,2,1) = 1.0d0
c      delrho(1,3,1) = 1.0d0
c      delrho(1,1,2) = 1.0d0
c      delrho(1,2,2) = 1.0d0
c      delrho(1,3,2) = 1.0d0
c Medium values
c      delrho(1,1,1) = 0.5d0
c      delrho(1,2,1) = 0.5d0
c      delrho(1,3,1) = 0.5d0
c      delrho(1,1,2) = 0.5d0
c      delrho(1,2,2) = 0.5d0
c      delrho(1,3,2) = 0.5d0
c Medium-small values
      delrho(1,1,1) = 0.1d0
      delrho(1,2,1) = 0.1d0
      delrho(1,3,1) = 0.1d0
      delrho(1,1,2) = 0.1d0
      delrho(1,2,2) = 0.1d0
      delrho(1,3,2) = 0.1d0
c Small values
c      delrho(1,1,1) = 0.01d0
c      delrho(1,2,1) = 0.01d0
c      delrho(1,3,1) = 0.01d0
c      delrho(1,1,2) = 0.01d0
c      delrho(1,2,2) = 0.01d0
c      delrho(1,3,2) = 0.01d0

c
c     Calculate the electron density gradient invariant, gamma.
c
c      gamma = delrho(1,1,1)*delrho(1,1,1) +
c     &        delrho(1,2,1)*delrho(1,2,1) +
c     &        delrho(1,3,1)*delrho(1,3,1)
c -----------------------------------------------------------------
c Real values from an N2/STO-3G calculation with PBEx + LYP.  These
c numbers are consistent with the density values above.
c -----------------------------------------------------------------
c Note that a large value of the density does not necessarily mean
c the density gradient will be large (although the data here might
c mislead you to think otherwise).
c Values are numbered from largest to smallest (1 is for the largest
c density), and each corresponds to the density given above.
c 1. Extremely large values (2e5)
c      gamma = 220378.38320076052 
c 2. Very large values (1e5) (Small a_sigma)
c      gamma = 138225.07971908242d0 
c 3. Large values (1e2)
c      gamma = 225.08843104207619d0
c 4. Medium-large values (7e0)
c      gamma = 6.9492129903820894d0
c 5. Medium values (2e-1)
c      gamma = 0.30531953669703915d0 
c 6. Medium-small values (1e-2)
c      gamma = 0.0233058215868636870d0
c 7. Small values (6e-3)
c      gamma = 5.78841884664339779d-003 
c 8. Small values (6e-5) (Intermediate a_sigma)
      gamma = 5.70962076927598785d-005 
c 9. Small values (3e-7)
c      gamma = 2.99784476069366267d-007 
c 10. Small values (2e-8)
c      gamma = 1.82801349440045699d-008
c 11. Small values (9e-10)
c      gamma = 9.02504183194566195d-010
c 12. Small values(2e-11) (Large a_sigma)
c      gamma = 1.96350879161959226d-011
c 13. Small values (2e-13)
c      gamma = 2.06998427591892019d-013 
c 14. Small values (7e-16)
c      gamma = 7.28717252544055553d-016 
c 15. Smallest values (6e-19)
c      gamma = 6.20927458004201392d-019 

c
c     Find the step sizes for the numerical derivatives.
c
      do i = 1, 3 
        drho(n,i) = 1d-3*rho(n,i)
        prho(n,i) = rho(n,i) + drho(n,i)
        mrho(n,i) = rho(n,i) - drho(n,i)
        p2rho(n,i) = rho(n,i) + 2d0*drho(n,i)
        m2rho(n,i) = rho(n,i) - 2d0*drho(n,i)
c        p3rho(n,i) = rho(n,i) + 3d0*drho(n,i)
c        m3rho(n,i) = rho(n,i) - 3d0*drho(n,i)
      enddo
c
      dgamma = 1d-2*gamma
      pgamma = gamma + dgamma
      mgamma = gamma - dgamma
      p2gamma = gamma + 2d0*dgamma
      m2gamma = gamma - 2d0*dgamma
c      p3gamma = gamma + 3d0*dgamma
c      m3gamma = gamma - 3d0*dgamma

      end
