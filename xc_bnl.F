c
c     -----------------------------------------------------------------------
c     Uniform electron gas exchange functional for the erfc(r)/r interaction
c     as implemented in the following paper:
c     "A well-tempered density functional theory of electrons in molecules"
c     Ester Livshits & Roi Baer, Phys. Chem. Chem. Phys., 9, 2932 (2007)
c     The other relevant publication is: 
c     R. Baer, D. Neuhauser, Phys. Rev. Lett., 94, 043002 (2005)
c     -----------------------------------------------------------------------
c
      subroutine xc_bnl(cam_omega, n, tol_rho, fac, lfac, nlfac, rho, 
     1              Amat, Amat2, Amat3, nq, ipol, Ex, qwght, ldew, 
     2              func)
c
      implicit none
c
c#include "errquit.fh"
c#include "stdio.fh"
cc
cccase...start
c#include "case.fh"
cccase...end
c
      integer nq, ipol, n
      double precision fac, Ex, total
      logical ldew, lfac, nlfac
      double precision func(*) ! value of the functional [output]
      double precision tol_rho
      double precision rho(n,3) ! charge density
      double precision qwght(nq)                 ! quadrature weights
      double precision Amat(n,2)             ! partial first derivatives  
      double precision F(n),RA(n),RB(n)
      double precision rhoA, rhoB, rhoTotal, rhoA1, rhoB1
      double precision gamma
c      double precision fA, fB, fpA, fpB, fppA, fppB
      double precision fA, fB, fpA, fpB
      double precision EpsX
      double precision EpsXprime
c      double precision EpsTwoXprime
c
      double precision Amat2(n,3)       ! partial second derivatives
      double precision fppA, fppB
      double precision EpsTwoXprime
      double precision Amat3(n,4)             
      double precision fpppA, fpppB
      double precision EpsThreeXprime
c
      double precision cam_omega
c
c     -----------------------------------------------------------------------
c     Preliminaries
c     -----------------------------------------------------------------------
c
      gamma = cam_omega
c
c      do n = 1,nq
         if (ipol.eq.1) then   ! spin-restricted
            rA(n) = rho(n,1)
            rB(n) = 0.d0
         else                  ! spin-unrestricted
            rA(n) = rho(n,2)
            rB(n) = rho(n,3)
         end if
c      end do
c
c     -----------------------------------------------------------------------
c     Calculate the first and second derivatives
c     -----------------------------------------------------------------------
c
      total = 0.d0
c      do n = 1,nq
         rhoA = rA(n)
         rhoB = rB(n)
         rhoTotal  = rhoA + rhoB   ! total density at point
         if (rhoTotal.gt.tol_rho) then

            if (ipol.eq.1) then    ! spin-restricted
              rhoA1 = rhoA
              rhoB1 = rhoB
            else                   ! spin-unrestricted
              rhoA1 = rhoA*2.0d0
              rhoB1 = rhoB*2.0d0
            end if

            fA   = EpsX(rhoA1,gamma)
            fB   = EpsX(rhoB1,gamma)
            fpA  = EpsXprime(rhoA1,gamma)
            fpB  = EpsXprime(rhoB1,gamma)

            f(n) = fA * rhoA + fB * rhoB
            Amat(n,1) = Amat(n,1) + (fpA*rhoA1+fA)*fac
            if (ipol.gt.1) then
              Amat(n,2) = Amat(n,2) + (fpB*rhoB1+fB)*fac
            end if

            fppA = EpsTwoXprime(rhoA1,gamma)
            fppB = EpsTwoXprime(rhoB1,gamma)

c Daniel (4-3-13): Factor of 4 here looks similar to the implementation
c of the PBE 96 exchange functional.  It should be a factor of 2,
c however, if this comes from spin-scaling...
c            Amat2(n,1) = Amat2(n,1) + 
c     &        ((fppA*rhoA+fpA)*4)*fac
c            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 
c     &        ((fppA*rhoA+fpA)*4.0d0)*fac
c            Amat2(n,1) = Amat2(n,1) + 
c     &        ((fppA*rhoA+fpA)*4.0d0)*fac
            Amat2(n,1) = Amat2(n,1) + 
     &        ((fppA*rhoA+2.0d0*fpA))*fac
c
c           Guard against case of no beta electrons, e.g. H atom
            if (ipol.gt.1) then
             if (rho(n,3).gt.tol_rho) then
c              Amat2(n,2) = Amat2(n,2) + 
c     &          ((fppB*rhoB+fpB)*4)*fac
c              Amat2(n,3) = Amat2(n,3) + 
c     &          ((fppB*rhoB+fpB)*4.0d0)*fac
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 
c     &          ((fppB*rhoB+fpB)*4.0d0)*fac
              Amat2(n,2) = Amat2(n,2) + 
     &          ((fppB*rhoB+2.0d0*fpB))*fac
             end if
            end if
c Daniel (4-3-13): Third derivatives
            fpppA = EpsThreeXprime(rhoA1,gamma)
            fpppB = EpsThreeXprime(rhoB1,gamma)
c
c            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
c     1        + ((fpppA*rhoa + 3.0d0*fppA))*fac 
            Amat3(n,1) = Amat3(n,1)
     1        + ((fpppA*rhoA + 3.0d0*fppA))*fac 
            if (ipol.gt.1) then
             if (rho(n,3).gt.tol_rho) then
              Amat3(n,4) = Amat3(n,4)
     1          + ((fpppB*rhoB + 3.0d0*fppB))*fac 
             end if
            end if

            if (ldew) func(n) = func(n) + f(n)*fac 
            total = total + f(n)*qwght(n)
         end if
c      end do

      Ex = Ex + total*fac

      return
      end
c
c
c
c     ---------------------------------------------------------------------------------------
c     Utility functions
c     ---------------------------------------------------------------------------------------
c
c     ---------------------------------------------------------------------------------------
c     Return the value of pi
c     ---------------------------------------------------------------------------------------
c
      double precision function ValueOfPi()
c 
      implicit none
c#include "xc_params.fh"
c      
c      ValueOfPi = pi          
      ValueOfPi = 3.1415926535897932385d0 

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Evaluates the actual function
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNL(q)

      implicit none
c#include "xc_params.fh"

      double precision q,TwoSqrtPi,OneOverQ,q2,DERF
      double precision pi

      pi = 3.1415926535897932385d0

      OneOverQ = 1.0d0/q
      TwoSqrtPi = 2.0d0*dsqrt(pi) 
      q2 = q**2.0d0

      if (q .lt. 1D-15) then
         HqBNL=1.d0
         return
      end if

      if (q .lt. 0.1d0) then
         HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi-q+q*(q2-2.0d0))
         return
      end if

      HqBNL=1.0d0-q*2.0d0/3.0d0*(TwoSqrtPi*DERF(OneOverQ)-q+
     $     q*(q2-2.0d0)*(1.0d0-exp(-OneOverQ*OneOverQ)))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the local Fermi vector for the provided density
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiK(den)

      implicit none
c#include "xc_params.fh"

      double precision F13, den
      double precision pi

      pi = 3.1415926535897932385d0

      F13 = 1.0D0 / 3.0D0
      FermiK = (3.d0*pi*pi*den)**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the function EpsX at the given density value and gamma
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsX(Rho,gamma)

      implicit none
c#include "xc_params.fh"

      double precision  kF,RHO,gamma,Cs
      double precision HqBNL
      double precision FermiK
      double precision pi

      pi = 3.1415926535897932385d0

      if (RHO.le.0D0) then
         EpsX = 0.0D0
         return
      end if

      kF = FermiK(Rho)
      Cs = -3.0D0/(4.0d0*pi)
      EpsX = Cs * kF * HqBNL(gamma/kF)

      return
      end      
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the function
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNLPrime(q)

      implicit none
c#include "xc_params.fh"

      double precision q,OneOverQ,q2,q3,DERF
      double precision pi

      pi = 3.1415926535897932385d0

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0
      q3 = q**3.0d0

      if (q .lt. 0.1d0) then
        HqBNLPrime = -4.0d0/3.0d0*(dsqrt(Pi)+2.0d0*q3-3.0d0*q)
        return
      end if

      HqBNLPrime = 4.0d0/3.0d0*(q*(exp(-OneOverQ*OneOverQ)*(2.0d0*q2
     $     -1.0d0)+(3.0d0-2.0d0*q2))-dsqrt(Pi)*DERF(OneOverQ))

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of the local Fermi vector (it depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiKPrime(den)

      implicit none
c#include "xc_params.fh"
   
      double precision F23, den
      double precision pi

      pi = 3.1415926535897932385d0

      F23 = 2.0D0 / 3.0D0
      FermiKPrime = (Pi/(3.0d0*den))**F23

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of q (q=gamma/kf) (it implicitly depends on the density)
c     ---------------------------------------------------------------------------------------
c
      double precision function QPrime(gamma,kF)

      implicit none

      double precision  kF, FermiK2, gamma

      FermiK2 = kF**2.0d0
      QPrime = -gamma/FermiK2

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the first derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsXprime(Rho,gamma)

      implicit none
c#include "xc_params.fh"

      double precision Rho,gamma
      double precision Cs,kF,CsPrime

      double precision HqBNL
      double precision HqBNLPrime
      double precision QPrime
      double precision FermiK
      double precision FermiKPrime
      double precision pi

      pi = 3.1415926535897932385d0

      kF = FermiK(Rho)
      CsPrime = -3.0D0/(4.0d0*Pi)
      Cs = CsPrime*kF

      if (Rho.le.0d0) then
         EpsXprime = 0.0d0
         return
      end if

      EpsXprime = FermiKPrime(Rho)*(CsPrime*HqBNL(gamma/kF)+
     $     QPrime(gamma,kF)*HqBNLPrime(gamma/kF)*Cs)

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the main function that consititutes the functional
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNLTwoPrime(q)

      implicit none
c#include "xc_params.fh"

      double precision q,OneOverQ,q2
      double precision pi

      pi = 3.1415926535897932385d0

      OneOverQ = 1.0d0/q
      q2 = q**2.0d0

      if (q .lt. 0.1d0) then
         HqBNLTwoPrime = 4.0d0-8.0d0*q2
         return
      end if

      HqBNLTwoPrime = exp(-OneOverQ*OneOverQ)*(4.0d0+8.0d0*q2)
     $     -8.0d0*q2+4.0d0

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of the local Fermi vector
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiKTwoPrime(den)

      implicit none
c#include "xc_params.fh"

      double precision F13, den
      double precision pi

      pi = 3.1415926535897932385d0

      F13 = 1.0D0/3.0D0
      FermiKTwoPrime =  -(8.0d0*Pi**2.0d0/(243.0d0*den**5.0d0))**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of q    
c     ---------------------------------------------------------------------------------------
c
      double precision function QTwoPrime(gamma,kF)

      implicit none

      double precision gamma, kF, FermiK3

      FermiK3 = kF**3.0d0
      QTwoPrime = (2.0d0*gamma)/FermiK3

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the second derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsTwoXprime(Rho,gamma)

      implicit none
c#include "xc_params.fh"

      double precision Rho,gamma
      double precision kF,kFPrim,kFprim2,kF2prim
      double precision q,qprim,qprim2,q2prim
      double precision g,gprim,g2prim
      double precision Cs,CsPrim

      double precision FermiK
      double precision FermiKPrime
      double precision FermiKTwoPrime
      double precision QPrime
      double precision QTwoPrime
      double precision HqBNL
      double precision HqBNLPrime
      double precision HqBNLTwoPrime
      double precision pi

      pi = 3.1415926535897932385d0

      if (Rho.le.0d0) then
         EpsTwoXprime = 0.0d0
         return
      end if

      kF = FermiK(Rho)
      kFPrim = FermiKPrime(Rho)
      kFPrim2=kFPrim**2.0d0
      kF2prim = FermiKTwoPrime(Rho)
      CsPrim = -3.0d0/(4.0d0*Pi)
      Cs = CsPrim * kF
      q = gamma / kF
      qprim = QPrime(gamma,kF)
      Qprim2=qprim**2.0d0
      q2prim = QTwoPrime(gamma,kF)
      g = HqBNL(q)
      gprim = HqBNLPrime(q)
      g2prim = HqBNLTwoPrime(q)

      EpsTwoXprime = 
     $     kFPrim2*(2.0d0*CsPrim*gprim*qprim
     $     +Cs*(QPrim2*g2prim+gprim*Q2Prim))
     $     +kF2Prim*(g*CsPrim+Cs*gprim*qprim)

      return
      end
c Daniel (4-3-13): Third derivative utility functions 
c
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of the main function that consititutes the functional
c     ---------------------------------------------------------------------------------------
c
      double precision function HqBNLThreePrime(q)

      implicit none
c#include "xc_params.fh"

      double precision q,OneOverQ
      double precision q2, q3, q4
      double precision pi

      pi = 3.1415926535897932385d0

      OneOverQ = 1.0d0/q
      q2 = q*q
      q3 = q2*q
      q4 = q3*q

      if (q .lt. 0.1d0) then
         HqBNLThreePrime = -16.0d0*q
         return
      end if

      HqBNLThreePrime = 8.0d0*( exp(-OneOverQ*OneOverQ)
     1                        + 2.0d0*q2*exp(-OneOverQ*OneOverQ)
     2                        - 2.0d0*q4
     3                        + 2.0d0*q4*exp(-OneOverQ*OneOverQ) ) / q3

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of the local Fermi vector
c     ---------------------------------------------------------------------------------------
c
      double precision function FermiKThreePrime(den)

      implicit none
c#include "xc_params.fh"

      double precision F13, den
      double precision pi

      pi = 3.1415926535897932385d0

      F13 = 1.0D0/3.0D0
      FermiKThreePrime =  (10.0d0/9.0d0)*
     1                  (Pi**2.0d0/(9.0d0*den**8.0d0))**F13

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of q    
c     ---------------------------------------------------------------------------------------
c
      double precision function QThreePrime(gamma,kF)

      implicit none

      double precision gamma, kF, FermiK4

      FermiK4 = kF**4.0d0
      QThreePrime = -(6.0d0*gamma)/FermiK4

      return
      end
c
c     ---------------------------------------------------------------------------------------
c     Calculate the third derivative of EpsX
c     ---------------------------------------------------------------------------------------
c
      double precision function EpsThreeXprime(Rho,gamma)

      implicit none
c#include "xc_params.fh"

      double precision Rho, gamma
c Fermi wavevector stuff
      double precision kF, kFPrim, kF2prim, kF3prim
      double precision kFPrim2, kFPrim3 
c q stuff
      double precision q, qprim, q2prim, q3prim
      double precision qprim2, qprim3 
c H(q) stuff
      double precision g, gprim, g2prim, g3prim
      double precision Cs

      double precision FermiK
      double precision FermiKPrime
      double precision FermiKTwoPrime
      double precision FermiKThreePrime
      double precision QPrime
      double precision QTwoPrime
      double precision QThreePrime
      double precision HqBNL
      double precision HqBNLPrime
      double precision HqBNLTwoPrime
      double precision HqBNLThreePrime
      double precision pi

      pi = 3.1415926535897932385d0

      if (Rho.le.0d0) then
         EpsThreeXprime = 0.0d0
         return
      end if

      kF = FermiK(Rho)
      kFPrim = FermiKPrime(Rho)
      kF2Prim = FermiKTwoPrime(Rho)
      kF3Prim = FermiKThreePrime(Rho)
c
      kFPrim2 = kFPrim**2.0d0
      kFPrim3 = kFPrim**3.0d0
c
      Cs = -3.0d0/(4.0d0*Pi)
c
      q = gamma / kF
      qprim = QPrime(gamma,kF)
      q2prim = QTwoPrime(gamma,kF)
      q3prim = QThreePrime(gamma,kF)
c
      qprim2 = qprim**2.0d0
      qprim3 = qprim**3.0d0
c
      g = HqBNL(q)
      gprim = HqBNLPrime(q)
      g2prim = HqBNLTwoPrime(q)
      g3prim = HqBNLThreePrime(q)

      EpsThreeXprime = Cs*kFPrim3*( 3.0d0*qprim2*g2prim
     1                            + 3.0d0*kF*qprim*g2prim*q2prim
     2                            + kF*qprim3*g3prim
     3                            + gprim*( 3.0d0*q2prim
     4                                    + kF*q3prim ) )
     5               + 3.0d0*Cs*kFPrim*kF2Prim*( kF*qprim2*g2prim
     6                                         + gprim*( 2.0d0*qprim
     7                                                 + kF*q2prim ) )
     8               + Cs*kF3Prim*( g + kF*gprim*qprim )

      return
      end
c
c $Id: xc_bnl.F 23959 2013-04-03 23:37:19Z dwsilverstein $
