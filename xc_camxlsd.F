c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
      Subroutine xc_camxlsd(cam_omega, n, tol_rho, fac, lfac, nlfac, 
     1                      rho, Amat, Amat2, Amat3, nq, ipol, Ex, 
     2                      qwght, ldew, func, a)
c
C$Id: xc_camxlsd.F 20247 2011-04-28 18:58:49Z d3y133 $
c
      Implicit none
c
      double precision a
      double precision cam_omega
      double precision cam_alpha, cam_beta
c
      integer nq, ipol
      double precision fac, Ex
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,3)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial First Derivatives of the Exchange Energy Functional
c
      double precision Amat(nq,2)
      double precision Etmp,Atmp,Ctmp,A2tmp,C2tmp,C3tmp
c
c
c     Partial Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,3)
c Daniel (4-10-13): Third derivatives
c
c     Partial Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,4)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
      double precision rhom23
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
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
      double precision rho13, rho32, rho33, one_third
      Parameter (one_third = 1.d0/3.d0)
c Daniel (4-10-13): For 3rd derivatives
      double precision two_ninth
      Parameter (two_ninth = 2.0d0/9.0d0)
      integer n
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
c         do 10 n = 1, nq
            if (rho(n,1).gt.tol_rho)then
             rho13=rho(n,1)**one_third
             Etmp = rho(n,1)*rho13*P2*fac
             if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
             Atmp = rho13*P1*fac
             Ctmp = 0.0d0
c
             A2tmp = (rho13/rho(n,1))*2.0d0*one_third*P1*fac
             C2tmp = 0.0d0
             C3tmp = 0.0d0
c Daniel (4-10-13): Thid derivatives
             rhom23 = rho13/rho(n,1)
             A3tmp = (rhom23/rho(n,1))*-4.0d0*two_ninth*P1*fac
             C4tmp = 0.0d0
             C5tmp = 0.0d0
             C6tmp = 0.0d0
c Daniel (4-10-13): Do the attenuation (I'm using the same alpha and
c beta from CAM-B3LYP for testing)
             cam_alpha = 0.19d0
             cam_beta = 0.46d0
             call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp,cam_omega,
     &           cam_alpha,cam_beta, a)
c
             Amat2(n,1) = Amat2(n,1) + A2tmp
c
             Amat3(n,1) = Amat3(n,1) + A3tmp
             Ex = Ex + qwght(n)*Etmp
             Amat(n,1) = Amat(n,1) + Atmp
            endif
c   10    continue
c
      else
cc
cc        ======> SPIN-UNRESTRICTED <======
cc
c         do 20 n = 1,nq
c           if (rho(n,1).gt.tol_rho)then
c             rho32=0.0d0
c             rho33=0.0d0
c             if (rho(n,2).gt.tol_rho) rho32=rho(n,2)**one_third
c             if (rho(n,3).gt.tol_rho) rho33=rho(n,3)**one_third
cc ---- alpha ----
c             Etmp = rho32*rho(n,2)*P4*fac
cc
c             Atmp = P3*rho32*fac
c             Ctmp = 0.0d0
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c             A2tmp = 0.0d0
c             C2tmp = 0.0d0
c             C3tmp = 0.0d0
cc
c             if (rho(n,2).gt.tol_rho) then
c               A2tmp = one_third*P3*rho32/rho(n,2)*fac
c             endif
c#endif
cc Daniel (4-10-13): Third derivatives
c#ifdef THIRD_DERIV
c             A3tmp = 0.0d0
c             C4tmp = 0.0d0
c             C5tmp = 0.0d0
c             C6tmp = 0.0d0
cc
c             if (rho(n,2).gt.tol_rho) then
c               A3tmp = -two_ninth*P3*rho32/(rho(n,2)**2)*fac 
c             endif
c#endif
cc Daniel (4-10-13): Do the attenuation
c#ifdef THIRD_DERIV
c             if (rho(n,2).gt.tol_rho) then
c               call xc_att_xc_d3(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c             endif
cc
c             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
cc
c             Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
c#elif defined(SECOND_DERIV)
c             if (rho(n,2).gt.tol_rho) then
c               call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp)
c             endif
cc
c             Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c#else
c             call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
c#endif
c             Amat(n,1) = Amat(n,1) + Atmp
c             Ex = Ex + qwght(n)*Etmp
cc ---- beta ----
c             Etmp = rho33*rho(n,3)*P4*fac
c             Atmp = P3*rho33*fac
c             Ctmp = 0.0d0
c#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c             A2tmp = 0.0d0
c             C2tmp = 0.0d0
c             C3tmp = 0.0d0
cc
c             if (rho(n,3).gt.tol_rho) then
c               A2tmp = one_third*P3*rho33/rho(n,3)*fac
c             end if
c#endif
cc Daniel (4-10-13): Third derivatives
c#ifdef THIRD_DERIV
c             A3tmp = 0.0d0
c             C4tmp = 0.0d0
c             C5tmp = 0.0d0
c             C6tmp = 0.0d0
cc
c             if (rho(n,3).gt.tol_rho) then
c               A3tmp = -two_ninth*P3*rho33/(rho(n,3)**2)*fac
c             endif
c#endif
cc Daniel (4-10-13): Do the attenuation
c#ifdef THIRD_DERIV
c             if (rho(n,3).gt.tol_rho) then
c               call xc_att_xc_d3(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c             endif
cc
c             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
cc
c             Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
c#elif defined(SECOND_DERIV)
c             if (rho(n,3).gt.tol_rho) then
c               call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
c     &           C2tmp,C3tmp)
c             end if
cc
c             Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c#else
c             call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
c#endif
c             Amat(n,2) = Amat(n,2) + Atmp
c             Ex = Ex + qwght(n)*Etmp
cc
c             if (ldew)func(n) = func(n) + ( rho32*rho(n,2) + 
c     &                                      rho33*rho(n,3)   )*P4*fac
c            endif
c   20    continue
cc 
      endif
c
      return
      end
