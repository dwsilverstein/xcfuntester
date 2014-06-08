      Subroutine xc_att_xc(rho,ipol,Ex,Amat,Cmat,Amat2,Cmat2,Cmat3,
     1                     Amat3,Cmat4,Cmat5,Cmat6, cam_omega,
     2                     cam_alpha, cam_beta, as)
c
C$Id: xc_att_xc.F 23977 2013-04-05 23:11:55Z dwsilverstein $
c
      implicit none
c
      double precision cam_omega
      double precision cam_alpha, cam_beta
c
      double precision rho, Ex, Amat, Cmat
      integer ipol
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2, Cmat2, Cmat3
c Daniel (4-4-13): Third derivatives
c
c     Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3, Cmat4, Cmat5, Cmat6
c
c
c References:
c
c
c***************************************************************************
c
      double precision as
      double precision a, b, c, btmp,bfactor
c
      double precision a_first,a2_first,btmp_first
c
      double precision sqrt_pi,t1,t2,t3,t4,t5,t6,t7
      double precision alpha,beta, DERF
      double precision f10, f01, b_first
c Daniel (4-8-13): New variables
      double precision a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
      double precision ta, ta2, ta3, ta4, ta5, ta6, ta7, ta8, ta9, 
     1                 ta10
      double precision f43, f23
      double precision expf, erff

      Parameter (sqrt_pi = 1.77245385090552d0)
      Parameter (t7 = 2.666666666666666667d0)
      Parameter (f43 = 4.0d0/3.0d0)
      Parameter (f23 = 2.0d0/3.0d0)
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision a_second, a2_second, f20
      double precision b_second, btmp_second, t8
      double precision a3_second
c Daniel (4-5-13): New variables
      double precision f11, f02
c Daniel (4-4-13): Third derivatives
      double precision a_third, a2_third, a3_third, a4_third
      double precision f30, f21, f12, f03, f02a
      double precision b_third, btmp_third
      double precision t9


c calculate the a_sigma parameter

c         write(luout,*) 'alpha',alpha
c         write(luout,*) 'beta',beta
c         write(luout,*) 'mu',mu
c           
c
          if (ipol.eq.1) then
             Ex = Ex/2d0
             rho = rho/2d0
          endif
          a = cam_omega*sqrt(-2d0*Ex)/(6d0*sqrt_pi*rho)
          as = a
c
          alpha = cam_alpha
          beta = cam_beta
c 
          f10 = Amat/(2d0*Ex) -1d0/rho
          a_first = f10*a 
          f01 = Cmat/(2d0*Ex)
          a2_first = f01*a 
c
          f20 = Amat2/(2d0*Ex) - Amat*Amat/(2d0*Ex*Ex) 
     &        + 1d0/(rho*rho)
c Daniel (4-5-13): New variables
          f11 = Cmat2/(2.0d0*Ex) - Amat*Cmat/(2.0d0*Ex*Ex) 
c
          f02 = Cmat3/(2.0d0*Ex) - Cmat*Cmat/(4.0d0*Ex*Ex)
c
          a_second = a*(f10*f10 + f20) 

c OLD CODE
          a2_second = a*(f10*f01 + Cmat2/(2d0*Ex)  
     &              - Amat*Cmat/(2d0*Ex*Ex))

          a3_second = a*(Cmat3/(2d0*Ex) - Cmat*Cmat/(4d0*Ex*Ex))
c OLD CODE
c
c NEW CODE
c          a2_second = a*( f10*f01 + f11 )  
c
c          a3_second = a*f02
c NEW CODE
c
c Daniel (4-5-13): Third derivatives
c
c Amat3 = drdrdr
c Cmat4 = drdrdg
c Cmat5 = drdgdg
c Cmat6 = dgdgdg
          f02a = Cmat3/(2.0d0*Ex) - Cmat*Cmat/(2.0d0*Ex*Ex)
c
          f30 = Amat3/(2.0d0*Ex) 
     1        - 3.0d0*Amat2*Amat/(2.0d0*Ex*Ex)
     2        + Amat*Amat*Amat/(Ex*Ex*Ex)
     3        - 2.0d0/(rho*rho*rho)
c
          f21 = Cmat4/(2.0d0*Ex)
     1        - Cmat2*Amat/(Ex*Ex)
     2        - Amat2*Cmat/(2.0d0*Ex*Ex)
     3        + Amat*Amat*Cmat/(Ex*Ex*Ex)
c
          f12 = Cmat5/(2.0d0*Ex)
     1        - Cmat2*Cmat/(Ex*Ex)
     2        - Amat*Cmat3/(2.0d0*Ex*Ex)
     3        + Amat*Cmat*Cmat/(Ex*Ex*Ex)
c
          f03 = Cmat6/(2.0d0*Ex)
     1        - Cmat3*Cmat/(Ex*Ex)
     2        + Cmat*Cmat*Cmat/(2.0d0*Ex*Ex*Ex)
c
          a_third = a*( f10*f10*f10 + 3.0d0*f10*f20 + f30 )
c
          a2_third = a*( f10*f10*f01 + f20*f01 + 2.0d0*f10*f11 + f21 )
c
          a3_third = a*( f10*f01*f01 + 2.0d0*f11*f01 + f10*f02a + f12 )
c
          a4_third = a*( f01*f02 + f03 )
c
c Daniel (4-8-13): New variables for handling exponents of a in the
c stable region. 
          a2 = a*a  
          a3 = a2*a  
          a4 = a3*a  
          a5 = a4*a  
          a6 = a5*a  
          a7 = a6*a  
          a8 = a7*a  
          a9 = a8*a  
          a10 = a9*a  
          a11 = a10*a  
c Daniel (4-8-13): New variables for the numerically unstable regions.
          ta = 2d0*a
          ta2 = ta*ta
          ta3 = ta2*ta
          ta4 = ta3*ta
          ta5 = ta4*ta
          ta6 = ta5*ta
          ta7 = ta6*ta
          ta8 = ta7*ta
          ta9 = ta8*ta
          ta10 = ta9*ta
c
          expf = exp(-1d0/(4d0*a2))
          erff = DERF(1d0/(2d0*a))
c
          if (a .lt. 0.14d0) then
c             write(6,*) 'First: a is small'
c OLD CODE
c              a = 2d0*a
c              btmp = 1d0-(4d0/3d0)*sqrt_pi*a + 2d0*a*a 
c     &             - (2d0/3d0)*a*a*a*a
c              btmp = 1d0-btmp
c
c              btmp_first = (4d0/3d0)*(-sqrt_pi + 3d0*a +
c     &                   (2d0*exp(-1/(a*a)) - 2d0)*a*a*a)  
c              btmp_first = 2d0*btmp_first
c              a = a /2d0
c OLD CODE
c
c NEW CODE
              btmp = 1.0d0 - f43*sqrt_pi*ta 
     1             + 2.0d0*ta2 - f23*ta4
              btmp = 1.0d0 - btmp

              btmp_first = f43*( -sqrt_pi + 3.0d0*ta +
     &                           (2.0d0*expf - 2.0d0)*ta3 )  
              btmp_first = 2.0d0*btmp_first
c NEW CODE
          else if (a .lt. 4.25d0) then
c            write(6,*) 'First: a is medium'
c             stop
c OLD CODE
c             b = exp(-1d0/(4d0*a*a)) - 1d0
             b = expf - 1d0
             c = 2d0*a*a*b + 0.5d0
c OLD CODE
c
c NEW CODE
c             b = expf - 1d0
c             c = 2.0d0*a2*b + 0.5d0
cc
cc             btmp = t7*a*( sqrt_pi*t2 + 2.0d0*a*(b-c) )
c             btmp = t7*a*( sqrt_pi*erff + 2.0d0*a*(b-c) )
c NEW CODE
c
c OLD CODE
c             btmp = (8d0/3d0)*a*(sqrt_pi*DERF(1/(2d0*a)) + 2d0*a*(b-c))
             btmp = (8d0/3d0)*a*(sqrt_pi*erff + 2d0*a*(b-c))
             t1 = 1/a
             t2 = a*a
             t3 = 1/t2
c             t4 = exp(-0.25d0*t3)
             t4 = expf
             t5 = t4 -1d0
             t6 = t4 -2d0*t2*t5 - 1.5d0
c             btmp_first = -t7*a *
c     &       (2*a*(t4/(2*a**3) - 4d0*a*t5 - t1*t4) + 2d0*t6 -t3*t4) -
c     &         t7*(2*a*t6 + sqrt_pi*DERF(0.5d0*t1))
c WORKS
c             btmp_first = -t7*a *
c     &       (2*a*(t4/(2*a**3) - 4d0*a*t5 - t1*t4) + 2d0*t6 -t3*t4) -
c     &         t7*(2*a*t6 + sqrt_pi*erff)
c WORKS
c
c This is a simplified form of what was written above.
             btmp_first = -t7*a *
     &       ( -16.0d0*a2*t5 + 2.0d0*expf - 6.0d0 ) - t7*sqrt_pi*erff
c 
c OLD CODE
c
c NEW CODE
c Daniel (4-8-13): t3 is what I call M in my notes
c             t3 = -2.0d0*expf + ( 4.0d0*expf - 6.0d0 )*a
c     1          + 10.0d0*( -expf + 1.0d0 )*a2 + sqrt_pi*erff
c     2          - expf/a2 + expf/a
cc t7 = 8/3
c             btmp_first = -t7*t3
c NEW CODE
          else
c            write(6,*) 'First: a is large'
c            stop
c OLD CODE
c             a = 2d0*a
c             btmp = 1d0 - 1d0/(9d0*a*a) + 1d0/(60d0*a**4d0) -
c     &           1d0/(420d0*a**6d0) + 1d0/(3240d0*a**8d0) -
c     &           1d0/(27720d0*a**10d0)
c
c             btmp_first = -1d0/(4.5d0*a**3) + 1d0/(15d0*a**5d0) -
c     &                  1d0/(70d0*a**7d0) + 1d0/(405d0*a**9d0) 
c             btmp_first = btmp_first*2d0
c             a = a /2d0
c OLD CODE
c
c NEW CODE
             btmp = 1.0d0 - 1.0d0/(9.0d0*ta2) + 1.0d0/(60.0d0*ta4) 
     1            - 1.0d0/(420.0d0*ta6) + 1.0d0/(3240.0d0*ta8)
     2            - 1.0d0/(27720.0d0*ta10)

             btmp_first = -1.0d0/(4.5d0*ta3) + 1.0d0/(15.0d0*ta5) 
     1                  - 1.0d0/(70.0d0*ta7) + 1.0d0/(405.0d0*ta9) 
             btmp_first = btmp_first*2.0d0
c NEW CODE
          end if
c
          if (abs(a) .lt. 1d-40) then 
c             write(6,*) 'Second: a is small'
c             btmp_second = 16d0
             btmp_second = 16.0d0
c Daniel (4-5-13): Should this be 4.25d0? (This probably doesn't
c matter since both exponentials and error functions vary slowly
c as the argument gets large).
          else if (a .ge. 5d0)  then
c             write(6,*) 'Second: a is large'
c OLD CODE
c             btmp_second = 1d0/(6d0*a**4d0) - 1d0/(48d0*a**6d0) +
c     &                  1d0/(640d0*a**8d0) - 1d0/(11520d0*a**10d0)
c OLD CODE
c
c NEW CODE
             btmp_second = 1.0d0/(6.0d0*a4) - 1.0d0/(48.0d0*a6) 
     1                   + 1.0d0/(640.0d0*a8) - 1.0d0/(11520.0d0*a10)
c NEW CODE
             
          else 
c             write(6,*) 'Second: a is medium'
c OLD CODE
             t1 = a*a
             t2 = 1d0/t1
c             t3 = exp(-0.25d0*t2)
             t3 = expf
             t4 = 1d0/(a*a*a)
             t5 = t3 - 1d0
             t6 = -t2*t3
             t8 = -t3/a + 0.5d0*t4*t3 - 4d0*a*t5
c
c             btmp_second = -(8d0*a*(2d0*a*(t3/(4*a**6d0) -
c     &       2d0*t3/(a**4d0) +t6 - 4d0*t5) -t3/(2*a**5d0) +
c     &       4d0*t8 + 2d0*t4*t3)/3d0 + 16d0*(2d0*a*t8 +
c     &       2d0*(t3 - 2d0*t1*t5-1.5d0) + t6)/3d0) 
c This is a simplified form of what was written above.
             btmp_second = 16.0d0 - 128.0d0*a2
     &                   + (16.0d0 + 128.0d0*a2)*expf
c 
c OLD CODE
c
c NEW CODE
c Daniel (4-8-13): We now rely on the definitions of t1 and t2 from
c above, so be careful if you change anything in this routine!
c             t4 = 1.0d0/a3
cc
cc Daniel (4-8-13): t5 is what I call N in my notes
c             t5 = -expf + expf*a - 6.0d0*expf*a2 
c     1          + ( 4.0d0*expf - 6.0d0 )*a3
c     2          + 24.0d0*( -expf + 1.0d0 )*a4
cc
c             btmp_second = -t7*t4*t5
c NEW CODE
          end if
c
c Daniel (4-4-13): Third derivatives
c
          if (abs(a) .lt. 1.0d-40) then
c             write(6,*) 'Third: a is small'
             btmp_third = 0.0d0
          else if (a .ge. 5.0d0) then
c             write(6,*) 'Third: a is large'
             btmp_third = -2.0d0/(3.0d0*a5)
     1                  + 1.0d0/(8.0d0*a7)
     2                  - 1.0d0/(80.0d0*a9)
     3                  + 1.0d0/(1152.0d0*a11)
          else
c             write(6,*) 'Third: a is medium'
c SAVE
cc Daniel (4-8-13): Again, be careful if you change anything in this
cc routine.  We still need t1 and t2.
c             t6 = t4*t4
cc Daniel (4-8-13): t8 is what I call O in my notes
c             t8 = -t1 + t1*a - 12.0d0*t1*a4
c     1          + 48.0d0*( -t1 + 1.0d0 )*a6
cc
c             btmp_third = -t7*t6*t8*0.50d0
c SAVE
c 
c Code here mimics Dalton
            t1 = 1d0/(a*a)
            t2 = expf
            t3 = t1*t1*t1
            t4 = t1*t1
            t5 = t3*a
            t6 = t2 - 1d0
            t8 = -t1*t2 - 2d0*t4*t2 + t3*t2/4d0 - 4d0*t6
            t9 = t4*a
c
c            btmp_third = -8d0*( a*( 2d0*a*( t2/(8*a**9d0)
c     1                                    - 5d0*t2/(2d0*a**7d0) 
c     2                                    + 15d0*t5*t2/2d0 )
c     3                            - t2/(4d0*a**8d0) + 6d0*t8 
c     4                            - 6d0*t4*t2 + 7d0*t3*t2/2d0 )/3d0
c     5                        + ( 4d0*( -t2/a + t9*t2/2d0 - 4d0*a*t6 )
c     6                          + 2d0*a*t8 + 2d0*t9*t2 
c     7                          - t5*t2/2d0 ) )
c This is a simplified form of what was written above.
             btmp_third = 8.0d0*( -32.0d0*a4
     1                          + ( 1.0d0 + 8.0d0*a2 
     2                            + 32.0d0*a4 )*expf )/a3
c
          endif
c
          bfactor = 1d0 - alpha - beta*btmp
          b_first = beta*btmp_first
c
          b_second = beta*btmp_second
c
c Daniel (4-4-13): Third derivatives
c
          b_third = beta*btmp_third
c
          Amat3 = bfactor*Amat3 
     1          + 3.0d0*Amat2*b_first*a_first
     2          + 3.0d0*Amat*( b_second*a_first*a_first
     3                       + b_first*a_second )
     4          + Ex*( b_third*a_first*a_first*a_first
     5               + 3.0d0*b_second*a_first*a_second
     6               + b_first*a_third )
c
          Cmat4 = bfactor*Cmat4
     1          + 2.0d0*Cmat2*b_first*a_first
     2          + Amat2*b_first*a2_first
     3          + 2.0d0*Amat*( b_second*a_first*a2_first
     4                       + b_first*a2_second )
     5          + Cmat*( b_second*a_first*a_first
     6                 + b_first*a_second )
     7          + Ex*( b_third*a_first*a_first*a2_first
     8               + b_second*( a2_first*a_second
     9                          + 2.0d0*a_first*a2_second )
     A               + b_first*a2_third )
c
          Cmat5 = bfactor*Cmat5
     1          + 2.0d0*Cmat2*b_first*a2_first
     2          + Amat*( b_second*a2_first*a2_first
     3                 + b_first*a3_second )
     4          + Cmat3*b_first*a_first
     5          + 2.0d0*Cmat*( b_second*a_first*a2_first
     6                       + b_first*a2_second )
     7          + Ex*( b_third*a_first*a2_first*a2_first
     8               + b_second*( a_first*a3_second
     9                          + 2.0d0*a2_first*a2_second )
     A               + b_first*a3_third )
c
          Cmat6 = bfactor*Cmat6
     1          + 3.0d0*Cmat3*b_first*a2_first
     2          + 3.0d0*Cmat*( b_second*a2_first*a2_first
     3                       + b_first*a3_second )
     4          + Ex*( b_third*a2_first*a2_first*a2_first
     5               + 3.0d0*b_second*a2_first*a3_second
     6               + b_first*a4_third )
c
c
c          b_second = beta*btmp_second
          Amat2 = bfactor*Amat2 + 2d0*Amat*b_first*a_first 
     &          + Ex*b_second*a_first*a_first  
     &          + Ex*b_first*a_second

          Cmat2 = bfactor*Cmat2 + Amat*b_first*a2_first
     &          + Cmat*b_first*a_first
     &          + Ex*b_second*a_first*a2_first
     &          + Ex*b_first*a2_second
          
          Cmat3 = bfactor*Cmat3 + 2d0*Cmat*b_first*a2_first
     &          + Ex*b_second*a2_first*a2_first
     &          + Ex*b_first*a3_second
c
          Amat = bfactor*Amat + Ex*b_first*a_first
          Cmat = bfactor*Cmat + Ex*b_first*a2_first
          Ex = Ex*bfactor

          if (ipol.eq.1) then
             Ex = 2d0*Ex
             rho = 2d0*rho
          endif

c         write(luout,*) 'bfactor',bfactor
c
c

c
      return
      end
