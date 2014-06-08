      subroutine array_init(Amat, Cmat, Amat2, Cmat2, Cmat11, 
     &  Amat3, Cmat3,
     &  Cmat21, Cmat12, d1E_rho, d1E_gamma, d2E_rho, d2E_gamma,
     &  d11E_rho_gamma, d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &  d12E_rho_gamma_gamma)

      implicit none

c Analytical derivatives
      double precision Amat(1,2),Amat2(1,3),Cmat(1,3),Cmat2(1,6)
      double precision Cmat11(1,6)
      double precision Amat3(1,4),Cmat21(1,9),Cmat12(1,12),Cmat3(1,10)
c Numerical derivatives
      double precision d1E_rho(1,2),d2E_rho(1,3)
      double precision d1E_gamma(1,3),d2E_gamma(1,6)
      double precision d11E_rho_gamma(1,6)
      double precision d3E_rho(1,4),d3E_gamma(1,10),
     &                 d21E_rho_rho_gamma(1,9)
      double precision d12E_rho_gamma_gamma(1,12)
c
      integer i

c ----------------------------------------------------------
c Subroutine array_init
c
c Purpose: Initializing all arrays used the XCfuntester code
c ----------------------------------------------------------

c Array indexing information:
c
c Amat   => 1 : dra,  2 : drb
c Cmat   => 1 : dgaa, 2 : dgab, 3 : dgbb
c
c Amat2  => 1 : dra2,     2 : dradrb,   3 : drbdrb 
c Cmat11 => 1 : dradgaa,  2 : dradgab,  3 : dradgbb 
c        => 4 : drbdgaa,  5 : drbdgab,  6 : drbdgbb 
c Cmat2  => 1 : dgaa2,    2 : dgaadgab, 3 : dgaadgbb 
c        => 4 : dgabdgab, 5 : dgabdgbb, 6 : dgbb2 
c
c Amat3  =>  1 : dra3,         2 : dra2drb,      3 : dradrb2,   4 : drb3
c Cmat21 =>  1 : dra2dgaa,     2 : dra2dgab,     3 : dra2dgbb 
c        =>  4 : dradrbdgaa,   5 : dradrbdgab,   6 : dradrbdgbb 
c        =>  7 : drb2dgaa,     8 : drb2dgab,     9 : drb2dgbb 
c Cmat12 =>  1 : dradgaa2,     2 : dradgaadgab,  3 : dradgaadgbb 
c        =>  4 : dradgabdgab,  5 : dradgabdgbb,  6 : dradgbb2
c        =>  7 : drbdgaa2,     8 : drbdgaadgab,  9 : drbdgaadgbb 
c        => 10 : drbdgabdgab, 11 : drbdgabdgbb, 12 : drbdgbb2
c Cmat3  =>  1 : dgaa3,        2 : dgaa2dgab,    3 : dgaa2dgbb 
c        =>  4 : dgaadgab2,    5 : dgaadgabdgbb, 6 : dgaadgbb2
c        =>  7 : dgab3,        8 : dgab2dgbb,    9 : dgabdgbb2
c        => 10 : dgbb3

c Initialize the arrays for the first derivatives wrt rho
      Amat(1,1) = 0.0d0      
      Amat(1,2) = 0.0d0      
      d1E_rho(1,1) = 0.0d0
      d1E_rho(1,2) = 0.0d0

c Initialize arrays with max index 3.
      do i = 1, 3
        Amat2(1,i) = 0.0d0
        Cmat(1,i) = 0.0d0
        d2E_rho(1,i) = 0.0d0
        d1E_gamma(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 4.
      do i = 1, 4
        Amat3(1,i) = 0.0d0
        d3E_rho(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 6.
      do i = 1, 6
        Cmat2(1,i) = 0.0d0
        Cmat11(1,i) = 0.0d0
        d2E_gamma(1,i) = 0.0d0
        d11E_rho_gamma(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 9.
      do i = 1, 9
        Cmat21(1,i) = 0.0d0
        d21E_rho_rho_gamma(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 10.
      do i = 1, 10
        Cmat3(1,i) = 0.0d0
        d3E_gamma(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 12.
      do i = 1, 12
        Cmat12(1,i) = 0.0d0
        d12E_rho_gamma_gamma(1,i) = 0.0d0
      enddo

      end
