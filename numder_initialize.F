      subroutine numder_initialize(func, E, Amat, Cmat, Amat2, 
     &  Cmat2, Cmat11, Amat3, Cmat3, Cmat21, Cmat12)

      implicit none

c Matrices storing derivatives derivatives
      double precision func(1), E
      double precision Amat(1,2),Amat2(1,3),Cmat(1,3),Cmat2(1,6)
      double precision Cmat11(1,6)
      double precision Amat3(1,4),Cmat21(1,9),Cmat12(1,12),Cmat3(1,10)
c
      integer i

c ----------------------------------------------------------
c Subroutine numder_initialize
c
c Purpose: Initializing all numerical derivative arrays used 
c          for the XCfuntester code
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

c Initialize the functional and energy
      func(1) = 0.0d0
      E = 0.0d0

c Initialize the arrays for the first derivatives wrt rho
      Amat(1,1) = 0.0d0      
      Amat(1,2) = 0.0d0      

c Initialize arrays with max index 3.
      do i = 1, 3
        Amat2(1,i) = 0.0d0
        Cmat(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 4.
      do i = 1, 4
        Amat3(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 6.
      do i = 1, 6
        Cmat2(1,i) = 0.0d0
        Cmat11(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 9.
      do i = 1, 9
        Cmat21(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 10.
      do i = 1, 10
        Cmat3(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 12.
      do i = 1, 12
        Cmat12(1,i) = 0.0d0
      enddo

      end

      subroutine numder_init_lda(func, E, Amat, Amat2, Amat3)

      implicit none

c Matrices storing derivatives derivatives
      double precision func(1), E
      double precision Amat(1,2), Amat2(1,3), Amat3(1,4)
c
      integer i

c ----------------------------------------------------------
c Subroutine numder_init_lda
c
c Purpose: Initializing all numerical derivative arrays used 
c          for LDA functionals in the XCfuntester code
c ----------------------------------------------------------

c Array indexing information:
c
c Amat   => 1 : dra,  2 : drb
c
c Amat2  => 1 : dra2,     2 : dradrb,   3 : drbdrb 
c
c Amat3  =>  1 : dra3,         2 : dra2drb,      3 : dradrb2,   4 : drb3

c Initialize the functional and energy
      func(1) = 0.0d0
      E = 0.0d0

c Initialize the arrays for the first derivatives wrt rho
      Amat(1,1) = 0.0d0      
      Amat(1,2) = 0.0d0      

c Initialize arrays with max index 3.
      do i = 1, 3
        Amat2(1,i) = 0.0d0
      enddo

c Initialize arrays with max index 4.
      do i = 1, 4
        Amat3(1,i) = 0.0d0
      enddo

      end
