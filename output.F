      subroutine output(Amat, Cmat, Amat2, Cmat2, Cmat11, Amat3, 
     &  Cmat3,
     &  Cmat21, Cmat12, d1E_rho, d1E_gamma, d2E_rho, d2E_gamma,
     &  d11E_rho_gamma, d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &  d12E_rho_gamma_gamma, functional, rho, gamma, errortype,
     &  asigma)      

      implicit none
c
      character*10 functional
      character*100 fname
      character*100 reference
      character*10 errortype
      character*68 text
c
      integer i      
c
      double precision Amat(1,2),Amat2(1,3),Cmat(1,3),Cmat2(1,6),
     &                 Cmat11(1,6)
      double precision Amat3(1,4),Cmat21(1,9),Cmat12(1,12),Cmat3(1,10)
c
      double precision d1E_rho(1,2),d2E_rho(1,3)
      double precision d1E_gamma(1,3),d2E_gamma(1,6)
c
      double precision d11E_rho_gamma(1,6)
c
      double precision d3E_rho(1,4),d3E_gamma(1,10),
     &                 d21E_rho_rho_gamma(1,9)
      double precision d12E_rho_gamma_gamma(1,12)
c
      character*255 frmt
c
      double precision rho(1,3), gamma
      double precision residual(9), maxres
      double precision ratio(9), maxratio
      double precision, parameter :: thresh = 1.0D-8 
      double precision, parameter :: ratiothresh = 1.0D-6
c 
      logical lpgamma
c
      logical lcam
c
      double precision asigma

c Logical flag for printing gamma values (false for LDA)
      lpgamma = .true.

c Logical flag for printing a_sigma for CAM-functionals
      lcam = .false.

c Initialize the residual vector.
      do i = 1, 9
        residual(i) = 0.0d0
      enddo

c Determine the functional name based on the user input.
      if (functional.eq.'xslater') then
        fname = 'Dirac (Slater) Exchange'
        reference = 'Phys. Rev. B, 5, 844 (1972)'
        lpgamma = .false.
      else if (functional.eq.'cvwn1') then
        fname = 'VWN 1 Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cvwn2') then
        fname = 'VWN 2 Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cvwn3') then
        fname = 'VWN 3 Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cvwn4') then
        fname = 'VWN 4 Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cvwn5') then
        fname = 'VWN 5 Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cvwn1rpa') then
        fname = 'VWN 1 RPA Correlation'
        reference = 'Can. J. Phys., 58, 1200 (1980)'
        lpgamma = .false.
      else if (functional.eq.'cpw91lda') then
        fname = 'Perdew-Wang 91 LDA Correlation'
        reference = 'Phys. Rev. B, 45, 13244 (1992)'
        lpgamma = .false.
      else if (functional.eq.'xbecke88') then
        fname = 'Becke 88 Exchange'
        reference = 'Phys. Rev. A, 88, 3098 (1988)'
      else if (functional.eq.'cperdew86') then
        fname = 'Perdew 86 Correlation'
        reference = 'Phys. Rev. B, 33, 8822 (1986)'
      else if (functional.eq.'clyp') then
        fname = 'LYP Correlation'
        reference = 'Phys. Rev. B, 37, 785 (1988)'
      else if (functional.eq.'xpbe96') then
        fname = 'PBE 96 Exchange'
        reference = 'Phys. Rev. Lett., 77, 3865 (1996); 78, 1396 (1997)'
      else if (functional.eq.'xrevpbe') then
        fname = 'revPBE Exchange'
        reference = 'Phys. Rev. Lett., 80, 890 (1998)'
      else if (functional.eq.'xrpbe') then
        fname = 'RPBE Exchange'
        reference = 'Phys. Rev. B, 59, 7413 (1999)'
      else if (functional.eq.'cpbe96') then
        fname = 'PBE 96 Correlation'
        reference = 'Phys. Rev. Lett., 77, 3865 (1996); 78, 1396 (1997)'
      else if (functional.eq.'xcamlsd') then
        fname = 'CAM-Slater/Dirac Exchange'
        reference = '(None, but this functional is similar to BNL)'
        lpgamma = .false.
        lcam = .true.
      else if (functional.eq.'xcamb88') then
        fname = 'CAM-Becke 88 Exchange'
        reference = 'Chem. Phys. Lett., 393, 51 (2004)'
        lcam = .true.
      else if (functional.eq.'xbnl07') then
        fname = 'BNL 07 Exchange'
        reference = 'Phys. Chem. Chem. Phys. 9, 2932 (2007)'
        lpgamma = .false.
      else if (functional.eq.'xcampbe96') then
        fname = 'CAM-PBE 96 Exchange'
        reference = 'J. Chem. Phys. 129, 034107 (2008)'
        lcam = .true.
      else if (functional.eq.'xwpbe') then
        fname = 'LC-wPBE Exchange'
        reference = 'J. Chem. Phys., 130, 054112 (2009)'
      end if

c Print the header.
      write(*,*)
      write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
      write(text,*) 'Functional Derivative Tester'
      call colorize(text, .false.)
      write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\'

      write(*,*) 'Results for functional: ', trim(fname)
      write(*,*) 'Reference: ', trim(reference)
      write(*,'(A20,E10.2E2)') ' Electron Density = ', rho(1,1)
      if (lpgamma) then
        write(*,'(A30,E10.2E2)') ' Density Gradient Invariant = ', gamma
      end if
      if (lcam) then
        write(*,'(A15,E10.2E2)') ' CAM a_sigma = ', asigma
        if ( asigma .lt. 0.14d0 ) then
          write(*,'(A33)') 'Small a_sigma equation was used.'
        else if ( asigma .gt. 4.25d0 ) then
          write(*,'(A33)') 'Large a_sigma equation was used.'
        else
          write(*,'(A40)') 'Intermediate a_sigma equation was used.'
        endif
      endif
      write(*,*)

      frmt = '(A,2X,E20.12,2X,A)'

      write(text,*) 'Analytical Derivative Results:'
      call colorize(text, .false.)
      write(*,*) '---------------------------------'
      write(*,frmt) '| Amat   |', Amat(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat   |', Cmat(1,1), '|'
      end if
      write(*,frmt) '| Amat2  |', Amat2(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat2  |', Cmat2(1,1), '|'
        write(*,frmt) '| Cmat11 |', Cmat11(1,1), '|'
      end if
      write(*,frmt) '| Amat3  |', Amat3(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat3  |', Cmat3(1,1), '|'
        write(*,frmt) '| Cmat21 |', Cmat21(1,1), '|'
        write(*,frmt) '| Cmat12 |', Cmat12(1,1), '|'
      end if
      write(*,*) '---------------------------------'
      write(*,*) 
      write(text,*) 'Numerical Differentiation Results:'
      call colorize(text, .false.)
      write(*,*) '---------------------------------'
      write(*,frmt) '| Amat   |', d1E_rho(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat   |', d1E_gamma(1,1), '|'
      end if
      write(*,frmt) '| Amat2  |', d2E_rho(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat2  |', d2E_gamma(1,1), '|'
        write(*,frmt) '| Cmat11 |', d11E_rho_gamma(1,1), '|'
      end if
      write(*,frmt) '| Amat3  |', d3E_rho(1,1), '|'
      if (lpgamma) then
        write(*,frmt) '| Cmat3  |', d3E_gamma(1,1), '|'
        write(*,frmt) '| Cmat21 |', d21E_rho_rho_gamma(1,1), '|'
        write(*,frmt) '| Cmat12 |', d12E_rho_gamma_gamma(1,1), '|'
      end if
      write(*,*) '---------------------------------'

c --------------
c Error Analysis
c --------------
c     Determine residuals for error between the analytical and 
c     numerical derivatives.  We arbitrarily define that
c     |Analytical - Numerical| is the residual.  If errortype is
c     ratio, then the error is given by Analytical/Numerical.

c There's sadly no systematic method for doing this besides redefining 
c all of the arrays, so each residual element is set individually.
      if (errortype.eq.'residual') then
        residual(1) = dabs(Amat(1,1)-d1E_rho(1,1))                ! dr
        residual(2) = dabs(Cmat(1,1)-d1E_gamma(1,1))              ! dg
        residual(3) = dabs(Amat2(1,1)-d2E_rho(1,1))               ! drdr
        residual(4) = dabs(Cmat2(1,1)-d2E_gamma(1,1))             ! dgdg
        residual(5) = dabs(Cmat11(1,1)-d11E_rho_gamma(1,1))       ! drdg
        residual(6) = dabs(Amat3(1,1)-d3E_rho(1,1))               ! drdrdr
        residual(7) = dabs(Cmat3(1,1)-d3E_gamma(1,1))             ! dgdgdg
        residual(8) = dabs(Cmat21(1,1)-d21E_rho_rho_gamma(1,1))   ! drdrdg
        residual(9) = dabs(Cmat12(1,1)-d12E_rho_gamma_gamma(1,1)) ! drdgdg
  
        write(*,*)
        write(text,*) 'Error Analysis (|Analytical - Numerical|):'
        call colorize(text, .false.)
        write(*,*) '---------------------------------'
        write(*,frmt) '| Amat   |', residual(1), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat   |', residual(2), '|'
        end if
        write(*,frmt) '| Amat2  |', residual(3), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat2  |', residual(4), '|'
          write(*,frmt) '| Cmat11 |', residual(5), '|'
        end if
        write(*,frmt) '| Amat3  |', residual(6), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat3  |', residual(7), '|'
          write(*,frmt) '| Cmat21 |', residual(8), '|'
          write(*,frmt) '| Cmat12 |', residual(9), '|'
        end if
        write(*,*) '---------------------------------'
      else
        ratio(1) = Amat(1,1)/d1E_rho(1,1)                ! dr
        ratio(2) = Cmat(1,1)/d1E_gamma(1,1)              ! dg
        ratio(3) = Amat2(1,1)/d2E_rho(1,1)               ! drdr
        ratio(4) = Cmat2(1,1)/d2E_gamma(1,1)             ! dgdg
        ratio(5) = Cmat11(1,1)/d11E_rho_gamma(1,1)       ! drdg
        ratio(6) = Amat3(1,1)/d3E_rho(1,1)               ! drdrdr
        ratio(7) = Cmat3(1,1)/d3E_gamma(1,1)             ! dgdgdg
        ratio(8) = Cmat21(1,1)/d21E_rho_rho_gamma(1,1)   ! drdrdg
        ratio(9) = Cmat12(1,1)/d12E_rho_gamma_gamma(1,1) ! drdgdg

        write(*,*)
        write(text,*) 'Error Analysis (Analytical/Numerical):'
        call colorize(text, .false.)
        write(*,*) '---------------------------------'
        write(*,frmt) '| Amat   |', ratio(1), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat   |', ratio(2), '|'
        end if
        write(*,frmt) '| Amat2  |', ratio(3), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat2  |', ratio(4), '|'
          write(*,frmt) '| Cmat11 |', ratio(5), '|'
        end if
        write(*,frmt) '| Amat3  |', ratio(6), '|'
        if (lpgamma) then
          write(*,frmt) '| Cmat3  |', ratio(7), '|'
          write(*,frmt) '| Cmat21 |', ratio(8), '|'
          write(*,frmt) '| Cmat12 |', ratio(9), '|'
        end if
        write(*,*) '---------------------------------'
      endif

c Print out term definitions
      write(*,*)
      write(text,*) 'VARIABLE DEFINITIONS'
      call colorize(text, .false.)
      write(*,*) '--------------------'
      write(*,*) 'r == rho'
      write(*,*) 'g == gamma (for GGAs)'
      write(*,*) 'E == exchange/correlation energy functional'
      write(*,*)
      write(*,*) '         dE                  dE '
      write(*,*) 'Amat == ----      Cmat ==   ----'
      write(*,*) '         dr                  dg '
      write(*,*)
      write(*,*) '          d^2E               d^2E              d^2E '
      write(*,*) 'Amat2 == ------   Cmat2 ==  ------  Cmat11 == ------'
      write(*,*) '          drdr               dgdg              drdg '
      write(*,*)
      write(*,*) '          d^3E               d^3E '
      write(*,*) 'Amat3 == ------   Cmat3 ==  ------'
      write(*,*) '         drdrdr             dgdgdg'
      write(*,*)
      write(*,*) '           d^3E              d^3E '
      write(*,*) 'Cmat21 == ------  Cmat12 == ------'
      write(*,*) '          drdrdg            drdgdg'
      write(*,*) 

c Check that the maximum residual is not above an arbitrary threshold
c (1E-8 for the moment).  Warn the user that it may be worthwhile to
c check for numerical instabilities in the implementation.
      if (errortype.eq.'residual') then
        maxres = 0.0d0
        do i = 1, 9
          maxres = max(maxres, residual(i))
        enddo

        if (maxres.gt.thresh) then
          write(*,*) 
          write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' //
     &      '~~~~~~~~~~~~~~~~~~~'
          write(text,*) 'WARNING: Numerical instabilities exist for' //
     &      ' this functional!!!'
          call colorize(text, .true.)
          write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' //
     &      '~~~~~~~~~~~~~~~~~~~'
          write(*,*) 
        endif
      else
        maxratio = 0.0d0
        do i = 1, 9
          maxratio = max(maxratio, ratio(i))
        enddo      
c       Check that the maximum ratio is not larger than the threshold
c       (1E-10) for ratios between analytical and numerical results.
        maxres = dabs(maxratio - 1.0d0)
        if (maxres.gt.ratiothresh) then
          write(*,*)
          write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' //
     &      '~~~~~~~~~~~~~~~~~~~'
          write(text,*) 'WARNING: Numerical instabilities exist for' //
     &      ' this functional!!!'
          call colorize(text, .true.)
          write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' //
     &      '~~~~~~~~~~~~~~~~~~~'
          write(*,*) 
        endif
      endif

      end
