      program xcfuntester

      implicit none

      character*10 functional
      character*10 errortype
      character*10 pol
      character*68 text

      integer ipol

c Electron density and gradient invariant
      double precision rho(1,3)
      double precision gamma
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
c CAM functional asigma variable
      double precision asigma

c --------------------------------------------------------------------
c Program xcfuntester
c
c Purpose: To test analytical formulas for the derivatives of exchange
c and correlation functionals in comparison to numerical derivatives,
c based on a 5-point differentiation formula.
c --------------------------------------------------------------------

c Read in the functional from the command line.
      call getarg(1,functional)

c Determine if the error analysis is to be done by a residual method, or
c by a ratio between the numerical and analytical derivatives.  Options
c are either: residual or ratio.
      call getarg(2,errortype)
      if (errortype.eq.' ') then
        errortype = 'residual'
      endif

c Get whether is this a restricted or unrestricted result.  Default is
c always restricted.  This isn't really used at the moment!!!
      call getarg(3, pol)
      if (pol.ne.' ') then
        if (pol.eq.'res') then
          ipol = 1
        else
          ipol = 2
          write(*,*) 'Unrestricted calculations are not implemented'
          stop
        endif
      else
        ipol = 1
      endif
 
c Initialize the values of all derivatives
      call array_init(Amat, Cmat, Amat2, Cmat2, Cmat11, Amat3, Cmat3,
     &  Cmat21, Cmat12, d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, 
     &  d11E_rho_gamma, d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &  d12E_rho_gamma_gamma)

c Initialize asigma
      asigma = 0.0d0

c ---------------------------------------------
c Evaluate Analytical and Numerical Derivatives
c ---------------------------------------------
c Determine what to do, based on the requested exchange or correlation
c functional.  Numerical derivatives use a 5-point differentiation 
c formula to obtain high accuracy.
      if (functional.eq.'xslater') then
        call dirac_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma)
      else if (functional.eq.'cvwn1') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cvwn2') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cvwn3') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cvwn4') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cvwn5') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cvwn1rpa') then
        call vwn_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, functional)
      else if (functional.eq.'cpw91lda') then
        call pw91lda_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma)
      else if (functional.eq.'xbecke88') then
        call becke88_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma)
      else if (functional.eq.'cperdew86') then
        call perdew86_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma)
      else if (functional.eq.'clyp') then
        call lyp88_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma)
      else if (functional.eq.'xpbe96') then
        call xpbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma, 'pbe')
      else if (functional.eq.'xrevpbe') then
        call xpbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma, 'revp')
      else if (functional.eq.'xrpbe') then
        call xpbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma, 'rpbe')
      else if (functional.eq.'cpbe96') then
        call cpbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma)
      else if (functional.eq.'xcamlsd') then
        call camlsd_wrapper(Amat, Amat2, Amat3,
     &    d1E_rho, d2E_rho, d3E_rho, rho, gamma, asigma)
      else if (functional.eq.'xcamb88') then
        call camb88_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma, asigma)
      else if (functional.eq.'xbnl07') then
        call bnl_wrapper(Amat, Amat2, Amat3,
     1    d1E_rho, d2E_rho, d3E_rho, rho, gamma)
      else if (functional.eq.'xcampbe96') then
        call campbe96_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma, asigma, 'pbe')
      else if (functional.eq.'xwpbe') then
        call wpbe_wrapper(Amat, Cmat, Amat2, Cmat2, Cmat11,
     &    Amat3, Cmat3, Cmat21, Cmat12,
     &    d1E_rho, d1E_gamma, d2E_rho, d2E_gamma, d11E_rho_gamma,
     &    d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &    d12E_rho_gamma_gamma, rho, gamma)
      else if (functional.eq.'') then
        write(*,*)
        write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
        write(*,*) 'WARNING!'
        write(text,*) 'No functional entered'
        call colorize(text, .true.)
        write(*,*)
        write(*,*) 'Try one of the following (key => meaning):'
        write(*,30) ' xslater   => Dirac Exchange',
     &              ' cvwn1     => VWN 1 Correlation',
     &              ' cvwn2     => VWN 2 Correlation',
     &              ' cvwn3     => VWN 3 Correlation',
     &              ' cvwn4     => VWN 4 Correlation',
     &              ' cvwn5     => VWN 5 Correlation',
     &              ' cvwn1rpa  => VWN 1 RPA Correlation',
     &              ' cpw91lda  => PW91 LDA Correlation',
     &              ' xbecke88  => Becke 88 Exchange',
     &              ' cperdew86 => Perdew 86 Correlation',
     &              ' clyp      => LYP Correlation',
     &              ' xpbe96    => PBE 96 Exchange',
     &              ' xrevpbe   => revPBE Exchange',
     &              ' xrpbe     => RPBE Exchange',
     &              ' cpbe96    => PBE 96 Correlation',
     &              ' xcamlsd   => CAM-Slater/Dirac Exchange',
     &              ' xcamb88   => CAM-Becke 88 Exchange',
     &              ' xbnl07    => BNL 07 Exchange',
     &              ' xcampbe96 => CAM-PBE 96 Exchange',
     &              ' xwpbe     => LC-wPBE Exchange'
        write(*,*)
        write(*,*) 'Please use one of the keys above'
        write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
        write(*,*)
        stop
      else
        write(*,*)
        write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
        write(*,*) 'WARNING!'
        write(text,*) 'Functional ' // trim(functional) //
     &    ' is not recognized.'
        call colorize(text, .true.)
        write(*,*)
        write(*,*) 'Try one of the following (key => meaning):'
        write(*,30) ' xslater   => Dirac Exchange',
     &              ' cvwn1     => VWN 1 Correlation',
     &              ' cvwn2     => VWN 2 Correlation',
     &              ' cvwn3     => VWN 3 Correlation',
     &              ' cvwn4     => VWN 4 Correlation',
     &              ' cvwn5     => VWN 5 Correlation',
     &              ' cvwn1rpa  => VWN 1 RPA Correlation',
     &              ' cpw91lda  => PW91 LDA Correlation',
     &              ' xbecke88  => Becke 88 Exchange',
     &              ' cperdew86 => Perdew 86 Correlation',
     &              ' clyp      => LYP Correlation', 
     &              ' xpbe96    => PBE 96 Exchange',
     &              ' xrevpbe   => revPBE Exchange',
     &              ' xrpbe     => RPBE Exchange',
     &              ' cpbe96    => PBE 96 Correlation',
     &              ' xcamlsd   => CAM-Slater/Dirac Exchange',
     &              ' xcamb88   => CAM-Becke 88 Exchange',
     &              ' xbnl07    => BNL 07 Exchange', 
     &              ' xcampbe96 => CAM-PBE 96 Exchange', 
     &              ' xwpbe     => LC-wPBE Exchange'
        write(*,*)
        write(*,*) 'Please use one of the keys above'
        write(*,*) '/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\'
        write(*,*)
        stop
      end if
  30  format(19(A,/),A)

c --------------
c Output Results
c --------------
      call output(Amat, Cmat, Amat2, Cmat2, Cmat11, Amat3, Cmat3,
     &  Cmat21, Cmat12, d1E_rho, d1E_gamma, d2E_rho, d2E_gamma,
     &  d11E_rho_gamma, d3E_rho, d3E_gamma, d21E_rho_rho_gamma,
     &  d12E_rho_gamma_gamma, functional, rho, gamma, errortype,
     &  asigma)

      end
