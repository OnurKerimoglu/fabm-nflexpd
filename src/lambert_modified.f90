! F90 module made from TOMS 743
! real gets 6 digits right, double precision 15
! modified by Onur Kerimoglu (OK): expressions like 1/0. were replaced by 1/zero
! where zero=0.0, in order to avoid compilation problems (search OK to see where)
MODULE lambert
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: dp=KIND(0D0)
  INTERFACE lambertw
     MODULE PROCEDURE wapd !, wapr
  END INTERFACE
CONTAINS
  ! lwm1: closed-form approximation of -1-branch of the Lambert-W function
  !       after Barry et al. (2000)
  !       the maximum relative error of this approximation is 0.025%
  ELEMENTAL FUNCTION lwm1 (x) RESULT (w)
    IMPLICIT NONE
    REAL(dp), PARAMETER :: m1=0.3361D0, m2=-0.0042D0, m3=-0.0201D0
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: w, s
    s = -1D0 - LOG(-x)
    w = -1D0 - s - 2D0/m1*(1D0 - 1D0/(1D0 + (m1*SQRT(0.5D0*s)&
                                             /(1D0 + m2*s*EXP(m3*SQRT(s))))))
  END FUNCTION lwm1
!     __________________________________________________________________
!
!     Approximating the W function
!     ____________________________
!
  ELEMENTAL FUNCTION WAPR (X, NB, L) RESULT (WAP)
!
!     WAPR - output
!     X - argument of W(X)
!     NB is the branch of the W function needed:
!        NB = 0 - upper branch
!        NB <> 0 - lower branch
!
!     NERROR is the output error flag:
!        NERROR = 0 -> routine completed successfully
!        NERROR = 1 -> X is out of range
!
!     Range: -exp(-1) <= X for the upper branch of the W function
!            -exp(-1) < X < 0 for the lower branch of the W function
!
!     L - determines how WAPR is to treat the argument X
!        L = 1 -> X is the offset from -exp(-1), so compute
!                 W(X-exp(-1))
!        L <> 1 -> X is the desired X, so compute W(X)
!
!     M - print messages from WAPR?
!         M = 1 -> Yes
!         M <> 1 -> No
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     NN is the output device number
!
!     NBITS is the number of bits (less 1) in the mantissa of the
!        floating point number number representation of your machine.
!        It is used to determine the level of accuracy to which the W
!        function should be calculated.
!
!        Most machines use a 24-bit matissa for single precision and
!        53-56 bits for double precision. The IEEE standard is 53
!        bits. The Fujitsu VP2200 uses 56 bits. Long word length
!        machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!        single precision.
!
    IMPLICIT NONE
    INTEGER, PARAMETER :: NN=6, NBITS=23, NITER=1
    REAL, PARAMETER ::EM=-0.367879441171442,&           ! -EXP(-1)
                      EM9=-1.234098040866796E-4,&       ! -EXP(-9)
                      C13=1.E0/3.E0,&
                      C23=2.E0*C13,&
                      EM2=2.E0/EM,&
                      E12=-EM2,&
                      TB=.5E0**NBITS,&
                      TB2=.5E0**(NBITS/2),&       ! SQRT(TB)
                      X0=0.0350769390096679055,&  ! TB**(1/6)*0.5E0
                      X1=0.302120119432784731,&   !(1 - 17*TB**(2/7))*EM
                      AN22=3.6E2/83.E0,&
                      AN11=135./83.E0,&
                      AN3=8.E0/3.E0,&
                      AN4=135.E0/83.E0,&
                      AN5=166.E0/39.E0,&
                      AN6=3167.E0/3549.E0,&
                      S2=1.41421356237310,& ! SQRT(2.E0)
                      S21=2.E0*S2-3.E0,&
                      S22=4.E0-3.E0*S2,&
                      S23=S2-2.E0
    REAL :: X, WAP, AN2, DELX, XX, RETA, ZL, T, TS, ETA, TEMP, TEMP2, ZN, zero
    INTEGER :: NB, L
    INTENT(IN) :: X, NB, L
!        Various mathematical constants
    
!
!     The following COMMON statement is needed only when testing this
!     function using BISECT, otherwise it can be removed.
!
!    COMMON/WAPCOM/NBITS
!    DATA INIT,NITER/0,1/
!     DATA NBITS/23/
!
!     IF(INIT.EQ.0) THEN
!        INIT=1
!
!        Code to calculate NBITS for the host machine. NBITS is the
!        mantissa length less one. This value is chosen to allow for
!        rounding error in the final bit. This need only be run once on
!        any particular machine. It can then be included in the above
!        DATA statement.
!
!        DO I=1,2000
!           B=2.E0**(-I)
!           V=1.E0+B
!           IF(V.EQ.1.E0)THEN
!              NBITS=I-1
!              J=-ALOG10(B)
!              IF(M.EQ.1) WRITE(NN,40)NBITS,J
!              EXIT
!           ENDIF
!        END DO
!
!        Remove to here after NBITS has been calculated once
!
!        The case of large NBITS
!
!        IF(NBITS.GE.56) NITER=2
!
!        Various mathematical constants
!
!        EM=-EXP(-1.E0)
!        EM9=-EXP(-9.E0)
!        C13=1.E0/3.E0
!        C23=2.E0*C13
!        EM2=2.E0/EM
!        E12=-EM2
!        TB=.5E0**NBITS
!        TB2=SQRT(TB)
!        X0=TB**(1.E0/6.E0)*.5E0
!        X1=(1.E0-17.E0*TB**(2.E0/7.E0))*EM
!        AN22=3.6E2/83.E0
!        AN11=135./83.E0
!        AN3=8.E0/3.E0
!        AN4=135.E0/83.E0
!        AN5=166.E0/39.E0
!        AN6=3167.E0/3549.E0
!        S2=SQRT(2.E0)
!        S21=2.E0*S2-3.E0
!        S22=4.E0-3.E0*S2
!        S23=S2-2.E0
!     ENDIF
    !OK: to avoid compilation problems due to division by zero, 
    ! set zero=0.0 and use it instead of 0.0
    zero=0.0
    IF(L.EQ.1) THEN
       DELX=X
       IF(DELX.LT.0.E0) THEN
          WAP = 1./zero
          RETURN
       END IF
       XX=X+EM
!        IF(E12*DELX.LT.TB**2.AND.M.EQ.1) WRITE(NN,60)DELX
    ELSE
       IF(X.LT.EM) THEN
          WAP = 1./zero
          RETURN
       END IF
       IF(X.EQ.EM) THEN
          WAP=-1.E0
          RETURN
       ENDIF
       XX=X
       DELX=XX-EM
!        IF(DELX.LT.TB2.AND.M.EQ.1) WRITE(NN,70)XX
    ENDIF
    IF(NB.EQ.0) THEN
!
!        Calculations for Wp
!
       IF(ABS(XX).LE.X0) THEN
          WAP=XX/(1.E0+XX/(1.E0+XX/(2.E0+XX/(.6E0+.34E0*XX))))
          RETURN
       ELSE IF(XX.LE.X1) THEN
          RETA=SQRT(E12*DELX)
          WAP=RETA/(1.E0+RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
               AN6+AN5))+AN3)))-1.E0
          RETURN
       ELSE IF(XX.LE.2.E1) THEN
          RETA=S2*SQRT(1.E0-XX/EM)
          AN2=4.612634277343749E0*SQRT(SQRT(RETA+&
               1.09556884765625E0))
          WAP=RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/&
               (S23*(AN2+RETA))))-1.E0
       ELSE
          ZL=ALOG(XX)
          WAP=ALOG(XX/ALOG(XX/ZL**EXP(-1.124491989777808E0/&
               (.4225028202459761E0+ZL))))
       ENDIF
    ELSE
!
!        Calculations for Wm
!
       IF(XX.GE.0.E0) THEN
          WAP = -1./zero
          RETURN
       END IF
       IF(XX.LE.X1) THEN
          RETA=SQRT(E12*DELX)
          WAP=RETA/(RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
               AN6-AN5))-AN3))-1.E0)-1.E0
          RETURN
       ELSE IF(XX.LE.EM9) THEN
          ZL=ALOG(-XX)
          T=-1.E0-ZL
          TS=SQRT(T)
          WAP=ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+&
               TS*127.0471381349219E0))*TS)
       ELSE
          ZL=ALOG(-XX)
          ETA=2.E0-EM2*XX
          WAP=ALOG(XX/ALOG(-XX/((1.E0-.5043921323068457E0*&
               (ZL+1.E0))*(SQRT(ETA)+ETA/3.E0)+1.E0)))
       ENDIF
    ENDIF
!     DO I=1,NITER
       ZN=ALOG(XX/WAP)-WAP
       TEMP=1.E0+WAP
       TEMP2=TEMP+C23*ZN
       TEMP2=2.E0*TEMP*TEMP2
       WAP=WAP*(1.E0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN))
!     END DO
    RETURN
! 40  FORMAT(/,' NBITS is',I4,'.',/,' Expect',I4,&
!          ' significant digits from WAPR.')
! 50  FORMAT(/,' Warning: the offset x is negative (it must be > 0)')
! 60  FORMAT(' Warning: For this offset (',E16.8,'),',/,&
!          ' W is negligibly different from -1')
! 70  FORMAT(' Warning: x (= ',E16.8,') is close to -exp(-1).',/,&
!          ' Enter x as an offset to -exp(-1) for greater accuracy')
!
!     END of WAPR
  END FUNCTION WAPR
!
!     __________________________________________________________________
!
!     Approximating the W function
!     ____________________________
!
  ELEMENTAL FUNCTION WAPD(X,NB,L) RESULT (WAP)
!
!     WAPR - output
!     X - argument of W(X)
!     NB is the branch of the W function needed:
!        NB = 0 - upper branch
!        NB <> 0 - lower branch
!
!     NERROR is the output error flag:
!        NERROR = 0 -> routine completed successfully
!        NERROR = 1 -> X is out of range
!
!     Range: -exp(-1) <= X for the upper branch of the W function
!            -exp(-1) < X < 0 for the lower branch of the W function
!
!     L - determines how WAPR is to treat the argument X
!        L = 1 -> X is the offset from -exp(-1), so compute
!                 W(X-exp(-1))
!        L <> 1 -> X is desired X, so compute W(X)
!
!     M - print messages from WAPR?
!         M = 1 -> Yes
!         M <> 1 -> No
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     NN is the output device number
!
!     NBITS is the number of bits (less 1) in the mantissa of the
!        floating point number number representation of your machine.
!        It is used to determine the level of accuracy to which the W
!        function should be calculated.
!
!        Most machines use a 24-bit matissa for single precision and
!        53-56 bits for double precision. The IEEE standard is 53
!        bits. The Fujitsu VP2200 uses 56 bits. Long word length
!        machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!        single precision.
!
    IMPLICIT NONE
    INTEGER, PARAMETER :: NN=6, NBITS=52
!        Various mathematical constants
    DOUBLE PRECISION, PARAMETER :: EM=-0.367879441171442D0,& ! -EXP(-1.D0)
         EM9=-1.234098040866796D-4,&  ! -EXP(-9.D0)
         C13=1.D0/3.D0,&
         C23=2.D0*C13,&
         EM2=2.D0/EM,&
         D12=-EM2,&
         TB=.5D0**NBITS,&
         TB2=.5D0**(NBITS/2),&        ! SQRT(TB)
         X0=1.23039165028796206D-3,&  ! TB**(1/6)*.5
         X1=0.367668719700312234D0,&  ! (1 - 17*TB**(2/7))*EM
         AN3=8.D0/3.D0,&
         AN4=135.D0/83.D0,&
         AN5=166.D0/39.D0,&
         AN6=3167.D0/3549.D0,&
         S2=1.41421356237310D0,&      ! SQRT(2.D0)
         S21=2.D0*S2 - 3.D0,&
         S22=4.D0 - 3.D0*S2,&
         S23=S2 - 2.D0
    DOUBLE PRECISION :: X, WAP, AN2, DELX, XX, RETA, TS, ZL, ZN, T, ETA, TEMP, TEMP2, zeroD0
    INTEGER :: NB, L
    INTENT(IN) :: X, NB, L
!
!     The following COMMON statement is needed only when testing this
!     function using BISECT, otherwise it can be removed.
!
! COMMON/WAPCOM/NBITS
! DATA INIT,NITER/0,1/
!     DATA NBITS/52/
!
!     IF(INIT.EQ.0) THEN
!        INIT=1
!
!        Code to calculate NBITS for the host machine. NBITS is the
!        mantissa length less one. This value is chosen to allow for
!        rounding error in the final bit. This need only be run once on
!        any particular machine. It can then be included in the above
!        DATA statement.
!
!        DO I=1,2000
!           B=2.D0**(-I)
!           V=1.D0+B
!           IF(V.EQ.1.D0)THEN
!              NBITS=I-1
!              J=-DLOG10(B)
!              IF(M.EQ.1) WRITE(NN,40)NBITS,J
!              EXIT
!           ENDIF
!        END DO
!
!        Remove to here after NBITS has been calculated once
!
!        The case of large NBITS
!
!        IF(NBITS.GE.56) NITER=2
!
!        Various mathematical constants
!
!        EM=-EXP(-1.D0)
!        EM9=-EXP(-9.D0)
!        C13=1.D0/3.D0
!        C23=2.D0*C13
!        EM2=2.D0/EM
!        D12=-EM2
!        TB=.5D0**NBITS
!        TB2=SQRT(TB)
!        X0=TB**(1.D0/6.D0)*.5D0
!        X1=(1.D0-17.D0*TB**(2.D0/7.D0))*EM
!        AN3=8.D0/3.D0
!        AN4=135.D0/83.D0
!        AN5=166.D0/39.D0
!        AN6=3167.D0/3549.D0
!        S2=SQRT(2.D0)
!        S21=2.D0*S2-3.D0
!        S22=4.D0-3.D0*S2
!        S23=S2-2.D0
!     ENDIF
    !OK: same as above, to avoid compilation problems due to division by zero,
    ! set zeroD0=0D0 and use it instead of 0D0
    zeroD0=0D0
    IF(L.EQ.1) THEN
       DELX=X
       IF(DELX.LT.0.D0) THEN
          WAP = 1D0/zeroD0
          RETURN
       END IF
       XX=X+EM
!        IF(D12*DELX.LT.TB**2.AND.M.EQ.1) WRITE(NN,60)DELX
    ELSE
       IF(X.LT.EM) THEN
          WAP = 1D0/zeroD0
          RETURN
       END IF
       IF(X.EQ.EM) THEN
          WAP=-1.D0
          RETURN
       ENDIF
       XX=X
       DELX=XX-EM
!        IF(DELX.LT.TB2.AND.M.EQ.1) WRITE(NN,70)XX
    ENDIF
    IF(NB.EQ.0) THEN
       !
       !        Calculations for Wp
       !
       IF(ABS(XX).LE.X0) THEN
          WAP=XX/(1.D0+XX/(1.D0+XX/(2.D0+XX/(.6D0+.34D0*XX))))
          RETURN
       ELSE IF(XX.LE.X1) THEN
          RETA=DSQRT(D12*DELX)
          WAP=RETA/(1.D0+RETA/(3.D0+RETA/(RETA/(AN4+RETA/(RETA*AN6+AN5))&
               +AN3)))-1.D0
          RETURN
       ELSE IF(XX.LE.2.D1) THEN
          RETA=S2*DSQRT(1.D0-XX/EM)
          AN2=4.612634277343749D0*DSQRT(DSQRT(RETA+1.09556884765625D0))
          WAP=RETA/(1.D0+RETA/(3.D0+(S21*AN2+S22)*RETA/(S23*(AN2+RETA))))-1.D0
       ELSE
          ZL=DLOG(XX)
          WAP=DLOG(XX/DLOG(XX/ZL**DEXP(-1.124491989777808D0/(.4225028202459761D0+ZL))))
       ENDIF
    ELSE
       !
       !        Calculations for Wm
       !
       IF(XX.GE.0.D0) THEN
          WAP = -1D0/zeroD0
          RETURN
       END IF
       IF(XX.LE.X1) THEN
          RETA=DSQRT(D12*DELX)
          WAP=RETA/(RETA/(3.D0+RETA/(RETA/(AN4+RETA/(RETA*AN6-AN5))-AN3))-1.D0)-1.D0
          RETURN
       ELSE IF(XX.LE.EM9) THEN
          ZL=DLOG(-XX)
          T=-1.D0-ZL
          TS=DSQRT(T)
          WAP=ZL-(2.D0*TS)/(S2+(C13-T/(2.7D2+TS*127.0471381349219D0))*TS)
       ELSE
          ZL=DLOG(-XX)
          ETA=2.D0-EM2*XX
          WAP=DLOG(XX/DLOG(-XX/((1.D0-.5043921323068457D0*(ZL+1.D0))&
               *(DSQRT(ETA)+ETA/3.D0)+1.D0)))
       ENDIF
    ENDIF
!     DO I=1,NITER
       ZN=DLOG(XX/WAP)-WAP
       TEMP=1.D0+WAP
       TEMP2=TEMP+C23*ZN
       TEMP2=2.D0*TEMP*TEMP2
       WAP=WAP*(1.D0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.D0*ZN))
!     END DO
    RETURN
! 40  FORMAT(/,' NBITS is',I4,'.',/,' Expect',I4,' significant digits from WAPD.')
! 50  FORMAT(/,' Warning: the offset x is negative (it must be > 0)')
! 60  FORMAT(' Warning: For this offset (',D16.8,'),',/,' W is negligibly different from -1')
! 70  FORMAT(' Warning: x (= ',D16.8,') is close to -exp(-1).',/,&
!          ' Enter x as an offset to -exp(-1) for greater accuracy')
!
!     END of WAPD
  END FUNCTION WAPD
END MODULE lambert
