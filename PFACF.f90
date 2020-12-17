!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    PRESSFCF - PRESS FLUCTUATION AUTO CORRELATION FUNCTION    **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                    **%%
!%%**  LICENSE    LGPL-V3                                                   **%%
!%%**  DATE       DECEMBER 17, 2020                                         **%%
!%%**                                                                       **%% 
!%%**  OBS        THIS PROGRAM DETERMINES THE CORRELATION FUNCTION OF THE   **%%
!%%**             PRESSURE FLUCTUATION AND DIAGONAL COMPONENTS USING THE    **%%
!%%**             TRACE AND DIAGONAL TENSOR PRESSURE COMPONENTS             **%%
!%%**                                                                       **%% 
!%%**             IT IS NECESSARY TO SPECIFIED THE NUMBER OF SAMPLES/OR     **%% 
!%%**             TIME INSTANTS WITHIN THE SAMPLES                          **%% 
!%%**                                                                       **%% 
!%%**             TMAXD MUST BE SET IN CONCORDANCE WITH THE TOBS TIME       **%%
!%%**             DESIRED FOR THE CORRELATION, ALSO, T0MAX MUST NOT BE      **%%
!%%**             GREATER THAN TMAXD                                        **%%
!%%**                                                                       **%%
!%%**             THIS PROGRAM IS BASED ON THE ALGORITHM IN FRENKEL'S       **%%
!%%**             BOOK FOR TIME CORRELATION FUNCTIONS                       **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE PFCVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D      = KIND(1.0D0)

 INTEGER, PARAMETER:: OUTNSAM= 1            !SAMPLES BETWEEN CORRELATION
 INTEGER, PARAMETER:: TMAXD  = 20000        !MAXIMUM CORRELATION TIME
 INTEGER, PARAMETER:: T0MAX  = 20000        !MAXIMUM NUMBER OF t=0
 INTEGER, PARAMETER:: IT0    = 1            !TAKE A NEW t=0 FOR CORRELATION

 REAL(D), PARAMETER:: H      = 0.001D0      !TIME STEP
 REAL(D), PARAMETER:: TOBS   = 2.0D0        !CORRELATION TIME
 
 INTEGER:: I,TMAX,NTEL,T0,TT0
 INTEGER:: TIME0(T0MAX),DELT,T
 INTEGER:: NTIME(TMAXD),NSAMP,K,IDUMM


 REAL(D):: PREXX,PREYY,PREZZ
 REAL(D):: PRESS,PREST,CTIME,DTIME
 REAL(D):: PF0(TMAXD),PFACF(TMAXD)
 REAL(D):: PFXX0(TMAXD),PFACFXX(TMAXD)
 REAL(D):: PFYY0(TMAXD),PFACFYY(TMAXD)
 REAL(D):: PFZZ0(TMAXD),PFACFZZ(TMAXD)
 REAL(D):: SXX,SXY,SXZ
 REAL(D):: SYX,SYY,SYZ
 REAL(D):: SZX,SZY,SZZ
END MODULE PFCVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM PRESSFCF
 USE PFCVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDTen.dat',STATUS='OLD')!STRESS TENSOR COMPONENTS
 OPEN(UNIT=11,FILE='PFACF.dat')             !AUTOCORRELATION PRESSURE FLUCTUATION 
 OPEN(UNIT=12,FILE='PFACFaa.dat')           !AUTOCORRELATION DIAGONAL PRESSURE FLUCTUATION 
 REWIND(10)
 
 NTEL=0
 DTIME=H*DBLE(OUTNSAM)
 TMAX=INT(TOBS/DTIME)
 NSAMP=1500001
 PRESS=0.0D0
 PREXX=0.0D0
 PREYY=0.0D0
 PREZZ=0.0D0

 DO I=1,NSAMP
    READ(10,*)IDUMM,SXX,SXY,SXZ,SYX,SYY,SYZ,SZX,SZY,SZZ
    PRESS=PRESS + (SXX + SYY + SZZ)/3.0D0
    PREXX=PREXX + SXX
    PREYY=PREYY + SYY
    PREZZ=PREZZ + SZZ
 ENDDO

 PRESS=PRESS/DBLE(NSAMP)                    !AVERAGE PRESSURE
 PREXX=PREXX/DBLE(NSAMP)                    !AVERAGE PRESSURE
 PREYY=PREYY/DBLE(NSAMP)                    !AVERAGE PRESSURE
 PREZZ=PREZZ/DBLE(NSAMP)                    !AVERAGE PRESSURE

 REWIND(10)

 DO I=1,TMAX
    NTIME(I)=0
    PFACF(I)=0.0D0

    PFACFXX(I)=0.0D0
    PFACFYY(I)=0.0D0
    PFACFXX(I)=0.0D0
 ENDDO

 DO K=1,NSAMP
    READ(10,*)IDUMM,SXX,SXY,SXZ,SYX,SYY,SYZ,SZX,SZY,SZZ
!
    NTEL=NTEL+1
    IF(MOD(NTEL,IT0) .EQ. 0)THEN
      T0=T0+1
      TT0=MOD(T0-1,T0MAX) + 1
      TIME0(TT0)=NTEL

      PF0(TT0)=(SXX + SYY + SZZ)/3.0D0 - PRESS
      PFXX0(TT0)=SXX - PREXX
      PFYY0(TT0)=SYY - PREYY
      PFZZ0(TT0)=SZZ - PREZZ
    ENDIF
 
    DO T=1,MIN(T0,T0MAX)
       DELT=NTEL -TIME0(T) + 1
       IF(DELT .LT. TMAX)THEN
          NTIME(DELT)=NTIME(DELT) + 1
          PREST=(SXX + SYY + SZZ)/3.0D0

          PFACF(DELT)=PFACF(DELT) + PF0(T)*(PREST - PRESS)

          PFACFXX(DELT)=PFACFXX(DELT) + PFXX0(T)*(SXX - PREXX)
          PFACFYY(DELT)=PFACFYY(DELT) + PFYY0(T)*(SYY - PREYY)
          PFACFZZ(DELT)=PFACFZZ(DELT) + PFZZ0(T)*(SZZ - PREZZ)
       ENDIF
    ENDDO
!
 ENDDO

 DO I=1,TMAX-1
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)
 
    PFACF(I)=PFACF(I)/DBLE(NTIME(I))

    PFACFXX(I)=PFACFXX(I)/DBLE(NTIME(I))
    PFACFYY(I)=PFACFYY(I)/DBLE(NTIME(I))
    PFACFZZ(I)=PFACFZZ(I)/DBLE(NTIME(I))
    WRITE(11,*)CTIME,PFACF(I)
    WRITE(12,*)CTIME,PFACFXX(I),PFACFYY(I),PFACFZZ(I)
 ENDDO

 CLOSE(UNIT=10)
 CLOSE(UNIT=11)
 CLOSE(UNIT=12)

 STOP
END PROGRAM PRESSFCF
