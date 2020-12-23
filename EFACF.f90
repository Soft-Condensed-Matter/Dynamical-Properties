!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    EFACF - ENERGY FLUX AUTO CORRELATION FUNCTION             **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                    **%%
!%%**  LICENSE    LGPL-V3                                                   **%%
!%%**  DATE       DECEMBER 17, 2020                                         **%%
!%%**                                                                       **%% 
!%%**  OBS        THIS PROGRAM DETERMINES THE CORRELATION FUNCTION OF THE   **%%
!%%**             ENERGY FLUX                                               **%%
!%%**                                                                       **%% 
!%%**             IT IS NECESSARY TO SPECIFIED THE NUMBER OF SAMPLES/OR     **%% 
!%%**             TIME INSTANTS WITHIN THE SAMPLE                           **%% 
!%%**                                                                       **%% 
!%%**             TMAXD MUST BE SET IN CONCORDANCE WITH THE TOBS TIME       **%%
!%%**             DESIRED FOR THE CORRELATION, ALSO, T0MAX MUST NOT BE      **%%
!%%**             GREATER THAN TMAXD                                        **%%
!%%**                                                                       **%%
!%%**             THIS PROGRAM IS BASED ON THE ALGORITHM IN FRENKEL'S       **%%
!%%**             BOOK FOR TIME CORRELATION FUNCTIONS                       **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE EFCVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D      = KIND(1.0D0)

 INTEGER, PARAMETER:: OUTNSAM= 1            !SAMPLES BETWEEN CORRELATION
 INTEGER, PARAMETER:: TMAXD  = 20000        !MAXIMUM CORRELATION TIME
 INTEGER, PARAMETER:: T0MAX  = 20000        !MAXIMUM NUMBER OF t=0
 INTEGER, PARAMETER:: IT0    = 1            !TAKE A NEW t=0 FOR CORRELATION

 REAL(D), PARAMETER:: H      = 0.001D0      !TIME STEP
 REAL(D), PARAMETER:: TOBS   = 5.0D0        !CORRELATION TIME
 
 INTEGER:: I,TMAX,NTEL,T0,TT0
 INTEGER:: TIME0(T0MAX),DELT,T
 INTEGER:: NTIME(TMAXD),NSAMP,K,IDUMM

 REAL(D):: EFX,EFY,EFZ,CTIME,DTIME,EFACT
 REAL(D):: EFX0(TMAXD),EFY0(TMAXD),EFZ0(TMAXD)
 REAL(D):: EFACX(TMAXD),EFACY(TMAXD),EFACZ(TMAXD)
END MODULE EFCVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM EFACF
 USE EFCVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDEfx.dat',STATUS='OLD')!ENERGY FLUX COMPONENTS
 OPEN(UNIT=11,FILE='EFACF.dat')             !ENERGY FLUX AUTOCORRELATION
 REWIND(10)
 
 NTEL=0
 DTIME=H*DBLE(OUTNSAM)
 TMAX=INT(TOBS/DTIME)
 NSAMP=1500001

 DO I=1,TMAX
    NTIME(I)=0

    EFACX(I)=0.0D0
    EFACY(I)=0.0D0
    EFACZ(I)=0.0D0
 ENDDO

 DO K=1,NSAMP
    READ(10,*)IDUMM,EFX,EFY,EFZ
!
    NTEL=NTEL+1
    IF(MOD(NTEL,IT0) .EQ. 0)THEN
      T0=T0+1
      TT0=MOD(T0-1,T0MAX) + 1
      TIME0(TT0)=NTEL

      EFX0(TT0)=EFX
      EFY0(TT0)=EFY
      EFZ0(TT0)=EFZ
    ENDIF
 
    DO T=1,MIN(T0,T0MAX)
       DELT=NTEL -TIME0(T) + 1
       IF(DELT .LT. TMAX)THEN
          NTIME(DELT)=NTIME(DELT) + 1

          EFACX(DELT)=EFACX(DELT) + EFX0(T)*EFX
          EFACY(DELT)=EFACY(DELT) + EFY0(T)*EFY
          EFACZ(DELT)=EFACZ(DELT) + EFZ0(T)*EFZ          
       ENDIF
    ENDDO
!
 ENDDO

 DO I=1,TMAX-1
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)

    EFACX(I)=EFACX(I)/DBLE(NTIME(I))
    EFACY(I)=EFACY(I)/DBLE(NTIME(I))
    EFACZ(I)=EFACZ(I)/DBLE(NTIME(I))
    EFACT=EFACX(I) + EFACY(I) + EFACZ(I)
    WRITE(11,*)CTIME,EFACX(I),EFACY(I),EFACZ(I),EFACT
 ENDDO

 CLOSE(UNIT=10)
 CLOSE(UNIT=11)

 STOP
END PROGRAM EFACF
