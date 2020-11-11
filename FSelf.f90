!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    FSELF - (SELF PART) INTERMEDIATE SCATTERING FUNCTION     **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                   **%%
!%%**  LICENSE    LGPL-V3                                                  **%%
!%%**  DATE       NOVEMBER 12, 2020                                        **%%
!%%**                                                                      **%%
!%%**  OBS        THIS PROGRAM DETERMINES THE SELF PART OF THE INTERMEDIATE**%%
!%%**             SCATTERING FUNCTION F_s(k,t) IN THREE DIMENSIONS.        **%%
!%%**                                                                      **%%
!%%**             THE PROGRAM REQUIERES AND INPUT FILE WITH THE MOLECULES  **%%
!%%**             CARTESIAN CORRDINATES MDFps.dat                          **%%
!%%**                                                                      **%%
!%%**             IT IS NECESSARY TO SPECIEFIED THE NUMBER OF MOLECULES    **%%
!%%**             AND THE TOTAL NUMBER OF SAMPLES IN THE COORDINATES FILE  **%%
!%%**                                                                      **%%
!%%**             TMAXD MUST BE SET IN CONCORDANCE WITH THE TOBS TIME      **%%
!%%**             DESIRED FOR THE CORRELATION, ALSO, T0MAX MUST NOT BE     **%%
!%%**             GREATER THAN TMAXD                                       **%%
!%%**                                                                      **%%
!%%**             THIS PROGRAM IS BASED ON THE ALGORITHM IN FRENKEL'S      **%%
!%%**             BOOK                                                     **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE VARFSE
 IMPLICIT NONE
 INTEGER, PARAMETER:: D       = KIND(1.0D0) !VARIABLES PRECISION
 INTEGER, PARAMETER:: NPARTX  = 3000        !MAXIMUM NUMBER OF PARTICLES

 INTEGER, PARAMETER:: OUTNSAM = 100         !SAMPLES BETWEEN CORRELATION 
 INTEGER, PARAMETER:: TMAXD   = 10000       !MAXIMUM CORRELATION TIME
 INTEGER, PARAMETER:: T0MAX   = 10000       !MAXIMUM NUMBER OF t=0
 INTEGER, PARAMETER:: IT0     = 100         !TAKE A NEW t=0 FOR CORRELATION

 REAL(D), PARAMETER:: H       = 0.001D0     !TIME STEP
 REAL(D), PARAMETER:: TOBS    = 1000.0D0    !CORRELATION TIME
 REAL(D), PARAMETER:: PI      =DACOS(-1.0D0)!PI NUMBER

 REAL(D), PARAMETER:: KMAX1   = 2.0d0*PI    !WAVE NUMBER OBSERVED

 INTEGER:: I,K,T,NPART,T0,TT0,DELT
 INTEGER:: NTEL,TMAX,NTIME(TMAXD)
 INTEGER:: TIME0(T0MAX),ISAM,NSAMP

 REAL(D):: RX(NPARTX),RX0(NPARTX,TMAXD)     !POSITIONS AT t_{i} AND t=0
 REAL(D):: RY(NPARTX),RY0(NPARTX,TMAXD)
 REAL(D):: RZ(NPARTX),RZ0(NPARTX,TMAXD)
 REAL(D):: X,DRX,DRY,DRZ,DR,DTIME,CTIME
 REAL(D):: FSELF1(TMAXD)
END MODULE VARFSE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM FSELF
 USE VARFSE
 IMPLICIT NONE 

 OPEN(UNIT=10,FILE='MDFps.dat',STATUS='OLD')!FILE WITH COORDINATES 
 OPEN(UNIT=11,FILE='FSelf.dat')             !FILE WITH FSELF AT WAVE NUMBER KMAX
 REWIND(10)
 
 NTEL=0                                     !NUMBER OF t=0
 DTIME=H*DBLE(OUTNSAM)                      !TIME CORRELATION
 TMAX=INT(TOBS/DTIME)                       !MAXIMUM TIME INSTANTS OBSERVED
 NSAMP=40001                                !NUMBER OF SAMPLES
 NPART=2048                                 !NUMBER OF PARTICLES
 
 DO I=1,TMAX                                !ARRAYS TIME INITIALIZATION
    NTIME(I)=0

    FSELF1(I)=0.0D0
 ENDDO
 
 DO K=1,NSAMP                               !STATISTICS OVER THE TOTAL SAMPLES
    READ(10,*)ISAM                          !NUMBER OF SAMPLES
    DO I=1,NPART
       READ(10,*)RX(I),RY(I),RZ(I)          !POSITION AT TIME t_{i}
    ENDDO
!
    NTEL=NTEL+1                             !UPDATE NUMBER OF t=0
    IF(MOD(NTEL,IT0) .EQ. 0)THEN            !SELECT A t=0
      T0=T0+1                               !n-th t=0
      TT0=MOD(T0-1,T0MAX) + 1               !RESET t=0 IF MAXIMUM NUMBER IS REACHED
      TIME0(TT0)=NTEL
 
      DO I=1,NPART                          !SET t=0 REFERENCE POSITIONS
         RX0(I,TT0)=RX(I)
         RY0(I,TT0)=RY(I)
         RZ0(I,TT0)=RZ(I)
      ENDDO
    ENDIF

    DO T=1,MIN(T0,T0MAX)                    !PERFORMS CORRELATION UP TO T0
       DELT=NTEL - TIME0(T) + 1             !OBSERVED t_{i}
       IF(DELT .LT. TMAX)THEN
         NTIME(DELT)=NTIME(DELT) + 1
         DO I=1,NPART                       !CORRELATION WITH POSITIONS AT TIME t_{i}
            DRX=(RX(I)-RX0(I,T))**2
            DRY=(RY(I)-RY0(I,T))**2
            DRZ=(RZ(I)-RZ0(I,T))**2

            DR=DSQRT(DRX + DRY + DRZ)

            X=KMAX1*DR
            FSELF1(DELT)=FSELF1(DELT) + DSIN(X)/X
         ENDDO
       ENDIF
    ENDDO
 ENDDO

 DO I=1,TMAX                                !AVERAGE WITH PARTICLE NUMBER AND t=0'S
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)
    FSELF1(I)=FSELF1(I)/(DBLE(NPART*NTIME(I)))

    WRITE(11,*)CTIME,FSELF1(I)
 ENDDO

 CLOSE(10)
 CLOSE(11)

 STOP
END PROGRAM FSELF
