!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    VACF - VELOCITY AUTOCORRELATION FUNCTION                 **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                   **%%
!%%**  LICENSE    LGPL-V3                                                  **%%
!%%**  DATE       NOVEMBER 27, 2020                                        **%%
!%%**                                                                      **%%
!%%**  OBS        THIS PROGRAM DETERMINES THE VACF IN THREE DIMENSIONS.    **%%
!%%**             THE PROGRAM REQUIERES AND INPUT FILE WITH THE MOLECULES  **%%
!%%**             CARTESIAN CORRDINATES MDVel.dat                          **%%
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
MODULE VACVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D      = KIND(1.0D0)  !VARIABLES PRECISION 
 INTEGER, PARAMETER:: NPARTX = 3000         !MAXIMUM NUMBER OF PARTICLES 

 INTEGER, PARAMETER:: OUTNSAM= 100          !SAMPLES BETWEEN CORRELATION
 INTEGER, PARAMETER:: TMAXD  = 20000        !MAXIMUM CORRELATION TIME
 INTEGER, PARAMETER:: T0MAX  = 20000        !MAXIMUM NUMBER OF t=0
 INTEGER, PARAMETER:: IT0    = 100          !TAKE A NEW t=0 FOR CORRELATION

 REAL(D), PARAMETER:: H      = 0.001D0      !TIME STEP
 REAL(D), PARAMETER:: TOBS   = 1000.0D0     !CORRELATION TIME
 
 INTEGER:: I,TMAX,NTEL,T0,TT0,NPART
 INTEGER:: TIME0(T0MAX),DELT,T,ISAM
 INTEGER:: NTIME(TMAXD),NSAMP,K

 REAL(D):: VX(NPARTX),VY(NPARTX),VZ(NPARTX) 
 REAL(D):: VX0(NPARTX,TMAXD),VY0(NPARTX,TMAXD)
 REAL(D):: VZ0(NPARTX,TMAXD),CTIME,DTIME,VACFT
 REAL(D):: VACFX(TMAXD),VACFY(TMAXD),VACFZ(TMAXD)
END MODULE VACVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM VACF
 USE VACVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDVel.dat',STATUS='OLD')!FILE WITH COORDINATES
 OPEN(UNIT=11,FILE='VACF.dat')              !FILE WITH THE MSD COMPONENTS
 REWIND(10)
 
 NTEL=0                                     !NUMBER OF t=0 
 DTIME=H*DBLE(OUTNSAM)                      !TIME CORRELATION
 TMAX=INT(TOBS/DTIME)                       !MAXIMUM TIME INSTANTS OBSERVED
 NSAMP=40001                                !NUMBER OF SAMPLES
 NPART=2048                                 !NUMBER OF PARTICLES

 DO I=1,TMAX                                !ARRAYS TIME INITIALIZATION 
    NTIME(I)=0

    VACFX(I)=0.0D0
    VACFY(I)=0.0D0
    VACFZ(I)=0.0D0
 ENDDO

 DO K=1,NSAMP                               !STATISTICS OVER THE TOTAL SAMPLES
    READ(10,*)ISAM                          !NUMBER OF SAMPLE
    DO I=1,NPART
       READ(10,*)VX(I),VY(I),VZ(I)          !POSITION AT TIME t_{i}
    ENDDO  
!
    NTEL=NTEL+1                             !UPDATE NUMBER OF t=0
    IF(MOD(NTEL,IT0) .EQ. 0)THEN
      T0=T0+1                               !n-th t=0 
      TT0=MOD(T0-1,T0MAX) + 1               !RESET t=0 IF MAXIMUM NUMBER IS REACHED
      TIME0(TT0)=NTEL

      DO I=1,NPART                          !SET t=0 REFERENCE POSITIONS
         VX0(I,TT0)=VX(I)
         VY0(I,TT0)=VY(I)
         VZ0(I,TT0)=VZ(I)
      ENDDO
    ENDIF
 
    DO T=1,MIN(T0,T0MAX)                    !PERFORMS CORRELATION UP TO T0
       DELT=NTEL -TIME0(T) + 1              !OBSERVED t_{i}
       IF(DELT .LT. TMAX)THEN
          NTIME(DELT)=NTIME(DELT) + 1
          DO I=1,NPART                      !CORRELATION WITH POSITIONS AT TIME t_{i}
             VACFX(DELT)=VACFX(DELT) + VX(I)*VX0(I,T)
             VACFY(DELT)=VACFY(DELT) + VY(I)*VY0(I,T)
             VACFZ(DELT)=VACFZ(DELT) + VZ(I)*VZ0(I,T)
          ENDDO
       ENDIF
    ENDDO
!
 ENDDO

 DO I=1,TMAX                                !AVERAGE WITH PARTICLE NUMBER AND t=0'S 
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)
    VACFX(I)=VACFX(I)/(DBLE(NPART*NTIME(I)))
    VACFY(I)=VACFY(I)/(DBLE(NPART*NTIME(I)))
    VACFZ(I)=VACFZ(I)/(DBLE(NPART*NTIME(I)))
    VACFT=VACFX(I) + VACFY(I) + VACFZ(I)
 
    WRITE(11,*)CTIME,VACFX(I),VACFY(I),VACFZ(I),VACFT
 ENDDO

 CLOSE(UNIT=10)
 CLOSE(UNIT=11)

 STOP
END PROGRAM VACF
