!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    MSD - MEAN SQUARE DISPLACEMENT                           **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                   **%%
!%%**  LICENSE    LGPL-V3                                                  **%%
!%%**  DATE       NOVEMBER 12, 2020                                        **%%
!%%**                                                                      **%%
!%%**  OBS        THIS PROGRAM DETERMINES THE MSD IN THREE DIMENSIONS. THE **%%
!%%**             PROGRAM REQUIERES AND INPUT FILE WITH THE MOLECULES      **%%
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
MODULE MSDVAR
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

 REAL(D):: RX(NPARTX),RY(NPARTX),RZ(NPARTX) 
 REAL(D):: RX0(NPARTX,TMAXD),RY0(NPARTX,TMAXD)
 REAL(D):: RZ0(NPARTX,TMAXD),CTIME,DTIME,MSDT
 REAL(D):: MSDX(TMAXD),MSDY(TMAXD),MSDZ(TMAXD)
END MODULE MSDVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM MSD
 USE MSDVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDFps.dat',STATUS='OLD')!FILE WITH COORDINATES
 OPEN(UNIT=11,FILE='MSD.dat')               !FILE WITH THE MSD COMPONENTS
 REWIND(10)
 
 NTEL=0                                     !NUMBER OF t=0 
 DTIME=H*DBLE(OUTNSAM)                      !TIME CORRELATION
 TMAX=INT(TOBS/DTIME)                       !MAXIMUM TIME INSTANTS OBSERVED
 NSAMP=40001                                !NUMBER OF SAMPLES
 NPART=2048                                 !NUMBER OF MOLECULES

 DO I=1,TMAX                                !ARRAYS TIME INITIALIZATION 
    NTIME(I)=0

    MSDX(I)=0.0D0
    MSDY(I)=0.0D0
    MSDZ(I)=0.0D0
 ENDDO

 DO K=1,NSAMP                               !STATISTICS AROUN THE TOTAL SAMPLES
    READ(10,*)ISAM                          !NUMBER OF SAMPLE
    DO I=1,NPART
       READ(10,*)RX(I),RY(I),RZ(I)          !POSITION AT TIME t_{i}
    ENDDO  
!
    NTEL=NTEL+1                             !UPDATE NUMBER OF t=0
    IF(MOD(NTEL,IT0) .EQ. 0)THEN
      T0=T0+1                               !n-TH t=0 
      TT0=MOD(T0-1,T0MAX) + 1               !RESET t=0 IF MAXIMUM NUMBER IS REACHED
      TIME0(TT0)=NTEL

      DO I=1,NPART                          !SET t=0 REFERENCE POSITIONS
         RX0(I,TT0)=RX(I)
         RY0(I,TT0)=RY(I)
         RZ0(I,TT0)=RZ(I)
      ENDDO
    ENDIF
 
    DO T=1,MIN(T0,T0MAX)                    !PERFORMS CORRELATION UP TO T0
       DELT=NTEL -TIME0(T) + 1              !OBSERVED t_{I}
       IF(DELT .LT. TMAX)THEN
          NTIME(DELT)=NTIME(DELT) + 1
          DO I=1,NPART                      !CORRELATION WITH POSITIONS AT TIME t_{i}
             MSDX(DELT)=MSDX(DELT) + (RX(I) - RX0(I,T))**2
             MSDY(DELT)=MSDY(DELT) + (RY(I) - RY0(I,T))**2
             MSDZ(DELT)=MSDZ(DELT) + (RZ(I) - RZ0(I,T))**2
          ENDDO
       ENDIF
    ENDDO
!
 ENDDO

 DO I=1,TMAX                                !AVERAGE WITH PARTICLE NUMBER AND t=0'S 
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)
    MSDX(I)=MSDX(I)/(DBLE(NPART*NTIME(I)))
    MSDY(I)=MSDY(I)/(DBLE(NPART*NTIME(I)))
    MSDZ(I)=MSDZ(I)/(DBLE(NPART*NTIME(I)))
    MSDT=MSDX(I) + MSDY(I) + MSDZ(I)
 
    WRITE(11,*)CTIME,MSDX(I),MSDY(I),MSDZ(I),MSDT
 ENDDO

 CLOSE(UNIT=10)
 CLOSE(UNIT=11)

 STOP
END PROGRAM MSD
