!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    BVISCOCITY - BULK VISCOCITY FROM THE GREEN-KUBO RELATION  **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                    **%%
!%%**  LICENSE    LGPL-V3                                                   **%%
!%%**  DATE       DECEMBER 17, 2020                                         **%%
!%%**                                                                       **%% 
!%%**  OBS        THIS PROGRAM DETERMINES THE BULK AND LONGITUDINAL         **%% 
!%%**             VISCOCITY FROM THE AUTOCORRELATION FUNCTION OF THE        **%%
!%%**             PRESSURE FLUCTUATION TROUGH THE GREEN-KUBO RELATION       **%%
!%%**             AUTOCORRELATION FUNCTION OF THE PRESSURE FLUCTUATION      **%%
!%%**             THE AUTOCORRELATION FUNCTIONS ARE READED FROM THE FILES   **%%
!%%**             PFACF.dat AND PFACFaa.dat, RESPECTIVELY                   **%%
!%%**                                                                       **%% 
!%%**             IT IS NECESSARY TO SPECIFIED THE NUMBER OF DATA (NDAT)    **%% 
!%%**             RELATED WITH THE TIME INSTANTS IN THE CORRELATION, THE    **%% 
!%%**             NUMBER OF PARTICLES AND THE VOLUME FRACTION               **%%
!%%**                                                                       **%%
!%%**             THE GREEN-KUBO INTEGRAL IS PERFORMED WITH THE SIMPLE      **%%
!%%**             SIMPSON RULE                                              **%%
!%%**                                                                       **%%
!%%**             THE AVERAGE IS ENHACED AVERAGING  ALL THE DIAGONAL        **%%
!%%**             CORRELATIONS                                              **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE VISVAR 
 IMPLICIT NONE
 INTEGER, PARAMETER:: D     = KIND(1.0D0)   !VARIABLES PRECISION
 INTEGER, PARAMETER:: NDATX = 5000          !MAXIMUM NUMBER OF DATA
 
 REAL(D), PARAMETER:: PI    = DACOS(-1.0D0) !PI NUMBER

 INTEGER:: NPART,I,NDAT

 REAL(D):: TIME(NDATX),PHI,RHOSTAR
 REAL(D):: DX,XI0,XI1,XI2,INTX,BVIS
 REAL(D):: TSTAR,PFACF(NDATX),LVIS
 REAL(D):: PFACFXX(NDATX),PFACFYY(NDATX)
 REAL(D):: PFACFZZ(NDATX)

 REAL(D):: XXI0,XXI1,XXI2,INTXX
 REAL(D):: YYI0,YYI1,YYI2,INTYY
 REAL(D):: ZZI0,ZZI1,ZZI2,INTZZ
END MODULE VISVAR 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM BVISCOCITY
 USE VISVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='PFACF.dat',STATUS='OLD')
 OPEN(UNIT=11,FILE='PFACFaa.dat',STATUS='OLD')
 REWIND(10)
 REWIND(11)

 NDAT=1999                                  !NUMBER OF TIME INSTANT 
 TSTAR=1.4737D0                             !SYSTEM TEMPERATURE
 NPART=500                                  !SYSTEM NUMBER OF PARTICLES
 PHI=0.10D0                                 !VOLUME FRACTION
 RHOSTAR=6.0D0*PHI/PI                       !SYSTEM NUMBER DENSITY

 DO I=1,NDAT                                !READ TIME CORRELATION COMPONENTS
    READ(10,*)TIME(I),PFACF(I)
    READ(11,*)TIME(I),PFACFXX(I),PFACFYY(I),PFACFZZ(I)
 ENDDO

 DX=TIME(2)-TIME(1)
 XI0=PFACF(1) + PFACF(NDAT)                 !FOR PRESSURE FLUCTUATION INTEGRAL
 XI1=0.0D0
 XI2=0.0D0


 XXI0=PFACFXX(1) + PFACFXX(NDAT)            !FOR XX PRESSURE FLUCTUATION INTEGRAL
 XXI1=0.0D0
 XXI2=0.0D0


 YYI0=PFACFYY(1) + PFACFYY(NDAT)            !FOR YY PRESSURE FLUCTUATION INTEGRAL
 YYI1=0.0D0
 YYI2=0.0D0


 ZZI0=PFACFZZ(1) + PFACFZZ(NDAT)            !FOR ZZ PRESSURE FLUCTUATION INTEGRAL
 ZZI1=0.0D0
 ZZI2=0.0D0

 DO I=2,NDAT-2
    IF(MOD(I,2) .EQ. 0)THEN
      XI1=XI1 + PFACF(I)

      XXI1=XXI1 + PFACFXX(I)
      YYI1=YYI1 + PFACFYY(I)
      ZZI1=ZZI1 + PFACFZZ(I)
    ELSE
      XI2=XI2 + PFACF(I)

      XXI2=XXI2 + PFACFXX(I)
      YYI2=YYI2 + PFACFYY(I)
      ZZI2=ZZI2 + PFACFZZ(I)
    ENDIF
 ENDDO

 INTX=DX*(XI0 + 2.0D0*XI1 + 4.0D0*XI2)/3.0D0

 INTXX=DX*(XXI0 + 2.0D0*XXI1 + 4.0D0*XXI2)/3.0D0
 INTYY=DX*(YYI0 + 2.0D0*YYI1 + 4.0D0*YYI2)/3.0D0
 INTZZ=DX*(ZZI0 + 2.0D0*ZZI1 + 4.0D0*ZZI2)/3.0D0


 BVIS=DBLE(NPART)*INTX/(RHOSTAR*TSTAR)
 LVIS=DBLE(NPART)*(INTXX + INTYY + INTZZ)/(3.0D0*RHOSTAR*TSTAR)
 WRITE(6,*)PHI,BVIS,LVIS

 CLOSE(UNIT=10)
 CLOSE(UNIT=11)

 STOP
END PROGRAM BVISCOCITY
