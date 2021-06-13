!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM     DYNAMICAL ARREST                                        **%%
!%%**  AUTHOR      ALEXIS TORRES CARBAJAL                                  **%%
!%%**  LICENSE     LGPL-V3                                                 **%%
!%%**  DATE        JUNE 13, 2021                                           **%%
!%%**                                                                      **%%
!%%**  OBS         THIS PROGRAM DETERMINES THE DYNAMIC ARREST OF A SYSTEM  **%%
!%%**              USING ONLY ITS STATIC STRUCTURE FACTOR, BASED IN        **%%
!%%**              PHYS. REV. E; 76, 062502 (2007)                         **%%
!%%**              (2007)                                                  **%%
!%%**                                                                      **%%
!%%**              THE PROGRAM REQUERIES THE FILE WITH S(k) AND THE NUMBER **%%
!%%**              OF DATA IN THE FILE AS WELL AS THE TEMPERATURE AND THE  **%%
!%%**              VOLUME FRACTION MUST BE SPECIFIED BEFORE TO COMPILE     **%%
!%%**                                                                      **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE SKVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D     = KIND(1.0D0)   !VARIABLES PRECISION 
 INTEGER, PARAMETER:: NDATX = 1000          !MAX NUMBER OF DATA

 REAL(D), PARAMETER:: PI    = DACOS(-1.0D0) !PI NUMBER
 REAL(D), PARAMETER:: KC    = 7.18D0        !CRITICAL WAVE NUMBER
 REAL(D), PARAMETER:: ERROR = 1.0D-6        !TOLERATION, CONVERGENCE 

 REAL(D), PARAMETER:: DGMMA = 0.000001D0    !GAMMA INCREMENT
 REAL(D), PARAMETER:: GAMMAM= 0.0D0         !MINIMUM GAMMA
 REAL(D), PARAMETER:: GAMMAX= 0.1D0         !MAXIMUM GAMMA
 REAL(D), PARAMETER:: DGAMMA= 0.00001D0     !GAMMA INCREMENT

 LOGICAL:: CONV
 INTEGER:: I,IT,NDAT,NG
 REAL(D):: K(NDATX),SK(NDATX)
 REAL(D):: LAMBDAK,GAMMA,DYNARR
 REAL(D):: GAMMAI,GAMMAF,TEM

 REAL(D):: XI0,XI1,XI2,FX,DX
 REAL(D):: NUM,DEN,INTE,PHI
END MODULE SKVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM ARREST
 USE SKVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDSk.dat',STATUS='OLD') !FILE WITH S(k)
 REWIND(10)

 NDAT=215                                   !NUMBER OF DATA IN FILE
 DO I=1,NDAT                                !READ DATA FROM FILE
    READ(10,*)K(I),SK(I)
 ENDDO

 GAMMA=DGAMMA                               !INITIAL GAMMA
 PHI=0.10D0                                 !VOLUME FRACTION
 TEM=0.45D0                                 !TEMPERATURE
 DYNARR=1.0D0
 GAMMAI=DGAMMA

 DO WHILE( DYNARR .GT. ERROR .AND. GAMMAI .LT. 1.D12)
    !COMPOSITE SIMPSON 1/3 INTEGRATION RULE
    LAMBDAK=1.0D0 + (K(1)/KC)**2            !INTERPOLATING FUNCTION
    LAMBDAK=1.0D0/LAMBDAK                   !INVERSE INTERPOLATING FUNCTION
    NUM=(LAMBDAK*(SK(1) - 1.0D0))**2        !NUMERATOR 
    DEN=(LAMBDAK + GAMMA*K(1)**2)*(LAMBDAK*SK(1) + GAMMA*K(1)**2)
    FX=(NUM*K(1)**4)/DEN
    XI0=FX

    LAMBDAK=1.0D0 + (K(NDAT)/KC)**2
    LAMBDAK=1.0D0/LAMBDAK
    NUM=(LAMBDAK*(SK(NDAT) - 1.0D0))**2     !NUMERATOR 
    DEN=(LAMBDAK + GAMMA*K(NDAT)**2)*(LAMBDAK*SK(NDAT) + GAMMA*K(NDAT)**2)
    FX=(NUM*K(NDAT)**4)/DEN
    XI0=XI0 + FX

    XI1=0.0D0
    XI2=0.0D0
    DX=K(2)-K(1)

    DO I=2,NDAT-1
       LAMBDAK=1.0D0 + (K(I)/KC)**2
       LAMBDAK=1.0D0/LAMBDAK
       NUM=((SK(I) - 1.0D0)*LAMBDAK)**2
       DEN=(LAMBDAK*SK(I) + GAMMA*K(I)**2)*(LAMBDAK + GAMMA*K(I)**2)
       FX=(NUM*K(I)**4)/DEN
       IF(MOD(I,2) .EQ. 0)THEN
         XI1=XI1 + FX
       ELSE
         XI2=XI2 + FX
       ENDIF
    ENDDO

    INTE=DX*(XI0 + 2.0D0*XI2 + 4.0D0*XI1)/3.0D0
    GAMMA=36.0D0*PI*PHI/INTE                !GAMMA NUMBER
    GAMMAF=GAMMA

    DYNARR=ABS((GAMMAF-GAMMAI)/GAMMAF)
    GAMMAI=GAMMA
 ENDDO

 IF(GAMMA .LE. 1)THEN                       !IF GAMMA^{-1} << 1 IS NOT ARRESTED
   IT=1                                     !IF THE SYSTEM IS ARRESTED
 ELSE
   IT=0                                     !IF THE SYSTEM IS NOT ARRESTED
 ENDIF

 WRITE(6,101)PHI,TEM,GAMMA,IT

 CLOSE(UNIT=10)
 101 FORMAT(2(F16.8),G16.4,I10)
 STOP
END PROGRAM ARREST
