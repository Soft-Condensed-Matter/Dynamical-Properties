!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%**  PROGRAM    SSCF - SHEAR STRESS AUTO CORRELATION FUNCTION             **%%
!%%**  AUTHOR     ALEXIS TORRES CARBAJAL                                    **%%
!%%**  LICENSE    LGPL-V3                                                   **%%
!%%**  DATE       DECEMBER 16, 2020                                         **%%
!%%**                                                                       **%% 
!%%**  OBS        THIS PROGRAM DETERMINES THE CORRELATION FUNCTIONS OF THE  **%%
!%%**             STRESS TENSOR COMPONENTS AND THEIR KYNETIC, POTENTIAL &   **%%
!%%**             KYNETIC-POTENTIAL CONTRIBUTIONS. THE PROGRAM REQUIERES    **%%
!%%**             THE INPUT FILES MDTen.dat, MDKtn.dat AND MDPtn.dat WITH   **%%
!%%**             THE TOTAL, KYNETIC AND POTENTIAL COMPONENTS OF THE STRESS **%%
!%%**             TENSOR, RESPECTIVELY                                      **%%
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
MODULE SSCFVAR
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

 REAL(D):: SXX,SXY,SXZ
 REAL(D):: SYX,SYY,SYZ
 REAL(D):: SZX,SZY,SZZ

 REAL(D):: KSXX,KSXY,KSXZ
 REAL(D):: KSYX,KSYY,KSYZ
 REAL(D):: KSZX,KSZY,KSZZ

 REAL(D):: PSXX,PSXY,PSXZ
 REAL(D):: PSYX,PSYY,PSYZ
 REAL(D):: PSZX,PSZY,PSZZ

 REAL(D):: SXX0(TMAXD),SXY0(TMAXD),SXZ0(TMAXD)
 REAL(D):: SYX0(TMAXD),SYY0(TMAXD),SYZ0(TMAXD)
 REAL(D):: SZX0(TMAXD),SZY0(TMAXD),SZZ0(TMAXD)

 REAL(D):: KSXX0(TMAXD),KSXY0(TMAXD),KSXZ0(TMAXD)
 REAL(D):: KSYX0(TMAXD),KSYY0(TMAXD),KSYZ0(TMAXD)
 REAL(D):: KSZX0(TMAXD),KSZY0(TMAXD),KSZZ0(TMAXD)

 REAL(D):: PSXX0(TMAXD),PSXY0(TMAXD),PSXZ0(TMAXD)
 REAL(D):: PSYX0(TMAXD),PSYY0(TMAXD),PSYZ0(TMAXD)
 REAL(D):: PSZX0(TMAXD),PSZY0(TMAXD),PSZZ0(TMAXD)

 REAL(D):: CTIME,DTIME,MSDT
 REAL(D):: SSCFXX(TMAXD),SSCFXY(TMAXD),SSCFXZ(TMAXD)
 REAL(D):: SSCFYX(TMAXD),SSCFYY(TMAXD),SSCFYZ(TMAXD)
 REAL(D):: SSCFZX(TMAXD),SSCFZY(TMAXD),SSCFZZ(TMAXD)

 REAL(D):: KSSCFXX(TMAXD),KSSCFXY(TMAXD),KSSCFXZ(TMAXD)
 REAL(D):: KSSCFYX(TMAXD),KSSCFYY(TMAXD),KSSCFYZ(TMAXD)
 REAL(D):: KSSCFZX(TMAXD),KSSCFZY(TMAXD),KSSCFZZ(TMAXD)

 REAL(D):: PSSCFXX(TMAXD),PSSCFXY(TMAXD),PSSCFXZ(TMAXD)
 REAL(D):: PSSCFYX(TMAXD),PSSCFYY(TMAXD),PSSCFYZ(TMAXD)
 REAL(D):: PSSCFZX(TMAXD),PSSCFZY(TMAXD),PSSCFZZ(TMAXD)

 REAL(D):: KPSSCFXX(TMAXD),KPSSCFXY(TMAXD),KPSSCFXZ(TMAXD)
 REAL(D):: KPSSCFYX(TMAXD),KPSSCFYY(TMAXD),KPSSCFYZ(TMAXD)
 REAL(D):: KPSSCFZX(TMAXD),KPSSCFZY(TMAXD),KPSSCFZZ(TMAXD)
END MODULE SSCFVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM SSCF
 USE SSCFVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE='MDTen.dat',STATUS='OLD')!STRESS TENSOR COMPONENTS
 OPEN(UNIT=11,FILE='MDKtn.dat',STATUS='OLD')!STRESS TENSOR KYNETIC COMPONENTS
 OPEN(UNIT=12,FILE='MDPtn.dat',STATUS='OLD')!STRESS TENSOR POTENTIAL COMPONENTS
 OPEN(UNIT=13,FILE='Phi010TSSCF.dat')       !AUTOCORRELATION STRESS TENSOR
 OPEN(UNIT=14,FILE='Phi010KSSCF.dat')       !AUTOCORRELATION KYNETIC STRESS TENSOR
 OPEN(UNIT=15,FILE='Phi010PSSCF.dat')       !AUTOCORRELATION POTENTIAL STRESS TENSOR
 OPEN(UNIT=16,FILE='Phi010KPSSCF.dat')      !AUTOCORRELATION KYNETIC-POTENTIAL STRESS TENSOR
 REWIND(10)
 REWIND(11)
 REWIND(12)
 
 NTEL=0
 DTIME=H*DBLE(OUTNSAM)
 TMAX=INT(TOBS/DTIME)
 NSAMP=1500001

 DO I=1,TMAX
    NTIME(I)=0

    SSCFXX(I)=0.0D0
    SSCFXY(I)=0.0D0
    SSCFXZ(I)=0.0D0

    SSCFYX(I)=0.0D0
    SSCFYY(I)=0.0D0
    SSCFYZ(I)=0.0D0

    SSCFZX(I)=0.0D0
    SSCFZY(I)=0.0D0
    SSCFZZ(I)=0.0D0


    KSSCFXX(I)=0.0D0
    KSSCFXY(I)=0.0D0
    KSSCFXZ(I)=0.0D0

    KSSCFYX(I)=0.0D0
    KSSCFYY(I)=0.0D0
    KSSCFYZ(I)=0.0D0

    KSSCFZX(I)=0.0D0
    KSSCFZY(I)=0.0D0
    KSSCFZZ(I)=0.0D0


    PSSCFXX(I)=0.0D0
    PSSCFXY(I)=0.0D0
    PSSCFXZ(I)=0.0D0

    PSSCFYX(I)=0.0D0
    PSSCFYY(I)=0.0D0
    PSSCFYZ(I)=0.0D0

    PSSCFZX(I)=0.0D0
    PSSCFZY(I)=0.0D0
    PSSCFZZ(I)=0.0D0


    KPSSCFXX(I)=0.0D0
    KPSSCFXY(I)=0.0D0
    KPSSCFXZ(I)=0.0D0

    KPSSCFYX(I)=0.0D0
    KPSSCFYY(I)=0.0D0
    KPSSCFYZ(I)=0.0D0

    KPSSCFZX(I)=0.0D0
    KPSSCFZY(I)=0.0D0
    KPSSCFZZ(I)=0.0D0
 ENDDO

 DO K=1,NSAMP
    READ(10,*)IDUMM,SXX,SXY,SXZ,SYX,SYY,SYZ,SZX,SZY,SZZ
    READ(11,*)IDUMM,KSXX,KSXY,KSXZ,KSYX,KSYY,KSYZ,KSZX,KSZY,KSZZ
    READ(12,*)IDUMM,PSXX,PSXY,PSXZ,PSYX,PSYY,PSYZ,PSZX,PSZY,PSZZ
!
    NTEL=NTEL+1
    IF(MOD(NTEL,IT0) .EQ. 0)THEN
      T0=T0+1
      TT0=MOD(T0-1,T0MAX) + 1
      TIME0(TT0)=NTEL

      SXX0(TT0)=SXX
      SXY0(TT0)=SXY
      SXZ0(TT0)=SXZ
      SYX0(TT0)=SYX
      SYY0(TT0)=SYY
      SYZ0(TT0)=SYZ
      SZX0(TT0)=SZX
      SZY0(TT0)=SZY
      SZZ0(TT0)=SZZ

      KSXX0(TT0)=KSXX
      KSXY0(TT0)=KSXY
      KSXZ0(TT0)=KSXZ
      KSYX0(TT0)=KSYX
      KSYY0(TT0)=KSYY
      KSYZ0(TT0)=KSYZ
      KSZX0(TT0)=KSZX
      KSZY0(TT0)=KSZY
      KSZZ0(TT0)=KSZZ

      PSXX0(TT0)=PSXX
      PSXY0(TT0)=PSXY
      PSXZ0(TT0)=PSXZ
      PSYX0(TT0)=PSYX
      PSYY0(TT0)=PSYY
      PSYZ0(TT0)=PSYZ
      PSZX0(TT0)=PSZX
      PSZY0(TT0)=PSZY
      PSZZ0(TT0)=PSZZ
    ENDIF
 
    DO T=1,MIN(T0,T0MAX)
       DELT=NTEL -TIME0(T) + 1
       IF(DELT .LT. TMAX)THEN
          NTIME(DELT)=NTIME(DELT) + 1
          SSCFXX(DELT)=SSCFXX(DELT) + SXX0(T)*SXX
          SSCFXY(DELT)=SSCFXY(DELT) + SXY0(T)*SXY
          SSCFXZ(DELT)=SSCFXZ(DELT) + SXZ0(T)*SXZ

          SSCFYX(DELT)=SSCFYX(DELT) + SYX0(T)*SYX
          SSCFYY(DELT)=SSCFYY(DELT) + SYY0(T)*SYY
          SSCFYZ(DELT)=SSCFYZ(DELT) + SYZ0(T)*SYZ

          SSCFZX(DELT)=SSCFZX(DELT) + SZX0(T)*SZX
          SSCFZY(DELT)=SSCFZY(DELT) + SZY0(T)*SZY
          SSCFZZ(DELT)=SSCFZZ(DELT) + SZZ0(T)*SZZ


          KSSCFXX(DELT)=KSSCFXX(DELT) + KSXX0(T)*KSXX
          KSSCFXY(DELT)=KSSCFXY(DELT) + KSXY0(T)*KSXY
          KSSCFXZ(DELT)=KSSCFXZ(DELT) + KSXZ0(T)*KSXZ

          KSSCFYX(DELT)=KSSCFYX(DELT) + KSYX0(T)*KSYX
          KSSCFYY(DELT)=KSSCFYY(DELT) + KSYY0(T)*KSYY
          KSSCFYZ(DELT)=KSSCFYZ(DELT) + KSYZ0(T)*KSYZ

          KSSCFZX(DELT)=KSSCFZX(DELT) + KSZX0(T)*KSZX
          KSSCFZY(DELT)=KSSCFZY(DELT) + KSZY0(T)*KSZY
          KSSCFZZ(DELT)=KSSCFZZ(DELT) + KSZZ0(T)*KSZZ


          PSSCFXX(DELT)=PSSCFXX(DELT) + PSXX0(T)*PSXX
          PSSCFXY(DELT)=PSSCFXY(DELT) + PSXY0(T)*PSXY
          PSSCFXZ(DELT)=PSSCFXZ(DELT) + PSXZ0(T)*PSXZ

          PSSCFYX(DELT)=PSSCFYX(DELT) + PSYX0(T)*PSYX
          PSSCFYY(DELT)=PSSCFYY(DELT) + PSYY0(T)*PSYY
          PSSCFYZ(DELT)=PSSCFYZ(DELT) + PSYZ0(T)*PSYZ

          PSSCFZX(DELT)=PSSCFZX(DELT) + PSZX0(T)*PSZX
          PSSCFZY(DELT)=PSSCFZY(DELT) + PSZY0(T)*PSZY
          PSSCFZZ(DELT)=PSSCFZZ(DELT) + PSZZ0(T)*PSZZ


          KPSSCFXX(DELT)=KPSSCFXX(DELT) + PSXX0(T)*KSXX + KSXX0(T)*PSXX
          KPSSCFXY(DELT)=KPSSCFXY(DELT) + PSXY0(T)*KSXY + KSXY0(T)*PSXY
          KPSSCFXZ(DELT)=KPSSCFXZ(DELT) + PSXZ0(T)*KSXZ + KSXZ0(T)*PSXZ

          KPSSCFYX(DELT)=KPSSCFYX(DELT) + PSYX0(T)*KSYX + KSYX0(T)*PSYX
          KPSSCFYY(DELT)=KPSSCFYY(DELT) + PSYY0(T)*KSYY + KSYY0(T)*PSYY
          KPSSCFYZ(DELT)=KPSSCFYZ(DELT) + PSYZ0(T)*KSYZ + KSYZ0(T)*PSYZ

          KPSSCFZX(DELT)=KPSSCFZX(DELT) + PSZX0(T)*KSZX + KSZX0(T)*PSZX
          KPSSCFZY(DELT)=KPSSCFZY(DELT) + PSZY0(T)*KSZY + KSZY0(T)*PSZY
          KPSSCFZZ(DELT)=KPSSCFZZ(DELT) + PSZZ0(T)*KSZZ + KSZZ0(T)*PSZZ
       ENDIF
    ENDDO
!
 ENDDO

 DO I=1,TMAX-1
    CTIME=H*DBLE(OUTNSAM)*DBLE(I-1)
    SSCFXX(I)=SSCFXX(I)/DBLE(NTIME(I))
    SSCFXY(I)=SSCFXY(I)/DBLE(NTIME(I))
    SSCFXZ(I)=SSCFXZ(I)/DBLE(NTIME(I))

    SSCFYX(I)=SSCFYX(I)/DBLE(NTIME(I))
    SSCFYY(I)=SSCFYY(I)/DBLE(NTIME(I))
    SSCFYZ(I)=SSCFYZ(I)/DBLE(NTIME(I))

    SSCFZX(I)=SSCFZX(I)/DBLE(NTIME(I))
    SSCFZY(I)=SSCFZY(I)/DBLE(NTIME(I))
    SSCFZZ(I)=SSCFZZ(I)/DBLE(NTIME(I))


    KSSCFXX(I)=KSSCFXX(I)/DBLE(NTIME(I))
    KSSCFXY(I)=KSSCFXY(I)/DBLE(NTIME(I))
    KSSCFXZ(I)=KSSCFXZ(I)/DBLE(NTIME(I))

    KSSCFYX(I)=KSSCFYX(I)/DBLE(NTIME(I))
    KSSCFYY(I)=KSSCFYY(I)/DBLE(NTIME(I))
    KSSCFYZ(I)=KSSCFYZ(I)/DBLE(NTIME(I))

    KSSCFZX(I)=KSSCFZX(I)/DBLE(NTIME(I))
    KSSCFZY(I)=KSSCFZY(I)/DBLE(NTIME(I))
    KSSCFZZ(I)=KSSCFZZ(I)/DBLE(NTIME(I))


    PSSCFXX(I)=PSSCFXX(I)/DBLE(NTIME(I))
    PSSCFXY(I)=PSSCFXY(I)/DBLE(NTIME(I))
    PSSCFXZ(I)=PSSCFXZ(I)/DBLE(NTIME(I))

    PSSCFYX(I)=PSSCFYX(I)/DBLE(NTIME(I))
    PSSCFYY(I)=PSSCFYY(I)/DBLE(NTIME(I))
    PSSCFYZ(I)=PSSCFYZ(I)/DBLE(NTIME(I))

    PSSCFZX(I)=PSSCFZX(I)/DBLE(NTIME(I))
    PSSCFZY(I)=PSSCFZY(I)/DBLE(NTIME(I))
    PSSCFZZ(I)=PSSCFZZ(I)/DBLE(NTIME(I))


    KPSSCFXX(I)=KPSSCFXX(I)/DBLE(NTIME(I))
    KPSSCFXY(I)=KPSSCFXY(I)/DBLE(NTIME(I))
    KPSSCFXZ(I)=KPSSCFXZ(I)/DBLE(NTIME(I))

    KPSSCFYX(I)=KPSSCFYX(I)/DBLE(NTIME(I))
    KPSSCFYY(I)=KPSSCFYY(I)/DBLE(NTIME(I))
    KPSSCFYZ(I)=KPSSCFYZ(I)/DBLE(NTIME(I))

    KPSSCFZX(I)=KPSSCFZX(I)/DBLE(NTIME(I))
    KPSSCFZY(I)=KPSSCFZY(I)/DBLE(NTIME(I))
    KPSSCFZZ(I)=KPSSCFZZ(I)/DBLE(NTIME(I))

    WRITE(13,*)CTIME,SSCFXX(I),SSCFXY(I),SSCFXZ(I),SSCFYX(I),SSCFYY(I),SSCFYZ(I), &
               SSCFZX(I),SSCFZY(I),SSCFZZ(I)

    WRITE(14,*)CTIME,KSSCFXX(I),KSSCFXY(I),KSSCFXZ(I),KSSCFYX(I),KSSCFYY(I), &
               KSSCFYZ(I),KSSCFZX(I),KSSCFZY(I),KSSCFZZ(I)

    WRITE(15,*)CTIME,PSSCFXX(I),PSSCFXY(I),PSSCFXZ(I),PSSCFYX(I),PSSCFYY(I), &
               PSSCFYZ(I),PSSCFZX(I),PSSCFZY(I),PSSCFZZ(I)

    WRITE(16,*)CTIME,KPSSCFXX(I),KPSSCFXY(I),KPSSCFXZ(I),KPSSCFYX(I),KPSSCFYY(I), &
               KPSSCFYZ(I),KPSSCFZX(I),KPSSCFZY(I),KPSSCFZZ(I)
 ENDDO

 DO I=10,16
    CLOSE(UNIT=I)
 ENDDO

 STOP
END PROGRAM SSCF