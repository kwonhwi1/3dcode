MODULE RHS_MODULE
  USE CONFIG_MODULE
  USE VARIABLE_MODULE
  USE GRID_MODULE
  USE EOS_MODULE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBORDINATE TO RHS MODULE  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  USE FLUX_MODULE
  USE VSFLUX_MODULE
  USE MUSCL_MODULE
  USE CAV_MODULE
  USE TURBSOURCE_MODULE
  USE UNSTEADY_MODULE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: T_RHS
  
  TYPE T_RHS
    PRIVATE
    INTEGER :: NPV,NDV,NTV,NGRD,IMAX,JMAX,KMAX
    LOGICAL :: L_FLUX,L_MUSCL,L_CAV,L_TURBSOURCE,L_VSFLUX,L_UNSTEADY
    REAL(8) :: PREF
    REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: RES,ICAV,ITT
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::OMEGA_CUT
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: EA,FA,GA,EVA,FVA,GVA
    CLASS(T_FLUX),       ALLOCATABLE :: FLUX
    CLASS(T_MUSCL),      ALLOCATABLE :: MUSCL
    CLASS(T_CAV),        ALLOCATABLE :: CAV
    CLASS(T_TURBSOURCE), ALLOCATABLE :: TURBSOURCE
    CLASS(T_VSFLUX),     ALLOCATABLE :: VSFLUX
    CLASS(T_UNSTEADY),   ALLOCATABLE :: UNSTEADY
    CONTAINS
      PROCEDURE :: CONSTRUCT
      PROCEDURE :: DESTRUCT
      PROCEDURE :: CALRHS
      PROCEDURE :: GETRES
      PROCEDURE :: GETICAV
      PROCEDURE :: GETITT
      PROCEDURE :: GETOMEGA_CUT
  END TYPE T_RHS
  
  CONTAINS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE CONSTRUCT(RHS,CONFIG,GRID,VARIABLE)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(OUT) :: RHS
      TYPE(T_CONFIG), INTENT(IN) :: CONFIG
      TYPE(T_GRID), INTENT(IN) :: GRID
      TYPE(T_VARIABLE), INTENT(IN) :: VARIABLE
      
      RHS%PREF = CONFIG%GETPREF()
      
      SELECT CASE(CONFIG%GETITURB())
      CASE(-3)

      CASE(-2)
        ALLOCATE(T_VSFLUX_LAMINAR::RHS%VSFLUX)

      CASE(-1)
        ALLOCATE(T_VSFLUX_TURBULENT::RHS%VSFLUX)
        ALLOCATE(T_KEPSILON::RHS%TURBSOURCE)

      CASE(0)
        ALLOCATE(T_VSFLUX_TURBULENT::RHS%VSFLUX)
        ALLOCATE(T_KWSST::RHS%TURBSOURCE)

      END SELECT
      
      SELECTCASE(CONFIG%GETNSCHEME())
      CASE(1)
        ALLOCATE(T_ROE::RHS%FLUX)
      CASE(2)
        ALLOCATE(T_ROEM::RHS%FLUX)
      CASE(3)
        ALLOCATE(T_AUSMPWP::RHS%FLUX)
      CASE(4)
        ALLOCATE(T_AUSMPUP::RHS%FLUX)
      END SELECT
            
      SELECT CASE(CONFIG%GETNMUSCL())
      CASE(0)
      CASE(1)
        ALLOCATE(T_TVD::RHS%MUSCL)
      CASE(2)
        ALLOCATE(T_MLP::RHS%MUSCL)
      END SELECT
                
      SELECT CASE(CONFIG%GETNCAV())
      CASE(0)
      CASE(1)
        ALLOCATE(T_MERKLE::RHS%CAV)
      CASE(2)
        ALLOCATE(T_KUNZ::RHS%CAV)
      CASE(3)
        ALLOCATE(T_SINGHAL::RHS%CAV)
      END SELECT

      SELECT CASE(CONFIG%GETNSTEADY())
      CASE(0)
      CASE(1)
        ALLOCATE(T_UNSTEADY::RHS%UNSTEADY)
      END SELECT

      RHS%L_FLUX       = .FALSE.
      RHS%L_MUSCL      = .FALSE.
      RHS%L_CAV        = .FALSE.
      RHS%L_TURBSOURCE = .FALSE.
      RHS%L_VSFLUX     = .FALSE.
      RHS%L_UNSTEADY   = .FALSE.
      
      IF(ALLOCATED(RHS%FLUX)) THEN
        CALL RHS%FLUX%CONSTRUCT(CONFIG,VARIABLE)
        RHS%L_FLUX = .TRUE.
      END IF
      
      IF(ALLOCATED(RHS%VSFLUX)) THEN
        CALL RHS%VSFLUX%CONSTRUCT(CONFIG,GRID,VARIABLE)
        RHS%L_VSFLUX = .TRUE.
      END IF
      
      IF(ALLOCATED(RHS%MUSCL)) THEN
        CALL RHS%MUSCL%CONSTRUCT(CONFIG,VARIABLE)
        RHS%L_MUSCL = .TRUE.
      END IF
   
      IF(ALLOCATED(RHS%CAV)) THEN
        CALL RHS%CAV%CONSTRUCT(CONFIG,GRID,VARIABLE)
        ALLOCATE(RHS%ICAV(0:4,2:GRID%GETIMAX(),2:GRID%GETJMAX(),2:GRID%GETKMAX()))
        RHS%L_CAV = .TRUE.
      END IF
      
      IF(ALLOCATED(RHS%TURBSOURCE)) THEN
        CALL RHS%TURBSOURCE%CONSTRUCT(CONFIG,GRID,VARIABLE)
        ALLOCATE(RHS%ITT(6,2:GRID%GETIMAX(),2:GRID%GETJMAX(),2:GRID%GETKMAX()))
        ALLOCATE(RHS%OMEGA_CUT(2:GRID%GETIMAX(),2:GRID%GETJMAX(),2:GRID%GETKMAX()))
        RHS%L_TURBSOURCE = .TRUE.
      END IF
      
      IF(ALLOCATED(RHS%UNSTEADY)) THEN
        CALL RHS%UNSTEADY%CONSTRUCT(CONFIG,GRID,VARIABLE)
        RHS%L_UNSTEADY   = .TRUE.
      END IF

      RHS%NGRD = GRID%GETNGRD()
      RHS%IMAX = GRID%GETIMAX()
      RHS%JMAX = GRID%GETJMAX()
      RHS%KMAX = GRID%GETKMAX()
      RHS%NPV = VARIABLE%GETNPV()
      RHS%NDV = VARIABLE%GETNDV()
      RHS%NTV = VARIABLE%GETNTV()

      ALLOCATE(RHS%RES(RHS%NPV,2:RHS%IMAX,2:RHS%JMAX,2:RHS%KMAX))
      ALLOCATE(RHS%EA(RHS%NPV,RHS%IMAX),RHS%FA(RHS%NPV,RHS%JMAX),RHS%GA(RHS%NPV,RHS%KMAX))
      ALLOCATE(RHS%EVA(RHS%NPV,RHS%IMAX),RHS%FVA(RHS%NPV,RHS%JMAX),RHS%GVA(RHS%NPV,RHS%KMAX))

    END SUBROUTINE CONSTRUCT
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE DESTRUCT(RHS)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(INOUT) :: RHS
      
      IF(RHS%L_FLUX) THEN
        CALL RHS%FLUX%DESTRUCT()  
        DEALLOCATE(RHS%FLUX)
      END IF
      IF(RHS%L_VSFLUX) THEN
        CALL RHS%VSFLUX%DESTRUCT()
        DEALLOCATE(RHS%VSFLUX)
      END IF
      IF(RHS%L_MUSCL) THEN
        CALL RHS%MUSCL%DESTRUCT()
        DEALLOCATE(RHS%MUSCL)
      END IF
      
      IF(RHS%L_CAV) THEN
        CALL RHS%CAV%DESTRUCT()
        DEALLOCATE(RHS%CAV,RHS%ICAV)
      END IF
      IF(RHS%L_TURBSOURCE) THEN
        CALL RHS%TURBSOURCE%DESTRUCT()
        DEALLOCATE(RHS%TURBSOURCE,RHS%ITT,RHS%OMEGA_CUT)
      END IF
      IF(RHS%L_UNSTEADY) THEN
        CALL RHS%UNSTEADY%DESTRUCT()
        DEALLOCATE(RHS%UNSTEADY)
      END IF
      
      DEALLOCATE(RHS%RES)
      DEALLOCATE(RHS%EA,RHS%FA,RHS%GA,RHS%EVA,RHS%FVA,RHS%GVA)
      
    END SUBROUTINE DESTRUCT
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE CALRHS(RHS,GRID,VARIABLE,EOS)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(INOUT) :: RHS
      TYPE(T_GRID), INTENT(IN) :: GRID
      TYPE(T_VARIABLE), INTENT(IN) :: VARIABLE
      TYPE(T_EOS), INTENT(IN) :: EOS
      INTEGER :: I,J,K
      INTEGER :: II,JJ,KK,LL
      REAL(8) :: NX(3)
      REAL(8) :: EX1(3),EX2(3),EX3(3),EX4(3)
      REAL(8) :: TX1(3),TX2(3),TX3(3),TX4(3)
      REAL(8) :: GRDL(RHS%NGRD),GRDR(RHS%NGRD)
      REAL(8) :: X(38,RHS%NPV)
      REAL(8) :: PVL(RHS%NPV),PVR(RHS%NPV),PV(RHS%NPV)
      REAL(8) :: DVL(RHS%NDV),DVR(RHS%NDV)
      REAL(8) :: TVL(RHS%NTV),TVR(RHS%NTV)
   
      RHS%EVA = 0.D0
      RHS%FVA = 0.D0
      RHS%GVA = 0.D0
      
      DO K=2,RHS%KMAX
        DO J=2,RHS%JMAX
          DO I=1,RHS%IMAX
            LL = 0
            DO KK = -1,1
              DO JJ = -1,1
                DO II = 0,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+JJ,K+KK)
                END DO
              END DO
            END DO
            IF(RHS%L_MUSCL) THEN
              DO KK=-1,1
                DO JJ= -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I-1,J+JJ,K+KK)
                END DO
              END DO
              DO KK=-1,1
                DO JJ= -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+2,J+JJ,K+KK)
                END DO
              END DO
              X(LL+1,:) = VARIABLE%GETPV(I-2,J,K)
              X(LL+2,:) = VARIABLE%GETPV(I+3,J,K)
              
              CALL RHS%MUSCL%SETPV(X)
              CALL RHS%MUSCL%INTERPOLATION(PVL,PVR)
              IF((PVL(6).LT.0.D0).OR.(PVL(6).GT.1.D0)) PVL(6) = X(9,6)
              IF((PVR(6).LT.0.D0).OR.(PVR(6).GT.1.D0)) PVR(6) = X(10,6)
              
              IF(RHS%L_TURBSOURCE) THEN
                IF((PVL(8).LT.0.D0).AND.(X(9,8).GT.0.D0))  PVL(8) = X(9,8)
                IF((PVR(8).LT.0.D0).AND.(X(10,8).GT.0.D0)) PVR(8) = X(10,8)
                IF((PVL(9).LT.0.D0).AND.(X(9,9).GT.0.D0))  PVL(9) = X(9,9)
                IF((PVR(9).LT.0.D0).AND.(X(10,9).GT.0.D0)) PVR(9) = X(10,9)
              END IF
              
              CALL EOS%DETEOS(PVL(1)+RHS%PREF,PVL(5),PVL(6),PVL(7),DVL)
              CALL EOS%DETEOS(PVR(1)+RHS%PREF,PVR(5),PVR(6),PVR(7),DVR)
            ELSE
              PVL = VARIABLE%GETPV(I,J,K) 
              PVR = VARIABLE%GETPV(I+1,J,K)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I+1,J,K)
            END IF
            
            NX = GRID%GETCX(I,J,K)
            
            CALL RHS%FLUX%SETNORM(NX)
            CALL RHS%FLUX%SETPV(PVL,PVR)
            CALL RHS%FLUX%SETDV(DVL,DVR)
            CALL RHS%FLUX%SETSDST(X(1:18,1))
            CALL RHS%FLUX%CALFLUX(EOS,RHS%EA(:,I))
            
            IF(RHS%L_VSFLUX) THEN
              EX1 = GRID%GETEX(I,J-1,K)
              EX2 = GRID%GETEX(I+1,J-1,K)
              EX3 = GRID%GETEX(I,J,K)
              EX4 = GRID%GETEX(I+1,J,K)
              TX1 = GRID%GETTX(I,J,K-1)
              TX2 = GRID%GETTX(I+1,J,K-1)
              TX3 = GRID%GETTX(I,J,K)
              TX4 = GRID%GETTX(I+1,J,K)              
              GRDL = GRID%GETGRD(I,J,K)
              GRDR = GRID%GETGRD(I+1,J,K)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I+1,J,K)
              TVL = VARIABLE%GETTV(I,J,K)
              TVR = VARIABLE%GETTV(I+1,J,K)
              CALL RHS%VSFLUX%SETNORM(NX,EX1,EX2,EX3,EX4,TX1,TX2,TX3,TX4)
              CALL RHS%VSFLUX%SETGRD(GRDL,GRDR)
              CALL RHS%VSFLUX%SETPV(X)
              CALL RHS%VSFLUX%SETDV(DVL,DVR)
              CALL RHS%VSFLUX%SETTV(TVL,TVR)
              CALL RHS%VSFLUX%CALFLUX(RHS%EVA(:,I))
            END IF
          END DO
          DO I=2,RHS%IMAX
            RHS%RES(:,I,J,K) = -( RHS%EA(:,I) - RHS%EA(:,I-1) ) + (RHS%EVA(:,I) - RHS%EVA(:,I-1) )
          END DO
        END DO
      END DO
      
      DO I=2,RHS%IMAX
        DO K=2,RHS%KMAX
          DO J=1,RHS%JMAX
            LL = 0
            DO II = -1,1
              DO KK = -1,1
                DO JJ = 0,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+JJ,K+KK)
                END DO
              END DO
            END DO
            
            IF(RHS%L_MUSCL) THEN
              DO II = -1,1
                DO KK = -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J-1,K+KK)
                END DO
              END DO
              DO II = -1,1
                DO KK = -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+2,K+KK)
                END DO
              END DO
              X(LL+1,:) = VARIABLE%GETPV(I,J-2,K)
              X(LL+2,:) = VARIABLE%GETPV(I,J+3,K)
              
              CALL RHS%MUSCL%SETPV(X)
              CALL RHS%MUSCL%INTERPOLATION(PVL,PVR)
              IF((PVL(6).LT.0.D0).OR.(PVL(6).GT.1.D0)) PVL(6) = X(9,6)
              IF((PVR(6).LT.0.D0).OR.(PVR(6).GT.1.D0)) PVR(6) = X(10,6)
              
              IF(RHS%L_TURBSOURCE) THEN
                IF((PVL(8).LT.0.D0).AND.(X(9,8).GT.0.D0))  PVL(8) = X(9,8)
                IF((PVR(8).LT.0.D0).AND.(X(10,8).GT.0.D0)) PVR(8) = X(10,8)
                IF((PVL(9).LT.0.D0).AND.(X(9,9).GT.0.D0))  PVL(9) = X(9,9)
                IF((PVR(9).LT.0.D0).AND.(X(10,9).GT.0.D0)) PVR(9) = X(10,9)
              END IF
              
              CALL EOS%DETEOS(PVL(1)+RHS%PREF,PVL(5),PVL(6),PVL(7),DVL)
              CALL EOS%DETEOS(PVR(1)+RHS%PREF,PVR(5),PVR(6),PVR(7),DVR)
            ELSE
              PVL = VARIABLE%GETPV(I,J,K) 
              PVR = VARIABLE%GETPV(I,J+1,K)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I,J+1,K)
            END IF

            NX = GRID%GETEX(I,J,K)
            
            CALL RHS%FLUX%SETNORM(NX)
            CALL RHS%FLUX%SETPV(PVL,PVR)
            CALL RHS%FLUX%SETDV(DVL,DVR)
            CALL RHS%FLUX%SETSDST(X(1:18,1))
            CALL RHS%FLUX%CALFLUX(EOS,RHS%FA(:,J))
            
            IF(RHS%L_VSFLUX) THEN
              EX1 = GRID%GETTX(I,J,K-1)
              EX2 = GRID%GETTX(I,J+1,K-1)
              EX3 = GRID%GETTX(I,J,K)
              EX4 = GRID%GETTX(I,J+1,K)
              TX1 = GRID%GETCX(I-1,J,K)
              TX2 = GRID%GETCX(I,J,K)
              TX3 = GRID%GETCX(I-1,J+1,K)
              TX4 = GRID%GETCX(I,J+1,K)
              GRDL = GRID%GETGRD(I,J,K)
              GRDR = GRID%GETGRD(I,J+1,K)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I,J+1,K)
              TVL = VARIABLE%GETTV(I,J,K)
              TVR = VARIABLE%GETTV(I,J+1,K)
              CALL RHS%VSFLUX%SETNORM(NX,EX1,EX2,EX3,EX4,TX1,TX2,TX3,TX4)
              CALL RHS%VSFLUX%SETGRD(GRDL,GRDR)
              CALL RHS%VSFLUX%SETPV(X)
              CALL RHS%VSFLUX%SETDV(DVL,DVR)
              CALL RHS%VSFLUX%SETTV(TVL,TVR)
              CALL RHS%VSFLUX%CALFLUX(RHS%FVA(:,J))
            END IF
          END DO
          DO J=2,RHS%JMAX
            RHS%RES(:,I,J,K) = RHS%RES(:,I,J,K) -( RHS%FA(:,J) - RHS%FA(:,J-1) ) + (RHS%FVA(:,J) - RHS%FVA(:,J-1) )
          END DO
        END DO
      END DO

      DO J=2,RHS%JMAX
        DO I=2,RHS%IMAX
          DO K=1,RHS%KMAX
            LL = 0
            DO JJ = -1,1
              DO II = -1,1
                DO KK = 0,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+JJ,K+KK)
                END DO
              END DO
            END DO
            
            IF(RHS%L_MUSCL) THEN
              DO JJ = -1,1
                DO II = -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+JJ,K-1)
                END DO
              END DO
              DO JJ = -1,1
                DO II = -1,1
                  LL = LL+1
                  X(LL,:) = VARIABLE%GETPV(I+II,J+JJ,K+2)
                END DO
              END DO
              X(LL+1,:) = VARIABLE%GETPV(I,J,K-2)
              X(LL+2,:) = VARIABLE%GETPV(I,J,K+3)
              
              CALL RHS%MUSCL%SETPV(X)
              CALL RHS%MUSCL%INTERPOLATION(PVL,PVR)
              IF((PVL(6).LT.0.D0).OR.(PVL(6).GT.1.D0)) PVL(6) = X(9,6)
              IF((PVR(6).LT.0.D0).OR.(PVR(6).GT.1.D0)) PVR(6) = X(10,6)
              
              IF(RHS%L_TURBSOURCE) THEN
                IF((PVL(8).LT.0.D0).AND.(X(9,8).GT.0.D0))  PVL(8) = X(9,8)
                IF((PVR(8).LT.0.D0).AND.(X(10,8).GT.0.D0)) PVR(8) = X(10,8)
                IF((PVL(9).LT.0.D0).AND.(X(9,9).GT.0.D0))  PVL(9) = X(9,9)
                IF((PVR(9).LT.0.D0).AND.(X(10,9).GT.0.D0)) PVR(9) = X(10,9)
              END IF
              
              CALL EOS%DETEOS(PVL(1)+RHS%PREF,PVL(5),PVL(6),PVL(7),DVL)
              CALL EOS%DETEOS(PVR(1)+RHS%PREF,PVR(5),PVR(6),PVR(7),DVR)
            ELSE
              PVL = VARIABLE%GETPV(I,J,K) 
              PVR = VARIABLE%GETPV(I,J,K+1)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I,J,K+1)
            END IF

            NX = GRID%GETTX(I,J,K)
            
            CALL RHS%FLUX%SETNORM(NX)
            CALL RHS%FLUX%SETPV(PVL,PVR)
            CALL RHS%FLUX%SETDV(DVL,DVR)
            CALL RHS%FLUX%SETSDST(X(1:18,1))
            CALL RHS%FLUX%CALFLUX(EOS,RHS%GA(:,K))
            
            IF(RHS%L_VSFLUX) THEN
              EX1 = GRID%GETCX(I-1,J,K)
              EX2 = GRID%GETCX(I-1,J,K+1)
              EX3 = GRID%GETCX(I,J,K)
              EX4 = GRID%GETCX(I,J,K+1)
              TX1 = GRID%GETEX(I,J-1,K)
              TX2 = GRID%GETEX(I,J,K)
              TX3 = GRID%GETEX(I,J-1,K+1)
              TX4 = GRID%GETEX(I,J,K+1)
              GRDL = GRID%GETGRD(I,J,K)
              GRDR = GRID%GETGRD(I,J,K+1)
              DVL = VARIABLE%GETDV(I,J,K)
              DVR = VARIABLE%GETDV(I,J,K+1)
              TVL = VARIABLE%GETTV(I,J,K)
              TVR = VARIABLE%GETTV(I,J,K+1)
              CALL RHS%VSFLUX%SETNORM(NX,EX1,EX2,EX3,EX4,TX1,TX2,TX3,TX4)
              CALL RHS%VSFLUX%SETGRD(GRDL,GRDR)
              CALL RHS%VSFLUX%SETPV(X)
              CALL RHS%VSFLUX%SETDV(DVL,DVR)
              CALL RHS%VSFLUX%SETTV(TVL,TVR)
              CALL RHS%VSFLUX%CALFLUX(RHS%GVA(:,K))
            END IF
          END DO
          DO K=2,RHS%KMAX
            RHS%RES(:,I,J,K) = RHS%RES(:,I,J,K) -( RHS%GA(:,K) - RHS%GA(:,K-1) ) + (RHS%GVA(:,K) - RHS%GVA(:,K-1) )
          END DO
        END DO
      END DO
      
      DO K=2,RHS%KMAX
        DO J=2,RHS%JMAX
          DO I=2,RHS%IMAX
            EX1 = GRID%GETCX(I-1,J,K)
            EX2 = GRID%GETCX(I,J,K)
            EX3 = GRID%GETEX(I,J-1,K)
            EX4 = GRID%GETEX(I,J,K)
            TX1 = GRID%GETTX(I,J,K-1)
            TX2 = GRID%GETTX(I,J,K)
            GRDL = GRID%GETGRD(I,J,K)
            X(1,:) = VARIABLE%GETPV(I-1,J,K)
            X(2,:) = VARIABLE%GETPV(I,J,K)
            X(3,:) = VARIABLE%GETPV(I+1,J,K)
            X(4,:) = VARIABLE%GETPV(I,J-1,K)
            X(5,:) = VARIABLE%GETPV(I,J+1,K)
            X(6,:) = VARIABLE%GETPV(I,J,K-1)
            X(7,:) = VARIABLE%GETPV(I,J,K+1)
            PV = X(2,:)
            DVL = VARIABLE%GETDV(I,J,K)
            TVL = VARIABLE%GETTV(I,J,K)
        
            IF(RHS%L_CAV) THEN
              CALL RHS%CAV%SETGRD(GRDL)
              CALL RHS%CAV%SETPV(PV)
              CALL RHS%CAV%SETDV(DVL)
              RHS%ICAV(:,I,J,K) = RHS%CAV%CAVSOURCE(EOS)
              RHS%RES(6,I,J,K) = RHS%RES(6,I,J,K) + RHS%ICAV(0,I,J,K)
            END IF
            
            IF(RHS%L_TURBSOURCE) THEN
              CALL RHS%TURBSOURCE%SETNORM(EX1,EX2,EX3,EX4,TX1,TX2)
              CALL RHS%TURBSOURCE%SETGRD(GRDL)
              CALL RHS%TURBSOURCE%SETPV(X)
              CALL RHS%TURBSOURCE%SETDV(DVL)
              CALL RHS%TURBSOURCE%SETTV(TVL)
              RHS%ITT(:,I,J,K) = RHS%TURBSOURCE%CALTURBSOURCE()
              RHS%OMEGA_CUT(I,J,K) = RHS%TURBSOURCE%GETOMEGA_CUT()
              RHS%RES(8:9,I,J,K) = RHS%RES(8:9,I,J,K) + RHS%ITT(1:2,I,J,K)
            END IF
            
            IF(RHS%L_UNSTEADY) THEN
              PVL = VARIABLE%GETQQ(1,I,J,K)
              PVR = VARIABLE%GETQQ(2,I,J,K)
              CALL RHS%UNSTEADY%SETGRD(GRDL)
              CALL RHS%UNSTEADY%SETPV(PV)
              CALL RHS%UNSTEADY%SETDV(DVL)
              CALL RHS%UNSTEADY%SETQQ(PVL,PVR)
              RHS%RES(:,I,J,K) = RHS%RES(:,I,J,K) - RHS%UNSTEADY%UNSTEADYSOURCE()
            END IF
            
          END DO
        END DO
      END DO
    END SUBROUTINE CALRHS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION GETRES(RHS,N,I,J,K)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(IN) :: RHS
      INTEGER, INTENT(IN) :: N,I,J,K
      REAL(8) :: GETRES
      
      GETRES = RHS%RES(N,I,J,K)
      
    END FUNCTION GETRES
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION GETICAV(RHS,I,J,K)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(IN) :: RHS
      INTEGER, INTENT(IN) :: I,J,K
      REAL(8) :: GETICAV(4)
      
      GETICAV = RHS%ICAV(1:4,I,J,K)
      
    END FUNCTION GETICAV
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION GETITT(RHS,I,J,K)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(IN) :: RHS
      INTEGER, INTENT(IN) :: I,J,K
      REAL(8) :: GETITT(4)
      
      GETITT= RHS%ITT(3:6,I,J,K)
      
    END FUNCTION GETITT
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION GETOMEGA_CUT(RHS,I,J,K)
      IMPLICIT NONE
      CLASS(T_RHS), INTENT(IN) :: RHS
      INTEGER, INTENT(IN) :: I,J,K
      REAL(8) :: GETOMEGA_CUT
      
      GETOMEGA_CUT= RHS%OMEGA_CUT(I,J,K)
      
    END FUNCTION GETOMEGA_CUT
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
END MODULE RHS_MODULE