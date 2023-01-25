!
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE computeSLOPE_drivingstress( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!

!******************************************************************************
!------------------------------------------------------------------------------
!******************************************************************************
  USE ParticleUtils
  USE DefUtils
  USE GeneralUtils
  USE SolverUtils
  USE ParticleUtils
  IMPLICIT NONE

!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  TYPE(ValueList_t), POINTER :: SolverParams
!
  CHARACTER(LEN=MAX_NAME_LEN) :: ObsFileName
  CHARACTER(LEN=MAX_NAME_LEN) :: ThiVarSolName,FluxDirVarSolName,SlopeVarSolName
  TYPE(Variable_t), POINTER :: ThiVariable,FluxDirectionVariable,SlopeVariable,Hsia
  REAL(KIND=dp), POINTER :: ThiValues(:),FluxDirectionValues(:),SlopeValues(:)
  REAL(KIND=dp),ALLOCATABLE :: Global_THIValues(:)
  INTEGER, POINTER :: ThiPerm(:),FluxDirectionPerm(:),SlopePerm(:)
!
  INTEGER :: ok
  INTEGER,SAVE :: nobs
!
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex
  real(kind=dp) :: xx,yy
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
!
  real(kind=dp) :: SLOPE,THI,SLOPEX,SLOPEY
  real(kind=dp) :: Arate,rhoice,grav,secy,alpha
  real(kind=dp) :: a1,a2,a3,b1,b2
  real(kind=dp) :: slope_threshold
  real(kind=dp) :: Average_Thi, AVE_THI, Sum_Thi, Global_Sum_Thi
  integer :: i,j,k,l,n,m,DIM,countit_S,countit,s
  integer :: nflow,ierr,Global_M
!
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
!
  INTEGER,PARAMETER :: IO=12
!
  save Firsttime,Parallel
  save SolverName,ThiVarSolName,FluxDirVarSolName,SlopeVarSolName


  SolverParams => GetSolverParams()

!-------------------------------------------------------------------------------
!
!                              MATERIAL CONSTANTS
!
!-------------------------------------------------------------------------------
     ! slope threshold ... if below SIA thickness is taken
     slope_threshold = GetConstReal( Model % Constants, 'slope threshold', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant slope threshold not found. &
                   &Setting to default (0.0)'
            CALL INFO(SolverName, Message, level=20)
            slope_threshold = 0.0
     End if

!-------------------------------------------------------------------------------
!
!                              CHECK IF INPUT VARIABLES ARE DEFINED
!
!-------------------------------------------------------------------------------
!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

  ThiVarSolName =  GetString( SolverParams,'Thickness Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Thickness Variable Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  FluxDirVarSolName =  GetString( SolverParams,'Flux Direction Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Flux Direction Variable Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  SlopeVarSolName =  GetString( SolverParams,'Slope Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Slope Variable Name< not found in section >Solver<')
      END IF

!-------------------------------------------------------------------------------
!
!                              FIRST TIME CHECKING
!
!-------------------------------------------------------------------------------

  If (Firsttime) then
    N = Solver%Mesh%NumberOfNodes

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'computeSLOPE_fluxdirection'

  !!! End of First visit
    Firsttime=.false.
  Endif

!-------------------------------------------------------------------------------
!
!                              MATERIAL CONSTANTS
!
!-------------------------------------------------------------------------------
     ! seconds per year
     secy = GetConstReal( Model % Constants, 'seconds per year', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Seconds Per Year not found. &
                   &Setting to (60*60*24*365.25)'
            CALL INFO(SolverName, Message, level=20)
            secy = 60.0*60.0*24.0*365.25
     End if

     ! ice density
     rhoice = GetConstReal( Model % Constants, 'ice density', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Ice Density not found. &
                   &Setting to 917 kg/m3'
            CALL INFO(SolverName, Message, level=20)
            rhoice = 917.0/(1.0e6*(secy)**2)
     End if

     ! gravitational acceleration
     grav = GetConstReal( Model % Constants, 'gravity constant', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Gravity Constant not found. &
                   &Setting to 9.81 m/s^2'
            CALL INFO(SolverName, Message, level=20)
            grav = -9.81*(secy)**2
     End if

!!! SOME INITIALISATION AT FIRST TIME
  If (Firsttime) then
    N = model % MaxElementNodes

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'computeSLOP_fluxdirection'

  !!! End of First visit
    Firsttime=.false.
  Endif

!-------------------------------------------------------------------------------
!
!                              GET VARIABLES
!
!-------------------------------------------------------------------------------

    ThiVariable => VariableGet( Solver % Mesh % Variables, ThiVarSolName  )
    IF ( ASSOCIATED( ThiVariable ) ) THEN
            ThiValues => ThiVariable % Values
            ThiPerm => ThiVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                    'No variable >',ThiVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    FluxDirectionVariable => VariableGet( Solver % Mesh % Variables, FluxDirVarSolName  )
    IF ( ASSOCIATED( FluxDirectionVariable ) ) THEN
            FluxDirectionValues => FluxDirectionVariable % Values
            FluxDirectionPerm => FluxDirectionVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',FluxDirVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    SlopeVariable => VariableGet( Solver % Mesh % Variables, SlopeVarSolName  )
    IF ( ASSOCIATED( SlopeVariable ) ) THEN
            SlopeValues => SlopeVariable % Values
            SlopePerm => SlopeVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',SlopeVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

!-------------------------------------------------------------------------------
!
!                              START LOOP
!
!-------------------------------------------------------------------------------

    N = Solver%Mesh%NumberOfNodes
    M = N
    K = Solver % Matrix % ParMatrix % ParEnv % MyPE

    IF (PARALLEL) THEN
        Sum_Thi = SUM(ThiValues)
      !IF (Solver % Matrix % ParMatrix % ParEnv % MyPE.EQ.0) THEN
      !CALL MPI_ALLREDUCE(ThiValues,Global_ThiValues,M,&
      !       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(M,Global_M,1,&
             MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      !print *,'THIAVE',N,M,Global_M

!      IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
        allocate(Global_ThiValues(Global_M))
        !CALL MPI_ALLGATHERV(ThiValues,M,MPI_DOUBLE_PRECISION,&
        !       Global_ThiValues,M,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)
        CALL MPI_ALLREDUCE(Sum_Thi,Global_Sum_Thi,1,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        Average_Thi = Global_Sum_Thi/Global_M
        !CALL MPI_GATHER(ThiValues,M,MPI_DOUBLE_PRECISION,&
        !       Global_ThiValues,M,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !print *,'THIAVE',Solver % Matrix % ParMatrix % ParEnv % MyPE , SUM(Global_ThiValues)/Global_M
        !Average_Thi = SUM(Global_ThiValues)/Global_M
      !  print *,'THIAVE',Solver % Matrix % ParMatrix % ParEnv % MyPE , Average_Thi
        !CALL MPI_SCATTER(Average_Thi,1,MPI_DOUBLE_PRECISION,AVE_THI,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        deallocate(Global_ThiValues)
!      ENDIF
       !print *,'THIAVE',Solver % Matrix % ParMatrix % ParEnv % MyPE ,Average_Thi,AVE_THI
    ELSE
       Average_Thi = SUM(ThiValues)/N
    ENDIF


    Do i=1,N


        THI  = ThiValues(ThiPerm(i))

        alpha = 3.0

        IF (THI.LT.1.0/alpha*Average_Thi) THEN
        THI = 1.0/alpha*Average_Thi
        ELSEIF (THI.GT.alpha*Average_Thi) THEN
        THI = alpha*Average_Thi
        ENDIF

        !!THI = Average_Thi+0.9*(THI-Average_Thi)
        !THI = Average_Thi

        SLOPEX = FluxDirectionValues(2*(FluxDirectionPerm(i)-1)+1)/abs(grav)/rhoice/thi*(-1.0)
        SLOPEY = FluxDirectionValues(2*(FluxDirectionPerm(i)-1)+2)/abs(grav)/rhoice/thi*(-1.0)
!
        SlopeValues(2*(SlopePerm(i)-1)+1) = SLOPEX
        SlopeValues(2*(SlopePerm(i)-1)+2) = SLOPEY

        !SLOPEX = FluxDirectionValues(2*(FluxDirectionPerm(i)-1)+1)/abs(grav)/rhoice/thi/(-1.0)
        !SLOPEY = FluxDirectionValues(2*(FluxDirectionPerm(i)-1)+2)/abs(grav)/rhoice/thi/(-1.0)
!
        !SlopeValues(2*(SlopePerm(i)-1)+1) = SLOPEX
        !SlopeValues(2*(SlopePerm(i)-1)+2) = SLOPEY

    End do

    !deallocate(xobs,THIobs,InElement)
    !deallocate(ElementNodes % x,ElementNodes % y,ElementNodes % z)

    Firsttime=.false.

END SUBROUTINE computeSLOPE_drivingstress

