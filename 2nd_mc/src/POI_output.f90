! GET OUTPUT EXACTLY AT POINT(S) OF INTEREST (POI)
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE POI_output( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE ParticleUtils
  USE GeneralUtils
  USE DefUtils
  USE Interpolation
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!
! default MAX_NAME_LEN=128
!  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'PosOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,PosFile,UsedDataFile
!  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: ThickName
  CHARACTER(LEN=MAX_NAME_LEN) :: ObsFileName,ResFileName,Fname
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar
!  TYPE(Variable_t), POINTER :: CostVar
  TYPE(Variable_t), POINTER :: ThickVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(Nodes_t) :: ElementNodes
  REAL(KIND=dp), POINTER :: ThickValues(:)
  INTEGER, POINTER :: NodeIndexes(:),NodeIndices(:)
  INTEGER, POINTER :: ThickPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  integer :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
!  real(kind=dp) :: Cost,Cost_S
!  REAL(KIND=dp),ALLOCATABLE :: ViscValues(:)
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  INTEGER,SAVE :: ThickDOFs
!  REAL(KIND=dp),ALLOCATABLE,SAVE :: KN_poi(:),KT_poi(:),Visc_poi(:)
!  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:),dirobs(:,:),maskobs(:)
  REAL(KIND=dp),ALLOCATABLE :: xobs(:),THI(:,:),xx(:,:),yy(:,:)
  real(KIND=dp) :: xll,yll,dx,noval
  INTEGER :: ii,jj,nx,ny
  INTEGER :: InElement
  INTEGER :: ElementIndex
  INTEGER,SAVE :: Nexec=0,SAVE_RES_INTERVAL=1
  INTEGER,SAVE :: NTOT=-1,NTOT_S
!  integer,SAVE :: nobs
  LOGICAL,SAVE :: FirstRound=.True.,SAVE_USED_DATA=.False.
  LOGICAL :: DOSAVE
  LOGICAL,SAVE :: SAVE_RESULTS=.False.
  CHARACTER*10 :: date,temps

  INTEGER,PARAMETER :: IO=12

  save Firsttime,Parallel 
  save SolverName
!,VariableName,ViscVariableName,FricVariableName
  save ThickName
  save ResFileName
!  save ElementNodes

   IF(Nexec.eq.0)THEN
 
!! Austfonna
!  nx=3230
!  ny=3726
!  xll=600490.644
!  yll=8773855.4
!  noval=-9999
!  dx=50.0

! Unteraar
  nx=501
  ny=367
  xll=431727.50
  yll=5153399.5
  dx=25.0
  noval=-9999
!  print *,'Jean Krug 0'
  Nexec=Nexec+1

  WRITE(Message,'(a,i1,a)')&
                'Call of VELOCITY OUTPUT SOLVER for the ',&
                nexec, ' . time'

  SolverParams => GetSolverParams()
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif


!------------------------------------------------------------------------------
!  Read Model constants
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  FIRST TIME
!------------------------------------------------------------------------------
! - get mesh details
! - get variable names from sif
! - get variable DOFs
! - get points of interest (POI)
!  print *,'Jean Krug 1'
!  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'POI_output'

!!!!!!!!!!! get Solver Variables

!!  thickness
    ThickName =  GetString( SolverParams,'Thickness Variable Name', Found)
    IF(.NOT.Found) THEN
        CALL FATAL(SolverName,'Keyword >Thickness Variable Name< not found in section >Solver<')
    END IF
    ThickVariable => VariableGet( Solver % Mesh % Variables, Trim(ThickName) )
    IF ( .NOT.ASSOCIATED( ThickVariable ) ) THEN
             WRITE(Message,'(A,A,A)') &
                                'No thickness variable >',Trim(ThickName),' < found'
             CALL FATAL(SolverName,Message)
    END IF
    ThickDOFs=ThickVariable%DOFs

   ResFileName =  GetString( SolverParams,'Result File Name', Found)
   SAVE_RESULTS = Found
   if (Parallel) write(ResFileName,'(A,A,A,I0)') trim(ResFileName),'par',ParEnv % MyPE
   SAVE_RES_INTERVAL = ListGetInteger(SolverParams,'Result save interval',Found)
   If (.NOT.Found) SAVE_RES_INTERVAL=1

   !SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)
   allocate(xobs(3))
   allocate(xx(nx,ny))
   allocate(THI(nx,ny))
   allocate(yy(nx,ny))
!   allocate(xobs(nobs,3),Vobs(nobs,VDOFs),InElement(nobs))
!   allocate(dirobs(nobs,2),maskobs(nobs))

   InElement=-1

!   Vobs=0.0_dp
!   xobs=0.0_dp
!   open(IO,file=trim(ObsFileName))
!   do i=1,nobs
!     read(IO,*) (xobs(i,j),j=1,DIM),(Vobs(i,j),j=1,VDOFs),(dirobs(i,j),j=1,VDOFs),maskobs(i)
!   end do
!   close(IO)

  !!! End of First visit
!    Firsttime=.false.
!  Endif

!  print *,'Jean Krug 2'
!------------------------------------------------------------------------------
!  READ VARIABLES
!------------------------------------------------------------------------------

!!  thickness
    ThickVariable => VariableGet( Solver % Mesh % Variables, Trim(ThickName) )
    IF ( ASSOCIATED( ThickVariable ) ) THEN
            ThickValues => ThickVariable % Values
            ThickPerm => ThickVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable 5 >',Trim(ThickName),' < found'
            CALL FATAL(SolverName,Message)
    END IF


    DOSAVE=SAVE_RESULTS.AND.(MOD(Nexec-1,SAVE_RES_INTERVAL).eq.0)

    IF (DOSAVE) then
        write(FName,'(A,A,A)') trim(ResFileName),'_init','.asc'
        open(IO,file=trim(FName))
    End if


    CALL StartAdvanceOutput(SolverName,'Locate POI and write output')
!  print *,'Jean Krug 3'
! Loop over observations
    do ii=1,nx
    do jj=1,ny
     xx(ii,jj) = (ii-1)*dx+xll
     yy(ii,jj) = (jj-1)*dx+yll

     xobs(1) = xx(ii,jj)
     xobs(2) = yy(ii,jj)
     xobs(3) = 0.0_dp

     THI(ii,jj) = noval
     ElementIndex=-1
!      print *,'Jean Krug 4' 
!     IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(1:DIM)
      CALL LocateParticleInMeshOctree( ElementIndex,Coord)
      If (ElementIndex.NE.0) InElement=ElementIndex

!      print *,InElement,Coord(1),Coord(2),DIM

!     ENDIF !End if FirstRound
!      print *,'Jean Krug 1'      

    ! Data Point has been found in one element
      IF (InElement>0) THEN
         Element => GetActiveElement(InElement)
         IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
         NodeIndices => Element % NodeIndexes
    ! set coords of highest occuring dimension to zero (to get correct path element)
          !-------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
!      print *,'Jean Krug 42'
         IF (DIM == 1) THEN !1D SSA
!      print *,'Jean Krug 421'
            ElementNodes % y(1:n) = 0.0_dp
            ElementNodes % z(1:n) = 0.0_dp
         ELSE IF (DIM == 2) THEN !2D SSA
!      print *,'Jean Krug 422'
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = 0.0_dp
         ELSE
!      print *,'Jean Krug 423'
            WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
            CALL Fatal( SolverName, Message)
         END IF
!  print *,'Jean Krug 5'
!         print *,InElement,xobs(1),xobs(2),ElementNodes % x(1),ElementNodes % y(1),ElementNodes % x(2),ElementNodes % y(2),ElementNodes % x(3),ElementNodes % y(3)
         IF (.NOT.PointInElement(Element,ElementNodes,xobs(1:3),UVW))  THEN
      !        CALL FATAL(SolverName,'Point was supposed to be found in this element')
         ELSE
            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )
!  print *,'Jean Krug 6'
            THI(ii,jj)  = SUM(ThickValues(ThickPerm(NodeIndexes(1:n)))*basis(1:n))

!jjf            if(maskobs(s).eq.1)then
!jjf            write(IO,'(7(e13.5,1x))') x,y,xobs(1),xobs(2),n1,n2,maskobs
!jjf            endif

!            IF (DOSAVE) Then
!              if ((DIM.eq.1)) Then
!                 write(IO,'(3(e15.8,1x))') THI(ii,jj),THI(ii,jj),THI(ii,jj)
!              Else if ((DIM.eq.2)) Then
!                 write(IO,'(3(e15.8,1x))') THI(ii,jj),THI(ii,jj),THI(ii,jj)
!              End if
!            END IF

          END IF

         ELSE

!            WRITE(Message,'(a,I0,I0,a)')&
!                'Data Point',ii,jj,'found in no element'
!            CALL Info( SolverName, Message,level=15)
         END IF

    END DO !do on jj
    END DO !do on ii

    WRITE(IO,'(A,T14,I0)') 'NCOLS',nx
    WRITE(IO,'(A,T14,I0)') 'NROWS',ny
    WRITE(IO,'(A,T14,F14.6)') 'XLLCORNER',xll
    WRITE(IO,'(A,T14,F14.6)') 'YLLCORNER',yll
    WRITE(IO,'(A,T14,F14.6)') 'CELLSIZE',dx
    WRITE(IO,'(A,T14,F7.0)') 'NODATA VLAUE',noval

    do jj=ny,1,-1

    WRITE(IO,'(501(F7.0,2x))') (THI(ii,jj),ii=1,nx) 
!(a(i,j), j=1,10)

    END DO !do on jj


    IF (DOSAVE) close(IO)

!   FirstRound=.False.
    ENDIF ! Nexec=0
   Return
!------------------------------------------------------------------------------
END SUBROUTINE POI_output
!------------------------------------------------------------------------------
! *****************************************************************************
