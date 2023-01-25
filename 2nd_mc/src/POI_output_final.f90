! GET OUTPUT EXACTLY AT POINT(S) OF INTEREST (POI)
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE POI_output_final( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: ThickName,FluxName,SlopeName,SurfaceName,VarName
  CHARACTER(LEN=MAX_NAME_LEN) :: ObsFileName,ResFileName,Fname
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar
!  TYPE(Variable_t), POINTER :: CostVar
  TYPE(Variable_t), POINTER :: ThickVariable,FluxVariable,SlopeVariable,SurfaceVariable,Variable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(Nodes_t) :: ElementNodes
  REAL(KIND=dp), POINTER :: ThickValues(:),FluxValues(:),SlopeValues(:),SurfaceValues(:),VarValues(:)
  INTEGER, POINTER :: NodeIndexes(:),NodeIndices(:)
  INTEGER, POINTER :: ThickPerm(:),FluxPerm(:),SlopePerm(:),SurfacePerm(:),VarPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  integer :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
!  real(kind=dp) :: Cost,Cost_S
!  REAL(KIND=dp),ALLOCATABLE :: ViscValues(:)
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  INTEGER,SAVE :: ThickDOFs,SlopeDOFs,VarDOFs
!  REAL(KIND=dp),ALLOCATABLE,SAVE :: KN_poi(:),KT_poi(:),Visc_poi(:)
!  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:),dirobs(:,:),maskobs(:)
  REAL(KIND=dp),ALLOCATABLE :: xobs(:,:),THI(:),THI_S(:),xx(:),yy(:),Vobs(:),VAR(:),VAR_S(:)
  REAL(KIND=dp):: SLOPE,SURF,FLUX,SLOPEX,SLOPEY
  real(KIND=dp) :: xll,yll,dx,noval
  INTEGER :: ii,jj,nobs
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
  save ThickName,FluxName,SlopeName,SurfaceName,VarName
  save ResFileName
!  save ElementNodes

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

!!  output variable
    VarName =  GetString( SolverParams,'Output Variable Name', Found)
    IF(.NOT.Found) THEN
        CALL FATAL(SolverName,'Keyword >Output Variable Name< not found in section >Solver<')
    END IF
    Variable => VariableGet( Solver % Mesh % Variables, Trim(VarName) )
    IF ( .NOT.ASSOCIATED( Variable ) ) THEN
             WRITE(Message,'(A,A,A)') &
                                'No variable >',Trim(VarName),' < found'
             CALL FATAL(SolverName,Message)
    END IF
    VarDOFs=Variable%DOFs

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

!!  ice flux
!    FluxName =  GetString( SolverParams,'Flux Variable Name', Found)
!    IF(.NOT.Found) THEN
!        CALL FATAL(SolverName,'Keyword >Flux Variable Name< not found in section >Solver<')
!    END IF
!    FluxVariable => VariableGet( Solver % Mesh % Variables, Trim(FluxName) )
!    IF ( .NOT.ASSOCIATED( FluxVariable ) ) THEN
!             WRITE(Message,'(A,A,A)') &
!                                'No thickness variable >',Trim(FluxName),' < found'
!             CALL FATAL(SolverName,Message)
!    END IF

!!  surface slope
    SlopeName =  GetString( SolverParams,'Slope Variable Name', Found)
    IF(.NOT.Found) THEN
        CALL FATAL(SolverName,'Keyword >Slope Variable Name< not found in section >Solver<')
    END IF
    SlopeVariable => VariableGet( Solver % Mesh % Variables, Trim(SlopeName) )
    IF ( .NOT.ASSOCIATED( SlopeVariable ) ) THEN
             WRITE(Message,'(A,A,A)') &
                                'No thickness variable >',Trim(SlopeName),' <found'
             CALL FATAL(SolverName,Message)
    END IF
    SlopeDOFs=SlopeVariable%DOFs

!!  surface elevation
    SurfaceName =  GetString( SolverParams,'Surface Variable Name', Found)
    IF(.NOT.Found) THEN
        CALL FATAL(SolverName,'Keyword >Surface Variable Name< not found in section >Solver<')
    END IF
    SurfaceVariable => VariableGet( Solver % Mesh % Variables, Trim(SurfaceName) )
    IF ( .NOT.ASSOCIATED( SurfaceVariable ) ) THEN
             WRITE(Message,'(A,A,A)') &
                                'No thickness variable >',Trim(SurfaceName),' <found'
             CALL FATAL(SolverName,Message)
    END IF

   ResFileName =  GetString( SolverParams,'Result File Name', Found)
   SAVE_RESULTS = Found
   !if (Parallel) write(ResFileName,'(A,A,A,I0)') trim(ResFileName),'_','par',ParEnv % MyPE
   SAVE_RES_INTERVAL = ListGetInteger(SolverParams,'Result save interval',Found)
   If (.NOT.Found) SAVE_RES_INTERVAL=1


 !! Get the obs
   ObsFileName =  GetString( SolverParams,'POI File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >POI File Name< not found in section >Solver<')
   END IF

   open(IO,file=trim(ObsFileName),status = 'old',iostat = ok)
   if(ok /= 0) then
       write(message,'(A,A)') 'Unable to open file ',TRIM(ObsFileName)
       CALL Fatal(Trim(SolverName),Trim(message))
   end if
   nobs=0
   do while(ok == 0)
     read(io,*,iostat = ok)
     if (ok == 0) nobs = nobs + 1
   end do
   close(IO)


   allocate(xobs(nobs,3),Vobs(nobs))!,InElement(nobs))
   !InElement(:)=-1
   InElement=-1

   Vobs = 0.0_dp
!   V    = 0.0_dp
   xobs = 0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),Vobs(i)
   end do
   close(IO)

   !SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)
   allocate(xx(nobs))
   allocate(THI(nobs),THI_S(nobs),VAR_S(nobs),VAR(nobs))
   allocate(yy(nobs))

  !!! End of First visit
!    Firsttime=.false.
!  Endif

!  print *,'Jean Krug 2'
!------------------------------------------------------------------------------
!  READ VARIABLES
!------------------------------------------------------------------------------

!!  output variable
    Variable => VariableGet( Solver % Mesh % Variables, Trim(VarName) )
    IF ( ASSOCIATED( Variable ) ) THEN
            VarValues => Variable % Values
            VarPerm => Variable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable 5 >',Trim(VarName),' < found'
            CALL FATAL(SolverName,Message)
    END IF

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

!!  ice flux
!    FluxVariable => VariableGet( Solver % Mesh % Variables, Trim(FluxName) )
!    IF ( ASSOCIATED( FluxVariable ) ) THEN
!            FluxValues => FluxVariable % Values
!            FluxPerm => FluxVariable % Perm
!    ELSE
!            WRITE(Message,'(A,A,A)') &
!                               'No variable 5 >',Trim(FluxName),' < found'
!            CALL FATAL(SolverName,Message)
!    END IF

!!  surface slope
    SlopeVariable => VariableGet( Solver % Mesh % Variables, Trim(SlopeName) )
    IF ( ASSOCIATED( SlopeVariable ) ) THEN
            SlopeValues => SlopeVariable % Values
            SlopePerm => SlopeVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable 5 >',Trim(SlopeName),' < found'
            CALL FATAL(SolverName,Message)
    END IF

!!  surface
    SurfaceVariable => VariableGet( Solver % Mesh % Variables, Trim(SurfaceName) )
    IF ( ASSOCIATED( SurfaceVariable ) ) THEN
            SurfaceValues => SurfaceVariable % Values
            SurfacePerm => SurfaceVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable 5 >',Trim(SurfaceName),' < found'
            CALL FATAL(SolverName,Message)
    END IF


    DOSAVE=SAVE_RESULTS.AND.(MOD(Nexec-1,SAVE_RES_INTERVAL).eq.0)

    IF (DOSAVE) then
        write(FName,'(A,A,A,A)') trim(ResFileName),'_','final','.asc'
        IF(.NOT.Parallel)THEN
        open(IO,file=trim(FName))
        ENDIF
    End if


    CALL StartAdvanceOutput(SolverName,'Locate POI and write output')
!  print *,'Jean Krug 3'
! Loop over observations
    !do ii=1,nx

    VAR_S = 0.0_dp
    THI_S = 0.0_dp

    do jj=1,nobs
     xx(jj) = xobs(jj,1)
     yy(jj) = xobs(jj,2)

     xobs(jj,3) = 0.0_dp

     noval = 0.0_dp
     !THI(jj) = noval
     ElementIndex=-1
!      print *,'Jean Krug 4' 
!     IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(jj,1:DIM)
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
         IF (.NOT.PointInElement(Element,ElementNodes,xobs(jj,1:3),UVW))  THEN
!          WRITE(IO,'(501(E16.8,3x))') xx(jj),yy(jj),-9999.
      !        CALL FATAL(SolverName,'Point was supposed to be found in this element')
          IF(.NOT.Parallel)THEN
          WRITE(IO,'(3(F16.8,1x))') xx(jj),yy(jj),-9999.
          ENDIF
         ELSE
            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )
!  print *,'Jean Krug 6'
            VAR_S(jj) = SUM(VarValues(VarPerm(NodeIndexes(1:n)))*basis(1:n))
            THI_S(jj)  = SUM(ThickValues(ThickPerm(NodeIndexes(1:n)))*basis(1:n))
            !FLUX = (SUM(FluxValues(FluxPerm(NodeIndexes(1:n)))*basis(1:n))**2.0)**0.5
            SLOPEX = SUM(SlopeValues(SlopeDOFs*(SlopePerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
            SLOPEY = SUM(SlopeValues(SlopeDOFs*(SlopePerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
            SLOPE = (SLOPEX*SLOPEX+SLOPEY*SLOPEY)**0.5
            SURF = SUM(SurfaceValues(SurfacePerm(NodeIndexes(1:n)))*basis(1:n))
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

          !IF(THI(jj).GT.0.0)THEN
!          IF(FLUX.lt.100.and.SLOPE.lt.1e-2.and.SURF.gt.100)then
!          ELSE
          IF(.NOT.Parallel)THEN
          WRITE(IO,'(3(F16.8,1x))') xx(jj),yy(jj),VAR_S(jj)
          ENDIF
!          ENDIF
          !ENDIF

          END IF

         ELSE

          IF(.NOT.Parallel)THEN
          WRITE(IO,'(3(F16.8,1x))') xx(jj),yy(jj),-9999.
          ENDIF

!          WRITE(IO,'(501(E16.8,3x))') xx(jj),yy(jj),-9999.

!            WRITE(Message,'(a,I0,I0,a)')&
!                'Data Point',ii,jj,'found in no element'
!            CALL Info( SolverName, Message,level=15)
         END IF

    END DO !do on jj
!    END DO !do on ii

!    WRITE(IO,'(A,T14,I0)') 'NCOLS',nx
!    WRITE(IO,'(A,T14,I0)') 'NROWS',ny
!    WRITE(IO,'(A,T14,F14.6)') 'XLLCORNER',xll
!    WRITE(IO,'(A,T14,F14.6)') 'YLLCORNER',yll
!    WRITE(IO,'(A,T14,F14.6)') 'CELLSIZE',dx
!    WRITE(IO,'(A,T14,F7.0)') 'NODATA VLAUE',noval

!    do jj=1,nobs,1
!
!    WRITE(IO,'(501(F16.8,3x))') xx(jj),yy(jj),THI(jj) 
!
!    END DO !do on jj

    IF (Parallel) THEN

       CALL MPI_ALLREDUCE(VAR_S,VAR,nobs,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(THI_S,THI,nobs,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
    IF(DOSAVE) open(IO,file=trim(FName))
    do jj=1,nobs
       if(THI(jj).gt.0.0)then
       WRITE(IO,'(3(E16.8,1x))') xx(jj),yy(jj),VAR(jj)
       else
       WRITE(IO,'(3(E16.8,1x))') xx(jj),yy(jj),-9999.
       endif
    enddo
    ENDIF

    ENDIF

    IF (DOSAVE) close(IO)

!   FirstRound=.False.

   Return
!------------------------------------------------------------------------------
END SUBROUTINE POI_output_final
!------------------------------------------------------------------------------
! *****************************************************************************
