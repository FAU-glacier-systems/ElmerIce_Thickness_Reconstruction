!
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE computeSIA_Hobs( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams,BC
!
  CHARACTER(LEN=MAX_NAME_LEN) :: ObsFileName
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,FluxVarSolName,SlopeVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: GradSlopeXSolName, GradSlopeYSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: SurfaceVarSolName!, ThiVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: Dist2MarVarSolName, Dist2ObsVarSolName
  CHARACTER(LEN=10) :: GlacierType
  TYPE(Variable_t), POINTER :: Variable,FluxVariable,SlopeVariable,Hsia
  TYPE(Variable_t), POINTER :: GradSlopeXVariable, GradSlopeYVariable
  TYPE(Variable_t), POINTER :: SurfaceVariable!,ThiVariable
  TYPE(Variable_t), POINTER :: Dist2MarVariable, Dist2ObsVariable
  REAL(KIND=dp), POINTER :: Values(:),FluxValues(:),SlopeValues(:)
  REAL(KIND=dp), POINTER :: GradSlopeXValues(:), GradSlopeYValues(:)
  REAL(KIND=dp), POINTER :: SurfaceValues(:)!,ThiValues(:)
  REAL(KIND=dp), POINTER :: Dist2MarValues(:), Dist2ObsValues(:)
  INTEGER, POINTER :: Perm(:),FluxPerm(:),SlopePerm(:)
  INTEGER, POINTER :: GradSlopeXPerm(:),GradSlopeYPerm(:)
  INTEGER, POINTER :: SurfacePerm(:)!,ThiPerm(:)
  INTEGER, POINTER :: Dist2MarPerm(:), Dist2ObsPerm(:)
  INTEGER :: ok
  INTEGER :: boundaryIDrange
  INTEGER,SAVE :: nobs,nbe,nbe_s
  INTEGER, ALLOCATABLE :: nbes(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),THIobs(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: AA(:), AA_S(:), BB_S(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: AA0(:), AA0_S(:), BB0_S(:)
  REAL(KIND=dp),ALLOCATABLE :: xx(:),yy(:),zz(:),xx_s(:),yy_s(:),zz_s(:),zz0(:),zz0_s(:)
  REAL(KIND=dp),ALLOCATABLE :: FLUX(:),SLOPE(:),FLUX_S(:),SLOPE_S(:)
  REAL(KIND=dp),ALLOCATABLE :: visc_scaling_margin(:),visc_scaling_margin_s(:)
  REAL(KIND=dp) :: SURF,DIST2MAR,DIST2OBS,VISC
  INTEGER,ALLOCATABLE :: COUNTER_S(:),COUNTER(:)
!
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t) :: ElementNodes
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex,bc_id
  !real(kind=dp) :: xx,yy
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
!
  real(kind=dp) :: DSDX,DSDY,THI,MFLUX,kappa,MFLUX_S,DSDXX,DSDYY
  real(kind=dp) :: FLUX_MIN_S,FLUX_MAX_S,FLUX_MIN,FLUX_MAX,ALPHA,BETA,OMEGA,MP
  real(kind=dp) :: HOBS,HMAR
  real(kind=dp) :: AR,AR_S,BR,BR_S
  real(kind=dp) :: AR0,AR0_S,BR0,BR0_S
  real(kind=dp) :: Arate,rhoice,grav,secy
  real(kind=dp) :: slope_threshold,slope_div_threshold,slope_div_trans,slope_factor
  real(kind=dp) :: a1,a2,a3,b1,b2
  real(kind=dp) :: visc_slope_thresh,visc_slope_grad
  real(kind=dp) :: visc_ref_el_z,visc_el_grad, el_min, el_max, el_z
  real(kind=dp) :: coupling_length,visc_threshold
!,rf_scaling,dist_scaling
  real(kind=dp) :: visc_slope_scaling, visc_el_scaling, visc_scaling
  real(kind=dp) :: visc_dist2mar_scaling, visc_dist2obs_scaling
  integer :: i,j,k,l,n,m,DIM,s,countit,countit_S,t,ii
  integer :: nflow,ierr,nPar,NNN,NNN_S
!
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  Logical :: ViscAverage,CalvingFlag

  INTEGER,PARAMETER :: IO=12,IO2=13,IO3=14

  save Firsttime,Parallel,nPar
  save SolverName,VarSolName,FluxVarSolName,SlopeVarSolName
  save GradSlopeXSolName,GradSlopeYSolName
!  save ElementNodes



  SolverParams => GetSolverParams()

!-------------------------------------------------------------------------------
!
!                              OPTION FLAGS
!
!-------------------------------------------------------------------------------



  CalvingFlag =  GetLogical( SolverParams,'Calving Flag', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Calving Flag< not found in section >Solver<')
              CALL WARN(SolverName,'Assuming no calving. CalvingFlag=.FALSE.')
              CalvingFlag = .FALSE.
      END IF

  ViscAverage =  GetLogical( SolverParams,'Viscosity Averaging', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Viscosity Averaging< not found in section >Solver<')
              CALL WARN(SolverName,'Assuming standard averaging of ice visocities.')
              ViscAverage = .TRUE.
      END IF


!-------------------------------------------------------------------------------
!
!                              CHECK IF INPUT VARIABLES ARE DEFINED
!
!-------------------------------------------------------------------------------
!! Dimension of the pb; ie with SSA it can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

!! Name if thickness variable
!  ThiVarSolName =  GetString( SolverParams,'Thickness Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Thickness Variable Name< not found in section >Solver<')
!      END IF

!! Name of the variable to regularise
  FluxVarSolName =  GetString( SolverParams,'Flux Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Flux Variable Name< not found in section >Solver<')
      END IF

!! Name if surface elevation variable
  SurfaceVarSolName =  GetString( SolverParams,'Surface Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Surface Variable Name< not found in section >Solver<')
      END IF

!! Name if surface elevation variable
  SurfaceVarSolName =  GetString( SolverParams,'Surface Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Surface Variable Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  SlopeVarSolName =  GetString( SolverParams,'Slope Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Slope Variable Name< not found in section >Solver<')
      END IF


!! Name of the variable to regularise
  GradSlopeXSolName =  GetString( SolverParams,'X Slope Gradient Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >X Slope Gradient Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  GradSlopeYSolName =  GetString( SolverParams,'Y Slope Gradient Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Y Slope Gradient Name< not found in section >Solver<')
      END IF

!! Name of variabe for distance to next location of thickness observation 
  Dist2ObsVarSolName =  GetString( SolverParams,'Obs Distance Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Obs Distance Variable Name< not found in section >Solver<')
      END IF

!! Name of variabe for distance to glacier margin/outline
  Dist2MarVarSolName =  GetString( SolverParams,'Margin Distance Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Margin Distance Variable Name< not found in section >Solver<')
      END IF


!-------------------------------------------------------------------------------
!
!                              FIRST TIME CHECKING
!
!-------------------------------------------------------------------------------

  If (Firsttime) then

    !N = Solver%Mesh%NumberOfNodes

    ! Check whether job is parallel
    Parallel = .FALSE.
    nPar     = 0
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
      IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
        Parallel = .TRUE.
        nPar = Solver %  Matrix % ParMatrix % ParEnv % PEs
        nbe_s = Solver % Mesh % NumberOfBoundaryElements
      END IF
    END IF

    WRITE(SolverName, '(A)') 'computeSIAthickness'

  Endif !!! End of First visit


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
     rhoice = rhoice*(1.0e6*(secy)**2)

     ! rate factor
     Arate = GetConstReal( Model % Constants, 'rate factor', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Rate Factor not found. &
                   &Setting to 100.0*1e-25 Pa^-3 yr-1'
            CALL INFO(SolverName, Message, level=20)
            Arate = 1.0e-16*(1.0e6**3)
     End if
     Arate = Arate/(1.0e6**3)

     ! flow exponent
     nflow = GetConstReal( Model % Constants, 'flow exponent', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Flow Exponent not found. &
                   &Setting to 3'
            CALL INFO(SolverName, Message, level=20)
            nflow = 3
     End if

     ! gravitational acceleration
     grav = GetConstReal( Model % Constants, 'gravity constant', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Gravity Constant not found. &
                   &Setting to 9.81 m/s^2'
            CALL INFO(SolverName, Message, level=20)
            grav = -9.81*(secy)**2
     End if
     grav = grav/(secy)**2

     ! SLOPE DEPENDENT VISCOSITY
     ! Define slope threshold (beyond which no scaling is applied)
     visc_slope_thresh = GetConstReal( Model % Constants, 'viscosity slope threshold', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<viscosity slope threshold> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            visc_slope_thresh = 0.0
     End if

     ! slope threshold ... if below SIA thickness is taken
     slope_threshold = GetConstReal( Model % Constants, 'slope threshold', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant slope threshold not found. &
                   &Setting to default (0.0)'
            CALL INFO(SolverName, Message, level=20)
            slope_threshold = 0.0
     End if

     ! slope gradient for viscosity computation
     visc_slope_grad = GetConstReal( Model % Constants, 'viscosity slope gradient', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<viscosity slope threshold> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            visc_slope_grad = 0.0
     End if


     ! MINIMUM ELEVATION
     el_min = GetConstReal( Model % Constants, 'minimum elevation', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<minimum elevation> not found. &
                   &Setting to default (0.0).'
            CALL INFO(SolverName, Message, level=20)
            el_min = 0.0
     End if

     ! MAXIMUM ELEVATION
     el_max = GetConstReal( Model % Constants, 'maximum elevation', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<maximum elevation> not found. &
                   &Setting to default (0.0).'
            CALL INFO(SolverName, Message, level=20)
            el_max = 0.0
     End if


     ! ELEVATION DEPENDENT VISCOSITY
     ! Define reference elevation where scaling is 1.0
     visc_ref_el_z = GetConstReal( Model % Constants, 'viscosity reference elevation', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<viscosity reference elevation> not found. &
                   &Setting to default (0.0).'
            CALL INFO(SolverName, Message, level=20)
            visc_ref_el_z = 0.0
     End if

     ! slope gradient for viscosity computation
     visc_el_grad = GetConstReal( Model % Constants, 'viscosity elevation gradient', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<viscosity elevation gradient> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            visc_el_grad = 0.0
     End if

     ! coupling length scale for glaciers and ice-sheets (effect of longitudinal stress coupling)
     coupling_length = GetConstReal( Model % Constants, 'coupling length', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant >Coupling Length< not found. &
                   &Setting to 10.0'
            CALL INFO(SolverName, Message, level=20)
            coupling_length=10.0_dp
     End if



!-------------------------------------------------------------------------------
!
!                              GET VARIABLES
!
!-------------------------------------------------------------------------------

!    Variable => VariableGet( Solver % Mesh % Variables, VarSolName  )
!    IF ( ASSOCIATED( Variable ) ) THEN
!            Values => Variable % Values
!            Perm => Variable % Perm
!    ELSE
!            WRITE(Message,'(A,A,A)') &
!                               'No variable >',VarSolName,' < found'
!            CALL FATAL(SolverName,Message)
!    END IF


    FluxVariable => VariableGet( Solver % Mesh % Variables, FluxVarSolName  )
    IF ( ASSOCIATED( FluxVariable ) ) THEN
            FluxValues => FluxVariable % Values
            FluxPerm => FluxVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',FluxVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

!    ThiVariable => VariableGet( Solver % Mesh % Variables, ThiVarSolName  )
!    IF ( ASSOCIATED( ThiVariable ) ) THEN
!            ThiValues => ThiVariable % Values
!            ThiPerm => ThiVariable % Perm
!    ELSE
!            WRITE(Message,'(A,A,A)') &
!                               'No variable >',ThiVarSolName,' < found'
!            CALL FATAL(SolverName,Message)
!    END IF

    SurfaceVariable => VariableGet( Solver % Mesh % Variables, SurfaceVarSolName  )
    IF ( ASSOCIATED( SurfaceVariable ) ) THEN
            SurfaceValues => SurfaceVariable % Values
            SurfacePerm => SurfaceVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',SurfaceVarSolName,' < found'
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

    GradSlopeXVariable => VariableGet( Solver % Mesh % Variables, GradSlopeXSolName  )
    IF ( ASSOCIATED( GradSlopeXVariable ) ) THEN
            GradSlopeXValues => GradSlopeXVariable % Values
            GradSlopeXPerm => GradSlopeXVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',GradSlopeXSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    GradSlopeYVariable => VariableGet( Solver % Mesh % Variables, GradSlopeYSolName  )
    IF ( ASSOCIATED( GradSlopeYVariable ) ) THEN
            GradSlopeYValues => GradSlopeYVariable % Values
            GradSlopeYPerm => GradSlopeYVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',GradSlopeYSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF


    Dist2MarVariable => VariableGet( Solver % Mesh % Variables, Dist2MarVarSolName)
    IF ( ASSOCIATED( Dist2MarVariable ) ) THEN
            Dist2MarValues => Dist2MarVariable % Values
            Dist2MarPerm => Dist2MarVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',Dist2MarVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF


    Dist2ObsVariable => VariableGet( Solver % Mesh % Variables, Dist2ObsVarSolName)
    IF ( ASSOCIATED( Dist2ObsVariable ) ) THEN
            Dist2ObsValues => Dist2ObsVariable % Values
            Dist2ObsPerm => Dist2ObsVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',Dist2ObsVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF


    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(nbes(nPar))


!! Get the obs
   ObsFileName =  GetString( SolverParams,'Observation File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Observation File Name< not found in section >Solver<')
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


   allocate(xobs(nobs,3),THIobs(nobs),InElement(nobs))
   allocate(SLOPE(nobs),SLOPE_S(nobs),FLUX(nobs),FLUX_S(nobs))
   allocate(AA(nobs),AA_S(nobs),BB_S(nobs))
   allocate(AA0(nobs),AA0_S(nobs),BB0_S(nobs))
   allocate(COUNTER_S(nobs),COUNTER(nobs))
   InElement(:)=-1

   THIobs=0.0_dp
   xobs=0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),(THIobs(i))
   end do
   close(IO)

!-------------------------------------------------------------------------------
!
!                              DETERMINE MEAN FLUX VALUE
!
!           NECESSARY FOR CORRECTING SMALL AND NEGATIVE FLUX VALUES BELOW
!
!-------------------------------------------------------------------------------

    visc_threshold = 1e+4**(-1.0/nflow) !1e-10**(-1.0/nflow)

    ! Initialisation
    MFLUX_S = 0.0_dp
    FLUX_MIN_S = 0.0_dp
    FLUX_MAX_S = 0.0_dp
    NNN_S   = 0

    Do i=1,Solver%Mesh%NumberOfNodes
      MFLUX_S = MFLUX_S + abs(FluxValues(FluxPerm(i)))
      NNN_S   = NNN_S+1
    Enddo

    FLUX_MIN_S = minval(FluxValues)
    FLUX_MAX_S = maxval(FluxValues)

    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(NNN_S,NNN,1,&
             MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(MFLUX_S,MFLUX,1,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(FLUX_MIN_S,FLUX_MIN,1,&
             MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(FLUX_MAX_S,FLUX_MAX,1,&
             MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      MFLUX = MFLUX/NNN
    ELSE
      MFLUX = MFLUX_S/NNN_S
      FLUX_MIN = FLUX_MIN_S
      FLUX_MAX = FLUX_MAX_S
    ENDIF

!-------------------------------------------------------------------------------
!
!                              START OBSERVATION LOOP
!
!-------------------------------------------------------------------------------

    IF(.NOT.Parallel)THEN
    if(FirstTime)then
    open(IO,file='./ice_viscosity.dat',status = 'unknown', position = 'rewind', iostat = ok)
    else
    open(IO,file='./ice_viscosity.dat',status = 'replace', iostat = ok)
    endif
    ENDIF

    CALL StartAdvanceOutput(SolverName,'Loop through thickness observations to determine best rate factor')

    ! Initialisation
    AR_S      = 0.0_dp
    AA_S      = 0.0_dp
    BB_S      = 0.0_dp
    AR0_S      = 0.0_dp
    AA0_S      = 0.0_dp
    BB0_S      = 0.0_dp

    countit_S = 0
    FLUX_S    = 0.0_dp
    SLOPE_S   = 0.0_dp
    COUNTER_S = 0
    
!    IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
    Do s=1,nobs
      AA_S(s) = 0.0_dp
      BB_S(s) = 0.0_dp
      CALL AdvanceOutput(s,nobs)

      !IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(s,1:DIM)
      CALL LocateParticleInMeshOctree( ElementIndex,Coord)

      If (ElementIndex.NE.0) InElement(s)=ElementIndex

        ! Data Point has been found in one element
        IF (InElement(s)>0) THEN
          Element => GetActiveElement(InElement(s))
          n = GetElementNOFNodes()
          NodeIndexes => Element % NodeIndexes

          ! set coords of highest occurring dimension to zero (to get correct path element)
          !-------------------------------------------------------------------------------
          ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
          IF (DIM == 1) THEN !1D SSA
            ElementNodes % y(1:n) = 0.0_dp
            ElementNodes % z(1:n) = 0.0_dp
          ELSE IF (DIM == 2) THEN !2D SSA
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = 0.0_dp
          ELSE
            WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
            CALL Fatal( SolverName, Message)
          END IF

          IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:3),UVW))  THEN
              CALL FATAL(SolverName,'Point was supposed to be found in this element')
          ELSE
            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )
            HOBS = THIobs(s)
            DSDX =  SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
            DSDY =  SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
            DSDXX = SUM(GradSlopeXValues(2*(GradSlopeXPerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
            DSDYY = SUM(GradSlopeYValues(2*(GradSlopeYPerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
            FLUX_S(s) = SUM(FluxValues(FluxPerm(NodeIndexes(1:n)))*basis(1:n))

            DIST2MAR = SUM(Dist2MarValues(Dist2MarPerm(NodeIndexes(1:n)))*basis(1:n))
            DIST2OBS = SUM(Dist2ObsValues(Dist2ObsPerm(NodeIndexes(1:n)))*basis(1:n))
            SURF = SUM(SurfaceValues(SurfacePerm(NodeIndexes(1:n)))*basis(1:n))


            kappa = 1.0-2.0/3.141592653589793*atan(FLUX_S(s)*FLUX_S(s)/(FLUX_MIN*FLUX_MIN))
            if(FLUX_S(s).le.0)then
              kappa = 1.0
            endif

            SLOPE_S(s) = (DSDX**2.0+DSDY**2.0)**0.5
            IF (SLOPE_S(s).LT.TAN(slope_threshold/180*3.14159)) THEN
              SLOPE_S(s)=TAN(slope_threshold/180*3.14159)
            ENDIF


            !! Correction strategy for divides and ridges
            !! proved not longer suitable on Vernagt and Hofsjokull
            !slope_div_threshold = 5e-4
            !slope_div_trans     = 2.0
            !slope_factor        = 10.0
            !IF (SLOPE_S(s).LT.TAN(slope_factor*slope_threshold/180*3.14159))THEN
            !  IF (DSDXX+DSDYY>slope_div_threshold) THEN
            !    OMEGA = slope_factor*slope_threshold
            !    SLOPE_S(s) = TAN(OMEGA/180*3.14159)
            !  ELSE IF ((DSDXX+DSDYY)>slope_div_threshold/slope_div_trans) THEN
            !    OMEGA = slope_factor*slope_threshold
            !    MP = ((DSDXX+DSDYY)-slope_div_threshold/slope_div_trans)/((1.0-1.0/slope_div_trans)*slope_div_threshold)
            !    SLOPE_S(s)=MP*TAN(OMEGA/180*3.14159)+(1.0-MP)*SLOPE_S(s)
            !  ENDIF
            !ENDIF


           ! SLOPE CORRECTION FACTOR 
           ! applicable to ice viscosities
           ! calibrated in the Alps by Christian Sommer
           if(180.0/3.14159*ATAN(SLOPE_S(s)).LT.visc_slope_thresh+1.0/visc_slope_grad)then
             !rf_scaling = (visc_slope_grad*(180.0/3.14159*ATAN(SLOPE_S(s))-visc_slope_thresh))
             visc_slope_scaling = (visc_slope_grad*(180.0/3.14159*ATAN(SLOPE_S(s))-visc_slope_thresh))
           else
             !rf_scaling = 1.0
             visc_slope_scaling = 1.0
           endif

           ! ELEVATION DEPENDENT SCALING
           if (SURF.gt.el_min.and.SURF.lt.el_max) then
             !rf_scaling = rf_scaling*(1.0+visc_el_grad*(SURF-visc_ref_elev))
             !visc_el_scaling = (1.0+visc_el_grad*(SURF-visc_ref_elev))
             el_z            = (SURF-el_min)/(el_max-el_min)
             visc_el_scaling = 1.0+visc_el_grad*(el_z-visc_ref_el_z)
           else
             visc_el_scaling = 1.0
           end if

           IF (SURF.LT.el_min) THEN
             el_z            = 0
             visc_el_scaling = 1.0+visc_el_grad*(el_z-visc_ref_el_z)
           END IF

           IF (SURF.GT.el_max) THEN
             el_z            = 1.0
             visc_el_scaling = 1.0+visc_el_grad*(el_z-visc_ref_el_z)
           END IF


           ! WEIGHTING ACCORDING TO DISTANCE TO NEXT THICKNESS OBSERVATIONS
           ! if there are no observations; distances are set high --> no effect
           IF (HOBS.GT.0)THEN
             !dist_scaling  = atan(DIST2OBS/(10.0*HOBS))*2.0/3.141592653589793
             visc_dist2obs_scaling  = atan(1.0*DIST2OBS/(10.0*HOBS))*2.0/3.141592653589793
           ELSE
             !dist_scaling  = 1.0
             visc_dist2obs_scaling = 1.0
           ENDIF

           ! WEIGHTING ACCORDING TO DISTANCE TO GLACIER MARGIN/OUTLINE
           ! account for marine termini
           IF (DIST2MAR.GE.0.AND.HOBS.GT.0) THEN
             !visc_dist2mar_scaling = (atan(1.0*DIST2MAR/(10.0*coupling_length*HOBS))*2.0/3.141592653589793)+1.0e-3
             visc_dist2mar_scaling = (atan(1.0*DIST2MAR/(5.0*HOBS))*2.0/3.141592653589793)+1.0e-3
           ELSE
             visc_dist2mar_scaling = 1.0e-3
           ENDIF

           !visc_dist2mar_scaling = 1.0
           visc_dist2obs_scaling = 1.0

           visc_scaling = (visc_dist2mar_scaling*(visc_el_scaling*visc_slope_scaling)**visc_dist2obs_scaling)


            a1 = abs(((1.0-kappa)*abs(FLUX_S(s))+kappa*abs(FLUX_MIN))/(SLOPE_S(s)**(nflow)))
            a2 = (1.0*nflow+2.0)/(2.0)*1.0/((rhoice*abs(grav))**nflow)
            a3 = 1.0/(HOBS)**(nflow+2.0)

            IF(a1*a2*a3.gt.0.and.visc_scaling.gt.0.and.visc_dist2mar_scaling.gt.0)THEN
              COUNTER_S(s) = COUNTER_S(s)+1
              !BB0_S(s) = (a1*a2*a3)**(-1.0/nflow)/(visc_scaling)
              BB0_S(s) = (a1*a2*a3)**(-1.0/nflow)!/visc_dist2mar_scaling
              AA0_S(s) = BB0_S(s)**(-nflow)
              BB_S(s)  = (a1*a2*a3)**(-1.0/nflow)/(visc_scaling)
              AA_S(s)  = BB_S(s)**(-nflow) 
              !AA_S(s)  = a1*a2*a3/(visc_dist2mar_scaling*rf_scaling**dist_scaling)**(-nflow)
              !print *, 'exit0',(visc_dist2mar_scaling*rf_scaling**dist_scaling), visc_scaling
              !print *, 'exit2', BB_S(s)**(-nflow), AA_S(s)
              AR_S = AR_S+AA_S(s)
              BR_S = BR_S+AA_S(s)**(-1.0/nflow)
              AR0_S = AR0_S+AA0_S(s)
              BR0_S = BR0_S+AA0_S(s)**(-1.0/nflow)
              countit_S = countit_S + 1

              IF(.NOT.Parallel)THEN
                if(AA_S(s).gt.0.and.AA_S(s).lt.visc_threshold**(-nflow))then
                WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xobs(s,1),xobs(s,2),AA_S(s)**(-1.0/nflow)
                endif
              ENDIF

            ENDIF

          ENDIF

        ELSE

        WRITE(Message,'(a,I0,a)')&
             'Data Point',s,'found in no element'
        CALL Info( SolverName, Message,level=15)
      END IF
    enddo ! s
!
!
! Collect information from parallel threads
    IF (Parallel) THEN
        CALL MPI_ALLREDUCE(AR_S,AR,1,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(AR0_S,AR0,1,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(countit_S,countit,1,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(AA_S,AA,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(AA0_S,AA0,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(SLOPE_S,SLOPE,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(COUNTER_S,COUNTER,nobs,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(FLUX_S,FLUX,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(nbe_s,nbe,1,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLGATHER(nbe_s,1,MPI_INTEGER,nbes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        allocate(xx_s(nbe),yy_s(nbe),zz_s(nbe),xx(nbe),yy(nbe),zz(nbe))
        allocate(zz0(nbe),zz0_s(nbe))
        allocate(visc_scaling_margin(nbe),visc_scaling_margin_s(nbe))

        xx_s = 0.0_dp
        yy_s = 0.0_dp
        zz_s = 0.0_dp
        zz0_s = 0.0_dp
    ELSE
        nbes    = Solver % Mesh % NumberOfBoundaryElements
        nbe     = Solver % Mesh % NumberOfBoundaryElements
        AR      = AR_S
        BR      = BR_S
        AR0     = AR0_S
        BR0     = BR0_S
        countit = countit_S
    END IF
!
!
! Summing up viscosity/ratefactor values
    IF(Parallel)THEN
    countit  = 0
    AR       = 0.0_dp
    BR       = 0.0_dp
    AR0      = 0.0_dp
    BR0      = 0.0_dp
      do s=1,nobs
        if(COUNTER(s).ge.1)THEN
          AA(s)    = AA(s)/COUNTER(s)
          AA0(s)   = AA0(s)/COUNTER(s)
          if(AA0(s).lt.visc_threshold**(-nflow))then!/1000.0)then
            SLOPE(s) = SLOPE(s)/COUNTER(s)
            FLUX(s)  = FLUX(s)/COUNTER(s)
            AR       = AR + AA(s)
            BR       = BR + AA(s)**(-1.0/nflow)
            AR0       = AR0 + AA0(s)
            BR0       = BR0 + AA0(s)**(-1.0/nflow)
            countit  = countit + 1
          endif
        endif
      enddo
    ENDIF
! 
!
! Averaging by total number
    if(countit.ge.1.0)then
      AR  = AR/countit
      BR  = BR/countit
      AR0 = AR0/countit
      BR0 = BR0/countit
    else
      AR = Arate
      BR = Arate**(-1.0/nflow)
      AR0 = Arate
      BR0 = Arate**(-1.0/nflow)
    endif

          !print *,'ladida',countit, AR**(-1.0/3.0),BR

!
! Write viscosity values for all observations
! (excluding unreasonable ratefactor values < 1e-10)
    IF (Parallel) THEN
      IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
        open(IO,file='./ice_viscosity.dat',status = 'replace', iostat = ok)
        do s=1, nobs
          if(AA(s).gt.0.and.AA(s).lt.visc_threshold**(-nflow))then
            WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xobs(s,1),xobs(s,2),AA(s)**(-1.0/nflow)
          endif
        enddo
        close(IO)
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!
!                              BOUNDARY ELEMENT LOOP
!
!-------------------------------------------------------------------------------

    ! Re-open data fiel
    IF(.NOT.Parallel)THEN
      if(FirstTime)then
        open(IO2,file='./ice_viscosity_mean.dat',status = 'unknown', position = 'rewind', iostat = ok)
      else
        open(IO2,file='./ice_viscosity_mean.dat',status = 'replace', iostat = ok)
      endif
    ENDIF

    ! Loop boundary elements
    DO t=1, Solver % Mesh % NumberOfBoundaryElements
      ! get element information
      Element => GetBoundaryElement(t)
      !print *,'Element'
      ! cycle non-active elements
      IF ( .NOT.ActiveBoundaryElement( Element ) ) CYCLE
      !print *,'Acitve Element'
      ! cycle halo elements
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      !print *,'Sort out halos'
      n = GetElementNOFNodes( Element )
      !print *,'Number of Nodes',n
      NodeIndexes => Element % NodeIndexes
      !print *,'Nodeindexes'

      ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      IF (DIM == 1) THEN !1D SSA
        ElementNodes % y(1:n) = 0.0_dp
        ElementNodes % z(1:n) = 0.0_dp
      ELSE IF (DIM == 2) THEN !2D SSA
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = 0.0_dp
      ELSE
        WRITE(Message,'(a,i1,a)')&
              'It is not possible to compute SSA problems with DOFs=',&
              DIM, ' . Aborting'
        CALL Fatal( SolverName, Message)
      END IF


      !print *, 'BE indices', t, Solver % Mesh % NumberOfBoundaryElements, NodeIndexes

      
      BC => GetBC()
      bc_id = GetBCId( Element )

      IF ( CalvingFlag ) THEN
        boundaryIDrange=2
      ELSE
        boundaryIDrange=1
      ENDIF

      IF(bc_id.LE.boundaryIDrange)THEN
          ! boundary ID 1 : always land terminated margin + holes
          ! boundary ID 2 : for land-terminating glaciers: nodes at the measurement location
          !                 for marine-terminating glaciers: calving front

            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )
           ! HMAR = 5.0
           ! DSDX =  SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
           ! DSDY =  SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
           ! !!DSDXX = SUM(GradSlopeXValues(2*(GradSlopeXPerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
           ! !!DSDYY = SUM(GradSlopeYValues(2*(GradSlopeYPerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
           ! !!FLUX_S(s) = SUM(FluxValues(FluxPerm(NodeIndexes(1:n)))*basis(1:n))

           ! DIST2MAR = SUM(Dist2MarValues(Dist2MarPerm(NodeIndexes(1:n)))*basis(1:n))
           ! !!DIST2OBS = SUM(Dist2ObsValues(Dist2ObsPerm(NodeIndexes(1:n)))*basis(1:n))
           ! !!SURF = SUM(SurfaceValues(SurfacePerm(NodeIndexes(1:n)))*basis(1:n))
           ! !HMAR = SUM(ThiValues(ThiPerm(NodeIndexes(1:n)))*basis(1:n))

           ! SLOPE_S(t) = (DSDX**2.0+DSDY**2.0)**0.5
           ! IF (SLOPE_S(t).LT.TAN(slope_threshold/180*3.14159)) THEN
           !   SLOPE_S(t)=TAN(slope_threshold/180*3.14159)
           ! ENDIF

           !! SLOPE CORRECTION FACTOR
           !! applicable to ice viscosities
           !! calibrated in the Alps by Christian Sommer
           !if(180.0/3.14159*ATAN(SLOPE_S(t)).LT.visc_slope_thresh+1.0/visc_slope_grad)then
           !  !rf_scaling = (visc_slope_grad*(180.0/3.14159*ATAN(SLOPE_S(t))-visc_slope_thresh))
           !  visc_slope_scaling = (visc_slope_grad*(180.0/3.14159*ATAN(SLOPE_S(t))-visc_slope_thresh))
           !else
           !  !rf_scaling = 1.0
           !  visc_slope_scaling = 1.0
           !endif


           !! WEIGHTING ACCORDING TO DISTANCE TO GLACIER MARGIN/OUTLINE
           !! account for marine termini
           !IF (DIST2MAR.GE.0.AND.HMAR.GT.0) THEN
           !  !visc_dist2mar_scaling = (atan(1.0*DIST2MAR/(10.0*coupling_length*HOBS))*2.0/3.141592653589793)+1.0e-3
           !  visc_dist2mar_scaling = (atan(10.0*DIST2MAR/(10.0*HMAR))*2.0/3.141592653589793)+1.0e-3
           !ELSE
           !  visc_dist2mar_scaling = 1.0e-3
           !ENDIF
           !
           !!visc_dist2mar_scaling = 1.0
           !visc_dist2obs_scaling = 1.0

        IF(PARALLEL)THEN
          if(Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
            ii = t
          else
            ii = sum(nbes(1:Solver % Matrix % ParMatrix % ParEnv % MyPE))+t
          endif
          !visc_scaling_margin_s(ii) = (visc_dist2mar_scaling*(visc_slope_scaling)**visc_dist2obs_scaling)
          !visc_scaling_margin_s(ii) = 1.0 !visc_dist2mar_scaling
          xx_s(ii) = ElementNodes % x(1)
          yy_s(ii) = ElementNodes % y(1)
          IF (ViscAverage) THEN
            !viscosity variant
            !print *,'whisky',countit, BR
            zz0_s(ii) = BR
            zz_s(ii) = BR0
            !/visc_scaling_margin_s(ii)
          ELSE
            !ratefactor variant
            !print *,'richy',countit, AR**(-1.0/3.0)
            zz0_s(ii) = AR**(-1.0/nflow)
            zz_s(ii) = AR0**(-1.0/nflow)
            !/visc_scaling_margin_s(ii)
          ENDIF
        ELSE
          ii = t
          !visc_scaling_margin_s(ii) = (visc_dist2mar_scaling*(visc_slope_scaling)**visc_dist2obs_scaling)
          !visc_scaling_margin_s(ii) = 1.0 !visc_dist2mar_scaling
          IF (ViscAverage) THEN
            !Viscosity variant
            zz0_s(ii) = BR
            zz_s(ii) = BR0
            !/visc_scaling_margin_s(ii)
          ELSE
            !Ratefactor variant
            zz0_s(ii) = AR**(-1.0/nflow)
            zz_s(ii) = AR0**(-1.0/nflow)
            !/visc_scaling_margin_s(ii)
          ENDIF
          if(zz_s(ii)**(-nflow).gt.0.and.zz_s(ii)**(-nflow).lt.visc_threshold**(-nflow))then
            WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') ElementNodes % x(1), ElementNodes % y(1),zz_s(ii)
            WRITE(IO2,'(2(E16.8,1x),1(E16.8,1x))') ElementNodes % x(1), ElementNodes % y(1),zz_s(ii)
          endif

        ENDIF
      ENDIF
    enddo

    IF (Parallel) THEN
      CALL MPI_ALLREDUCE(xx_s,xx,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(yy_s,yy,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(zz_s,zz,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(zz0_s,zz0,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
        open(IO,file='./ice_viscosity.dat',status = 'old', position = 'append', iostat = ok)
        open(IO2,file='./ice_viscosity_mean.dat',status = 'replace', iostat = ok)
        do ii=1, nbe
          !print *,'ZZ',zz(ii),zz0(ii)
          if(zz(ii)**(-nflow).gt.0.and.zz(ii)**(-nflow).lt.visc_threshold**(-nflow))then
            !print *,'final', zz(ii)
            WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xx(ii),yy(ii),zz(ii)
            WRITE(IO2,'(2(E16.8,1x),1(E16.8,1x))') xx(ii),yy(ii),zz(ii)
          endif
        enddo
      deallocate(xx,yy,zz,zz0)
      ENDIF
    deallocate(xx_s,yy_s,zz_s,zz0_s)
    ENDIF
    close(IO)
    close(IO2)


    deallocate(xobs,THIobs,InElement)
    deallocate(AA,AA_S,BB_S)
    deallocate(AA0,AA0_S,BB0_S)
    deallocate(ElementNodes % x,ElementNodes % y,ElementNodes % z)
    deallocate(nbes)
    deallocate(visc_scaling_margin,visc_scaling_margin_s)

    Firsttime=.false.

END SUBROUTINE computeSIA_Hobs

