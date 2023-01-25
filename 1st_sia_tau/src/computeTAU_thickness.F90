!
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE computeTAU_thickness( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,SlopeVarSolName
!  CHARACTER(LEN=MAX_NAME_LEN) :: TauVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: RefTauVarSolName,EffTauVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: Dist2ObsVarSolName,SurfaceVarSolName,Dist2MarVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: GradSlopeXSolName, GradSlopeYSolName
  TYPE(Variable_t), POINTER :: RefTauVariable, EffTauVariable
  TYPE(Variable_t), POINTER :: Variable,SlopeVariable,Hsi !,TauVariable
  TYPE(Variable_t), POINTER :: GradSlopeXVariable, GradSlopeYVariable
  TYPE(Variable_t), POINTER :: Dist2ObsVariable,SurfaceVariable,Dist2MarVariable
  REAL(KIND=dp), POINTER :: RefTauValues(:), EffTauValues(:)
  REAL(KIND=dp), POINTER :: Values(:),SlopeValues(:) !,TauValues(:)
  REAL(KIND=dp), POINTER :: GradSlopeXValues(:), GradSlopeYValues(:)
  REAL(KIND=dp), POINTER :: Dist2ObsValues(:),SurfaceValues(:),Dist2MarValues(:)
  INTEGER, POINTER :: RefTauPerm(:), EffTauPerm(:)
  INTEGER, POINTER :: Perm(:),SlopePerm(:) !,TauPerm(:)
  INTEGER, POINTER :: GradSlopeXPerm(:),GradSlopeYPerm(:)
  INTEGER, POINTER :: Dist2ObsPerm(:),SurfacePerm(:),Dist2MarPerm(:)
!
  INTEGER :: ok
  INTEGER,SAVE :: nobs
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),THIobs(:)
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
  real(kind=dp) :: DSDX,DSDY,THI,SLOPE,HOBS,TAU_REF,TAU_EFF
  real(kind=dp) :: DSDXX,DSDYY,ALPHA,BETA,OMEGA,MP
  real(kind=dp) :: SURF,THI0,DIST2OBS,DIST2MAR
  real(kind=dp) :: rhoice,grav,secy
  real(kind=dp) :: a1
  real(kind=dp) :: slope_threshold,slope_div_threshold,slope_div_trans,slope_factor
  !real(kind=dp) :: ys_slope_thresh, ys_slope_grad, rf_scaling, dist_scaling
  !real(kind=dp) :: ys_ref_elev, ys_el_grad
  real(kind=dp) :: tau_slope_thresh,tau_slope_grad
  real(kind=dp) :: tau_ref_el_z,tau_el_grad, el_min, el_max, el_z
  real(kind=dp) :: coupling_length
  real(kind=dp) :: tau_slope_scaling, tau_el_scaling, tau_scaling
  real(kind=dp) :: tau_dist2mar_scaling, tau_dist2obs_scaling
  integer :: i,j,k,l,n,m,DIM,countit_S,countit,s
  integer :: nflow,ierr
!
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
!
  INTEGER,PARAMETER :: IO=12
!
  save Firsttime,Parallel
  save SolverName,VarSolName,SlopeVarSolName
  save GradSlopeXSolName,GradSlopeYSolName,Dist2ObsVarSolName,SurfaceVarSolName,Dist2MarVarSolName


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

!! Central Variable
  VarSolName =  GetString( SolverParams,'Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Variable Name< not found in section >Solver<')
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

!! Name of referene yield stress variable
  RefTauVarSolName =  GetString( SolverParams,'Reference Yield stress Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Reference Yield Stress Variable Name< not found in section >Solver<')
      END IF

!! Name of effective yield stress variable
  EffTauVarSolName =  GetString( SolverParams,'Effective Yield stress Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Effective Yield Stress Variable Name< not found in section >Solver<')
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

!! Name of variabe of margin distance
  Dist2MarVarSolName =  GetString( SolverParams,'Margin Distance Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Margin Distance Variable Name< not found in section >Solver<')
      END IF

!! Name of variabe for distance to next observational point
  Dist2ObsVarSolName =  GetString( SolverParams,'Obs Distance Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Obs Distance Variable Name< not found in section >Solver<')
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

     WRITE(SolverName, '(A)') 'computeTAU_thickness'

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
     rhoice = rhoice*(1.0e6*(secy)**2)

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

     ! SLOPE DEPENDENT YIELD STRESS
     ! Define slope threshold (beyond which no scaling is applied)
     tau_slope_thresh = GetConstReal( Model % Constants, 'yield stress slope threshold', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<yield stress slope threshold> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            tau_slope_thresh = 0.0
     End if

     ! slope gradient for yield stress computation
     tau_slope_grad = GetConstReal( Model % Constants, 'yield stress slope gradient', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<yield stress slope threshold> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            tau_slope_grad = 0.0
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


     ! ELEVATION DEPENDENT YIELD STRESS
     ! Define reference elevation where scaling is 1.0
     tau_ref_el_z = GetConstReal( Model % Constants, 'yield stress reference elevation', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<yield stress reference elevation> not found. &
                   &Setting to default (0.0).'
            CALL INFO(SolverName, Message, level=20)
            tau_ref_el_z = 0.0
     End if

     ! slope gradient for yield stress computation
     tau_el_grad = GetConstReal( Model % Constants, 'yield stress elevation gradient', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '<yield stress elevation gradient> not found. &
                   &Setting to default (0.0). No slope dependence.'
            CALL INFO(SolverName, Message, level=20)
            tau_el_grad = 0.0
     End if

     ! coupling length scale for glaciers and ice-sheets (effect of longitudinal stress coupling)
     coupling_length = GetConstReal( Model % Constants, 'coupling length', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant >Coupling Length< not found. &
                   &Setting to 10.0'
            CALL INFO(SolverName, Message, level=20)
            coupling_length=10.0_dp
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

     WRITE(SolverName, '(A)') 'computeSIAthickness'

  !!! End of First visit
    Firsttime=.false.
  Endif

!-------------------------------------------------------------------------------
!
!                              GET VARIABLES
!
!-------------------------------------------------------------------------------

    Variable => VariableGet( Solver % Mesh % Variables, VarSolName  )
    IF ( ASSOCIATED( Variable ) ) THEN
            Values => Variable % Values
            Perm => Variable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                    'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

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

    RefTauVariable => VariableGet( Solver % Mesh % Variables, RefTauVarSolName  )
    IF ( ASSOCIATED( RefTauVariable ) ) THEN
            RefTauValues => RefTauVariable % Values
            RefTauPerm => RefTauVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',RefTauVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    EffTauVariable => VariableGet( Solver % Mesh % Variables, EffTauVarSolName  )
    IF ( ASSOCIATED( EffTauVariable ) ) THEN
            EffTauValues => EffTauVariable % Values
            EffTauPerm => EffTauVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',EffTauVarSolName,' < found'
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

!    N = model % MaxElementNodes
!    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
!
!! Get the obs
!   ObsFileName =  GetString( SolverParams,'Observation File Name', Found)
!   IF(.NOT.Found) THEN
!       CALL FATAL(SolverName,'Keyword >Observation File Name< not found in section >Solver<')
!   END IF
!
!   open(IO,file=trim(ObsFileName),status = 'old',iostat = ok)
!   if(ok /= 0) then
!       write(message,'(A,A)') 'Unable to open file ',TRIM(ObsFileName)
!       CALL Fatal(Trim(SolverName),Trim(message))
!   end if
!   nobs=0
!   do while(ok == 0)
!     read(io,*,iostat = ok)
!     if (ok == 0) nobs = nobs + 1
!   end do
!   close(IO)
!
!
!   allocate(xobs(nobs,3),THIobs(nobs),InElement(nobs))
!   InElement(:)=-1
!
!   THIobs=0.0_dp
!   xobs=0.0_dp
!   open(IO,file=trim(ObsFileName))
!   do i=1,nobs
!     read(IO,*) (xobs(i,j),j=1,DIM),(THIobs(i))
!   end do
!   close(IO)

!   open(IO,file='./ratefactor.csv',status = 'replace',iostat = ok)
!   close(IO)

!-------------------------------------------------------------------------------
!
!                              START OBSERVATION LOOP
!
!-------------------------------------------------------------------------------

!   CALL StartAdvanceOutput(SolverName,'Loop through thickness observations to determine best rate factor')
!    AR_S = 0.0_dp
!    countit_S = 0
!    Do s=1,nobs
!
!     CALL AdvanceOutput(s,nobs)
!
!     !IF (FirstRound) then
!      !Need to find in which Element the data point resides
!      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
!      Coord=0._dp
!      Coord(1:DIM)=xobs(s,1:DIM)
!      CALL LocateParticleInMeshOctree( ElementIndex,Coord)
!      If (ElementIndex.NE.0) InElement(s)=ElementIndex
!     !ENDIF !End if FirstRound
!
!    ! Data Point has been found in one element
!      IF (InElement(s)>0) THEN
!         Element => GetActiveElement(InElement(s))
!         n = GetElementNOFNodes()
!         NodeIndexes => Element % NodeIndexes
!    ! set coords of highest occurring dimension to zero (to get correct path
!    ! element)
!
!          !-------------------------------------------------------------------------------
!         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
!         IF (DIM == 1) THEN !1D SSA
!            ElementNodes % y(1:n) = 0.0_dp
!            ElementNodes % z(1:n) = 0.0_dp
!         ELSE IF (DIM == 2) THEN !2D SSA
!            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
!            ElementNodes % z(1:n) = 0.0_dp
!         ELSE
!            WRITE(Message,'(a,i1,a)')&
!                'It is not possible to compute SSA problems with DOFs=',&
!                DIM, ' . Aborting'
!            CALL Fatal( SolverName, Message)
!         END IF
!
!         IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:3),UVW))  THEN
!              CALL FATAL(SolverName,'Point was supposed to be found in this element')
!         ELSE
!            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
!                              Basis,dBasisdx )
!            HOBS = THIobs(s)
!            DSDX = SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+1)*basis(1:n))
!            DSDY = SUM(SlopeValues(2*(SlopePerm(NodeIndexes(1:n))-1)+2)*basis(1:n))
!
!            FLUX =  SUM(FluxValues(FluxPerm(NodeIndexes(1:n)))*basis(1:n))
!
!            SLOPE = (DSDX**2.0+DSDY**2.0)**0.5
!
!            IF (SLOPE.LT.TAN(slope_threshold/180*3.14159)) THEN
!            SLOPE=TAN(slope_threshold/180*3.14159)
!            ENDIF
!
!!
!            a1 = abs(FLUX/(SLOPE**(nflow)))
!            a2 = (1.0*nflow+2.0)/(2.0)*1.0/((rhoice*abs(grav))**nflow)
!            a3 = 1.0/(HOBS)**(nflow+2.0)
!
!            if(a1*a2*a3.lt.1e-13.and.a1*a2*a3.gt.0)then
!            countit_S = countit_S + 1
!            AR_S = AR_S+a1*a2*a3
!            endif
!!
!          END IF
!
!       ELSE
!
!            WRITE(Message,'(a,I0,a)')&
!                'Data Point',s,'found in no element'
!            CALL Info( SolverName, Message,level=15)
!       END IF
!    enddo ! s
!
!    IF (Parallel) THEN
!        CALL MPI_ALLREDUCE(AR_S,AR,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!        CALL MPI_ALLREDUCE(countit_S,countit,1,&
!               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!    ELSE
!        AR      = AR_S
!        countit = countit_S
!    END IF
!
!    if(countit.ge.1.0)then
!    AR = AR/countit
!    else
!    AR = Arate
!    endif
!    !print *,'Part rate factor : ', AR_S,countit_S,countit
!    print *,'Ratefactor: ',Arate,AR,countit
!    !print *,'Rate factor :',AR,' Pa^3 yr'



!-------------------------------------------------------------------------------
!
!                              START LOOP
!
!-------------------------------------------------------------------------------


    N = Solver%Mesh%NumberOfNodes
    Do i=1,N

        DSDX  = SlopeValues(2*(SlopePerm(i)-1)+1)
        DSDY  = SlopeValues(2*(SlopePerm(i)-1)+2)
        TAU_REF = RefTauValues(RefTauPerm(i))
        TAU_EFF = EffTauValues(EffTauPerm(i))
        DSDXX = GradSlopeXValues(2*(GradSlopeXPerm(i)-1)+1)
        DSDYY = GradSlopeYValues(2*(GradSlopeYPerm(i)-1)+2)
        DIST2MAR = Dist2MarValues(Dist2MarPerm(i))
        DIST2OBS = Dist2ObsValues(Dist2ObsPerm(i))
        SURF     = SurfaceValues(SurfacePerm(i))

        THI0 = Values(Perm(i))

        SLOPE = (DSDX**2.0+DSDY**2.0)**0.5

        IF (SLOPE.LT.TAN(slope_threshold/180*3.14159)) THEN
        SLOPE=TAN(slope_threshold/180*3.14159)
        ENDIF

        !! Correction strategy for divides and ridges
        !! proved not longer suitable on Vernagt and Hofsjokull
        !slope_div_threshold = 5e-4
        !slope_div_trans     = 2.0
        !slope_factor        = 10.0
        !IF (SLOPE.LT.TAN(slope_factor*slope_threshold/180*3.14159))THEN
        !  IF (DSDXX+DSDYY>slope_div_threshold) THEN
        !    OMEGA = slope_factor*slope_threshold
        !    SLOPE = TAN(OMEGA/180*3.14159)
        !  ELSE IF ((DSDXX+DSDYY)>slope_div_threshold/slope_div_trans) THEN
        !    OMEGA = slope_factor*slope_threshold
        !    MP = ((DSDXX+DSDYY)-slope_div_threshold/slope_div_trans)/((1.0-1.0/slope_div_trans)*slope_div_threshold)
        !    !MP = ATAN(5e-5/(DSDXX+DSDYY)-1.0)*2.0/3.14159
        !    SLOPE=MP*TAN(OMEGA/180*3.14159)+(1.0-MP)*SLOPE
        !  ENDIF
        !ENDIF

        THI = 0.0

        IF (TAU_REF.GT.0.0) THEN

           ! SLOPE CORRECTION FACTOR 
           ! applicable to yield stress
           IF(180.0/3.14159*ATAN(SLOPE).LT.tau_slope_thresh+1.0/tau_slope_grad)THEN
             tau_slope_scaling = (tau_slope_grad*(180.0/3.14159*ATAN(SLOPE)-tau_slope_thresh))
           ELSE
             tau_slope_scaling = 1.0
           ENDIF

           ! ELEVATION DEPENDENT SCALING
           IF (SURF.gt.el_min.and.SURF.lt.el_max) THEN
             !rf_scaling = rf_scaling*(1.0+tau_el_grad*(SURF-tau_ref_elev))
             !tau_el_scaling = (1.0+tau_el_grad*(SURF-tau_ref_elev))
             el_z            = (SURF-el_min)/(el_max-el_min)
             tau_el_scaling = 1.0+tau_el_grad*(el_z-tau_ref_el_z)
           ELSE
             tau_el_scaling = 1.0
           ENDIF

           IF (SURF.LT.el_min) THEN
             el_z            = 0
             tau_el_scaling = 1.0+tau_el_grad*(el_z-tau_ref_el_z)
           END IF

           IF (SURF.GT.el_max) THEN
             el_z            = 1.0
             tau_el_scaling = 1.0+tau_el_grad*(el_z-tau_ref_el_z)
           END IF


           ! WEIGHTING ACCORDING TO DISTANCE TO NEXT THICKNESS OBSERVATIONS
           ! if there are no observations; distances are set high --> no effect
           IF (DIST2OBS.GT.0)THEN
             tau_dist2obs_scaling = atan(1.0*DIST2OBS/(10.0*THI0))*2.0/3.141592653589793
           ELSE
             tau_dist2obs_scaling = 1.0
           ENDIF

           ! WEIGHTING ACCORDING TO DISTANCE TO GLACIER MARGIN/OUTLINE
           ! account for marine termini
           IF (DIST2MAR.GE.0.AND.THI0.GT.0) THEN
             tau_dist2mar_scaling = (atan(1.0*DIST2MAR/(5.0*THI0))*2.0/3.141592653589793)+1.0e-3
             !tau_dist2mar_scaling = (atan(1.0*DIST2MAR/(0.5*THI0))*2.0/3.141592653589793)+1.0e-3
           ELSE
             tau_dist2mar_scaling = 1.0e-3
           ENDIF

           !tau_dist2mar_scaling = 1.0
           tau_dist2obs_scaling = 1.0

           tau_scaling = tau_dist2mar_scaling*(tau_el_scaling*tau_slope_scaling)**tau_dist2obs_scaling
           !tau_scaling = (tau_el_scaling*tau_slope_scaling)**tau_dist2obs_scaling


           !! SLOPE CORRECTION FACTOR 
           !if(180.0/3.14159*ATAN(SLOPE).LT.tau_slope_thresh+1.0/tau_slope_grad)then
           !  rf_scaling = (tau_slope_grad*(180.0/3.14159*ATAN(SLOPE)-tau_slope_thresh))
           !else
           !  rf_scaling = 1.0
           !endif
           !
           !! ELEVATION DEPENDENT SCALING
           !rf_scaling = rf_scaling*(1.0+tau_el_grad*(SURF-tau_ref_elev))
           !
           !
           !! WEIGHTING ACCORDING TO DISTANCE TO NEXT THICKNESS OBSERVATIONS
           !! if there are no observations; distances are set high --> no effect
           !IF (DIST2OBS.GT.0)THEN
           !dist_scaling  = atan(DIST2OBS/(coupling_length*THI0))*2.0/3.141592653589793
           !ELSE
           !dist_scaling = 1.0
           !ENDIF
           !
           !dist_scaling = 1.0

           !if (DIST.GT.0.0) then
           !dist_scaling  = atan(DIST/(1000.0))*2.0/3.141592653589793
           !TAU_REF = rf_scaling**dist_scaling*TAU
           TAU_EFF = tau_scaling*TAU_REF
           !print *, 'SCALING',rf_scaling, tau_scaling

           IF(SLOPE.GT.0.AND.DIST2MAR.GT.0)THEN
           !!a1  = atan(DIST2MAR/50.0)*2.0/3.141592653589793
           !!THI = TAU_EFF/(rhoice*abs(grav)*SLOPE)*a1
           THI = TAU_EFF/(rhoice*abs(grav)*SLOPE)
           !a1  = atan(DIST2MAR/(0.5*THI0))*2.0/3.141592653589793+1.e-3
           !print *,'DIST2MAR',a1,tau_dist2mar_scaling
           !THI = THI*a1
           ELSE
           THI = 0.0
           ENDIF

        ENDIF

!
        Values(Perm(i)) = THI

        EffTauValues(EffTauPerm(i)) = TAU_EFF

    End do

    !deallocate(xobs,THIobs,InElement)
    !deallocate(ElementNodes % x,ElementNodes % y,ElementNodes % z)

    Firsttime=.false.

END SUBROUTINE computeTAU_thickness

