!
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE computeTAU_Hobs( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,SlopeVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: GradSlopeXSolName, GradSlopeYSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: SurfaceVarSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: Dist2MarVarSolName, Dist2ObsVarSolName
  CHARACTER(LEN=10) :: GlacierType
  TYPE(Variable_t), POINTER :: Variable,SlopeVariable,Hsia
  TYPE(Variable_t), POINTER :: GradSlopeXVariable, GradSlopeYVariable
  TYPE(Variable_t), POINTER :: SurfaceVariable
  TYPE(Variable_t), POINTER :: Dist2MarVariable, Dist2ObsVariable
  REAL(KIND=dp), POINTER :: Values(:),SlopeValues(:)
  REAL(KIND=dp), POINTER :: GradSlopeXValues(:), GradSlopeYValues(:)
  REAL(KIND=dp), POINTER :: SurfaceValues(:)
  REAL(KIND=dp), POINTER :: Dist2MarValues(:), Dist2ObsValues(:)
  INTEGER, POINTER :: Perm(:),SlopePerm(:)
  INTEGER, POINTER :: GradSlopeXPerm(:),GradSlopeYPerm(:)
  INTEGER, POINTER :: SurfacePerm(:)
  INTEGER, POINTER :: Dist2MarPerm(:), Dist2ObsPerm(:)
  INTEGER :: ok
  INTEGER :: boundaryIDrange
  INTEGER,SAVE :: nobs,nbe,nbe_s
  INTEGER, ALLOCATABLE :: nbes(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),THIobs(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: TAU0_S(:), TAU0(:),TAU_S(:), TAU(:)
  REAL(KIND=dp),ALLOCATABLE :: xx(:),yy(:),zz(:),xx_s(:),yy_s(:),zz_s(:),SLOPE(:),SLOPE_S(:)
  REAL(KIND=dp) :: SURF,DIST2MAR,DIST2OBS
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
  real(kind=dp) :: DSDX,DSDY,THI,DSDXX,DSDYY
  real(kind=dp) :: ALPHA,BETA,OMEGA,MP
  real(kind=dp) :: HOBS,meanTAU,sumTAU,meanTAU0,sumTAU0
  real(kind=dp) :: rhoice,grav,secy,tau_ref
  real(kind=dp) :: slope_threshold,slope_div_threshold,slope_div_trans,slope_factor
  real(kind=dp) :: tau_slope_thresh,tau_slope_grad
  real(kind=dp) :: tau_ref_el_z,tau_el_grad, el_min, el_max, el_z
  real(kind=dp) :: coupling_length
  real(kind=dp) :: tau_slope_scaling, tau_el_scaling, tau_scaling
  real(kind=dp) :: tau_dist2mar_scaling, tau_dist2obs_scaling
  real(kind=dp) :: a1,a2,a3,b1,b2
  integer :: i,j,k,l,n,m,DIM,s,countit,countit_S,t,ii
  integer :: nflow,ierr,nPar,NNN,NNN_S
!
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  Logical :: ViscAverage,CalvingFlag

  INTEGER,PARAMETER :: IO=12,IO2=13,IO3=14

  save Firsttime,Parallel,nPar
  save SolverName,VarSolName,SlopeVarSolName
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
!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

!  VarSolName =  GetString( SolverParams,'Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Variable< not found in section >Solver<')
!      END IF

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
    N = Solver%Mesh%NumberOfNodes

!!!!!!! Check for parallel run 
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

  !!! End of First visit
!    Firsttime=.false.
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

     ! slope threshold ... if below SIA thickness is taken
     slope_threshold = GetConstReal( Model % Constants, 'slope threshold', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant slope threshold not found. &
                   &Setting to default (0.0)'
            CALL INFO(SolverName, Message, level=20)
            slope_threshold = 0.0
     End if

     ! yield stress
     tau_ref = GetConstReal( Model % Constants, 'yield stress', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Yield Stress not found. &
                   &Setting to 1.0e+5 Pa'
            CALL INFO(SolverName, Message, level=20)
            tau_ref = 1.0e5
     End if

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


     WRITE(SolverName, '(A)') 'computeSIAthickness'

  !!! End of First visit
!    Firsttime=.false.
  Endif


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


   allocate(xobs(nobs,3),THIobs(nobs),InElement(nobs),SLOPE(nobs),SLOPE_S(nobs))
   allocate(TAU(nobs),TAU_S(nobs))
   allocate(TAU0(nobs),TAU0_S(nobs))
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
!                              START OBSERVATION LOOP
!
!-------------------------------------------------------------------------------

    IF(.NOT.Parallel)THEN
    if(FirstTime)then
    open(IO,file='./yield_stress.dat',status = 'unknown', position = 'rewind', iostat = ok)
    !open(IO3,file='./ratefactor_init.csv',status = 'unknown', position = 'rewind', iostat = ok)
    else
    open(IO,file='./yield_stress.dat',status = 'replace', iostat = ok)
    endif
    ENDIF

   CALL StartAdvanceOutput(SolverName,'Loop through thickness observations to determine best rate factor')
    countit_S = 0
    SLOPE_S   = 0.0_dp
    COUNTER_S = 0
    TAU_S     = 0.0_dp
    TAU0_S    = 0.0_dp
    sumTAU    = 0.0_dp
    !TAU       
    
!    IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
!    open(IO,file='./ratefactor.csv',status = 'old', position = 'append', iostat = ok)
    Do s=1,nobs
     TAU_S(s)  = 0.0_dp
     TAU0_S(s) = 0.0_dp
     CALL AdvanceOutput(s,nobs)

     !IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(s,1:DIM)
      CALL LocateParticleInMeshOctree( ElementIndex,Coord)
      If (ElementIndex.NE.0) InElement(s)=ElementIndex
     !ENDIF !End if FirstRound

    ! Data Point has been found in one element
      IF (InElement(s)>0) THEN
         Element => GetActiveElement(InElement(s))
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
    ! set coords of highest occurring dimension to zero (to get correct path
    ! element)

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

            DIST2MAR = SUM(Dist2MarValues(Dist2MarPerm(NodeIndexes(1:n)))*basis(1:n))
            DIST2OBS = SUM(Dist2ObsValues(Dist2ObsPerm(NodeIndexes(1:n)))*basis(1:n))
            SURF = SUM(SurfaceValues(SurfacePerm(NodeIndexes(1:n)))*basis(1:n))

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
         !    !MP = ATAN(5e-5/(DSDXX+DSDYY)-1.0)*2.0/3.14159
         !    SLOPE_S(s)=MP*TAN(OMEGA/180*3.14159)+(1.0-MP)*SLOPE_S(s)
         !  ENDIF
         !ENDIF

           ! SLOPE CORRECTION FACTOR 
           ! applicable to yield stress
           if(180.0/3.14159*ATAN(SLOPE_S(s)).LT.tau_slope_thresh+1.0/tau_slope_grad)then
             tau_slope_scaling = (tau_slope_grad*(180.0/3.14159*ATAN(SLOPE_S(s))-tau_slope_thresh))
           else
             tau_slope_scaling = 1.0
           endif

           ! ELEVATION DEPENDENT SCALING
           if (SURF.gt.el_min.and.SURF.lt.el_max) then
             !rf_scaling = rf_scaling*(1.0+tau_el_grad*(SURF-tau_ref_elev))
             !tau_el_scaling = (1.0+tau_el_grad*(SURF-tau_ref_elev))
             el_z            = (SURF-el_min)/(el_max-el_min)
             tau_el_scaling = 1.0+tau_el_grad*(el_z-tau_ref_el_z)
           else
             tau_el_scaling = 1.0
           end if

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
           IF (HOBS.GT.0)THEN
             tau_dist2obs_scaling  = atan(1.0*DIST2OBS/(10.0*HOBS))*2.0/3.141592653589793
           ELSE
             tau_dist2obs_scaling = 1.0
           ENDIF

           ! WEIGHTING ACCORDING TO DISTANCE TO GLACIER MARGIN/OUTLINE
           ! account for marine termini
           IF (DIST2MAR.GE.0.AND.HOBS.GT.0) THEN
             tau_dist2mar_scaling = (atan(1.0*DIST2MAR/(5.0*HOBS))*2.0/3.141592653589793)+1.0e-3
             !tau_dist2mar_scaling = (atan(1.0*DIST2MAR/(0.5*HOBS))*2.0/3.141592653589793)+1.0e-3
           ELSE
             tau_dist2mar_scaling = 1.0e-3
           ENDIF

           !tau_dist2mar_scaling = 1.0
           tau_dist2obs_scaling = 1.0

           tau_scaling = (tau_dist2mar_scaling*(tau_el_scaling*tau_slope_scaling)**tau_dist2obs_scaling)

               COUNTER_S(s) = COUNTER_S(s)+1

               IF(HOBS.GT.0.AND.DIST2MAR.GT.0)THEN
                 !a1        = atan(DIST2MAR/50.0)*2.0/3.141592653589793
                 a1        = atan(DIST2MAR/(0.5*HOBS))*2.0/3.141592653589793+1.0e-3
                 !TAU_S(s)  = rhoice*abs(grav)*HOBS*abs(SLOPE_S(s))/a1
                 TAU_S(s)  = rhoice*abs(grav)*HOBS*abs(SLOPE_S(s))/tau_scaling
                 TAU0_S(s) = rhoice*abs(grav)*HOBS*abs(SLOPE_S(s))/tau_dist2mar_scaling
                 countit_S = countit_S + 1
                 sumTAU    = sumTAU+TAU_S(s)
                 sumTAU0    = sumTAU0+TAU0_S(s)
               ELSE
                 TAU_S(s)  = 0.0
                 TAU0_S(s) = 0.0
               ENDIF
               
               IF(.NOT.Parallel)THEN
               WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xobs(s,1),xobs(s,2),TAU_S(s)
               ENDIF
!
            END IF

         ELSE

            WRITE(Message,'(a,I0,a)')&
                'Data Point',s,'found in no element'
            CALL Info( SolverName, Message,level=15)
         END IF
    enddo ! s


    IF (Parallel) THEN
        CALL MPI_ALLREDUCE(TAU_S,TAU,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(TAU0_S,TAU0,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(countit_S,countit,1,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(SLOPE_S,SLOPE,nobs,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(COUNTER_S,COUNTER,nobs,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(nbe_s,nbe,1,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLGATHER(nbe_s,1,MPI_INTEGER,nbes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        !CALL MPI_ALLGATHER(nbe_s,1,MPI_INTEGER,nbes,1,MPI_INTEGER,MPI_COMM_WORLD)
        allocate(xx_s(nbe),yy_s(nbe),zz_s(nbe),xx(nbe),yy(nbe),zz(nbe))
        xx_s = 0.0_dp
        yy_s = 0.0_dp
        zz_s = 0.0_dp
    ELSE
        nbes    = Solver % Mesh % NumberOfBoundaryElements
        nbe     = Solver % Mesh % NumberOfBoundaryElements
        countit = countit_S
    END IF

    IF(Parallel)THEN
    countit = 0
    sumTAU  = 0.0_dp
    sumTAU0 = 0.0_dp
    do s=1,nobs
        if(COUNTER(s).ge.1)THEN
            TAU(s)  = TAU(s)/COUNTER(s)
            TAU0(s) = TAU0(s)/COUNTER(s)            
            sumTAU  = sumTAU+TAU(s)
            sumTAU0 = sumTAU0+TAU0(s)
            countit = countit + 1
        endif
    enddo
    ENDIF
! 
!
    if(countit.ge.1.0)then
    meanTAU  = sumTAU/countit
    meanTAU0 = sumTAU0/countit
    else
    meanTAU  = tau_ref
    meanTAU0 = tau_ref
    endif

    IF (Parallel) THEN
    IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
    open(IO,file='./yield_stress.dat',status = 'replace', iostat = ok)
    do s=1, nobs
      WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xobs(s,1),xobs(s,2),TAU(s)
    enddo
    close(IO)
    ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!
!                              BOUNDARY ELEMENT LOOP
!
!-------------------------------------------------------------------------------

    IF(.NOT.Parallel)THEN
    if(FirstTime)then
    open(IO2,file='./yield_stress_mean.dat',status = 'unknown', position = 'rewind', iostat = ok)
    else
    open(IO2,file='./yield_stress_mean.dat',status = 'replace', iostat = ok)
    endif
    ENDIF

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

        BC => GetBC()
        bc_id = GetBCId( Element )


        IF ( CalvingFlag ) THEN
           boundaryIDrange=2
        ELSE
           boundaryIDrange=1
        ENDIF

        !IF(bc_id.EQ.1.or.bc_id.eq.2)THEN
        IF(bc_id.LE.boundaryIDrange)THEN
          ! boundary ID 1 : always land terminated margin + holes
          ! boundary ID 2 : for land-terminating glaciers: nodes at the measurement location
          !                 for marine-terminating glaciers: calving front

          IF(PARALLEL)THEN
            if(Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
            ii = t
            else
            ii = sum(nbes(1:Solver % Matrix % ParMatrix % ParEnv % MyPE))+t
            endif
            xx_s(ii) = ElementNodes % x(1)
            yy_s(ii) = ElementNodes % y(1)
            !zz_s(ii) = meanTAU
            zz_s(ii) = meanTAU0
          ELSE
            WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') ElementNodes % x(1), ElementNodes % y(1),meanTAU
            WRITE(IO2,'(2(E16.8,1x),1(E16.8,1x))') ElementNodes % x(1), ElementNodes % y(1),meanTAU

          ENDIF
        ELSE
        ENDIF

    enddo

    IF (Parallel) THEN
        CALL MPI_ALLREDUCE(xx_s,xx,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(yy_s,yy,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(zz_s,zz,nbe,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) THEN
      open(IO,file='./yield_stress.dat',status = 'old', position = 'append', iostat = ok)
      open(IO2,file='./yield_stress_mean.dat',status = 'replace', iostat = ok)
    do ii=1, nbe
        if(zz(ii).gt.0)then
            WRITE(IO,'(2(E16.8,1x),1(E16.8,1x))') xx(ii),yy(ii),zz(ii)
            WRITE(IO2,'(2(E16.8,1x),1(E16.8,1x))') xx(ii),yy(ii),zz(ii)
        endif
    enddo
    deallocate(xx,yy,zz)
    ENDIF
    deallocate(xx_s,yy_s,zz_s)
    ENDIF
    close(IO)
    close(IO2)


    deallocate(xobs,THIobs,InElement)
    deallocate(TAU,TAU_S,TAU0,TAU0_S)
    deallocate(ElementNodes % x,ElementNodes % y,ElementNodes % z)
    !deallocate(nbes,xx_s,yy_s,zz_s,xx,yy,zz)
    deallocate(nbes)

    Firsttime=.false.

END SUBROUTINE computeTAU_Hobs

