!
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE computeSIA_deltaflux( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarSolName,FluxVarSolName,deltaFluxVarSolName,SlopeVarSolName,ViscVarSolName,ThiUncVarSolName
  TYPE(Variable_t), POINTER :: Variable,FluxVariable,deltaFluxVariable,SlopeVariable,Hsia,ViscVariable,ThiUncVariable
  REAL(KIND=dp), POINTER :: Values(:),FluxValues(:),deltaFluxValues(:),SlopeValues(:),ViscValues(:),ThiUncValues(:)
  INTEGER, POINTER :: Perm(:),FluxPerm(:),deltaFluxPerm(:),SlopePerm(:),ViscPerm(:),ThiUncPerm(:)
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
  real(kind=dp) :: DSDX,DSDY,FLUX,THI,SLOPE,RF,VISC,HOBS,dFLUX,dTHI,THIunc
  real(kind=dp) :: MFLUX,MFLUX_S,kappa
  real(kind=dp) :: FLUX_MIN_S,FLUX_MIN
  real(kind=dp) :: AR,AR_S
  real(kind=dp) :: Arate,rhoice,grav,secy
  real(kind=dp) :: a1,a2,a3,b1,b2
  real(kind=dp) :: slope_threshold
  integer :: i,j,k,l,n,m,DIM,countit_S,countit,s
  integer :: nflow,ierr,NNN,NNN_S
!
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
!
  INTEGER,PARAMETER :: IO=12
!
  save Firsttime,Parallel
  save SolverName,VarSolName,FluxVarSolName,deltaFluxVarSolName,SlopeVarSolName,ViscVarSolName,ThiUncVarSolName


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

  VarSolName =  GetString( SolverParams,'Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Variable Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  FluxVarSolName =  GetString( SolverParams,'Flux Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Flux Variable Name< not found in section >Solver<')
      END IF

!! Name of thickness uncertainty variable
  ThiUncVarSolName =  GetString( SolverParams,'Thickness Uncertainty Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Thickness Uncertainty Name< not found in section >Solver<')
      END IF

!! Name of the variable to regularise
  SlopeVarSolName =  GetString( SolverParams,'Slope Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Slope Variable Name< not found in section >Solver<')
      END IF

!! Name of viscosity variable
  ViscVarSolName =  GetString( SolverParams,'Viscosity Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Viscosity Variable Name< not found in section >Solver<')
      END IF

!!! Name of rate factor variable
!  MRFVarSolName =  GetString( SolverParams,'Mean Ratefactor Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Mean Ratefactor Variable Name< not found in section >Solver<')
!      END IF

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

     WRITE(SolverName, '(A)') 'computeSIAthickness'

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


     ! observational error on ice thickness data
     dTHI = GetConstReal( Model % Constants, 'deltahobs', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') '> deltahobs< not found. &
                   &Setting to 0.0 m'
            CALL INFO(SolverName, Message, level=20)
            dTHI = 0.0
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

    FluxVariable => VariableGet( Solver % Mesh % Variables, FluxVarSolName  )
    IF ( ASSOCIATED( FluxVariable ) ) THEN
            FluxValues => FluxVariable % Values
            FluxPerm => FluxVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',FluxVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    ThiUncVariable => VariableGet( Solver % Mesh % Variables, ThiUncVarSolName  )
    IF ( ASSOCIATED( ThiUncVariable ) ) THEN
            ThiUncValues => ThiUncVariable % Values
            ThiUncPerm => ThiUncVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',ThiUncVarSolName,' < found'
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

    ViscVariable => VariableGet( Solver % Mesh % Variables, ViscVarSolName  )
    IF ( ASSOCIATED( ViscVariable ) ) THEN
            ViscValues => ViscVariable % Values
            ViscPerm => ViscVariable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                   'No variable >',ViscVarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF
!
!    MRFVariable => VariableGet( Solver % Mesh % Variables, MRFVarSolName  )
!    IF ( ASSOCIATED( MRFVariable ) ) THEN
!            MRFValues => MRFVariable % Values
!            MRFPerm => MRFVariable % Perm
!    ELSE
!            WRITE(Message,'(A,A,A)') &
!                   'No variable >',MRFVarSolName,' < found'
!            CALL FATAL(SolverName,Message)
!    END IF


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
!                              DETERMINE MEAN FLUX VALUE
!
!           NECESSARY FOR CORRECTING SMALL AND NEGATIVE FLUX VALUES BELOW
!
!-------------------------------------------------------------------------------

    N = Solver%Mesh%NumberOfNodes
    MFLUX_S = 0.0_dp
    FLUX_MIN_S = 0.0_dp
    NNN_S   = 0
    Do i=1,N

        MFLUX_S = MFLUX_S + abs(FluxValues(FluxPerm(i)))
        NNN_S   = NNN_S+1

    Enddo

    FLUX_MIN_S = minval(FluxValues)

    IF (Parallel) THEN
        CALL MPI_ALLREDUCE(NNN_S,NNN,1,&
               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(MFLUX_S,MFLUX,1,&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(FLUX_MIN_S,FLUX_MIN,1,&
               MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
        MFLUX = MFLUX/NNN
    ELSE
        MFLUX = MFLUX_S/NNN_S
        FLUX_MIN = FLUX_MIN_S
    ENDIF


!-------------------------------------------------------------------------------
!
!                              START LOOP
!
!-------------------------------------------------------------------------------

!    MRF  = MRFValues(MRFPerm(1))**(-3.0)

    N = Solver%Mesh%NumberOfNodes
    Do i=1,N

        !dFLUX = deltaFluxValues(deltaFluxPerm(i))
        FLUX = FluxValues(FluxPerm(i))
        DSDX = SlopeValues(2*(SlopePerm(i)-1)+1)
        DSDY = SlopeValues(2*(SlopePerm(i)-1)+2)
        VISC = ViscValues(ViscPerm(i))
        RF   = VISC**(-nflow)
        THIunc = ThiUncValues(ThiUncPerm(i))

        SLOPE = (DSDX**2.0+DSDY**2.0)**0.5

        IF (SLOPE.LT.TAN(slope_threshold/180*3.14159)) THEN
        SLOPE=TAN(slope_threshold/180*3.14159)
        ENDIF

        !IF(abs(FLUX).lt.0.01)THEN
        !    FLUX = 0.01
        !ENDIF

        kappa = 1.0-2.0/3.141592653589793*atan(FLUX*FLUX/(FLUX_MIN*FLUX_MIN))
        if(FLUX.le.0)then
          kappa = 1.0
        endif

        !THI = -9999.0
!
        !if(RF.lt.1.0e-20.or.RF.gt.1.0e-12)then
        !   RF=MRF
        !endif
        if(RF.lt.1.0e-20)then
          RF=1.0e-20
        endif



           !a1 = abs(FLUX)
           a1 = abs(((1.0-kappa)*abs(FLUX)+kappa*abs(FLUX_MIN)))
           a2 = (1.0/SLOPE**(nflow)*nflow+2.0)/(2.0*RF)*1.0/((rhoice*abs(grav))**nflow)

           !a1 = abs(FLUX)
           !a2 = (SLOPE**(nflow)/(nflow+2.0))*(2.0*RF)*((rhoice*abs(grav))**nflow)

           a3 = 1.0/(nflow+2.0)
           !b1 = (nflow+2.0)*a1**((-1.0)*(a3-1.0))
           b1 = (1.0/(nflow+2.0))*a1**((+1.0)*(a3-1.0))
           b2 = a2**a3

           !THI = b1*b2*dFLUX
           dFLUX = THIunc/b1/b2

!
         Values(Perm(i)) = dFLUX

    End do

    !deallocate(xobs,THIobs,InElement)
    !deallocate(ElementNodes % x,ElementNodes % y,ElementNodes % z)

    Firsttime=.false.

END SUBROUTINE computeSIA_deltaflux

