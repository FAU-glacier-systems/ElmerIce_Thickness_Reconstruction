!
!  CALL to quadrant tree will be problematic if we work on a lower dimension
!  than the mesh.....
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE CostSolver_Hobs( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!
!     OUTPUT are : J and DJDvar
!                      
!
!    !!!!! BE carefull it will reset Cost and DJ to 0 by default !!!!
!      !!! If other cost and gradient are computed before (i.e. from the adjoint pb, 
!       use "<Reset Cost Value> = False" to add cost and gradient to previously
!       computed values !!!
!
!
!       Required Sif parameters are:
!
!          In the solver section:
!               Problem Dimension=Integer (default:coordinate sytem dimension),
!               Cost Filename=File (default: CostOfIT.dat),
!               Optimized Variable Name= String (default='Beta'),
!               Gradient Variable Name= String (default = 'DJDBeta'),
!               Cost Variable Name= String (default='Cost Value'),
!               Lambda= Real (default 1.0),
!               Reset Cost Value= Logical (default = True),
!
!          In Body Force section:
!               <VarSolName> a priori value = Real (defualt =0.0),
!               <VarSolName> variance = real (defualt=1.0)
!
!
!******************************************************************************
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfIT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile,UsedDataFile,ObsFileName
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: Variable,DJDVariable
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: ElementNodes
  REAL(KIND=dp), POINTER :: DJDValues(:),Values(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: DJDPerm(:),Perm(:)
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  integer :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
  real(kind=dp) :: Cost,Cost_S,scaler,Lambda
  real(kind=dp) :: xx,yy
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: THI,WEIGHT,aa,bb,ff,z0,z1,ww
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),THIobs(:), REFELobs(:), obsFLAG(:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex
  INTEGER,SAVE :: NTOT=-1,NTOT_S
  integer,SAVE :: nobs
  LOGICAL,SAVE :: FirstRound=.True.,SAVE_USED_DATA=.False.
  LOGICAL :: Reset
  CHARACTER*10 :: date,temps

  INTEGER,PARAMETER :: IO=12

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile,VarSolName,GradSolName
  save ElementNodes


  SolverParams => GetSolverParams()

!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

! get some needed solver parameters
!! Cost File for Output
  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile

!! Name of the variable to regularise
  VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >Beta<')
              WRITE(VarSolName,'(A)') 'Beta'
      END IF

!! Name of the variable to regularise
  GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
      IF(.NOT.Found) THEN
              CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found in section >Solver<')
              CALL WARN(SolverName,'Taking default value >DJDBeta<')
              WRITE(GradSolName,'(A)') 'DJDBeta'
      END IF


!! Name of the variable with the cost function
   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF

!! Optional weighting term
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=1.0')
           Lambda = 1.0
   End if

!! Do we need to reset cost and DJDVar to 0? Default YES
   Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
            IF(.NOT.Found) Reset=.True.

!!! SOME INITIALISATION AT FIRST TIME
  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'CostSolver_Hobs'


!!!!!!!!!!!  initiate Cost File
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (IO, FILE=CostFile)
                    write(IO,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(IO)
         End if
    Else
           OPEN (IO, FILE=CostFile)
                    write(IO,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(IO)
    End if

 !! Get the obs
   ObsFileName =  GetString( SolverParams,'Observation File Name', Found)
   IF(.NOT.Found) THEN
       CALL FATAL(SolverName,'Keyword >Observation File Name< not found in section >Solver<')
   END IF
   ParallelFile = .False.
   ParallelFile = GetLogical(SolverParams,'Parallel Observation Files', Found)
   if (Parallel.AND.ParallelFile) &
    write(ObsFileName,'(A,A,I0)') trim(ObsFileName),'.',ParEnv % MyPE  
               

   SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)

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


   allocate(xobs(nobs,3),THIobs(nobs),InElement(nobs),REFELobs(nobs),obsFLAG(nobs))
   InElement(:)=-1

   THIobs=0.0_dp
   xobs=0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),(THIobs(i)),(REFELobs(i)),(obsFLAG(i))
   end do
   close(IO)

  !!! End of First visit
    Firsttime=.false.
  Endif
!!!! INITIALISATION DONE

    Variable => VariableGet( Solver % Mesh % Variables, VarSolName  )
    IF ( ASSOCIATED( Variable ) ) THEN
            Values => Variable % Values
            Perm => Variable % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF
    DJDVariable => VariableGet( Solver % Mesh % Variables, GradSolName  )
    IF ( ASSOCIATED( DJDVariable ) ) THEN
            DJDValues => DJDVariable % Values
            DJDPerm => DJDVariable % Perm
!            Vb => VelocitybSol % Values
!            VbPerm => VelocitybSol % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF
    IF (Reset) DJDValues=0.0_dp


    IF (DJDVariable%DOFs.NE.Variable%DOFs) & 
        CALL FATAL(Trim(SolverName),'DOFs do not correspond')



!!!! Start looping
    Cost=0.0_dp




    CALL StartAdvanceOutput(SolverName,'Compute cost')
    Do s=1,nobs

     CALL AdvanceOutput(s,nobs)
     
     IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(s,1:DIM)
      CALL LocateParticleInMeshOctree( ElementIndex,Coord)
      If (ElementIndex.NE.0) InElement(s)=ElementIndex
     ENDIF !End if FirstRound
      

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


            xx = xobs(s,1)
            yy = xobs(s,2)

!            ! Variable at obs point
!            Do i=1,VDOFs
!              Vobs_mag = Vobs_mag + Vobs(s,i)*Vobs(s,i)
!              V(i)=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
!            End do
!            Vobs_mag = sqrt(Vobs_mag)

!            ! Update cost
!            scaler=Vobs_mag
!            scaler=1.0

            WEIGHT = 1.0
            ww     = 0.01
            z0     = 0.0
            z1     = 1200.0
            ff     = 1.0/1000.0
            IF (obsFLAG(s).eq.0) THEN
              IF(REFELobs(s).lt.z1)THEN
                 !aa = (1.0-1.0/100.0)/(1200.0**2.0-0.0**2.0)
                 !bb = 1-aa*1200.0**2.0
                 !WEIGHT = aa*REFELobs(s)**2.0+bb
                 !print *,'REFEL',REFELobs(s),WEIGHT
                 z0=0.0
                 z1=1200.0
                 ff=1.0/1000.0
                 bb=LOG(ff)/(z0-z1)
                 aa=exp(-bb*z1)
                 WEIGHT=ww*aa*exp(bb*REFELobs(s))
              ELSE
                 WEIGHT=ww
              ENDIF
              !WEIGHT=1.0
              !bb=LOG(ff)/(z0-z1)
              !aa=exp(-bb*z1)
              !WEIGHT=ww*aa
            ENDIF
            !print *,'REFEL',REFELobs(s),WEIGHT,obsFLAG(s)
 
            THI =  SUM(Values(Perm(NodeIndexes(1:n)))*basis(1:n))

            Cost=Cost+0.5*WEIGHT*(THI-THIobs(s))*(THI-THIobs(s))*Lambda

            !PRINT *,V,Vobs

            !Derivative of Cost at nodal Point
            Do j=1,n
               k=DJDPerm(NodeIndexes(j))
               DJDValues(k)=DJDValues(k)+WEIGHT*(THI-THIobs(s))*basis(j)*Lambda
            End Do

          END IF

         ELSE

            WRITE(Message,'(a,I0,a)')&
                'Data Point',s,'found in no element'
            CALL Info( SolverName, Message,level=15)
         END IF

    END DO !Do on s


    if (NTOT < 0) NTOT=COUNT(InElement(:)>0)
!! Save Cost as a function of time

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(NTOT,NTOT_S,1,&
                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr) 
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
               IF (Reset) then
                 CostVar % Values(1)=Cost_S
               Else
                 CostVar % Values(1)=CostVar % Values(1)+Cost_S
               Endif
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (IO, FILE=CostFile,POSITION='APPEND')
                 write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,sqrt(1.0*Cost_S/NTOT_S)
                 CLOSE(IO)
                 write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT_S
                 call INFO(SolverName,Message,level=3)
         End if
   ELSE
            CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                 IF (Reset) then
                    CostVar % Values(1)=Cost
                 Else
                    CostVar % Values(1)=CostVar % Values(1)+Cost
                 Endif
            END IF
            OPEN (IO, FILE=CostFile,POSITION='APPEND')
            write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,sqrt(1.0*Cost/NTOT)
            close(IO)
            write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT
            call INFO(SolverName,Message,level=3)
   END IF
   
   Return
!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Hobs
!------------------------------------------------------------------------------
! *****************************************************************************
