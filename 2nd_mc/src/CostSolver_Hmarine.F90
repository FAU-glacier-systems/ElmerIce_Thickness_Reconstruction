!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: f. Gillet-Chaulet (LGGE, Grenoble,France)
! *  Email:   gillet-chaulet@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!
! *****************************************************************************
SUBROUTINE CostSolver_Hmarine( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!               Lambda= Real (default 1.0)
!       use "<Reset Cost Value> = False" to add cost and gradient to previously
!       computed values !!!
!
!       Required Sif parameters are:
!
!          In the solver section:
!               Problem Dimension=Integer (default:coordinate sytem dimension),
!               Cost Filename=File (default: CostOfIt.dat),
!               Cost Variable Name= String (default='Cost Value'),
!               Lambda= Real (default 1.0),
!               Reset Cost Value= Logical (default = True),
!
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfIt.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VelocitybSol,HSol,minBathySol,surfSol
  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce,BC
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: Vb(:),HVal(:),minBathyVal(:),surfVal(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: VbPerm(:),HPerm(:),minBathyPerm(:),surfPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c,ii,jj,m,kk
  integer :: bc_id,maskit(1:2)
  INTEGER, ALLOCATABLE :: NodeIndexMemory(:)
  real(kind=dp) :: Cost,Cost_S,Lambda
  real(kind=dp) :: THI(1:2),BATHY(1:2),SURF(1:2),HMIN(1:2),HMAX(1:2)
  real(kind=dp) :: SLVL,rhowater,rhoice
  real(kind=dp) :: u,v,w,s,coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: NodeCost(Model % MaxElementNodes)
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: NodeCostb(Model % MaxElementNodes),NodeCost_der(3,Model %MaxElementNodes)
  CHARACTER*10 :: date,temps

  LOGICAL :: Reset

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile
  save ElementNodes


     WRITE(SolverName, '(A)') 'CostSolver_Hmarine'

  SolverParams => GetSolverParams()

!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

!-------------------------------------------------------------------------------
!
!                              MATERIAL CONSTANTS
!
!-------------------------------------------------------------------------------
     rhowater = GetConstReal( Model % Constants, 'water density', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Water Density not found. &
                   &Setting to 1028 kg/m3'
            CALL INFO(SolverName, Message, level=20)
            rhowater = 1028.0/(1.0e6*(60.0*60.0*24.0*365.25)**2)
     End if  

     rhoice = GetConstReal( Model % Constants, 'ice density', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Ice Density not found. &
                   &Setting to 917 kg/m3'
            CALL INFO(SolverName, Message, level=20)
            rhoice = 917.0/(1.0e6*(60.0*60.0*24.0*365.25)**2)
     End if 


!! Optional weighting term
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=1.0')
           Lambda = 1.0
   End if
   !print *,'Lambda Hmarine',Lambda

!! Sea-level
   SLVL =  GetConstReal( Model % Constants,'sea level', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Sea Level< not found  in section >Constants<')
           CALL WARN(SolverName,'Taking default value SLVL=0.0')
           SLVL = 0.0
   End if

!! Do we need to reset cost and DJDVar to 0? Default YES
   Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
            IF(.NOT.Found) Reset=.True.

  !print *,'FIRST TIME'
  If (Firsttime) then
!-------------------------------------------------------------------------------
!
!                              MESH DETAILS / PARALLEL
!
!-------------------------------------------------------------------------------
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

!-------------------------------------------------------------------------------
!
!                              READ INPUT FIELDS
!
!-------------------------------------------------------------------------------

!!!!!!!!!!! get Solver Variables

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
    IF (.NOT. Found) CostFile = DefaultCostFile
    CALL DATE_AND_TIME(date,temps)
    If (Parallel) then
        if (ParEnv % MyPe.EQ.0) then
           OPEN (12, FILE=CostFile)
                    write(12,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(12)
         End if
    Else
           OPEN (12, FILE=CostFile)
                    write(12,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') '#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
           CLOSE(12)
    End if

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF
  
  !!! End of First visit
    Firsttime=.false.
  Endif

  !print *,'END FIRST TIME'
    HSol => VariableGet( Solver % Mesh % Variables, 'H'  )
    IF ( ASSOCIATED( HSol ) ) THEN
            HVal => HSol % Values
            HPerm => HSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > H < found'
            CALL FATAL(SolverName,Message)
    END IF  

    VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Hb'  )
    IF ( ASSOCIATED( VelocitybSol ) ) THEN
            Vb => VelocitybSol % Values
            VbPerm => VelocitybSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > Velocityb < found'
            CALL FATAL(SolverName,Message)
    END IF  

! Surface elevation
    surfSol => VariableGet( Solver % Mesh % Variables, 'surface'  )
    IF ( ASSOCIATED( surfSol ) ) THEN
            surfVal => surfSol % Values
            surfPerm => surfSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > surface < found'
            CALL FATAL(SolverName,Message)
    END IF

! Maximum ice thickness from minimum bathymetry in a certain radius
! (based on IBCAO which comes at a 500m resolution)
    minBathySol => VariableGet( Solver % Mesh % Variables, 'minel'  )
    IF ( ASSOCIATED( minBathySol ) ) THEN
            minBathyVal => minBathySol % Values
            minBathyPerm => minBathySol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > minel < found'
            CALL FATAL(SolverName,Message)
    END IF


!-------------------------------------------------------------------------------
!
!                              INITIALISE FIELDS
!
!-------------------------------------------------------------------------------

    !print *,'INITIALISE'
    allocate(NodeIndexMemory(Solver % Mesh % NumberOfBoundaryElements))

    IF (Reset)  Vb=0.0_dp


    Cost=0._dp


!-------------------------------------------------------------------------------
!
!                              BOUNDARY ELEMENT LOOP
!
!-------------------------------------------------------------------------------

     NodeIndexMemory = -1
     kk = 0

     !print *,'BEFORE LOOP', Solver % Mesh % NumberOfBoundaryElements, Solver % Matrix % ParMatrix % ParEnv % MyPE
     DO t=1, Solver % Mesh % NumberOfBoundaryElements
        if (t.eq.Solver % Mesh % NumberOfBoundaryElements) CYCLE
        ! get element information
        Element => GetBoundaryElement(t)
        !print *,'Element'
        ! cycle non-active elements
!        IF ( .NOT.ActiveBoundaryElement( Element ) ) print *,'NOT ACTIVE'
        IF ( .NOT.ActiveBoundaryElement( Element ) ) CYCLE
        !print *,'Acitve Element'
        ! cycle halo elements
!        IF (ParEnv % myPe .NE. Element % partIndex) print *,'HALO'
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
!        !print *,'Sort out halos'
        n = GetElementNOFNodes( Element )
!        !print *,'Number of Nodes',n
        NodeIndexes => Element % NodeIndexes
!        !print *,'Nodeindexes'
        m = kk
!        !print *,'BOUNDARY'
!        !print *,'BEFORE BOUNDARY CONDITION'
        BC => GetBC()
        bc_id = GetBCId( Element )
!        !print *,'BCBCBC',bc_id,n
!        !print *,t,n,m,bc_id
        IF(bc_id.EQ.2.0)THEN
!
        THI = HVal(Hperm(NodeIndexes(1:n)))
        SURF = surfVal(surfPerm(NodeIndexes(1:n)))
        BATHY = minBathyVal(minBathyPerm(NodeIndexes(1:n)))
!        if (Solver % Matrix % ParMatrix % ParEnv % MyPE == 23) then
!        print *,'BCTHI',THI(1),THI(2)
!        print *,'BCSUR',SURF(1),SURF(2)
!        print *, 'BCBATHY',BATHY(1),BATHY(2)
!        endif
!
        do ii=1,n
            maskit(ii) = 1
            if (SURF(ii).gt.0)then
            HMIN(ii) = rhowater/(rhowater-rhoice)*(SURF(ii)-SLVL)
            else
            HMIN(ii) = -1.0
            maskit(ii) = 0.0
            endif
 
            if (BATHY(ii).lt.0)then
            HMAX(ii) = rhowater/rhoice*(SLVL-BATHY(ii))
            else
            HMAX(ii) = -1.0
            maskit(ii) = 0.0
            endif
            if(kk.gt.0)then 
            do jj=1,m
               if(NodeIndexes(ii).eq.NodeIndexMemory(jj))then
                 maskit(ii) = 0.0
               endif
            enddo
            endif
        enddo
!        if (Solver % Matrix % ParMatrix % ParEnv % MyPE == 23) then
!        print *,'BCMASK',maskit(1),maskit(2)
!        print *,'BCHMIN',HMIN(1),HMIN(2)
!        print *,'BCHMAX',HMAX(1),HMAX(2)
!        endif
        do ii=1,n
           if(THI(ii).lt.Hmin(ii).and.maskit(ii).gt.0.0)then
           kk = kk+1
           NodeIndexMemory(kk)=NodeIndexes(ii)
!           !print *,'SMALLER HMIN',SURF(ii),HMIN(ii),THI(ii)
           Cost = Cost + maskit(ii)*0.5*Lambda*(THI(ii)-HMIN(ii))*(THI(ii)-HMIN(ii))
           Vb(VbPerm(NodeIndexes(ii)))=Vb(VbPerm(NodeIndexes(ii)))+maskit(ii)*(-THI(ii)+HMIN(ii))*Lambda
           else
!           !print *,'ELSE1'
           Cost = Cost
           Vb(VbPerm(NodeIndexes(ii)))=Vb(VbPerm(NodeIndexes(ii)))
           endif
           if(THI(ii).gt.Hmax(ii).and.maskit(ii).gt.0.0)then
           kk = kk+1
           NodeIndexMemory(kk)=NodeIndexes(ii)
!           !print *,'LARGER HMAX',BATHY(ii),HMAX(ii),THI(ii)
           Cost = Cost + maskit(ii)*0.5*Lambda*(THI(ii)-HMAX(ii))*(THI(ii)-HMAX(ii))
           Vb(VbPerm(NodeIndexes(ii)))=Vb(VbPerm(NodeIndexes(ii)))+maskit(ii)*(THI(ii)-HMAX(ii))*Lambda
           else
!           !print *,'ELSE2'
           Cost = Cost
           Vb(VbPerm(NodeIndexes(ii)))=Vb(VbPerm(NodeIndexes(ii)))
           endif
        enddo 
        ENDIF !bc_id
!      
!        if (Solver % Matrix % ParMatrix % ParEnv % MyPE == 23) then
!        print *,'LALA',t,n,bc_id,Solver % Matrix % ParMatrix % ParEnv % MyPE
!        endif
!
        !do ii=1,n
        !    NodeIndexMemory(kk+ii)=NodeIndexes(ii)
        !enddo
        !kk = kk+ii
!        if (Solver % Matrix % ParMatrix % ParEnv % MyPE == 23) then
!        print *,'LALA2',t,n,bc_id,Solver % Matrix % ParMatrix % ParEnv % MyPE
!        endif
     ENDDO
!     print *,'BE loop done',Solver % Matrix % ParMatrix % ParEnv % MyPE
     deallocate(NodeIndexMemory)

!-------------------------------------------------------------------------------
!
!                              WRITE OUTPUT FILES
!
!-------------------------------------------------------------------------------

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

!    print *,Parallel,TimeVar % Values(1)
   IF (Parallel) THEN
           !print *,'BEFORE MPI_ALLREDUCE'
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
           !print *,'AFTER MPI_ALLREDUCE'
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
               IF (Reset) then
                 CostVar % Values(1)=Cost_S
               Else
                 CostVar % Values(1)=CostVar % Values(1)+Cost_S
               Endif
               print *,'RESET OPTION CLEAR'
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (12, FILE=CostFile,POSITION='APPEND')
                 write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost_S
                 CLOSE(12)
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
                    OPEN (12, FILE=CostFile,POSITION='APPEND')
                       write(12,'(e13.5,2x,e15.8)') TimeVar % Values(1),Cost
                    close(12)
   END IF
   
   Return
!------------------------------------------------------------------------------
END SUBROUTINE CostSolver_Hmarine
!------------------------------------------------------------------------------
! *****************************************************************************
