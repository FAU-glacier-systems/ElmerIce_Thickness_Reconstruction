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
SUBROUTINE CostSolver_FLUXpos( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!               Lambda= Real (default 1.0)
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,MaskVarSolName,SpeedVarSolName,HsiaVarSolName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VelocitybSol,HSol,MaskSol,SpeedSol,HsiaSol
  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  REAL(KIND=dp), POINTER :: Vb(:),HVal(:),MaskVal(:),SpeedVal(:),HsiaVal(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: VbPerm(:),HPerm(:),MaskPerm(:),SpeedPerm(:),HsiaPerm(:)
  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c
  real(kind=dp) :: Cost,Cost_surf,Cost_S,H,Lambda
  real(kind=dp) :: u,v,w,s,coeff,SqrtElementMetric,x
  real(kind=dp) :: speed_threshold,MASK,SPEED,HSIA,factor
  REAL(KIND=dp) :: NodeCost(Model % MaxElementNodes)
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: NodeCostb(Model % MaxElementNodes),NodeCost_der(3,Model %MaxElementNodes)
  REAL(KIND=dp) :: Hmin
  LOGICAL :: Reset
  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile
  save ElementNodes


     WRITE(SolverName, '(A)') 'CostSolver_FLUXpos'

  SolverParams => GetSolverParams()

!! Dimension of the pb; ie with SSA we can be 1D or 2D on a 2D mesh, or 2D on a 3D mesh
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

!  SolverParams => GetSolverParams()
!  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
!  If (.NOT.Found) then
!     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
!     DIM = CoordinateSystemDimension()
!  Endif

!-------------------------------------------------------------------------------
!
!                              MATERIAL CONSTANTS
!
!-------------------------------------------------------------------------------
!     ! speed threshold ... if below SIA thickness is taken
!     speed_threshold = GetConstReal( Model % Constants, 'speed threshold', Found )
!     If (.NOT.Found) Then
!            WRITE(Message,'(A)') 'Constant speed threshold not found. &
!                   &Setting to default (0.0)'
!            CALL INFO(SolverName, Message, level=20)
!            speed_threshold = 0.0
!     End if


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

    HSol => VariableGet( Solver % Mesh % Variables, 'Flux'  )
    IF ( ASSOCIATED( HSol ) ) THEN
            HVal => HSol % Values
            HPerm => HSol % Perm
    ELSE
            WRITE(Message,'(A)') &
                               'No variable > Flux < found'
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


!    MaskVarSolName =  GetString( SolverParams,'Mask Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Mask Variable Name< not found in section >Solver<')
!      END IF
!    MaskSol => VariableGet( Solver % Mesh % Variables, MaskVarSolName  )
!    IF ( ASSOCIATED( MaskSol ) ) THEN
!            MaskVal => MaskSol % Values
!            MaskPerm => MaskSol % Perm
!    ELSE
!            WRITE(Message,'(A)') &
!                               'No variable > Mask Variable Name < found'
!            CALL FATAL(SolverName,Message)
!    END IF

!    SpeedVarSolName =  GetString( SolverParams,'Speed Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Speed Variable Name< not found in section >Solver<')
!      END IF
!    SpeedSol => VariableGet( Solver % Mesh % Variables, SpeedVarSolName  )
!    IF ( ASSOCIATED( SpeedSol ) ) THEN
!            SpeedVal => SpeedSol % Values
!            SpeedPerm => SpeedSol % Perm
!    ELSE
!            WRITE(Message,'(A)') &
!                               'No variable > Speed Variable Name < found'
!            CALL FATAL(SolverName,Message)
!    END IF

!    HsiaVarSolName =  GetString( SolverParams,'Hsia Variable Name', Found)
!      IF(.NOT.Found) THEN
!              CALL WARN(SolverName,'Keyword >Hsia Variable Name< not found in section >Solver<')
!      END IF
!    HsiaSol => VariableGet( Solver % Mesh % Variables, HsiaVarSolName  )
!    IF ( ASSOCIATED( SpeedSol ) ) THEN
!            HsiaVal => HsiaSol % Values
!            HsiaPerm => HsiaSol % Perm
!    ELSE
!            WRITE(Message,'(A)') &
!                               'No variable > Hsia Variable Name < found'
!            CALL FATAL(SolverName,Message)
!    END IF

    IF (Reset)  Vb=0.0_dp


    Cost=0._dp
    Cost_surf=0._dp

    Hmin = 0.0_dp

    Do t=1,Solver % Mesh % NumberOfNodes
       H     = HVal(Hperm(t))
       !MASK  = MaskVal(MaskPerm(t))
       !SPEED = SpeedVal(SpeedPerm(t))
       !HSIA  = HsiaVal(HsiaPerm(t))
       IF (H.LT.0.0) THEN
         !IF (SPEED.LT.speed_threshold.OR.MASK.LT.1)THEN
         Cost=Cost+0.5*(H-Hmin)*(H-Hmin)*Lambda
         Vb(VbPerm(t))=Vb(VbPerm(t))+(H-Hmin)*Lambda
         !ENDIF
       ENDIF
    END DO
      



    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
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
END SUBROUTINE CostSolver_FLUXpos
!------------------------------------------------------------------------------
! *****************************************************************************
