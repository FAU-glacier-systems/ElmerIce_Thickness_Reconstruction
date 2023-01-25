!
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE DJDp_Adjoint_IceFluxSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE Differentials
  USE MaterialModels
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    external variables
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------

  LOGICAL ::&
       firstTime=.TRUE., Found, AllocationsDone = .FALSE., stat, &
       LimitDisp,  Bubbles = .False.,&
       SubstantialSurface = .TRUE.,&
       UseBodyForce = .TRUE., ApplyDirichlet=.FALSE.,  ALEFormulation=.FALSE. , &
       ConvectionVar,Compute_dhdt
  LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:), ActiveNode(:,:)

  INTEGER :: & 
       i,j,K,L, p, q, R, t,N,NMAX,MMAX,nfamily, deg, Nmatrix,&
       edge, bf_id,DIM,istat,LocalNodes,nocorr,&
       NSDOFs,NonlinearIter,iter, numberofsurfacenodes
  INTEGER, POINTER ::&
       ThickPerm(:), DHDTPrem(:),FlowPerm(:), NodeIndexes(:), EdgeMap(:,:), SlopePerm(:)

!adjoint
  INTEGER, POINTER :: ThickbPerm(:)

  REAL(KIND=dp) :: &
       at,st,totat,totst,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh,LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss

  REAL(KIND=dp), POINTER :: ForceVector(:), Thick(:),DHDT(:),PreH(:,:), &
       FlowSolution(:),  PointerToResidualVector(:), SlopeSolution(:)

!!adjoint
  REAL(KIND=dp), POINTER :: Thickb(:)

  REAL(KIND=dp), ALLOCATABLE :: ResidualVector(:), &
       STIFF(:,:),SourceFunc(:),FORCE(:), TimeForce(:), LOAD(:),&
       MASS(:,:), Velo(:,:), Flux(:,:), LowerLimit(:), UpperLimit(:), &
       OldValues(:), OldRHS(:),StiffVector(:),MeshVelocity(:,:)
!!adjoint
  REAL(KIND=dp), ALLOCATABLE ,SAVE:: LOADb(:),Velob(:,:),Dirichlet(:),DirichletC(:)

  LOGICAL :: Reset, ComputeDJDsmbTop,ComputeDJDsmbBot,ComputeDJDUV
  TYPE(Variable_t), POINTER :: DJDsmbTopSol,DJDsmbBotSol,DJDUVSol
  INTEGER,POINTER :: DJDsmbTopPerm(:),DJDsmbBotPerm(:),DJDUVPerm(:)
  REAL(KIND=dp), POINTER :: DJDsmbTop(:),DJDsmbBot(:),DJDUV(:)

  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName,  EquationName, FlowSolName, StabilizeFlag, SlopeSolName
  CHARACTER(LEN=MAX_NAME_LEN) :: ThickSolName,ThickbSolName

  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: FlowSol, ThickSol,ThickbSol, SlopeSol
  TYPE(ValueList_t), POINTER :: BodyForce, SolverParams, Material, Equation,BC
  TYPE(Matrix_t), POINTER :: Systemmatrix
  !-----------------------------------------------------------------------------
  !      remember these variables
  !----------------------------------------------------------------------------- 
  SAVE STIFF, MASS, FORCE, &
       LOAD, &
       ElementNodes, AllocationsDone, Velo,  TimeForce, &
       UseBodyForce, LimitedSolution, LowerLimit, UpperLimit, ActiveNode, OldValues, OldRHS, &
       ResidualVector, StiffVector, MeshVelocity

  !------------------------------------------------------------------------------
  !    Get variabel/solver name
  !------------------------------------------------------------------------------
  SolverName = 'DJDp_Adjoint_ThicknessSolver'

  !------------------------------------------------------------------------------
  !    Get constants and solver params
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  smallestpossiblenumber = TINY(smallestpossiblenumber)
  SolverParams => GetSolverParams()


  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
     CALL Info(SolverName, 'No keyword > Apply Dirichlet < found. No limitation of solution',Level=6 )
  ELSE
     IF (ApplyDirichlet) THEN
        CALL Info(SolverName, 'Using Dirichlet method for limitation',Level=6 )
     ELSE
        CALL Info(SolverName, 'No limitation of solution',Level=6 )
     END IF
  END IF

  ALEFormulation = GetLogical( SolverParams, &
       'ALE Formulation', Found)
  IF ( .NOT.Found ) THEN
     ALEFormulation = .FALSE.
  END IF
  IF (ALEFormulation) THEN 
     CALL Info(SolverName, 'Using horizontal ALE Formulation',Level=6 )
  ELSE
     CALL Info(SolverName, 'Using horizontal Eulerian Formulation',Level=6 )
  END IF

  StabilizeFlag = GetString( SolverParams, &
       'Stabilization Method',Found )
  SELECT CASE(StabilizeFlag)
     CASE('stabilized')
        Bubbles = .FALSE.
     CASE('bubbles')
        Bubbles = .TRUE.
     CASE DEFAULT
        Bubbles = .FALSE.
  END SELECT
  IF (Bubbles) THEN
     CALL Info(SolverName, 'Using residual free bubble stabilization',Level=6 )
  ELSE
     CALL Info(SolverName, &
          'Using residual squared-stabilized formulation.',Level=6 )
  END IF


  WRITE(Message,'(A,I0)') 'Mesh dimension: ', DIM
  CALL Info( SolverName, Message, Level=8 )

  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------

  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     NMAX = Model % MaxElementNodes
     MMAX = Model % Mesh % NumberOfNodes 

     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,    &
             ElementNodes % y,    &
             ElementNodes % z,    &
             TimeForce,        &
             FORCE,    &
             STIFF, &
             MASS,  &
             LOAD,&
             LOADb, &
             Velo,  &
             Velob,&
             MeshVelocity, &
             Flux, &
             SourceFunc, &
             LowerLimit,                      &
             UpperLimit, &
             LimitedSolution,  &
             ActiveNode,                      & 
             ResidualVector, &
             OldValues,&
             OldRHS,&
             StiffVector)

     END IF


     IF (Bubbles) THEN
        Nmatrix = 2*NMAX
     ELSE
        Nmatrix = NMAX
     END IF

     ALLOCATE( ElementNodes % x( NMAX ),    &
          ElementNodes % y( NMAX ),    &
          ElementNodes % z( NMAX ),    &
          LOAD(NMAX) , LOADb(NMAX),&
          Dirichlet(NMAX),DirichletC(NMAX),&
          Velo( 3, NMAX ), Velob( 3, NMAX ),&
          MeshVelocity( 3,NMAX ), &
          Flux( 3, NMAX), &
          SourceFunc( NMAX ), &
          LowerLimit( MMAX ), &
          UpperLimit( MMAX ), &
          LimitedSolution( MMAX, 2 ),  &
          ActiveNode( MMAX, 2 ),       &  
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error 1, Aborting.')
     END IF

     CALL Info(SolverName,'Memory allocations done' )
     AllocationsDone = .TRUE.
  END IF


  !------------------------------------------------------------------------------
  !    Get variables 
  !------------------------------------------------------------------------------
   
   ThickSolName =  GetString(SolverParams ,'Thickness Solution Name', Found)
   IF(.NOT.Found) THEN
     CALL WARN(SolverName,'<Thickness Solution Name> not Found; assume default H')
     ThickSolName = 'H'
   ENDIF
   ThickSol => VariableGet( Solver % Mesh % Variables, ThickSolName)
   IF (ASSOCIATED(ThickSol)) THEN
     Thick => ThickSol % Values
     ThickPerm => ThickSol % Perm
   ELSE
     WRITE(Message,'(A,A,A)') &
     'No variable >',ThickSolName,'< found'
     CALL FATAL(SolverName,Message)
   ENDIF

   ThickbSolName = GetString(SolverParams ,'Adjoint Solution Name', Found)
   IF(.NOT.Found) THEN
     CALL WARN(SolverName,'<Adjoint Solution Name> not Found; assume default Adjoint')
     ThickSolName = 'Adjoint'
   ENDIF
   ThickbSol => VariableGet( Solver % Mesh % Variables, ThickbSolName)
   IF (ASSOCIATED(ThickbSol)) THEN
     Thickb => ThickbSol % Values
     ThickbPerm => ThickbSol % Perm
   ELSE
     WRITE(Message,'(A,A,A)') &
     'No variable >',ThickbSolName,'< found'
     CALL FATAL(SolverName,Message)
   ENDIF

   ComputeDJDsmbTop = GetLogical(SolverParams ,'ComputeDJDsmbTop' , Found)
   IF (ComputeDJDsmbTop) THEN
     DJDsmbTopSol => VariableGet( Solver % Mesh % Variables, 'DJDsmbTop_IceFlux' )
     IF (ASSOCIATED(DJDsmbTopSol)) THEN
       DJDsmbTop => DJDsmbTopSol % Values
       DJDsmbTopPerm => DJDsmbTopSol % Perm
     ELSE
       CALL FATAL(SolverName,'No variable <DJDsmbTop_IceFlux> Found')
     ENDIF
     Reset =  GetLogical( SolverParams,'Reset DJDsmbTop', Found)
     if (Reset.OR.(.NOT.Found)) DJDsmbTop = 0.0
   ENDIF

   ComputeDJDsmbBot = GetLogical(SolverParams ,'ComputeDJDsmbBot' , Found)
   IF (ComputeDJDsmbBot) THEN
     DJDsmbBotSol => VariableGet( Solver % Mesh % Variables, 'DJDsmbBot' )
     IF (ASSOCIATED(DJDsmbBotSol)) THEN
       DJDsmbBot => DJDsmbBotSol % Values
       DJDsmbBotPerm => DJDsmbBotSol % Perm
     ELSE
       CALL FATAL(SolverName,'No variable <DJDsmBot> Found')
     ENDIF
     Reset =  GetLogical( SolverParams,'Reset DJDsmbBot', Found)
     if (Reset.OR.(.NOT.Found)) DJDsmbBot = 0.0
   ENDIF


  !------------------------------------------------------------------------------
  !    Get Flow solution
  !------------------------------------------------------------------------------
   ConvectionVar=.True.
   FlowSolName =  GetString(SolverParams ,'Flow Solution Name', Found)
   SlopeSolName =  GetString(SolverParams ,'Slope Variable Name', Found)
   IF(.NOT.Found) THEN        
        WRITE(Message,'(A)') &
           '<Flow Solution Name> Not Found; will look for <convection velocity> in body forces'
        CALL Info(SolverName,Message,level=10)
        ConvectionVar=.False.
        NSDOFS=GetInteger(SolverParams ,'Convection Dimension',Found)
        IF(.NOT.Found) &
            CALL Fatal(SolverName,'if <Flow Solution Name> not given prescribe <Convection Dimension>')
   ELSE
        FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
        IF ( ASSOCIATED( FlowSol ) ) THEN
             FlowPerm     => FlowSol % Perm
             NSDOFs     =  FlowSol % DOFs
             FlowSolution => FlowSol % Values
        ELSE
             WRITE(Message,'(A,A,A)') &
                    'No variable >',FlowSolName,'< found'
             CALL Fatal(SolverName,Message)              
        END IF
        SlopeSol => VariableGet( Solver % Mesh % Variables, SlopeSolName )
        IF ( ASSOCIATED( SlopeSol ) ) THEN
             SlopePerm     => SlopeSol % Perm
             NSDOFs     =  SlopeSol % DOFs
             SlopeSolution => SlopeSol % Values
        ELSE
             WRITE(Message,'(A,A,A)') &
                    'No variable >',SlopeSolName,'< found'
             CALL Fatal(SolverName,Message)
        END IF
   END IF

   ComputeDJDUV = GetLogical(SolverParams ,'ComputeDJDUV' , Found)
   IF (ComputeDJDUV) THEN
     DJDUVSol => VariableGet( Solver % Mesh % Variables, 'DJDuv_IceFlux' )
     IF (ASSOCIATED(DJDUVSol)) THEN
       DJDUV => DJDUVSol % Values
       DJDUVPerm => DJDUVSol % Perm
       IF (DJDUVSol % DOFs.NE.NSDOFs) &
          CALL FATAL(SolverName,'DJDUV DOFs is different from the velocity DOFs')
     ELSE
       CALL FATAL(SolverName,'No variable <DJDuv_IceFlux> Found')
     ENDIF
     Reset =  GetLogical( SolverParams,'Reset DJDUV', Found)
     if (Reset.OR.(.NOT.Found)) DJDUV = 0.0
   ENDIF



  DO t=1,GetNOFBoundaryElements()
     CurrentElement => GetBoundaryElement(t)
     BC => GetBC()
     IF (.NOT.ASSOCIATED(BC)) CYCLE
     n = GetElementNOFNodes()
     NodeIndexes => CurrentElement % NodeIndexes
     Dirichlet(1:n)=ListGetReal(BC,'H',n, NodeIndexes, Found)
     IF (Found) Then
       DirichletC(1:n)=ListGetReal(BC,'H Condition',n, NodeIndexes, Found)
       Do i=1,n
          IF (Found.AND.(DirichletC(i).LT.0._dp)) then
              CYCLE
          ELSE
                !PRINT *,Thickb(ThickbPerm(NodeIndexes(i)))
              Thickb(ThickbPerm(NodeIndexes(i)))=0._dp
          ENDIF
        END DO
       END IF
    END DO



     !------------------------------------------------------------------------------
     !    Do the assembly
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes

        ! set coords of highest occuring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (NSDOFs == 1) THEN
           ElementNodes % y(1:n) = 0.0
           ElementNodes % z(1:n) = 0.0
        ELSE IF (NSDOFs == 2) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i0,a)')&
                'It is not possible to compute Thickness evolution if Flow Sol DOFs=',&
                NSDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message) 
           STOP   
        END IF


        ! get pointers on Equation, Material and body-Force section input
        !----------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()
        BodyForce => GetBodyForce()

        Dirichlet(1:n)=ListGetReal(BodyForce,'H',n, NodeIndexes, Found)
        IF (Found) Then
          DirichletC(1:n)=ListGetReal(BodyForce,'H Condition',n, NodeIndexes, Found)
          Do i=1,n
             IF (Found.AND.(DirichletC(i).LT.0._dp)) then
                CYCLE
             ELSE
                !PRINT *,Thickb(ThickbPerm(NodeIndexes(i)))
                Thickb(ThickbPerm(NodeIndexes(i)))=0._dp
             ENDIF
          END DO
       END IF


        ! get lower limit for solution 
        !-----------------------------
        LowerLimit( Nodeindexes(1:N)) = &
             ListGetReal(Material,'Min ' // TRIM(ThickSolName),n, NodeIndexes, Found) 
        LimitedSolution( Nodeindexes(1:N), 1) = Found
        ! get upper limit for solution 
        !-----------------------------
        UpperLimit(Nodeindexes(1:N)) = &
             ListGetReal(Material,'Max ' // TRIM(ThickSolName),n, NodeIndexes, Found)              
        LimitedSolution(Nodeindexes(1:N), 2) = Found

  ! Adjoint
  !! IF solutions has been limited then derivative is 0
       IF (ApplyDirichlet) THEN
         Do i=1,N
              l= ThickPerm(Nodeindexes(i)) 
              IF ((LimitedSolution(Nodeindexes(i),1)) &
                 .AND.(Thick(l)-LowerLimit(Nodeindexes(i)).LE.0.0_dp+smallestpossiblenumber )) THEN
                Thickb(ThickbPerm(Nodeindexes(i)))=0._dp
              END IF
              IF ((LimitedSolution(Nodeindexes(i),2)) &
                 .AND.(Thick(l)-UpperLimit(Nodeindexes(i)).GE.0.0_dp-smallestpossiblenumber )) THEN
                Thickb(ThickbPerm( Nodeindexes(i)))=0._dp
              END IF
          END DO
       END IF

        ! get flow soulution and velocity field from it
        !----------------------------------------------
        Velo = 0.0_dp
        !----------------------------------------------------

        ! get velocity profile
        IF (ConvectionVar) Then
          DO i=1,n
             j = NSDOFs*SlopePerm(NodeIndexes(i))
              !2D problem - 1D Thickness evolution
              IF((DIM == 2) .AND. (NSDOFs == 1)) THEN
                 Velo(1,i) = (-1.0)*SlopeSolution( j )
                 Velo(2,i) = 0.0_dp
              !2D problem - 2D Thickness evolution (plane view pb)
              ELSE IF ((DIM == 2) .AND. (NSDOFs == 2)) THEN
                 Velo(1,i) = (-1.0)*SlopeSolution( j-1 )
                 Velo(2,i) = (-1.0)*SlopeSolution( j )
              !3D problem - 2D Thickness evolution 
              ELSE IF ((DIM == 3) .AND. (NSDOFs == 2)) THEN
                 Velo(1,i) = (-1.0)*SlopeSolution( j-1 )
                 Velo(2,i) = (-1.0)*SlopeSolution( j )
              ELSE
!             j = NSDOFs*FlowPerm(NodeIndexes(i))
!              !2D problem - 1D Thickness evolution
!              IF((DIM == 2) .AND. (NSDOFs == 1)) THEN 
!                 Velo(1,i) = FlowSolution( j ) 
!                 Velo(2,i) = 0.0_dp
!              !2D problem - 2D Thickness evolution (plane view pb)
!              ELSE IF ((DIM == 2) .AND. (NSDOFs == 2)) THEN
!                 Velo(1,i) = FlowSolution( j-1 ) 
!                 Velo(2,i) = FlowSolution( j ) 
!              !3D problem - 2D Thickness evolution 
!              ELSE IF ((DIM == 3) .AND. (NSDOFs == 2)) THEN
!                 Velo(1,i) = FlowSolution( j-1 ) 
!                 Velo(2,i) = FlowSolution( j ) 
!              ELSE
                 WRITE(Message,'(a,i0,a,i0,a)')&
                      'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                 CALL Fatal( SolverName, Message)
              END IF
           END DO
       ELSE
          IF (ASSOCIATED( BodyForce ) ) THEN
               Velo(1,1:n) = GetReal( BodyForce, 'Convection Velocity 1',Found )
               if (NSDOFs.eq.2) Velo(2,1:n) = GetReal( BodyForce, 'Convection Velocity 2',Found )
          END IF
        END IF
        !------------------------------------------------------------------------------
        ! Get mesh velocity
        !------------------------------------------------------------------------------
        MeshVelocity = 0.0_dp
        CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity',CurrentElement)
        !

        !------------------------------------------------------------------------------
        !      get the accumulation/ablation rate (i.e. normal surface flux)
        !      from the body force section
        !------------------------------------------------------------------------------
        LOAD=0.0_dp
        IF (ASSOCIATED( BodyForce ) ) THEN
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Top Surface Accumulation', Found )
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Bottom Surface Accumulation', Found )
        END IF


        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix( STIFF, MASS, FORCE,&
             LOAD, LOADb,  Velo, Velob, NSDOFs, MeshVelocity, &
             CurrentElement, n, ElementNodes, NodeIndexes, &
             TransientSimulation,&
              ALEFormulation)

        If (ComputeDJDsmbTop) &
           DJDsmbTop(DJDsmbTopPerm(NodeIndexes(1:n)))=DJDsmbTop(DJDsmbTopPerm(NodeIndexes(1:n)))+LOADb(1:n)
        If (ComputeDJDsmbBot) &
           DJDsmbBot(DJDsmbBotPerm(NodeIndexes(1:n)))=DJDsmbBot(DJDsmbBotPerm(NodeIndexes(1:n)))+LOADb(1:n)
        If (ComputeDJDUV) then

        Do j=1,n
          Do i=1,NSDOFs
             DJDUV(NSDOFs*(DJDUVPerm(NodeIndexes(j))-1)+i)=DJDUV(NSDOFs*(DJDUVPerm(NodeIndexes(j))-1)+i)+&
               Velob(i,j)
           End do
          End Do
        End if

     END DO ! End loop bulk elements

     !------------------------------------------------------------------------------
   CONTAINS

     !------------------------------------------------------------------------------
     !==============================================================================
     SUBROUTINE LocalMatrix( STIFF, MASS, FORCE,&
          LOAD,LOADb,  Velo, Velob,  NSDOFs, MeshVelo, &
          Element, nCoord, Nodes, NodeIndexes, &
          TransientSimulation,&
          ALEFormulation)
       !------------------------------------------------------------------------------
       !------------------------------------------------------------------------------
       !      external variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            STIFF(:,:), MASS(:,:), FORCE(:), LOAD(:), LOADb(:), &
            Velo(:,:), Velob(:,:),MeshVelo(:,:)

       INTEGER :: nCoord, NodeIndexes(:), NSDOFs
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: TransientSimulation,ALEFormulation

       !------------------------------------------------------------------------------
       !      internal variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            Basis(2*nCoord),dBasisdx(2*nCoord,3), &
            Vgauss(3),  Source, &
            X,Y,Z,U,V,W,S,SqrtElementMetric, SU(2*nCoord),SW(2*nCoord),Tau,hK,UNorm,divu


       REAL(KIND=dp) :: Sourceb,divub,Vgaussb(3),SUb(2*nCoord),SWb(2*nCoord),Taub,Unormb

       TYPE(ElementType_t), POINTER :: SaveElementType
       INTEGER :: LinType(2:4) = [202,303,404]

       LOGICAL :: Stat, UseLinear
       INTEGER :: i,j,t,p,q, n
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       !------------------------------------------------------------------------------


       LOADb = 0.0_dp
       Velob = 0.0_dp

       IF (Bubbles) THEN
          n = nCoord * 2
       ELSE
          n = nCoord
       END IF


       hK = ElementDiameter( Element, Nodes )

       !
       !      Numerical integration:
       !      ----------------------
       IF (Bubbles) THEN
          IntegStuff = GaussPoints( Element, Element % TYPE % gausspoints2)
       ELSE
          IntegStuff = GaussPoints( Element )
       END IF

       SU = 0.0_dp
       SW = 0.0_dp

       DO t = 1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
          !
          !        Basis function values & derivatives at the integration point:
          !        -------------------------------------------------------------
          stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
               Basis,dBasisdx, Bubbles=Bubbles )

          !        Correction from metric
          !        ----------------------
          S = S * SqrtElementMetric

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             X = SUM( Nodes % x(1:nCoord) * Basis(1:nCoord) )
             Y = SUM( Nodes % y(1:nCoord) * Basis(1:nCoord) )
             Z = SUM( Nodes % z(1:nCoord) * Basis(1:nCoord) )
             S = S * X
          END IF
          !
          !        Velocities and (norm of) gradient of free surface and source function 
          !        at Gauss point
          !        ---------------------------------------------------------------------

          Vgauss=0.0_dp

          IF (.NOT.ALEFormulation) THEN
             DO i=1,NSDOFs
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord)))
                norm = norm + Vgauss(i)**2.0
             END DO
          ELSE
             DO i=1,NSDOFs
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord) - MeshVelo(i,1:nCoord)))
                norm = norm + Vgauss(i)**2.0
             END DO
          END IF
          norm = norm**0.5

          IF (norm.gt.0.0) then
            DO i=1,NSDOFs
              Vgauss(i) = Vgauss(i)/norm
            ENDDO
          ENDIF

          divu = 0.0_dp
          DO i=1,NSDOFs
             divu = divu +  SUM( dBasisdx(1:nCoord,i)*(Velo(i,1:nCoord)))
          END DO

          UNorm = SQRT( SUM( Vgauss(1:NSDOFs)**2 ) )
          IF (UNorm .NE. 0.0_dp) THEN
             Tau = hK / ( 2*Unorm )
          ELSE
             Tau = 0.0_dp
          END IF

          IF ( .NOT. Bubbles ) THEN
             DO p=1,n
                SU(p) = 0.0_dp
                DO i=1,NSDOFs
                   SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
                END DO

                SW(p) = 0.0_dp
                DO i=1,NSDOFs
                   SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
                END DO
             END DO
          END IF

          !        Stiffness matrix:
          !        -----------------
          !DO p=1,n
          !   DO q=1,n
          !      DO i=1,NSDOFs
          !         STIFF(p,q) = STIFF(p,q) + &
          !              s * Vgauss(i) * dBasisdx(q,i) * Basis(p)
          !      END DO
          !      STIFF(p,q) =  STIFF(p,q) + s * Tau * SU(q) * SW(p)
          !      STIFF(p,q) =  STIFF(p,q) + s * divu * Basis(q) * (Basis(p) + Tau*SW(p))
          !   END DO
          !END DO

          !! Adjoint part of the stiffness matrix
          divub = 0._dp
          Vgaussb = 0._dp
          SUb = 0._dp
          SWb = 0._dp
          Taub =0._dp
          Do p=1,n
            Do q=1,n
               SUb(q) = SUb(q) + s * Tau * SW(p) * &
                      (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               SWb(p) = SWb(p) + s * Tau * SU(q) * &
                      (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               SWb(p) = SWb(p) + s * divu * Basis(q) * Tau * &
                      (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
            !!!! semble moins bon si j'ajoute dependance de Tau à V??? erreur
               Taub = Taub + s *  SU(q) * SW(p) * &
                      (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               Taub = Taub + s * divu * Basis(q) * SW(p) * &
                     (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               divub = divub + s  * Basis(q) * (Basis(p) + Tau*SW(p)) * &
                      (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               Do i=1,n
                 Vgaussb(i)=Vgaussb(i) + s  * dBasisdx(q,i) * Basis(p) * &
                     (-Thick(ThickPerm(NodeIndexes(q)))*Thickb(ThickbPerm(NodeIndexes(p))))
               End do
            End do
          End do

          !        Get accumulation/ablation function 
          !        --------------------------------------------------------- 
          Source = 0.0_dp
          Source=SUM(Basis(1:nCoord)*LOAD(1:nCoord))

          !        Assemble force vector:
          !        ---------------------
          !FORCE(1:n) = FORCE(1:n) &
          !     + Source * (Basis(1:n) + Tau*SW(1:n)) * s
          
          !! Adjoint Part of the accumulation/ablation function
          Sourceb=0._dp
          Do p=1,n
            Sourceb = Sourceb + (Basis(p) + Tau*SW(p)) * s * Thickb(ThickbPerm(NodeIndexes(p)))
            SWb(p) = SWb(p) + Source * Tau * s * Thickb(ThickbPerm(NodeIndexes(p)))
          End Do
          LOADb(1:n) = LOADb(1:n) + Sourceb * Basis(1:n)
          !!

          !IF ( .NOT. Bubbles ) THEN
          !   DO p=1,n
          !      SU(p) = 0.0_dp
          !      DO i=1,NSDOFs
          !         SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
          !      END DO

          !      SW(p) = 0.0_dp
          !      DO i=1,NSDOFs
          !         SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
          !      END DO
          !   END DO
          !ENDIF

          IF ( .NOT. Bubbles ) THEN
            DO p=1,n
              DO i=1,NSDOFs
                Vgaussb(i)=Vgaussb(i)+SUb(p)*dBasisdx(p,i)
                Vgaussb(i)=Vgaussb(i)+SWb(p)*dBasisdx(p,i)
              END DO
            END DO
          ENDIF

        !!!! semble moins bon si j'ajoute dependance de Tau à V??? erreur
          !UNorm = SQRT( SUM( Vgauss(1:NSDOFs)**2 ) )
          !IF (UNorm .NE. 0.0_dp) THEN
          !   Tau = hK / ( 2*Unorm )
          !ELSE
          !   Tau = 0.0_dp
          !END IF

          IF (UNorm .NE. 0.0_dp) THEN
             Unormb = -0.5*hK*(Unorm**(-2)) * Taub
          ELSE
             Unormb = 0._dp
          END IF

          !Do i=1,NSDOFs
          !  IF (UNorm .NE. 0.0_dp) &
          !     Vgaussb(i) = Vgaussb(i) + Unormb * 0.5*(SUM( Vgauss(1:NSDOFs)**2)**(-0.5))*2.0*Vgauss(i)
          !END DO

          Do p=1,n
            Do i=1,NSDOFs
              Velob(i,p)=Velob(i,p)+Vgaussb(i)*Basis(p)+divub*dBasisdx(p,i)
            End do
          End do


       END DO

       !------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix

     !------------------------------------------------------------------------------
   END SUBROUTINE DJDp_Adjoint_IceFluxSolver
!------------------------------------------------------------------------------
