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
! *  Authors: Olivier Gagliardini             
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 30. April 2010
! * 
! *****************************************************************************
!> SSolver to inquire the velocity from the SSA solution            
SUBROUTINE SMOOTHdrivingSTRESS_tau( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the in-plane basal velocity with the SSA solution !
!  To be computed only at the base. Use then the SSASolver to export verticaly 
!  the basal velocity and compute the vertical velocity and pressure (if needed)
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, HSol, &
                               VeloSol

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, CalvingFront 

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs
  INTEGER :: NonlinearIter, NewtonIter, iter, other_body_id
          
  INTEGER, POINTER :: Permutation(:), &
       ZsPerm(:), HPerm(:), &
       NodeIndexes(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Zs(:), HH(:)
                            
  REAL(KIND=dp) :: UNorm, cn, dd, NonlinearTol, rhow, sealevel, &
                   PrevUNorm, relativeChange,minv

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalGravity(:), NodalDensity(:), &
           NodalZs(:), NodalH(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!  REAL(KIND=dp) :: at, at0, CPUTime, RealTime
  REAL(KIND=dp) :: at, at0, coupling_length
       
  SAVE rhow,sealevel,coupling_length
  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, ElementNodes
  SAVE NodalGravity, NodalDensity, &
           NodalZs, NodalH,   &
           NodeIndexes

!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs 
  WRITE(SolverName, '(A)') 'FLUXdirection SOLVER'

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
        DIM = CoordinateSystemDimension()


        HSol => VariableGet( Solver % Mesh % Variables, 'htau' )
        IF (ASSOCIATED(HSol)) THEN
           HH => HSol % Values
           HPerm => HSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >tau<')
        END IF

        ZsSol => VariableGet( Solver % Mesh % Variables, 'surface' )
        IF (ASSOCIATED(ZsSol)) THEN
           Zs => ZsSol % Values
           ZsPerm => ZsSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >surface<')
        END IF
  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

     ! Get some constants
     rhow = GetConstReal( Model % Constants, 'Water Density', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant Water Density not found. &
                   &Setting to 1.03225e-18'
            CALL INFO(SolverName, Message, level=20)
            rhow = 1.03225e-18_dp
     End if

     sealevel = GetConstReal( Model % Constants, 'Sea Level', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant >Sea Level< not found. &
                   &Setting to 0.0'
            CALL INFO(SolverName, Message, level=20)
            sealevel=0.0_dp
     End if


     coupling_length = GetConstReal( Model % Constants, 'coupling length', Found )
     If (.NOT.Found) Then
            WRITE(Message,'(A)') 'Constant >Coupling Length< not found. &
                   &Setting to 10.0'
            CALL INFO(SolverName, Message, level=20)
            coupling_length=10.0_dp
     End if

     ! Allocate

     N = Model % MaxElementNodes
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
                       NodalDensity,  &
                       NodalH, NodalZs, &
                       ElementNodes % x, &
                       ElementNodes % y, ElementNodes % z )

     ALLOCATE( FORCE(STDOFs*N), LOAD(N), STIFF(STDOFs*N,STDOFs*N), &
          NodalGravity(N), NodalDensity(N), &
          NodalH(N), NodalZs(N) ,&
          ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
           STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      NonlinearTol = GetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

      NonlinearIter = GetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

      IF ( .NOT.GotIt ) NonlinearIter = 1


!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, &
                   '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'SSA BASAL VELOCITY NON-LINEAR ITERATION', iter
       CALL Info( SolverName, Message, Level=4 )
       CALL Info( SolverName, ' ', Level=4 )
       CALL Info( SolverName, &
                   '-------------------------------------',Level=4 )
       CALL Info( SolverName, ' ', Level=4 )


  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()

  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

 ! set coords of highest occuring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN !1D SSA
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN !2D SSA
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                STDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

     ! Read the gravity in the Body Force Section 
     BodyForce => GetBodyForce()
     NodalGravity = 0.0_dp
     IF ( ASSOCIATED( BodyForce ) ) THEN
           IF (STDOFs==1) THEN 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
           ELSE 
           NodalGravity(1:n) = ListGetReal( &
                   BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
           END IF
     END IF

     ! Read the Viscosity eta, density, and exponent m in MMaterial Section
     ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
     Material => GetMaterial(Element)

     
     NodalDensity=0.0_dp
     NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n,NodeIndexes,Found)
     IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop.  >SSA Mean Density<')


     ! Get the Nodal value of H and Zs
     NodalH(1:n) = HH(HPerm(NodeIndexes(1:n)))
     NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))


     CALL LocalMatrixTAUS (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
        NodalDensity, NodalH, NodalZs, STDOFs, coupling_length)

     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO
  CALL DefaultFinishBulkAssembly()
  
  CALL DefaultFinishAssembly()

  ! Dirichlet 
  CALL DefaultDirichletBCs()
  

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()


      RelativeChange = Solver % Variable % NonlinChange
      !IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
      !   RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      !ELSE
      !   RelativeChange = 0.0d0
      !END IF

      WRITE( Message, * ) 'Result Norm   : ', UNorm, PrevUNorm
      CALL Info(SolverName, Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ', RelativeChange
      CALL Info(SolverName, Message, Level=4 )


!------------------------------------------------------------------------------
      IF ( RelativeChange < NonLinearTol ) EXIT
!------------------------------------------------------------------------------

  END DO ! Loop Non-Linear Iterations

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixTAUS(  STIFF, FORCE, Element, n, Nodes, gravity, &
           Density, LocalH, LocalZs, STDOFs , cl)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
                     LocalH(:), LocalZs(:)
    INTEGER :: n, cp , STDOFs
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    REAL(KIND=dp) :: g, rho, h , kappa
    REAL(KIND=dp) :: gradS(2),tauD(2),cl
    LOGICAL :: Stat
    INTEGER :: i, j, t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    STIFF = 0.0d0
    FORCE = 0.0d0


    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

! Needed Intergration Point value

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       h = SUM( LocalH(1:n) * Basis(1:n) )
       IF(h.lt.10) h = 10
       kappa = (cl*h)
       kappa = kappa*kappa
       gradS = 0._dp
       gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
       tauD(1)  = (-1.0)*rho*g*h*gradS(1)
       if (STDOFs == 2) then
           gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
           tauD(2)  = (-1.0)*rho*g*h*gradS(2)
       endif

       DO p=1,n
         DO q=1,n

           DO i=1,STDOFs
             STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                  Basis(q) * Basis(p) * IP % S(t) * detJ + kappa * SUM(dBasisdx(p,1:STDOFs) * dBasisdx(q,1:STDOFs)) * IP % S(t) * detJ
           END DO

         END DO

         DO i=1,STDOFs
         FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) - &   
            tauD(i) * IP % s(t) * detJ * Basis(p) 
         END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixTAUS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SMOOTHdrivingSTRESS_tau
!------------------------------------------------------------------------------


