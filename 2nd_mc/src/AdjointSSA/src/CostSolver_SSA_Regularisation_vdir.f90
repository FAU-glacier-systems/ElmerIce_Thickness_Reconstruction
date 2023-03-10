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
! *****************************************************************************
SUBROUTINE CostSolver_SSA_Regularisation_vdir( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Compute a regularisation term for SSA inverse problems and update the
!  gradient of the cost function with respect to the regularized variable.
!
!   Regularisation by default is: int_{Pb dimension} 0.5 * (d(var)/dx)**2 
!   A priori regularisation can also be used ( A priori Regularisation=True) :
!                                 int_{Pb dimension} 0.5 *(1/sigma**2)*(var-var{a_priori})**2
!
!     OUTPUT are : J and DJDvar
!                      
!
!    !!!!! BE carefull it will reset Cost and DJ to 0 by default !!!!
!      !!! If other cost and gradient are computed before (i.e. from the adjoint pb, 
!       use "<Reset Cost Value> = False" to add cost and gradient to previously computed values !!!
!
!
!       Required Sif parameters are:
!
!          In the solver section:
!               Problem Dimension=Integer (default:coordinate sytem dimension),
!               Cost Filename=File (default: CostOfT.dat),
!               Optimized Variable Name= String (default='Beta'),
!               Gradient Variable Name= String (default = 'DJDBeta'),
!               Cost Variable Name= String (default='Cost Value'),
!               Lambda= Real (default 1.0),
!               Reset Cost Value= Logical (default = True),
!               A priori Regularisation= Logical (defualt = .False.),
!
!          In Body Force section:
!               <VarSolName> a priori value = Real (defualt =0.0),
!               <VarSolName> variance = real (defualt=1.0)
!
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
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,varname,DirectionSolName1,DirectionSolName2
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar

  TYPE(Variable_t), POINTER :: Variable,DJDVariable,DirectionVariable1,DirectionVariable2
  REAL(KIND=dp), POINTER :: Values(:),DJDValues(:),DirectionValues1(:),DirectionValues2(:)
  INTEGER, POINTER :: Perm(:),DJDPerm(:),DirectionPerm1(:),DirectionPerm2(:)

  TYPE(ValueList_t), POINTER :: SolverParams,BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
  INTEGER, POINTER :: NodeIndexes(:)

  Logical :: Firsttime=.true.,Found,Parallel,stat,Gotit
  integer :: i,j,k,l,t,n,NMAX,DIM,ierr,c

  real(kind=dp) :: Cost,Cost_S,Lambda,normit,scaler
  real(kind=dp) :: u,v,w,s,coeff_reg,SqrtElementMetric,x
  real(kind=dp) :: coeff_reg_per,coeff_reg_par,nn1,nn2

  REAL(KIND=dp),dimension(:),allocatable,SAVE :: coeff_reg_par2,coeff_reg_per2
  REAL(KIND=dp),dimension(:),allocatable,SAVE :: NodeAp,NodeRMS,NodeValues,NodalRegb,NodalRegb2
  REAL(KIND=dp),dimension(:),allocatable,SAVE :: NodeDirectionValues1,NodeDirectionValues2
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: IPerr,IPvar

  LOGICAL :: Apriori,Reset

  CHARACTER*10 :: date,temps

  save Firsttime,Parallel 
  save SolverName,CostSolName,VarSolName,Lambda,CostFile,DirectionSolName1,DirectionSolName2
  save ElementNodes

   WRITE(SolverName, '(A)') 'CostSolver_Regular'

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


!! Name of the variable with teh cost function
   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
          END IF

!! Direction for weighting factors - X-xomponent
   DirectionSolName1 =  GetString( SolverParams,'Direction Component X', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Direction Component X< not found in section >Solver<')
          END IF

!! Direction for weighting factors - Y-component
   DirectionSolName2 =  GetString( SolverParams,'Direction Component Y', Found)
          IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Direction Component Y< not found in section >Solver<')
          END IF

!! Optional weighting term
   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value Lambda=1.0')
           Lambda = 1.0
   End if

!! Optional direction weighting term
   scaler =  GetConstReal( SolverParams,'DirectionScaler', Found)
   IF(.NOT.Found) THEN
           CALL WARN(SolverName,'Keyword >DirectionScaler< not found  in section >Solver<')
           CALL WARN(SolverName,'Taking default value DirectionScaler=1.0')
           Scaler = 1.0
   End if

!! Do we need to reset cost and DJDVar to 0? Default YES
   Reset =  GetLogical( SolverParams,'Reset Cost Value', Found)
            IF(.NOT.Found) Reset=.True.

!! What type of regularistaion ? Default penalise 1st derivative
   Apriori =  GetLogical( SolverParams,'A priori Regularisation', Found)
            IF(.NOT.Found) Apriori=.False.

   !print *,'SCALER',scaler
   !print *,'DirectionSolName1',DirectionSolName1
   !print *,'DirectionSolName2',DirectionSolName2

!!! SOME INITIALISATION AT FIRST TIME
  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(NodeAp(N),NodeRMS(N),NodeValues(N),NodalRegb(N),NodeDirectionValues1(N),NodeDirectionValues2(N))
    allocate(NodalRegb2(N))
   allocate(coeff_reg_par2(N),coeff_reg_per2(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

!!!!!!!!!!!  initiate Cost File

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
  
  !!! End of First visit
    Firsttime=.false.
  Endif
!!!! INITIALISATION DONE

    DirectionVariable1 => VariableGet( Solver % Mesh % Variables, DirectionSolName1  )
    IF ( ASSOCIATED( DirectionVariable1 ) ) THEN
            DirectionValues1 => DirectionVariable1 % Values
            DirectionPerm1 => DirectionVariable1 % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable DSN1>',DirectionSolName1,' < found'
            CALL FATAL(SolverName,Message)
    END IF

    DirectionVariable2 => VariableGet( Solver % Mesh % Variables, DirectionSolName2  )
    IF ( ASSOCIATED( DirectionVariable2 ) ) THEN
            DirectionValues2 => DirectionVariable2 % Values
            DirectionPerm2 => DirectionVariable2 % Perm
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable DSN2>',DirectionSolName2,' < found'
            CALL FATAL(SolverName,Message)
    END IF
    !print *,'READ IN DIRECTIONS'

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
    ELSE
            WRITE(Message,'(A,A,A)') &
                               'No variable >',VarSolName,' < found'
            CALL FATAL(SolverName,Message)
    END IF  
    IF (Reset) DJDValues=0.0_dp


    Cost=0._dp

    DO t=1,Solver % NumberOfActiveElements
       Element => GetActiveElement(t)
       IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
       n = GetElementNOFNodes()

       NodeIndexes => Element % NodeIndexes

 ! set coords of highest occuring dimension to zero (to get correct path element)
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
           STOP
        END IF

 ! Compute inetgrated cost


      IF (Apriori) then
         BodyForce => GetBodyForce()
         write(varname,'(A,A)') trim(VarSolName),' a priori value'
         NodeAp(1:n) = 0._dp
         NodeAp(1:n) = ListGetReal( BodyForce, trim(varname), n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A,A,A)') &
                     'No variable >',trim(varname),'< found in "Body Forces" default is 0'
                  CALL Info(SolverName,Message,level=6)
          END IF 
          write(varname,'(A,A)') trim(VarSolName),' variance'
          NodeRMS(1:n)=ListGetReal( BodyForce, trim(varname), n, NodeIndexes, GotIt)
          IF (.NOT.GotIt) Then
                  WRITE(Message,'(A,A,A)') &
                     'No variable >',trim(varname),'< found in "Body Forces" default is 1'
                  CALL Info(SolverName,Message,level=6)
                  NodeRMS=1.0_dp
          END IF 
     END IF
      !print *,'NodeDirectionValues 0'
 ! Nodal values of the variable        
      NodeValues(1:n)=Values(Perm(NodeIndexes(1:n)))
      NodeDirectionValues1(1:n)=DirectionValues1(DirectionPerm1(NodeIndexes(1:n)))
      NodeDirectionValues2(1:n)=DirectionValues2(DirectionPerm2(NodeIndexes(1:n)))
      !print *,'NodeDirectionValues 1'

      normit = 0.0_dp
      do i=1,n
      normit = (NodeDirectionValues1(i)**2.0+NodeDirectionValues2(i)**2.0)**0.5
      if(normit.gt.0)then
      NodeDirectionValues1(i)=NodeDirectionValues1(i)/normit
      NodeDirectionValues2(i)=NodeDirectionValues2(i)/normit
      endif
      enddo

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------

        NodalRegb = 0.0_dp

        IntegStuff = GaussPoints( Element )

        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )

          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          s = 1.0d0

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             s = 2.0d0 * PI * x 
          END IF
          s = s * SqrtElementMetric * IntegStuff % s(i)
         
          IF (Apriori) then
             IPerr = SUM((NodeValues(1:n)-NodeAp(1:n))*Basis(1:n))
             IPvar = SUM(NodeRMS(1:n)*Basis(1:n))
             coeff_reg=IPerr/IPvar
             coeff_reg =  coeff_reg*coeff_reg

             !Now compute the derivative
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*IPerr*Basis(1:n)/(IPVar**2.0)
!          Else
!             coeff_reg = SUM(NodeValues(1:n) * dBasisdx(1:n,1))
!             coeff_reg =  coeff_reg*coeff_reg
!             IF (DIM.eq.2) then
!                  coeff_reg=coeff_reg+ &
!                  SUM(NodeValues(1:n)*dBasisdx(1:n,2))*SUM(NodeValues(1:n) * dBasisdx(1:n,2))
!             END IF
!
!             !Now compute the derivative
!               NodalRegb(1:n)=NodalRegb(1:n)+&
!                    s*Lambda*SUM(NodeValues(1:n)* dBasisdx(1:n,1))*dBasisdx(1:n,1)
!               IF (DIM.eq.2) then
!                  NodalRegb(1:n)=NodalRegb(1:n)+&
!                          s*Lambda*SUM(NodeValues(1:n)*dBasisdx(1:n,2))*dBasisdx(1:n,2)
!               End if
! 
          Else
             coeff_reg = SUM(NodeValues(1:n) * dBasisdx(1:n,1))  !!!!
             coeff_reg =  coeff_reg*coeff_reg                    !!!!
             IF (DIM.eq.1) THEN
             coeff_reg = SUM(NodeValues(1:n) * dBasisdx(1:n,1))
             coeff_reg = coeff_reg*coeff_reg
             ELSEIF (DIM.eq.2) then
             nn1 = SUM(NodeDirectionValues1(1:n)*Basis(1:n))
             nn2 = SUM(NodeDirectionValues2(1:n)*Basis(1:n))
             normit = nn1*nn1+nn2*nn2
             nn1 = nn1/sqrt(normit)
             nn2 = nn2/sqrt(normit)
             coeff_reg=coeff_reg+ &                              !!!!
                   SUM(NodeValues(1:n)*dBasisdx(1:n,2))*SUM(NodeValues(1:n) * dBasisdx(1:n,2)) !!!!
             coeff_reg_par =                 SUM(NodeValues(1:n) * dBasisdx(1:n,1))*nn1
             coeff_reg_par = coeff_reg_par + SUM(NodeValues(1:n) * dBasisdx(1:n,2))*nn2
             coeff_reg_par = coeff_reg_par*coeff_reg_par
             coeff_reg_per =                 SUM(NodeValues(1:n) * dBasisdx(1:n,1))*nn2*(-1.0)
             coeff_reg_per = coeff_reg_per + SUM(NodeValues(1:n) * dBasisdx(1:n,2))*nn1
             coeff_reg_per = coeff_reg_per*coeff_reg_per
             !print *,'SCALERY',scaler,coeff_reg,scaler*coeff_reg_par + coeff_reg_per,nn1*nn1+nn2*nn2
 
             coeff_reg = coeff_reg_par + scaler*coeff_reg_per
             END IF
       
             !Now compute the derivative
               NodalRegb2(1:n)=NodalRegb(1:n)+&                                           !!!!
                    s*Lambda*SUM(NodeValues(1:n) * dBasisdx(1:n,1))*dBasisdx(1:n,1)       !!!!
               IF (DIM.eq.1) THEN
               NodalRegb(1:n)=NodalRegb(1:n)+&
                    s*Lambda*SUM(NodeValues(1:n)*dBasisdx(1:n,1))*dBasisdx(1:n,1)
               ELSEIF (DIM.eq.2) then
                  NodalRegb2(1:n)=NodalRegb2(1:n)+&                                       !!!!
                          s*Lambda*SUM(NodeValues(1:n)*dBasisdx(1:n,2))*dBasisdx(1:n,2)   !!!!
                  coeff_reg_par2(1:n) = 0.0
                  coeff_reg_par2(1:n) = (       SUM(NodeValues(1:n) * dBasisdx(1:n,1))*nn1 + SUM(NodeValues(1:n) * dBasisdx(1:n,2))*nn2)*(       dBasisdx(1:n,1)*nn1+dBasisdx(1:n,2)*nn2)
                  !coeff_reg_par2(1:n) = dBasisdx(1:n,1)
                  coeff_reg_per2(1:n) = 0.0
                  coeff_reg_per2(1:n) = ((-1.0)*SUM(NodeValues(1:n) * dBasisdx(1:n,1))*nn2 + SUM(NodeValues(1:n) * dBasisdx(1:n,2))*nn1)*((-1.0)*dBasisdx(1:n,1)*nn2+dBasisdx(1:n,2)*nn1)
                  !coeff_reg_per2(1:n) = coeff_reg_per2(1:n) + SUM(NodeValues(1:n) * dBasisdx(1:n,2))*dBasisdx(1:n,2)*nn1

                  NodalRegb(1:n)=NodalRegb(1:n) +        s*Lambda*coeff_reg_par2(1:n)
                  NodalRegb(1:n)=NodalRegb(1:n) + scaler*s*Lambda*coeff_reg_per2(1:n)
                  !print *,'SCALERY',NodalRegb2(2),NodalRegb(2),nn1*nn1+nn2*nn2

               End if
          Endif


          Cost=Cost+0.5*Lambda*coeff_reg*s

        End do !IP

        DJDValues(DJDPerm(NodeIndexes(1:n)))=DJDValues(DJDPerm(NodeIndexes(1:n))) + NodalRegb(1:n)

    End do !Elements

! 30/01/2014 Found this error
!   Cost=0.5*Lambda*Cost

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
END SUBROUTINE CostSolver_SSA_Regularisation_vdir
!------------------------------------------------------------------------------
! *****************************************************************************
