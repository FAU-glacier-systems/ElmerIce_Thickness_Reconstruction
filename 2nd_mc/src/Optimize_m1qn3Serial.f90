!  Optimize a cost function 
!    using quasi-Newton M1QN3 Routine in Reverse Communication
!    Using Euclidian inner product
!
!  Serial Only   2D/3D
!
!  Need:
!  - Value of the Cost function
!  - Value of the variable to optimize
!  - Value of gradient of cost function with respect to Variable
!      (sum the contribution of each partition shared node)
!  - Optimisation Mask Variable (Optional): name of a mask variable. If
!  mask.lt.0 the variable is considered fixed
!
! => Update the new value of variable to optimize
!
! *****************************************************************************
SUBROUTINE Optimize_m1qn3Serial( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  TYPE(Element_t),POINTER ::  Element
! Variables Beta,DJDbeta and Cost  
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: BetaVar,CostVar,GradVar,MaskVar,TimeVar
  REAL(KIND=dp), POINTER :: BetaValues(:),CostValues(:),GradValues(:),MaskValues(:)
  INTEGER, POINTER :: BetaPerm(:),GradPerm(:),MaskPerm(:),NodeIndexes(:)

  REAL(KIND=dp),allocatable :: x(:),g(:)
  REAL(KIND=dp) :: f,NormG

  integer :: i,j,t,n,NMAX,NActiveNodes
  integer,allocatable :: ActiveNodes(:)
  integer,allocatable :: NewNode(:)
  integer :: VarDOFs

  Logical :: FirstVisit=.true.,Found,UseMask,ComputeNormG=.False.
  logical,allocatable :: VisitedNode(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VarSolName,GradSolName,NormM1QN3,MaskVarName,NormFile
  CHARACTER*10 :: date,temps

!Variables for m1qn3
  external simul_rc,euclid,ctonbe,ctcabe
  character*3 normtype
  REAL(KIND=dp) :: dxmin,df1,epsrel,dzs(1)
  real(kind=dp), allocatable :: dz(:)
  REAL :: rzs(1)
  integer :: imp,io=20,imode(3),omode=-1,niter,nsim,iz(5),ndz,reverse,indic,izs(1)
  integer :: ierr,Npes,ntot
  CHARACTER(LEN=MAX_NAME_LEN) :: IOM1QN3
  logical :: DISbool
!
  save NActiveNodes
  save x,g
  save ActiveNodes
  save VarDOFs
  save ComputeNormG
  save normtype,dxmin,df1,epsrel,dz,dzs,rzs,imp,io,imode,omode,niter,nsim,iz,ndz,reverse,indic,izs
  save FirstVisit
  save SolverName
  save CostSolName,VarSolName,GradSolName,IOM1QN3,NormFile


!  Read Constant from sif solver section
      IF(FirstVisit) Then
            FirstVisit=.FALSE.
            WRITE(SolverName, '(A)') 'Optimize_m1qn3Serial'

           ! Check we have a parallel run
           IF(ASSOCIATED(Solver %  Matrix % ParMatrix)) Then
             CALL FATAL(SolverName,'ParMatrix associated! Ths solver for serial only!!')
           End if
           !!


          SolverParams => GetSolverParams()
          MaskVarName = GetString( SolverParams,'Optimisation Mask Variable',UseMask)
            IF (UseMask) Then
                MaskVar => VariableGet( Solver % Mesh % Variables, MaskVarName ) 
                IF (ASSOCIATED(MaskVar)) THEN 
                   MaskValues => MaskVar % Values 
                   MaskPerm => MaskVar % Perm 
               ELSE 
                   WRITE(Message,'(A,A,A)') 'No variable >',MaskVarName,'< found' 
                   CALL FATAL(SolverName,Message) 
               ENDIF
            ENDIF

            VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >Beta<')
                    WRITE(VarSolName,'(A)') 'Beta'
                END IF
            GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >DJDB<')
                    WRITE(GradSolName,'(A)') 'DJDB'
                END IF
            BetaVar => VariableGet( Solver % Mesh % Variables, VarSolName ) 
            IF (ASSOCIATED(BetaVar)) THEN 
              VarDOFs = BetaVar % DOFs
            ELSE 
               WRITE(Message,'(A,A,A)') 'No variable >',VarSolName,'< found' 
               CALL FATAL(SolverName,Message) 
            ENDIF

           GradVar => VariableGet( Solver % Mesh % Variables, GradSolName) 
           IF (ASSOCIATED(GradVar)) THEN 
              if (GradVar % DOFs.NE.VarDOFs) then
                  WRITE(Message,'(A,A,A,A,A)') &
                     'variables',GradSolName,'and',VarSolName,'do not have the same DOFs'
                  CALL FATAL(SolverName,Message)
              endif
           ELSE 
             WRITE(Message,'(A,A,A)') 'No variable >',GradSolName,'< found' 
             CALL FATAL(SolverName,Message)    
          END IF
         !!
          NMAX=Solver % Mesh % NumberOfNodes
          allocate(VisitedNode(NMAX),NewNode(NMAX))

             

!!!!!!!!!!!!find active nodes 
           VisitedNode=.false.  
           NewNode=-1

           NActiveNodes=0 
           DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t)
              n = GetElementNOFNodes()
              NodeIndexes => Element % NodeIndexes
              Do i=1,n
                 if (VisitedNode(NodeIndexes(i))) then
                     cycle
                 else
                     VisitedNode(NodeIndexes(i))=.true.
                     IF (UseMask) Then
                             IF (MaskValues(MaskPerm(NodeIndexes(i))).lt.0) cycle
                     END IF
                     NActiveNodes=NActiveNodes+1
                     NewNode(NActiveNodes)=NodeIndexes(i)
                 endif
             End do
           End do

           if (NActiveNodes.eq.0) THEN
              WRITE(Message,'(A)') 'NActiveNodes = 0 !!!'
              CALL FATAL(SolverName,Message)
           End if


           allocate(ActiveNodes(NActiveNodes),x(VarDOFs*NActiveNodes),g(VarDOFs*NActiveNodes))
           ActiveNodes(1:NActiveNodes)=NewNode(1:NActiveNodes)

           deallocate(VisitedNode,NewNode)

!!!!!!!  Solver Params

            CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >CostValue<')
                    WRITE(CostSolName,'(A)') 'CostValue'
                END IF


             NormFile=GetString( SolverParams,'gradient Norm File',Found)
             IF(Found)  Then
                 ComputeNormG=.True.
                 open(io,file=trim(NormFile))
                    CALL DATE_AND_TIME(date,temps)
                    write(io,'(a1,a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)')'#',date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                 close(io)
             END IF

!!  initialization of m1qn3 variables
            dxmin=GetConstReal( SolverParams,'M1QN3 dxmin', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 dxmin< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-10<')
                    dxmin=1.e-10
                END IF
            epsrel=GetConstReal( SolverParams,'M1QN3 epsg', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 epsg< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >1.e-06<')
                    epsrel=1.e-6
                END IF
            niter=GetInteger(SolverParams,'M1QN3 niter', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 niter< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    niter=200
                END IF
            nsim=GetInteger(SolverParams,'M1QN3 nsim', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 nsim< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >200<')
                    nsim=200
                END IF
            imp=GetInteger(SolverParams,'M1QN3 impres', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 impres< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >5<')
                    imp=5
                END IF
              ndz=GetInteger( SolverParams,'M1QN3 ndz', Found)
                  IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 ndz< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >5< update')
                       ndz=5
                   END IF
            DISbool=GetLogical( SolverParams, 'M1QN3 DIS Mode', Found)
                IF(.NOT.Found) THEN
                    CALL WARN(SolverName,'Keyword >M1QN3 DIS Mode< not found  in section >Solver<')
                    CALL WARN(SolverName,'Taking default value >FALSE<')
                    DISbool=.False.
                END IF
                if(DISbool) then
                    imode(1)=0 !DIS Mode
                else
                    imode(1)=1 !SIS Mode
                End if
                IF (DISbool) then
                   ndz=4*NActiveNodes+ndz*(2*NActiveNodes+1)+10
                else
                   ndz=3*NActiveNodes+ndz*(2*NActiveNodes+1)+10
               end if
           allocate(dz(ndz))
            df1=GetConstReal( SolverParams,'M1QN3 df1', Found)
                IF(.NOT.Found) THEN
                   CALL WARN(SolverName,'Keyword >M1QN3 df1< not found  in section >Solver<')
                   CALL WARN(SolverName,'Taking default value >0.2<')
                   df1=0.2
                End if
                CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
                IF (ASSOCIATED(CostVar)) THEN
                    CostValues => CostVar % Values
                 ELSE
                     WRITE(Message,'(A,A,A)') 'No variable >',CostSolName,'< found'
                     CALL FATAL(SolverName,Message)
                 ENDIF
                 df1=CostValues(1)*df1
             NormM1QN3 = GetString( SolverParams,'M1QN3 normtype', Found)
                 IF((.NOT.Found).AND.((NormM1QN3(1:3).ne.'dfn').OR.(NormM1QN3(1:3).ne.'sup') &
                     .OR.(NormM1QN3(1:3).ne.'two'))) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 normtype< not good in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >dfn<')
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = 'dfn'
                  ELSE
                       PRINT *,'M1QN3 normtype  ',NormM1QN3(1:3)
                       normtype = NormM1QN3(1:3)
                  END IF
              IOM1QN3 = GetString( SolverParams,'M1QN3 OutputFile', Found)
                 IF(.NOT.Found) THEN
                       CALL WARN(SolverName,'Keyword >M1QN3 OutputFile< not found  in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >M1QN3.out<')
                       WRITE(IOM1QN3,'(A)') 'M1QN3.out'
                 END IF
                 open(io,file=trim(IOM1QN3))
                    CALL DATE_AND_TIME(date,temps)
                    write(io,*) '******** M1QN3 Output file ************'
                    write(io,'(a2,a1,a2,a1,a4,5x,a2,a1,a2,a1,a2)') date(5:6),'/',date(7:8),'/',date(1:4), &
                                 temps(1:2),':',temps(3:4),':',temps(5:6)
                    write(io,*) '*****************************************'
                 close(io)

                    imode(2)=0 
                    imode(3)=0 
                    reverse=1 
                    omode=-1 
                    dzs=0.0 
                    rzs=0.0
                    izs=0

        End if


! Omode from previous iter; if > 0 m1qn3 has terminated => return 
     IF (omode.gt.0) then 
             WRITE(Message,'(a,I3)') 'm1qn3 finished; omode=',omode 
             CALL Info(SolverName, Message, Level=1) 
             return  
     End if

!  Get Variables CostValue, Beta and DJDBeta
     CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
     IF (ASSOCIATED(CostVar)) THEN 
             CostValues => CostVar % Values 
     ELSE
            WRITE(Message,'(A,A,A)') 'No variable >',CostSolName,'< found' 
            CALL FATAL(SolverName,Message) 
    ENDIF 
    f=CostValues(1)

     BetaVar => VariableGet( Solver % Mesh % Variables, VarSolName ) 
     IF (ASSOCIATED(BetaVar)) THEN 
             BetaValues => BetaVar % Values 
             BetaPerm => BetaVar % Perm 
     ELSE 
             WRITE(Message,'(A,A,A)') 'No variable >',VarSolName,'< found' 
             CALL FATAL(SolverName,Message) 
     ENDIF

     GradVar => VariableGet( Solver % Mesh % Variables, GradSolName) 
     IF (ASSOCIATED(GradVar)) THEN 
             GradValues   => GradVar % Values 
             GradPerm => GradVar % Perm 
     ELSE 
             WRITE(Message,'(A,A,A)') 'No variable >',GradSolName,'< found' 
             CALL FATAL(SolverName,Message)    
     END IF

     Do i=1,NActiveNodes
       Do j=1,VarDOFs
       x(VarDOFs*(i-1)+j)=BetaValues(VarDOFs*(BetaPerm(ActiveNodes(i))-1)+j)
       g(VarDOFs*(i-1)+j)=GradValues(VarDOFs*(GradPerm(ActiveNodes(i))-1)+j)
       End DO
     End Do

     If (ComputeNormG) then
             TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
             NormG=0.0_dp
             Do i=1,NActiveNodes
                NormG=NormG+g(i)*g(i)
             End do
             open(io,file=trim(NormFile),position='append')
                write(io,'(e13.5,2x,e15.8)') TimeVar % Values(1),sqrt(NormG)
             close(io)
     End if        

     ! go to minimization
      open(io,file=trim(IOM1QN3),position='append')
       call m1qn3 (simul_rc,Euclid,ctonbe,ctcabe,VarDOFs*NActiveNodes,x,f,g,dxmin,df1, &
                        epsrel,normtype,imp,io,imode,omode,niter,nsim,iz, &
                        dz,ndz,reverse,indic,izs,rzs,dzs)

     close(io)
     WRITE(Message,'(a,E15.8,x,I2)') 'm1qn3: Cost,omode= ',f,omode
     CALL Info(SolverName, Message, Level=3)

     ! Update Beta Values 
     Do i=1,NActiveNodes
       Do j=1,VarDOFs
         BetaValues(VarDOFs*(BetaPerm(ActiveNodes(i))-1)+j)=x(VarDOFs*(i-1)+j)
       End DO
     End Do

   Return
!------------------------------------------------------------------------------
END SUBROUTINE Optimize_m1qn3Serial
!------------------------------------------------------------------------------


