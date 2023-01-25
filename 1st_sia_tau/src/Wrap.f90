
      SUBROUTINE Wrap( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      TYPE(Variable_t), POINTER :: DJDU,DJDSMBT,DJVar
      INTEGER :: i


      !Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )
      !DJDU => VariableGet( Solver % Mesh % Variables, 'DJDUV'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')
      DJDSMBT => VariableGet( Solver % Mesh % Variables, 'DJDsmbTop')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')
      DJVar => VariableGet( Solver % Mesh % Variables, 'DJDVar')

      Do i=1,Solver%Mesh%NumberOfNodes
        !Var%Values(3*(Var%Perm(i)-1)+1)=Velocity%Values(2*(Velocity%Perm(i)-1)+1)
        !Var%Values(3*(Var%Perm(i)-1)+2)=Velocity%Values(2*(Velocity%Perm(i)-1)+2)
        !Var%Values(3*(Var%Perm(i)-1)+3)=TOPACCU%Values(TOPACCU%Perm(i))
        Var%Values(Var%Perm(i))=TOPACCU%Values(TOPACCU%Perm(i))/1

        !DJVar%Values(3*(DJVar%Perm(i)-1)+1)=DJDU%Values(2*(DJDU%Perm(i)-1)+1)
        !DJVar%Values(3*(DJVar%Perm(i)-1)+2)=DJDU%Values(2*(DJDU%Perm(i)-1)+2)
        !DJVar%Values(3*(DJVar%Perm(i)-1)+3)=DJDSMBT%Values(DJDSMBT%Perm(i))
        DJVar%Values(DJVar%Perm(i))=DJDSMBT%Values(DJDSMBT%Perm(i))/1

      End do

      END SUBROUTINE Wrap

      SUBROUTINE UnWrap( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      INTEGER :: i

      !Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')

      Do i=1,Solver%Mesh%NumberOfNodes
        !Velocity%Values(2*(Velocity%Perm(i)-1)+1)=Var%Values(3*(Var%Perm(i)-1)+1)
        !Velocity%Values(2*(Velocity%Perm(i)-1)+2)=Var%Values(3*(Var%Perm(i)-1)+2)
        !TOPACCU%Values(TOPACCU%Perm(i))=Var%Values(3*(Var%Perm(i)-1)+3)
        TOPACCU%Values(TOPACCU%Perm(i))=Var%Values(Var%Perm(i))*1
      End do

      END SUBROUTINE UnWrap

!!!!!!!!      
!!!! changement de variable U_opt=U/UAdim SMB_opt=SMB/SMB_Adim
      SUBROUTINE WrapC( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      TYPE(Variable_t), POINTER :: DJDU,DJDSMBT,DJVar
      REAL(KIND=dp), PARAMETER :: UAdim=50.0_dp,SMBAdim=1.0_dp
      INTEGER :: i


      Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )
      DJDU => VariableGet( Solver % Mesh % Variables, 'DJDUV'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')
      DJDSMBT => VariableGet( Solver % Mesh % Variables, 'DJDsmbTop')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')
      DJVar => VariableGet( Solver % Mesh % Variables, 'DJDVar')

      Do i=1,Solver%Mesh%NumberOfNodes
        Var%Values(3*(Var%Perm(i)-1)+1)=Velocity%Values(2*(Velocity%Perm(i)-1)+1)/UAdim
        Var%Values(3*(Var%Perm(i)-1)+2)=Velocity%Values(2*(Velocity%Perm(i)-1)+2)/UAdim
        Var%Values(3*(Var%Perm(i)-1)+3)=TOPACCU%Values(TOPACCU%Perm(i))/SMBAdim

        DJVar%Values(3*(DJVar%Perm(i)-1)+1)=DJDU%Values(2*(DJDU%Perm(i)-1)+1)*UAdim
        DJVar%Values(3*(DJVar%Perm(i)-1)+2)=DJDU%Values(2*(DJDU%Perm(i)-1)+2)*UAdim
        DJVar%Values(3*(DJVar%Perm(i)-1)+3)=DJDSMBT%Values(DJDSMBT%Perm(i))*SMBAdim

      End do

      END SUBROUTINE WrapC

      SUBROUTINE UnWrapC( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      REAL(KIND=dp), PARAMETER :: UAdim=50.0_dp,SMBAdim=1.0_dp
      INTEGER :: i

      Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')

      Do i=1,Solver%Mesh%NumberOfNodes
        Velocity%Values(2*(Velocity%Perm(i)-1)+1)=Var%Values(3*(Var%Perm(i)-1)+1)*UAdim
        Velocity%Values(2*(Velocity%Perm(i)-1)+2)=Var%Values(3*(Var%Perm(i)-1)+2)*UAdim
        TOPACCU%Values(TOPACCU%Perm(i))=Var%Values(3*(Var%Perm(i)-1)+3)*SMBAdim
      End do

      END SUBROUTINE UnWrapC
!!!!!!!!      
!!!! changement de variable U_opt=U/UAdim SMB_opt=SMB/SMB_Adim from sif file
      SUBROUTINE WrapC_param( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      TYPE(Variable_t), POINTER :: DJDU,DJDSMBT,DJVar
      REAL(KIND=dp):: UAdim,SMBAdim
      INTEGER :: i

      TYPE(ValueList_t), POINTER :: Material
      CHARACTER(LEN=MAX_NAME_LEN) ::readstring
      LOGICAL :: FirstTime = .True.
      LOGICAL :: GotIt

      SAVE FirstTime,UAdim,SMBAdim

       IF (FirstTime) THEN
          FirstTime = .FALSE.
          Material => GetMaterial(Model % CurrentElement)
          IF (.NOT.ASSOCIATED(Material)) CALL FATAL("MyRoutineName","Material not found. Exit")

      UAdim =  ListGetCReal( Material, 'UAdim', GotIt )
      IF (.NOT.GotIt) CALL FATAL("MyRoutineName","UAdim not found! Exit")
      SMBAdim =  ListGetCReal( Material, 'SMBAdim', GotIt )
      IF (.NOT.GotIt) CALL FATAL("MyRoutineName","SMBAdim not found! Exit")
     ENDIF

      Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )
      DJDU => VariableGet( Solver % Mesh % Variables, 'DJDUV'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')
      DJDSMBT => VariableGet( Solver % Mesh % Variables, 'DJDsmbTop')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')
      DJVar => VariableGet( Solver % Mesh % Variables, 'DJDVar')

      Do i=1,Solver%Mesh%NumberOfNodes
        Var%Values(3*(Var%Perm(i)-1)+1)=Velocity%Values(2*(Velocity%Perm(i)-1)+1)/UAdim
        Var%Values(3*(Var%Perm(i)-1)+2)=Velocity%Values(2*(Velocity%Perm(i)-1)+2)/UAdim
        Var%Values(3*(Var%Perm(i)-1)+3)=TOPACCU%Values(TOPACCU%Perm(i))/SMBAdim

        DJVar%Values(3*(DJVar%Perm(i)-1)+1)=DJDU%Values(2*(DJDU%Perm(i)-1)+1)*UAdim
        DJVar%Values(3*(DJVar%Perm(i)-1)+2)=DJDU%Values(2*(DJDU%Perm(i)-1)+2)*UAdim
        DJVar%Values(3*(DJVar%Perm(i)-1)+3)=DJDSMBT%Values(DJDSMBT%Perm(i))*SMBAdim

      End do

      END SUBROUTINE WrapC_param

      SUBROUTINE UnWrapC_param( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE(Variable_t), POINTER :: Velocity,TOPACCU,Var
      REAL(KIND=dp):: UAdim,SMBAdim
      INTEGER :: i

      TYPE(ValueList_t), POINTER :: Material
      LOGICAL :: FirstTime = .True.
      LOGICAL :: GotIt

      SAVE FirstTime,UAdim,SMBAdim

       IF (FirstTime) THEN
          FirstTime = .FALSE.
          Material => GetMaterial(Model % CurrentElement)
          IF (.NOT.ASSOCIATED(Material)) CALL FATAL("MyRoutineName","Material not found. Exit")

      UAdim =  ListGetCReal( Material, 'UAdim', GotIt )
      IF (.NOT.GotIt) CALL FATAL("MyRoutineName","UAdim not found! Exit")
      SMBAdim =  ListGetCReal( Material, 'SMBAdim', GotIt )
      IF (.NOT.GotIt) CALL FATAL("MyRoutineName","SMBAdim not found! Exit")
      ENDIF

      Velocity => VariableGet( Solver % Mesh % Variables, 'Velocity'  )

      TOPACCU => VariableGet( Solver % Mesh % Variables, 'TopAccumulation')

      Var => VariableGet( Solver % Mesh % Variables, 'Var')

      Do i=1,Solver%Mesh%NumberOfNodes
        Velocity%Values(2*(Velocity%Perm(i)-1)+1)=Var%Values(3*(Var%Perm(i)-1)+1)*UAdim
        Velocity%Values(2*(Velocity%Perm(i)-1)+2)=Var%Values(3*(Var%Perm(i)-1)+2)*UAdim
        TOPACCU%Values(TOPACCU%Perm(i))=Var%Values(3*(Var%Perm(i)-1)+3)*SMBAdim
      End do

      END SUBROUTINE UnWrapC_param
