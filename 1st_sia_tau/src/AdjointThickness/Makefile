
PLATFORM := $(shell uname -s)
EXECUTABLES := AdjointThicknessSolvers.$(PLATFORM)

SRC=$(SRC_H)/src

_SOURCES := ThicknessSolver.f90\
            ThicknessSolver_mask.f90\
            ThicknessUncertaintySolver.f90\
            IceFluxSolver.f90\
            IceFluxUncertaintySolver.f90\
            inverseIceFluxUncertaintySolver.f90\
	    CostSolver_H_DISC_Quad.f90\
	    DJDp_Adjoint_ThicknessSolver.f90\
            DJDp_Adjoint_IceFluxSolver.f90\
	    AdjointSolverH.f90\
            AdjointSolverIceFlux.f90

SOURCES=$(patsubst $(SRC)/%, % ,$(_SOURCES))

_SOLVERS := $(patsubst %.f90, %.o ,$(SOURCES))
SOLVERS = $(patsubst %, $(OBJ)/% ,$(_SOLVERS))

compile: $(EXECUTABLES)

AdjointThicknessSolvers.$(PLATFORM): $(SOLVERS)
	@if test '$(PLATFORM)' = 'Darwin'; then \
	   echo "Dynamic linking solvers to $(PLATFORM)"; \
	   elmerf90 $(SOLVERS) -o AdjointThicknessSolvers.dylib ; \
	else \
	   if test '$(PLATFORM)' = 'Linux'; then \
	      echo "Dynamic linking solvers to $(PLATFORM)"; \
	      elmerf90  $(SOLVERS) -o AdjointThicknessSolvers.so ; \
	   fi \
	 fi

$(SOLVERS): $(OBJ)/%.o: $(SRC)/%.f90
	elmerf90 -c $< -o $@


clean: 
	rm -f $(EXECUTABLES)
	rm -f $(SOLVERS)
