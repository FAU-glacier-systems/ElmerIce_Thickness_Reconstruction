# SOURCES and EXECUTABLES DIRECTORIES
EXEDIR=./Executables
SRCDIR=./src

export OBJ=./objects
export SRC_H=$(PWD)/src/AdjointThickness
export SRC_SSA=./src/AdjointSSA

# SOLVERS TO COMPILE
LIBS= $(EXEDIR)/AdjointThicknessSolvers\
     $(EXEDIR)/AdjointSSASolvers\
     $(EXEDIR)/Optimize_m1qn3Serial\
     $(EXEDIR)/Optimize_m1qn3Parallel_Multi.so\
     $(EXEDIR)/UpdateExport \
     $(EXEDIR)/POI_output \
     $(EXEDIR)/POI_output_final \
     $(EXEDIR)/POI_output_thickness \
     $(EXEDIR)/USFs\
     $(EXEDIR)/Wrap\
     $(EXEDIR)/EffectiveMB \
     $(EXEDIR)/NodalSlopeSolver \
     $(EXEDIR)/CostSolver_Hpos \
     $(EXEDIR)/CostSolver_Hobs \
     $(EXEDIR)/CostSolver_Hmarine \
     $(EXEDIR)/direction_coupling


all: $(LIBS) 


$(EXEDIR)/AdjointThicknessSolvers:
	make -f $(SRC_H)/Makefile
	mv AdjointThicknessSolvers.so $(EXEDIR)/AdjointThicknessSolvers.so

$(EXEDIR)/direction_coupling: $(SRCDIR)/direction_coupling.f90
	elmerf90 $^ -o $@

$(EXEDIR)/CostSolver_Hpos: $(SRCDIR)/CostSolver_Hpos.F90
	elmerf90 $^ -o $@

$(EXEDIR)/CostSolver_Hobs: $(SRCDIR)/CostSolver_Hobs.f90
	elmerf90 $^ -o $@

$(EXEDIR)/CostSolver_Hmarine: $(SRCDIR)/CostSolver_Hmarine.F90
	elmerf90 $^ -o $@

$(EXEDIR)/POI_output: $(SRCDIR)/POI_output.f90
	elmerf90 $^ -o $@

$(EXEDIR)/POI_output_final: $(SRCDIR)/POI_output_final.f90
	elmerf90 $^ -o $@

$(EXEDIR)/POI_output_thickness: $(SRCDIR)/POI_output_thickness.f90
	elmerf90 $^ -o $@

$(EXEDIR)/NodalSlopeSolver : $(SRCDIR)/NodalSlopeSolver.F90
	elmerf90 $^ -o $@
	cp $(EXEDIR)/NodalSlopeSolver $(EXEDIR)/NodalSlopeSolver2
 
$(EXEDIR)/EffectiveMB : $(SRCDIR)/EffectiveMB.F90
	elmerf90 $^ -o $@

$(EXEDIR)/AdjointSSASolvers:
	make -f $(SRC_SSA)/Makefile
	mv AdjointSSASolvers.so $(EXEDIR)/AdjointSSASolvers.so

$(EXEDIR)/Wrap:$(SRCDIR)/Wrap.f90
	elmerf90 $^ -o $@

$(EXEDIR)/USFs: $(SRCDIR)/USFs.f90
	elmerf90 $^ -o $@

$(EXEDIR)/UpdateExport:$(SRCDIR)/UpdateExport.f90
	elmerf90 $^ -o $@

$(EXEDIR)/Optimize_m1qn3Serial: $(SRCDIR)/Optimize_m1qn3Serial.f90 $(SRCDIR)/m1qn3.f
	elmerf90 $^  -o $@ -lblas -L/home/hpc/gwgi/gwgi17/libraries/BLAS/

$(EXEDIR)/Optimize_m1qn3Parallel_Multi.so: $(SRCDIR)/Optimize_m1qn3Parallel.Multi.f90 $(SRCDIR)/m1qn3.f
	elmerf90 $(SRCDIR)/Optimize_m1qn3Parallel.Multi.f90 $(SRCDIR)/m1qn3.f -o $(EXEDIR)/Optimize_m1qn3Parallel_Multi.so 


clean:
	rm -f $(LIBS)
	rm -f objects/*
