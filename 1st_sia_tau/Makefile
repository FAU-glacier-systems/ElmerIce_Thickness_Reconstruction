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
     $(EXEDIR)/POI_output_final \
     $(EXEDIR)/POI_output_thickness \
     $(EXEDIR)/USFs\
     $(EXEDIR)/Wrap\
     $(EXEDIR)/EffectiveMB \
     $(EXEDIR)/NodalSlopeSolver \
     $(EXEDIR)/CostSolver_FLUXpos \
     $(EXEDIR)/computeSIA_deltathickness \
     $(EXEDIR)/computeSIA_deltaflux \
     $(EXEDIR)/computeTAU_Hobs \
     $(EXEDIR)/computeTAU_thickness \
     ${EXEDIR}/computeSIA_Hobs \
     ${EXEDIR}/computeSIA_thickness \
     $(EXEDIR)/computeSLOPE_drivingSTRESS \
     $(EXEDIR)/SMOOTHdrivingSTRESS_sia \
     $(EXEDIR)/SMOOTHdrivingSTRESS_tau 


all: $(LIBS) 


$(EXEDIR)/AdjointThicknessSolvers:
	make -f $(SRC_H)/Makefile
	mv AdjointThicknessSolvers.so $(EXEDIR)/AdjointThicknessSolvers.so

$(EXEDIR)/computeSLOPE_drivingSTRESS: $(SRCDIR)/computeSLOPE_drivingSTRESS.F90
	elmerf90 $^ -o $@

$(EXEDIR)/computeSIA_deltathickness: $(SRCDIR)/computeSIA_deltathickness.F90
	elmerf90 $^ -o $@

$(EXEDIR)/computeSIA_deltaflux: $(SRCDIR)/computeSIA_deltaflux.F90
	elmerf90 $^ -o $@

$(EXEDIR)/computeTAU_Hobs: $(SRCDIR)/computeTAU_Hobs.F90
	elmerf90 $^ -o $@
	cp $(EXEDIR)/computeTAU_Hobs $(EXEDIR)/computeTAU_Hobs_init

$(EXEDIR)/computeTAU_thickness: $(SRCDIR)/computeTAU_thickness.F90
	elmerf90 $^ -o $@
	cp $(EXEDIR)/computeTAU_thickness $(EXEDIR)/computeTAU_thickness_init

$(EXEDIR)/computeSIA_Hobs: $(SRCDIR)/computeSIA_Hobs.F90
	elmerf90 $^ -o $@
	cp $(EXEDIR)/computeSIA_Hobs $(EXEDIR)/computeSIA_Hobs_init

$(EXEDIR)/computeSIA_thickness: $(SRCDIR)/computeSIA_thickness.F90
	elmerf90 $^ -o $@
	cp $(EXEDIR)/computeSIA_thickness $(EXEDIR)/computeSIA_thickness_init

$(EXEDIR)/CostSolver_FLUXpos: $(SRCDIR)/CostSolver_FLUXpos.F90
	elmerf90 $^ -o $@

$(EXEDIR)/POI_output_thickness: $(SRCDIR)/POI_output_thickness.f90
	elmerf90 $^ -o $@

$(EXEDIR)/POI_output_final: $(SRCDIR)/POI_output_final.f90
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

$(EXEDIR)/SMOOTHdrivingSTRESS_sia:$(SRCDIR)/SMOOTHdrivingSTRESS_sia.f90
	elmerf90 $^ -o $@

$(EXEDIR)/SMOOTHdrivingSTRESS_tau:$(SRCDIR)/SMOOTHdrivingSTRESS_tau.f90
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
