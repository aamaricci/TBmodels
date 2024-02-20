##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu


##$ SET THE TARGET DIRECTORY WHERE TO PUT THE EXECUTABLE (default if $HOME/.bin in the PATH)
DIREXE=$(HOME)/.bin


##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
EXE=hm_2d


##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin


##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)
#-ldmft_tools -lscifor

ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp -D_MPI
endif 

ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime
#-O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPFLAG =-cpp -D_MPI
endif




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90


all:
	@echo ""
	@echo "compiling $(EXE).f90 .... "
	@echo ""
	$(FC) $(FFLAG) $(FPPFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCARGS) $(LIBARGS)
	@echo "Done"

debug:
	@echo ""
	@echo "compiling $(EXE).f90 .... "
	@echo ""
	$(FC) $(DFLAG) $(FPPFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCARGS) $(LIBARGS)
	@echo "Done"
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~

version:
	@echo $(VER)


#########################################################################
