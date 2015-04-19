############################
## Makefile Chenming       #
############################
## make ou make all : compiling
## make clean       : cleaning objects, modules et executables

## Create gversion file for version control
GV=gversion.f90
##echo "#define GIT_REF \"`git show-ref refs/heads/master | cut -d " " -f 1 | cut -c 31-40`\"" >
## Compiler 

#COMPILER=g95_64
#COMPILER=gfortran -fopenmp
#COMPILER=ifort
#COMPILER=ifort -openmp -parallel
COMPILER=gfortran

## compiling options

#OPTIONS= -O2 -IPF-flt-eval-method0
#OPTIONS= -O2 -g -traceback -fp-stack-check
#OPTIONS= -O2 -ffree-line-length-none
#OPTIONS= -O2 -openmp -parallel -fpp
#OPTIONS= -O0 -ftrace=full -fbounds-check -ffree-line-length-none
#OPTIONS= -O0 -fbounds-check -ffree-line-length-none
#OPTIONS= -fbacktrace 
#OPTIONS= -ffpe-trap=zero
OPTIONS= -O3  # -g -traceback -check bounds

## executable file
EXECUTABLE=sutraset_gf

## Sources files
#SOURCES=subroutine1.f90           \
#        subroutine2.f90           \
#        subroutine3.f90           \
#        subroutine3.f90	      \
#        main.f90
# 1. $(gv) has to have the bracket
# 2. no space or anything else after the slash file !!!!!!!!!!!!!!
# 3. it is ok to use $(GV) in SOURCES
SOURCES=    fmods_2_2.f \
	indatet.f90 \
	unsat.f \
	ssubs_2_2.f \
	usubs_2_2.f \
	ft03.f90\
	SinkareaRegular.f\
	Others.f\
	surfrsis.f90\
	$(GV)    \
	sutra_2_2.f 
## object file from source.f90 (.f90 -> .o)               
OBJECTS_1=$(SOURCES:.f90=.o)
OBJECTS=$(OBJECTS_1:.f=.o)

OBJECTS=$(SOURCES:.f=.o)
## Libraries
#LIBS = -L/usr/lib -llapack -lblas
#LIBS = -L/usr/lib -llapack.so.3 -lblas.so.3
#LIBS = /usr/lib64/atlas/liblapack.so.3 /usr/lib64/libblas.so.3
LIBS  = 

#target: prerequisites


.PHONY: do_script

do_script: 
	@echo 'SUBROUTINE GVERSION (K3)' >$(GV) 
	@echo "  WRITE (K3,*) \" GIT VERSION: `git show-ref refs/heads/master | cut -d " " -f 1 `\" ">> $(GV)
	@echo '      RETURN '  >>$(GV)
	@echo 'END SUBROUTINE' >>$(GV)


#target: prerequisites 


## make all
all: $(EXECUTABLE) 
#all:
## 1) make: links and executable creation
## this line has to start without indentation
$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(OPTIONS) $^ -o $@ $(LIBS)
## 2) compiling separated objects
%.o: %.f
	$(COMPILER) -c $(OPTIONS) $^ 

prerequisites: do_script
## clean files
clean:
	@rm -f *.o *.mod *~ $(EXECUTABLE)
	@rm $(GV)
	@echo " "
	@echo "cleaning OK."
	@echo "-------------"

