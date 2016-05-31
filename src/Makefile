############################
## Makefile Chenming       #
############################
## make ou make all : compiling
## make clean       : cleaning objects, modules et executables
## THE MERIT OF USING MAKEFILE IS
## 1. many complex make commands can be convinently put into one file
## 2. one does not need to recompile the object file of which the source code
## didn't change.
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
OPTIONDB=-ggdb
## executable file
EXECUTABLE=sutraset_gf
DEBUGEXE  =sutra_debug_gf
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
	perfilm.f \
	ssubs_2_2.f \
	usubs_2_2.f \
	ft03.f90\
	SinkareaRegular.f\
	others.f\
	evaporation.f90\
	surfrsis.f90\
	$(GV)    \
	sutra_2_2.f 
## object file from source.f90 (.f90 -> .o)               
OBJECTS_1=$(SOURCES:.f90=.o)
OBJECTS=$(OBJECTS_1:.f=.o)

#OBJECTS=$(SOURCES:.f=.o)
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
#	@echo $(OBJECTS) >aaa.txt 

#target: prerequisites 


## make all
all: $(EXECUTABLE) 
#all:
## 1) make: links and executable creation
## this line has to start without indentation
$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(OPTIONS) $^ -o $@ $(LIBS)
## 2) compiling separated objects
# '%.f' is all files in '$(OBJECTS)' with exetension name as f 
# '%.o'?????????????????
%.o: %.f
	$(COMPILER) -c $(OPTIONS) $^ 
%.o: %.f90
	$(COMPILER) -c $(OPTIONS) $^ 
#	@echo %.o >b.txt
# question: why it seems that the program first executes process 
#    2) --i.e., make objects then process 1) --i.e., create 
#    binary files?


prerequisites: do_script
## clean files
clean:
	@rm -f *.o *.mod *~ $(EXECUTABLE) $(DEBUGEXE)
	@rm $(GV)
	@echo " "
	@echo "cleaning OK."
	@echo "-------------"
# debug files
# 'debug: $(DEBUGEXE)' 
# 'debug'      -- the name that this subroutine will be called 
#                 for example 'make debug'
# '$(DEBUGEXE)'-- the binary file created after this subroutine 
#
debug: $(DEBUGEXE)
# '${DEBUGEXE}:${SOURCES}'
# '${DEBUGEXE}'  -- the binary that is going to be created. it 
#                 is referred to as '$@' in the following line
# '${SOURCES}'   -- the source file that is going to be used
#                 it is referred to as '$^' in the following line
${DEBUGEXE}: ${SOURCES}
	$(COMPILER) $(OPTIONDB) $^ -o $@
