#========== Definicion de variables ==========
# Actualizado el 21/09/2021
# Actualizado el 15/07/2021
# Actualizado el 3/11/2021
#	compiler
FC = ifort

#OBJ_DIR = ./obj	#FOR OBJECTS FILES
#BIN_DIR = ./bin  #FOR EXECUTABLE
#SRC_DIR = .

#	compiler flags
#	standard
#CFLAGS = -std=f2008ts
#	debugger option
CFLAGS += -g
#	warning flags
CFLAGS += -warn all 
#	optimization flags 
CFLAGS += -O0 -heap-arrays
#	error finding options
CFLAGS += -check all -traceback -fp-stack-check -CB #-fp-stack-check -check noarg_temp_created 
# at the end of the tests return to -check all option

#	source files
SRCS = mod_param mod_biunit mod_library main_CDR

OBJS = $(SRCS:=.o)

#	executable 
MAIN = CDR3d.x
#========== Fin variables ===========

#	compile project
all : $(MAIN)
#	@echo Compiling files . . . . .
#	@echo Making objects  . . . . . 
#	@echo Building an executable . . . . . 
	@echo ' '
	@echo '======================'
	@echo Compilation completed . . . . .
	@echo '======================'

$(MAIN) : $(OBJS)
	@$(FC) $(CFLAGS) -g -mkl -O0 $(OBJS) -o $(MAIN) 
#next line is for debugging (GDB does not work if the executable was compiled with the -mkl flag)
#@$(FC) $(CFLAGS) -g $(OBJS) -o $(MAIN) 
# @$(FC) $(CFLAGS) $(OBJS) -o $(MAIN)
#These three flags are needed to correct execution of the LAPACK library
#If the libray is comented (not used) in the code, must them desactivated the flegs

.SUFFIXES : .f90 .o
#.o.f90 :Dos opciones, cual sera la correcta?

.f90.o :
	@$(FC) $(CFLAGS) -c $<

#	Regla ficticia, es decir que no tiene dependencias (phony rules)
clean :
	@$(RM) *.o *.mod *.exe $(MAIN)
	@$(RM) -rf Res/*.txt
#@$(RM) Res/ *.txt
#	clean no tiene dependencias pero si objetivos
	@echo '* * * * * '
	@echo ' - Everything is clean -'
	@echo '* * * * * '
