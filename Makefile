# weather-buoys Makefile
FC = ifort
# FC = ifx
SRCDIR = SRC
OBJDIR = OBJ
OBJS = $(addprefix $(OBJDIR)/, mod_param.o mod_geometry.o mod_inputInfo.o mod_biunit.o mod_BVs.o mod_source.o mod_library.o mod_exacSol.o mod_timeInt.o)


#	Para eliminar mensaje para actualizar compilador de ifort a ifx
CFLAGS = -diag-disable=10448
#	compiler flags
CFLAGS += -stand f08#f03 #f90 #f08
#	debugger option
CFLAGS += -g -debug all
#	warning flags
CFLAGS += -warn all
#	optimization flags
CFLAGS += -O0 -heap-arrays
#	error finding options
CFLAGS +=  -check all -traceback -mcmodel large -fp-stack-check -check noarg_temp_created
#	mkl library
CFLAGS += -qmkl
  

all: main_CDR.exe #weather_stats_parallel
#	@echo Compiling files . . . . .
#	@echo Making objects  . . . . .
#	@echo Building an executable . . . . .
	@echo ' '
	@echo '======================'
	@echo Compilation completed . . . . .
	@echo '======================'

.SUFFIXES: .f90 .o

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@$(FC) $(CFLAGS) -c $< -o $@ -module $(OBJDIR)

main_CDR.exe: $(SRCDIR)/main_CDR.f90 $(OBJS)
	@$(FC) $(CFLAGS) $< $(OBJS) -o $@ -module $(OBJDIR)

%.o: %.mod

$(OBJS): | $(OBJDIR)


$(OBJDIR):
	mkdir -p $(OBJDIR)


.PHONY: run

run: main_CDR.exe
	@./main_CDR.exe


.PHONY: clean

clean:
	@$(RM) -r $(OBJDIR)/*.o $(OBJDIR)/*.mod *.dat main_CDR.exe
	@echo ' '
	@echo '* * * * * '
	@echo ' - Everything is clean -'
	@echo '* * * * * '
