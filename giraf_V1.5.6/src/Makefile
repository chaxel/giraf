.IGNORE:

.SUFFIXES: .o .f90

target: giraf.exe

VERSION=1.56

OBJS=	params.o getdimensions.o giraf.o check.o  wcal2d_proj.o tools.o alloc.o fonctions.o chimie.o physique.o io_netcdf.o io_ascii.o io.o proj_xy.o

FC=$(F90)

LIBS=-L$(NETCDFHOME)/lib -lnetcdf -lnetcdff -L$(PROJECTIONLIB) -lproj

FCFLAGS=-O2 -I$(NETCDFHOME)/include $(FCOPT)

RM=rm

RM_LIST = $(OBJS) *.mod *.exe 

######################################################################

giraf.exe: $(OBJS)
	($(FC) -o $@ $(OBJS) $(FCFLAGS)  $(LIBS))

#-----------------------------------------------------------------------------
# Don't touch anything below this line
#-----------------------------------------------------------------------------

params.o:
	$(FC) -o params.o -c $(FCFLAGS) params.f90
	$(FC) -c $(FCFLAGS) params.f90	

.f90.o:
#	$(RM) $@
	$(FC) -c $(FCFLAGS) $*.f90

clean:
	$(RM) $(RM_LIST)
