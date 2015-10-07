# Makefile

# execution file name
EXE_POST = post
EXE_MAIN = main

# src / obj / lib / bin dir
SRC_DIR     =   src
OBJ_DIR     =   obj
BIN_DIR     =   bin
LIB_DIR     =   lib
INC_DIR     =   inc

ifeq ($(os),mac)
CGNS_INC    =   /usr/local/include
CGNS_LIB    =   /usr/local/lib
HDF_LIB     =   /usr/local/lib
else ifeq ($(os),tachyon2)
CGNS_INC    =   /applic/wa/cgnslib_3.2.1/intel-2015/include
CGNS_LIB    =   /applic/wa/cgnslib_3.2.1/intel-2015/lib
HDF_LIB     =   /applic/compilers/intel/2015/applib1/HDF5/1.8.13/lib
else
CGNS_INC    =   /usr/local/include
CGNS_LIB    =   /usr/local/lib
HDF_LIB     =   /usr/lib/x86_64-linux-gnu/hdf5/serial
endif

# f90 compiler
ifeq ($(mpi),ompi)
FC          =   mpif90
else
FC          =   mpiifort
endif

# flags for maximum performance
ifeq ($(compiler),gcc)
ifeq ($(mode),debug)
FFLAGS = -O0 -cpp -ffpe-trap=invalid -fcheck=all -Wall -pedantic -ffree-line-length-none -J $(OBJ_DIR)
else
FFLAGS = -O3 -cpp -ffpe-trap=invalid -ffree-line-length-none -J $(OBJ_DIR)
endif
else
ifeq ($(mode),debug)
FFLAGS = -O0 -fpp -debug full -warn all -ftrapuv -fp-model precise -traceback -check all -fpe0 -module $(OBJ_DIR)
else
FFLAGS = -O3 -fpp -fpe0 -fp-model precise -module $(OBJ_DIR)
endif
endif

# include needed for compile
IFLAGS = -I $(INC_DIR) -I $(CGNS_INC)

# library needed for linking
LFLAGS = -L $(LIB_DIR) -L $(CGNS_LIB) -lcgns -L $(HDF_LIB) -lhdf5 -lstdc++

# defin needed for preprocessor
DFLAGS = -D $(def)
DFLAGS_R = -D $(defr)

# whar archiving to use
AR = xiar rcs

# target obj list
OBJS_POST =   $(OBJ_DIR)/eos_module.o\
              $(OBJ_DIR)/prop_module.o\
              $(OBJ_DIR)/config_module.o\
              $(OBJ_DIR)/postgrid_module.o\
              $(OBJ_DIR)/postvariable_module.o\
              $(OBJ_DIR)/postprocess.o

OBJS_MAIN =   $(OBJ_DIR)/eos_module.o\
              $(OBJ_DIR)/prop_module.o\
              $(OBJ_DIR)/config_module.o\
              $(OBJ_DIR)/grid_module.o\
              $(OBJ_DIR)/variable_module.o\
              $(OBJ_DIR)/cav_module.o\
              $(OBJ_DIR)/turbsource_module.o\
              $(OBJ_DIR)/flux_module.o\
              $(OBJ_DIR)/muscl_module.o\
              $(OBJ_DIR)/unsteady_module.o\
              $(OBJ_DIR)/vsflux_module.o\
              $(OBJ_DIR)/bc_module.o\
              $(OBJ_DIR)/jacobian_module.o\
              $(OBJ_DIR)/lhs_module.o\
              $(OBJ_DIR)/timestep_module.o\
              $(OBJ_DIR)/residual_module.o\
              $(OBJ_DIR)/eddy_module.o\
              $(OBJ_DIR)/rhs_module.o\
              $(OBJ_DIR)/update_module.o\
              $(OBJ_DIR)/initial_module.o\
              $(OBJ_DIR)/solve_module.o\
              $(OBJ_DIR)/main.o

.PHONY: all clean distclean

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IFLAGS) $(DFLAGS) $(DFLAGS_R)

all: main post

main: $(OBJS_MAIN)
	$(FC) -o $(BIN_DIR)/$(EXE_MAIN) $(OBJS_MAIN) $(FFLAGS) $(LFLAGS)

post: $(OBJS_POST)
	$(FC) -o $(BIN_DIR)/$(EXE_POST) $(OBJS_POST) $(FFLAGS) $(LFLAGS)

clean:
	find $(OBJ_DIR) -name '*~' -exec rm {} \;
	find $(OBJ_DIR) -name '*.o' -exec rm {} \;
	find $(OBJ_DIR) -name '*.mod' -exec rm {} \;
	find $(OBJ_DIR) -name '*__genmod.f90' -exec rm {} \;

distclean: clean
	rm -f $(BIN_DIR)/$(EXE_POST) $(BIN_DIR)/$(EXE_MAIN)
