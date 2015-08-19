# Makefile

# execution file name
EXE_POST = post
EXE_MAIN = main

# src / obj / lib / bin dir
SRC_DIR     =   src
OBJ_DIR     =   obj
LIB_DIR     =   /usr/local/lib
HDF_LIB     =   /usr/lib/x86_64-linux-gnu/hdf5/serial
HDF_INC     =   /usr/include/hdf5/serial
BIN_DIR     =   bin
INC_DIR     =   /usr/local/include

# f90 compiler
ifeq ($(mpi),open)
FC          =   mpif90
else
FC          =   mpiifort
endif

# flags for maximum performance
ifeq ($(compiler),gcc)
ifeq ($(mode),debug)
FFLAGS = -O0 -fbounds-check -Wall -pedantic -ffree-line-length-none $(OBJ_DIR)
else
FFLAGS = -O3 -ffree-line-length-none $(OBJ_DIR)
endif
else
ifeq ($(mode),debug)
FFLAGS = -O0 -fpp -debug full -warn all -ftrapuv -fp-model strict -traceback -check all -fpe0 -module $(OBJ_DIR)
else
FFLAGS = -O3 -fpp -fpe0 -fp-model strict -module $(OBJ_DIR)
endif
endif

# include needed for compile
IFLAGS = -I $(INC_DIR) 

# library needed for linking
LFLAGS = -L $(LIB_DIR) -lcgns -L $(HDF_LIB) -lhdf5_serial -lstdc++

# defin needed for preprocessor
DFLAGS = -D $(def)

# whar archiving to use
AR = xiar rcs

# target obj list
OBJS_POST =   $(OBJ_DIR)/eos_module.o\
              $(OBJ_DIR)/prop_module.o\
              $(OBJ_DIR)/config_module.o\
              $(OBJ_DIR)/postgrid_module.o\
              $(OBJ_DIR)/postvariable_module.o\
              $(OBJ_DIR)/datawriting_module.o\
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
              $(OBJ_DIR)/post_module.o\
              $(OBJ_DIR)/residual_module.o\
              $(OBJ_DIR)/eddy_module.o\
              $(OBJ_DIR)/rhs_module.o\
              $(OBJ_DIR)/update_module.o\
              $(OBJ_DIR)/initial_module.o\
              $(OBJ_DIR)/solve_module.o\
              $(OBJ_DIR)/main.o

.PHONY: all clean distclean

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IFLAGS) $(DFLAGS)

post: $(OBJS_POST) 
	$(FC) -o $(BIN_DIR)/$(EXE_POST) $(OBJS_POST) $(FFLAGS) $(LFLAGS)

main: $(OBJS_MAIN)
	$(FC) -o $(BIN_DIR)/$(EXE_MAIN) $(OBJS_MAIN) $(FFLAGS) $(LFLAGS)

clean:
	find $(OBJ_DIR) -name '*~' -exec rm {} \;
	find $(OBJ_DIR) -name '*.o' -exec rm {} \;
	find $(OBJ_DIR) -name '*.mod' -exec rm {} \;
	find $(OBJ_DIR) -name '*__genmod.f90' -exec rm {} \;

distclean: clean
	rm -f $(BIN_DIR)/$(EXE_POST) $(BIN_DIR)/$(EXE_MAIN)
