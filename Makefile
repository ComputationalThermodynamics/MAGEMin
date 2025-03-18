# CC=gcc
CC=clang 

UNAME_S := $(shell uname -s)

EXE_NAME:= MAGEMin

# Check if USE_MPI is set (default: 1)
USE_MPI ?= 1

CCFLAGS = -Wall -O3 -g -fPIC -Wno-unused-variable -Wno-unused-but-set-variable -march=native -funroll-loops -flto
ifeq ($(UNAME_S),Darwin)
	INC      = -I/opt/homebrew/include 
	LIBS     = -lm -framework Accelerate /opt/homebrew/lib/libnlopt.dylib
	ifeq ($(USE_MPI),1)
		CCFLAGS += -DUSE_MPI
		LIBS    += /opt/homebrew/lib/libmpi.dylib
	endif
endif
ifeq ($(UNAME_S),Linux)
	LIBS     = -lm -llapacke -lnlopt -g -L/usr/lib 
	ifeq ($(USE_MPI),1)
		CCFLAGS += -DUSE_MPI
		LIBS    += -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
		INC      = -I/usr/lib/x86_64-linux-gnu/openmpi/include/
	endif
endif
	EXE_NAME = MAGEMin

SOURCES=src/MAGEMin.c 							\
		src/initialize.c 						\
		src/TC_database/TC_init_database.c		\
		src/TC_database/TC_endmembers.c			\
		src/TC_database/TC_gem_function.c		\
		src/SB_database/SB_init_database.c		\
		src/SB_database/SB_endmembers.c			\
		src/SB_database/SB_gem_function.c		\
		src/toolkit.c							\
		src/io_function.c						\
		src/gem_function.c 						\
		src/TC_database/tc_gss_init_function.c	\
		src/TC_database/tc_gss_function.c		\
		src/SB_database/sb_gss_init_function.c	\
		src/SB_database/sb_gss_function.c		\
		src/TC_database/NLopt_opt_function.c 	\
		src/SB_database/SB_NLopt_opt_function.c \
		src/TC_database/objective_functions.c	\
		src/SB_database/sb_objective_functions.c\
		src/TC_database/SS_xeos_PC_mp.c			\
		src/TC_database/SS_xeos_PC_mb.c			\
		src/TC_database/SS_xeos_PC_ig.c			\
		src/TC_database/SS_xeos_PC_igad.c		\
		src/TC_database/SS_xeos_PC_um.c			\
		src/TC_database/SS_xeos_PC_mtl.c		\
		src/TC_database/SS_xeos_PC_mpe.c		\
		src/SB_database/SS_xeos_PC_sb11.c		\
		src/SB_database/SS_xeos_PC_sb21.c		\
		src/pp_min_function.c 					\
		src/ss_min_function.c 					\
		src/simplex_levelling.c 				\
		src/PGE_function.c 						\
		src/phase_update_function.c				\
		src/dump_function.c

OBJECTS=$(SOURCES:.c=.o)



.c.o:
	$(CC) $(CCFLAGS) -c $< -o $@ $(INC)
 
all: $(OBJECTS)
	$(CC)  -o $(EXE_NAME) $(OBJECTS) $(INC) $(LIBS)  -flto
	rm src/*.o src/TC_database/*.o

lib: $(OBJECTS)
	$(CC) -shared -fPIC -o libMAGEMin.dylib $(OBJECTS) $(INC) $(LIBS) -flto
 
clean:
	rm -f src/*.o  src/TC_database/*.o src/SB_database/*.o *.dylib MAGEMin
