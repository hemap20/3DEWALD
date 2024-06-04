CC = /usr/bin/g++	
DEBUGFLAGS = -Wall
OPTFLAGS = -O3 -fopenmp -pthread
FFTFLAGS = -lfftw3_threads -lfftw3 -lm -lfftw3_omp

RUN_DIR=./run
RAT_OUTPUT=$(RUN_DIR)/coulomb.x

# Include folders

INC_LIST= -I ./inc /usr/include/
	#   -I/home/prateek/eigen3/

# Source Folders
SRC_DIR=./src

#For if ./obj is not present: Types of Prerequisites (https://www.gnu.org/software/make/manual/html_node/Prerequisite-Types.html) 
# Object folders
OBJ_DIR=./obj

# Library folders
LIB_DIR=./lib

# Object list
OBJ_FILES=$(OBJ_DIR)/main.o \
	  $(OBJ_DIR)/dSFMT.o \
	  $(OBJ_DIR)/print.o \
	  $(OBJ_DIR)/self_e.o \
	  $(OBJ_DIR)/real_e.o \
	  $(OBJ_DIR)/dist.o \
	  $(OBJ_DIR)/reciprocal.o \
	  $(OBJ_DIR)/bspline_reci.o \
	  $(OBJ_DIR)/error.o

# Make Targets
all:$(OBJ_FILES) output

output:$(RAT_OUTPUT)

# Build object files
$(OBJ_DIR)/main.o:$(SRC_DIR)/main.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o $(OBJ_DIR)/main.o $(INC_LIST)
$(OBJ_DIR)/dSFMT.o:$(SRC_DIR)/dSFMT.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/dSFMT.o $(INC_LIST)
$(OBJ_DIR)/print.o:$(SRC_DIR)/print.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/print.o $(INC_LIST)
$(OBJ_DIR)/self_e.o:$(SRC_DIR)/self_e.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/self_e.o $(INC_LIST)
$(OBJ_DIR)/real_e.o:$(SRC_DIR)/real_e.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/real_e.o $(INC_LIST)
$(OBJ_DIR)/dist.o:$(SRC_DIR)/dist.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/dist.o $(INC_LIST)
$(OBJ_DIR)/reciprocal.o:$(SRC_DIR)/reciprocal.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reciprocal.o $(INC_LIST)
$(OBJ_DIR)/bspline_reci.o:$(SRC_DIR)/bspline_reci.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/bspline_reci.o $(INC_LIST)
$(OBJ_DIR)/error.o:$(SRC_DIR)/error.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/error.o $(INC_LIST)

$(RAT_OUTPUT):$(OBJ_FILES)
	$(CC) $(OPTFLAGS) $(INC_LIST) -o $(RAT_OUTPUT) $(OBJ_FILES) $(FFTFLAGS)

# Clean objects and library
clean:
	$(RM) $(OBJ_FILES)
