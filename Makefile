
# Specifying some environmental variables for Linux, note this can be done in the shell prompt

COMP	= GCC
TECIO	= NO
CODE	= DEBUG
OPENMP	= NO

# Specifing Standard Variables:
CXX	= g++ -std=gnu++11 #-pedantic-errors # c++ gcc compiler
CXXFLAGS=       # C++ compiler flags
LDLFLAGS=	# linker flags
CPPFLAGS=	# c/c++ preprocessor flags

OPTS	= 	# optimization flags and other options

# Includes

OPTS	+= -I include
#OPTS    += -I /home/mhawwary/Libraries/Eigen3.3.2/Eigen
#OPTS += -I /home/mhawwary/work/hpMusic/contrib/eigen/Eigen

ifeq ($(TECIO),YES)
	OPTS += -I $(TECIO_DIR)/tecsrc
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif

ifeq ($(OPENMP),YES)	
	OPTS	+=  -lgomp -fopenmp 
endif

ifeq ($(CODE),DEBUG)
	ifeq ($(COMP),GCC)
		OPTS	+= -fbuiltin -g -Wall #-Werror
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif


# Source

SRC	= src/
OBJ	= obj/
BIN	= bin/
INC	= include/      

vpath %.cpp src
vpath %.c src
vpath %.o   obj
vpath %.h include src
vpath %.hpp include src

# Objects
OBJS	= $(OBJ)FD1D.o  # objects 
INCLS	= 

# Compile

.PHONY: default help clean


default: all
help:	
	@echo 'help'

all: FD1DFlow.exe

FD1DFlow.exe: $(OBJS)
	$(CXX) $(OPTS) -o $(BIN)$@ $+


$(OBJ)%.o : %.cpp 
	$(CXX) $(OPTS) -c -o $@ $<


$(OBJ)FD1DFlow.o:   FD1D.cpp 
#$(OBJ)SimData.o:   SimData.hpp SimData.cpp
#$(OBJ)general_tools.o:   general_tools.h general_tools.cpp
#$(OBJ)GridData.o:   GridData.hpp GridData.cpp

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe 
	#rm -rf ./output/ 
	@echo  removing all object and executable files
