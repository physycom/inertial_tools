OPT = -O3 -std=c++11
SRC_FOLDER = src
OBJ_FOLDER = obj
EXE_FOLDER = bin
FILES = data_calibrator.cpp \
data_filter.cpp \
data_gen.cpp \
data_purge.cpp \
data_reco.cpp \
data_rotator.cpp \
inertial_stats.cpp \
mutual_entropy.cpp \
speed_compare.cpp

DEPFILES = math_func.cpp \
math_lib.c

SRC = $(addprefix $(SRC_FOLDER)/, $(FILES))
DEP_SRC = $(addprefix $(SRC_FOLDER)/, $(DEPFILES))

EXES = $(addprefix $(EXE_FOLDER)/,$(addsuffix .exe, $(basename $(FILES))))
DEPS = $(addprefix $(OBJ_FOLDER)/,$(addsuffix .o, $(basename $(DEPFILES))))

all : foldertree
all : $(DEPS)
all : $(EXES)

openmp : OPT += -fopenmp 
openmp : $(DEPS)
openmp : $(EXES)

debug : OPT = -O0 -g -std=c++0x 
debug : $(DEPS)
debug : $(EXES)

foldertree:
	mkdir -p bin
	mkdir -p obj

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.c*
	$(CXX) $(OPT) -c -o $@ $<

$(EXE_FOLDER)/%.exe: $(SRC_FOLDER)/%.cpp
	$(CXX) $(OPT) -I$(SRC_FOLDER)/jsoncons/src -o $@ $(DEPS) $< $(LIB) 

clean:
	rm -f $(DEPS)

cleanall: clean
	rm -f $(EXES) 


