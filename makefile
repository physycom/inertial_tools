SRC_FOLDER = src
OPT_CXX = -O3 -std=c++11
OPT_CC = -O3
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

DEPFILES = math_func.cpp

LIBRARIES = bbmutils.c

SRC = $(addprefix $(SRC_FOLDER)/, $(FILES))
DEP_SRC = $(addprefix $(SRC_FOLDER)/, $(DEPFILES))

EXES = $(addprefix $(EXE_FOLDER)/,$(addsuffix .exe, $(basename $(FILES))))
DEPS = $(addprefix $(OBJ_FOLDER)/,$(addsuffix .o, $(basename $(DEPFILES))))
LIBS = $(addprefix $(OBJ_FOLDER)/,$(addsuffix .o, $(basename $(LIBRARIES))))

all : foldertree
all : $(LIBS)
all : $(DEPS)
all : $(EXES)

openmp : foldertree
openmp : OPT_CXX += -fopenmp 
openmp : OPT_CC += -fopenmp 
openmp : $(LIBS)
openmp : $(DEPS)
openmp : $(EXES)

debug : foldertree
debug : OPT_CXX = -O0 -g -std=c++0x 
debug : OPT_CC = -O0 -g
debug : $(LIBS)
debug : $(DEPS)
debug : $(EXES)

foldertree:
	mkdir -p bin
	mkdir -p obj

$(OBJ_FOLDER)/math_func.o: $(SRC_FOLDER)/math_func.cpp
	$(CXX) $(OPT_CXX) -c -o $@ $<

$(OBJ_FOLDER)/bbmutils.o: $(SRC_FOLDER)/libbbmutils/src/bbmutils.c
	$(CC) $(OPT_CC) -I$(SRC_FOLDER) -c -o $@ $<

$(EXE_FOLDER)/%.exe: $(SRC_FOLDER)/%.cpp
	$(CXX) $(OPT_CXX) -I$(SRC_FOLDER)/jsoncons/src -o $@ $(DEPS) $(LIBS) $< 

clean:
	rm -f $(DEPS)

cleanall: 
	rm -rf $(EXE_FOLDER) $(OBJ_FOLDER)


