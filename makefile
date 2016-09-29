OPT = -O3 -std=c++11
SRC_FOLDER = src
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

EXES = $(addsuffix .exe, $(basename $(FILES)))
DEPS = $(addsuffix .o, $(basename $(DEPFILES)))

all : $(DEPS)
all : $(EXES)

debug : OPT = -O0 -g -std=c++0x 
debug : $(DEPS)
debug : $(EXES)

%.o: $(SRC_FOLDER)/%.c*
	$(CXX) $(OPT) -c -o $@ $<

%.exe: $(SRC_FOLDER)/%.cpp
	$(CXX) $(OPT) -I$(SRC_FOLDER) -o $@ $(DEPS) $< $(LIB) 

clean:
	rm -f $(DEPS)

cleanall: clean
	rm -f $(EXES) 


