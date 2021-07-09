#MAI USATI
VargenoUML_EXEC ?= genolight
GBF_EXEC ?= gbf

LIBS = -lm

#Non possiamo modificare CC con gcc. Questo perchè gcc è il compilatore per c. Noi abbiamo sia c(bwa) che cpp(vargeno).
#La mia soluzione (provvisoria?) è stata quella di cambiare estensione a tutti i file .c di bwa che generavano 
#eccezzioni compilandoni con g++.
CC = g++

#SENZA -fpermissive da errore a bwt

WARNINGS = -fpermissive
CFLAGS = -std=c++14 -march=native -O3 -flto -fstrict-aliasing $(WARNINGS) #ERA 11 in origine
LFLAGS = -march=native -O3 -flto

SRC_DIR = src
BUILD_DIR = build
INCDIR = .

SRCS := $(shell find $(SRC_DIR) ! -samefile "src/lib_aln_inexact_matching/src/example.c" \( -name *.cc -or -name *.c -or -name *.cpp \) )
$(info $$SRCS is [${SRCS}])
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
INC_DIRS := $(shell find $(SRC_DIRS) -type d)

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

CXX=g++ 
SDSL_PATH=$(PREFIX)/lib
CXXFLAGS= -Wall -g -std=c++11 -O3 -I$(PREFIX)/include

LD_LIB=-L $(SDSL_PATH)
LD_FLAG= -lpthread -lz -lrt 					    
.PHONY: all clean
all: genolight

HEADERS = $(wildcard $(INCDIR)/*.h)

$(BUILD_DIR)/src/%.c.o: $(SRC_DIR)/%.c $(HEADERS) src/*.h
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -fpermissive $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ #IN ORIGINE ERA CXX


$(BUILD_DIR)/src/%.cc.o: $(SRC_DIR)/%.cc 
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) -fpermissive -c -o $@ $^ $(CXXFLAGS) -I$(INCDIR)

$(BUILD_DIR)/src/%.cpp.o: $(SRC_DIR)/%.cpp 
	$(MKDIR_P) $(dir $@)
	$(CXX) -std=c++14 $(CPPFLAGS) -c $< -o $@ 

.PRECIOUS: $(TARGET) $(OBJECTS)
genolight: $(OBJS)
	$(CC) $(LFLAGS) $(CXXFLAGS) $(OBJS) $(LIBS) -std=c++11 -lstdc++fs -Wodr -o $@ $(LD_LIB) $(LD_FLAG) 

clean:
	-rm -f $(BUILD_DIR)

MKDIR_P ?= mkdir -p

