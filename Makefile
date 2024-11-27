CXX = g++

CPPFLAGS = -g -std=c++11 -O3 -ffast-math -DTIMING 

TARGET = rabbitbam

DIR_SRC := .
DIR_OBJ := .

LIB =
INCLUDE =
HTSLIB_INSTALL_PATH = ./ext/htslib-1.20
LIBDEFLATE_INSTALL_PATH = ./ext/libdeflate-1.20

LIB += $(HTSLIB_INSTALL_PATH)
LIB += $(LIBDEFLATE_INSTALL_PATH)/build

INCLUDE += $(HTSLIB_INSTALL_PATH)
INCLUDE += $(LIBDEFLATE_INSTALL_PATH)


INC_FLAGS = $(addprefix -I,$(INCLUDE)) 
LIB_FLAGS = $(addprefix -L,$(LIB)) -fopenmp -lpthread -lz -lhts -ldeflate

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))


${TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LIB_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) -c $< -o $@


.PHONY: clean

clean:
	rm -rf $(DIR_OBJ)/*.o
	rm -rf $(HTSLIB_INSTALL_PATH)/*.o
	rm -rf $(LIBDEFLATE_INSTALL_PATH)/build
	rm -rf $(TARGET)

