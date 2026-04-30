CXX = clang++
CXXFLAGS = -Wall -O2 -Iinclude 

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

COMMON_SRC = myntcpfunctions.cpp
COMMON_OBJ = $(COMMON_SRC:.cpp=.o)

SRC_SINGLE = myntcpanalysis.cpp
OBJ_SINGLE = $(SRC_SINGLE:.cpp=.o)

TARGET_SINGLE = myntcpanalysis

all: $(TARGET_SINGLE)

$(TARGET_SINGLE): $(OBJ_SINGLE) $(COMMON_OBJ)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(OBJ_SINGLE) $(COMMON_OBJ) $(ROOTLIBS) -lz -o $(TARGET_SINGLE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_SINGLE)  $(COMMON_OBJ) $(TARGET_SINGLE)