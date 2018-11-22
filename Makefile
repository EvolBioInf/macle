CXX ?= g++
BUILD64BIT ?= 1
USE_SDSL ?= 0
LOCAL_LIBDIVSUFSORT ?= 1
PARALLEL_DIVSUFSORT ?= 0

CXXFLAGS := -std=c++11 -Isrc -Wall -Wextra -O3 -g -ggdb -Wshadow # -pg
LDFLAGS := -lm -ldivsufsort
ifeq ($(PARALLEL_DIVSUFSORT), 0)
LDFLAGS += -ldivsufsort64
else
CXXFLAGS += -DPARALLEL
endif

TARGET := macle

SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

TEST_SRC=$(wildcard tests/*.cpp)
TEST_OBJ=$(TEST_SRC:.cpp=.o)
TESTS=$(patsubst tests/%.cpp,build/%,$(TEST_SRC))

###############################################################################
#### Add configuration dependent compiler flags

ifeq ($(BUILD64BIT), 1)
  CXXFLAGS += -DU64
endif
ifeq ($(USE_SDSL), 1)
  CXXFLAGS += -DUSE_SDSL -Isdsl/include -msse4.2
  LDFLAGS += -lsdsl -Lsdsl/lib
endif
ifeq ($(LOCAL_LIBDIVSUFSORT), 1)
ifeq ($(PARALLEL_DIVSUFSORT), 1)
  CXXFLAGS += -Iparallel-divsufsort/include -fopenmp
  LDFLAGS += -Lparallel-divsufsort/lib -Lparallel-divsufsort/external/libprange/lib -lgomp -llibprange
else
  CXXFLAGS += -Ilibdivsufsort/include
  LDFLAGS += -Llibdivsufsort/lib
endif
endif

###############################################################################
#### Build Rules

all: build tests

build: build/$(TARGET) $(TESTS)

build/$(TARGET): $(OBJECTS)
	mkdir -p build
	$(CXX) $^ -o build/$(TARGET) $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

tests/%.o: tests/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(TESTS): $(OBJECTS) $(TEST_OBJ)
	$(CXX) $(filter-out src/$(TARGET).o,$(OBJECTS)) $(@:build/%=tests/%).o -o $@ $(LDFLAGS)

clean:
	$(RM) -r build
	$(RM) tests/*.o src/*.o

tests: $(TESTS)
	bash ./tests/runtests.sh

valgrind:
	VALGRIND="valgrind --leak-check=full" $(MAKE)

format:
	clang-format -i src/*.cpp src/*.h tests/*.cpp tests/*.h

divsufsort:
	-git clone https://github.com/y-256/libdivsufsort.git
	cd libdivsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 && make

parallel-divsufsort:
	-git clone https://github.com/jlabeit/parallel-divsufsort.git
	sed -i '/(examples)/d' parallel-divsufsort/CMakeLists.txt
	-cd parallel-divsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 -DOPENMP=1 -DCILKP=0 && make
	sed -i 's/cilkrts//' parallel-divsufsort/external/libprange/demo/CMakeLists.txt
	cd parallel-divsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 -DOPENMP=1 -DCILKP=0 && make

sdsl:
	-git clone https://github.com/simongog/sdsl-lite.git
	sdsl-lite/install.sh sdsl

# for neomake/syntastic in vim
show_cxxflags:
	@echo $(CXXFLAGS)

.PHONY: all build clean tests valgrind format divsufsort parallel-divsufsort sdsl show_cxxflags
