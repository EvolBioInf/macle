CXX ?= g++
USE_SDSL ?= 1
LOCAL_LIBDIVSUFSORT ?= 1

CXXFLAGS := -std=c++11 -Isrc -Wall -Wextra -O3 -g -ggdb -Wshadow # -pg
LDFLAGS := -lm -lgsl -lblas -ldivsufsort -ldivsufsort64 # -pg
TARGET := dnalc

SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

TEST_SRC=$(wildcard tests/*.cpp)
TEST_OBJ=$(TEST_SRC:.cpp=.o)
TESTS=$(patsubst tests/%.cpp,build/%,$(TEST_SRC))

###############################################################################
#### Add configuration dependent compiler flags

ifeq ($(USE_SDSL), 1)
  CXXFLAGS += -DUSE_SDSL -Isdsl/include -msse4.2
  LDFLAGS += -lsdsl -Lsdsl/lib
endif
ifeq ($(LOCAL_LIBDIVSUFSORT), 1)
  CXXFLAGS += -Ilibdivsufsort/include
  LDFLAGS += -Llibdivsufsort/lib
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

libdivsufsort:
	git clone https://github.com/y-256/libdivsufsort.git
	cd libdivsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 && make

sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	sdsl-lite/install.sh sdsl

# for neomake/syntastic in vim
show_cxxflags:
	@echo $(CXXFLAGS)

.PHONY: all build clean tests valgrind format libdivsufsort sdsl show_cxxflags
