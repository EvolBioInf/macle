CXX := g++
CFLAGS := -DUSE_SDSL -std=c++11 -Isrc -Ilibdivsufsort/include -Igsl/include -Isdsl/include -Wall -Wextra -O2 -g -ggdb -Wshadow # -pg -Wshadow
LDFLAGS := -lm -lgsl -lgslcblas -lblas -ldivsufsort -ldivsufsort64 -lsdsl -Llibdivsufsort/lib -Lgsl/lib -Lsdsl/lib # -pg
TARGET := dnalc

SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

TEST_SRC=$(wildcard tests/*.cpp)
TEST_OBJ=$(TEST_SRC:.cpp=.o)
TESTS=$(patsubst tests/%.cpp,build/%,$(TEST_SRC))

all: build tests

build: build/$(TARGET) $(TESTS)

build/$(TARGET): $(OBJECTS)
	@echo " mkdir -p build"; mkdir -p build
	@echo " $(CXX) $^ -o $(TARGET) $(LIB)"; $(CXX) $^ -o build/$(TARGET) $(LDFLAGS)

src/%.o: src/%.cpp
	@echo " $(CXX) $(CFLAGS) -c -o $@ $<"; $(CXX) $(CFLAGS) -c -o $@ $<

tests/%.o: tests/%.cpp
	@echo " $(CXX) $(CFLAGS) -c -o $@ $<"; $(CXX) $(CFLAGS) -c -o $@ $<

$(TESTS): $(OBJECTS) $(TEST_OBJ)
	@echo " $(CXX) $(filter-out src/$(TARGET).o,$(OBJECTS)) $@.o -o $@ $(LDFLAGS)"
	$(CXX) $(filter-out src/$(TARGET).o,$(OBJECTS)) $(@:build/%=tests/%).o -o $@ $(LDFLAGS)

clean:
	@echo " $(RM) -r build"; $(RM) -r build
	@echo " $(RM) $(OBJECTS) $(TEST_OBJ) $(TESTS)"; $(RM) $(OBJECTS) $(TEST_OBJ) $(TESTS)

tests: $(TESTS)
	bash ./tests/runtests.sh

valgrind:
	VALGRIND="valgrind --leak-check=full" $(MAKE)

format:
	clang-format -i src/*.cpp src/*.h tests/*.cpp tests/*.h

libdivsufsort:
	git clone https://github.com/y-256/libdivsufsort.git
	cd libdivsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 && make

gsl:
	wget http://ftp.halifax.rwth-aachen.de/gnu/gsl/gsl-2.1.tar.gz
	tar xf gsl-2.1.tar.gz && rm gsl-2.1.tar.gz
	cd gsl-2.1 && ./configure --prefix=$(shell pwd)/gsl && make && make install
	find gsl -name "*.so.*" -exec rm {} \;

sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	sdsl-lite/install.sh sdsl

.PHONY: all build clean tests valgrind format libdivsufsort gsl sdsl
