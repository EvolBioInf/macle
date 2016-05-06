CXX := g++
CFLAGS := -std=c++11 -Isrc -Ilibdivsufsort/include -Wall -Wextra -Wshadow -O2 -g -ggdb # -pg
LDFLAGS := -lm -lgsl -lgslcblas -lblas -ldivsufsort64 -Llibdivsufsort/lib # -pg
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
	git clone git@github.com:y-256/libdivsufsort.git
	cd libdivsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 && make

.PHONY: all build clean lint tests valgrind format libdivsufsort
