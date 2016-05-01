CC := g++
CFLAGS := -std=c++11 -Isrc -Ilibdivsufsort/include -Wall -Wextra -Wshadow -O2 -g # -Isdsl/include -msse4.2 -pg
LDFLAGS := -lm -lgsl -lgslcblas -lblas -ldivsufsort64 -Llibdivsufsort/lib # -Lsdsl/lib -lsdsl -pg
TARGET := dnalc

SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

TEST_SRC=$(wildcard tests/*.cpp)
TEST_OBJ=$(TEST_SRC:.cpp=.o)
TESTS=$(patsubst tests/%.cpp,build/%,$(TEST_SRC))

all: build/$(TARGET) tests

build/$(TARGET): $(OBJECTS)
	@echo " mkdir -p build"; mkdir -p build
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o build/$(TARGET) $(LDFLAGS)

src/%.o: src/%.cpp
	@echo " $(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<

tests/%.o: tests/%.cpp
	@echo " $(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<

$(TESTS): $(OBJECTS) $(TEST_OBJ)
	@echo " $(CC) $(filter-out src/$(TARGET).o,$(OBJECTS)) $@.o -o $@ $(LDFLAGS)"
	$(CC) $(filter-out src/$(TARGET).o,$(OBJECTS)) $(@:build/%=tests/%).o -o $@ $(LDFLAGS)

clean:
	@echo " $(RM) -r build"; $(RM) -r build
	@echo " $(RM) $(OBJECTS) $(TEST_OBJ) $(TESTS)"; $(RM) $(OBJECTS) $(TEST_OBJ) $(TESTS)

tests: $(TESTS)
	bash ./tests/runtests.sh

valgrind:
	VALGRIND="valgrind --leak-check=full" $(MAKE)

format:
	clang-format -i src/*.cpp src/*.h tests/*.cpp tests/*.h

sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	sdsl-lite/install.sh sdsl

libdivsufsort:
	git clone git@github.com:y-256/libdivsufsort.git
	cd libdivsufsort && cmake -DBUILD_SHARED_LIBS=0 -DBUILD_EXAMPLES=0 -DBUILD_DIVSUFSORT64=1 && make

.PHONY: all clean lint tests valgrind format sdsl libdivsufsort
