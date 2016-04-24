CC := g++
CFLAGS := -std=c++11 -Isrc -Isdsl/include -Wall -Wextra -O2 -msse4.2 -g # -pg -Wshadow 
LDFLAGS := -lm -lz -lgsl -lgslcblas -lblas -Lsdsl/lib -ldivsufsort # -ldivsufsort64 -lsdsl -pg
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
	sh ./tests/runtests.sh

valgrind:
	VALGRIND="valgrind --leak-check=full" $(MAKE)

format:
	clang-format -i src/*.cpp src/*.h tests/*.cpp tests/*.h

sdsl:
	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite
	sdsl-lite/install.sh sdsl

.PHONY: all clean lint tests valgrind
