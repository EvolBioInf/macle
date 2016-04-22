CC ?= g++
CFLAGS := -std=c++11 -Isrc -Wall -Wextra -Wshadow -O2 -g # -pg
LDFLAGS := -lm -lz -lgsl -lgslcblas -lblas -ldivsufsort # -pg
TARGET := dnalc
VERSION=0.1

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

.PHONY: all clean lint tests valgrind
