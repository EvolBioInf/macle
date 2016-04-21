CC ?= g++
CFLAGS := -std=c++11 -Isrc -Wall -Wextra -Wshadow -O2 -g # -pg
LDFLAGS := -lm -lz -lgsl -lgslcblas -lblas -ldivsufsort # -pg
TARGET := dnalc
VERSION=0.1

SOURCES := $(wildcard src/*.c)
OBJECTS := $(SOURCES:.c=.o)

TEST_SRC=$(wildcard tests/*.c)
TEST_OBJ=$(TEST_SRC:.c=.o)
TESTS=$(patsubst tests/%.c,build/%,$(TEST_SRC))

all: build/$(TARGET) tests

build/$(TARGET): $(OBJECTS)
	@echo " mkdir -p build"; mkdir -p build
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o build/$(TARGET) $(LDFLAGS)

src/%.o: src/%.c
	@echo " $(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<

tests/%.o: tests/%.c
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
	clang-format -i src/*.[ch] tests/*.[ch]

.PHONY: all clean lint tests valgrind
