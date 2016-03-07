CC := clang
CFLAGS := -std=gnu99 -Wall -Wextra -Wshadow -pedantic -O0 -g
LDFLAGS := -lm -lz -lgsl -lgslcblas -lblas -ldivsufsort
TARGET := dnalc
VERSION=0.1

SOURCES := $(wildcard src/*.c)
OBJECTS := $(SOURCES:.c=.o)

all: build/$(TARGET)

build/$(TARGET): $(OBJECTS)
	@echo " mkdir -p build"; mkdir -p build
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o build/$(TARGET) $(LDFLAGS)

src/%.o: src/%.c
	@echo " $(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo " $(RM) -r build"; $(RM) -r build
	@echo " $(RM) src/*.o"; $(RM) src/*.o

format:
	clang-format -i src/*.[ch]

.PHONY: all clean lint
