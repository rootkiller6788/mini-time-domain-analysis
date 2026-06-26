# Mini Steady-State Error Analysis - Makefile
# Reference: Nise Ch.7, Ogata Ch.5, Franklin Ch.4
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
LDFLAGS = -lm
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
TEST_SRC = $(wildcard tests/*.c)
TEST_BIN = $(TEST_SRC:.c=.exe)
EXAMPLE_SRC = $(wildcard examples/*.c)
EXAMPLE_BIN = $(EXAMPLE_SRC:.c=.exe)

.PHONY: all clean test examples

all: libsse.a

libsse.a: $(OBJ)
	ar rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: libsse.a $(TEST_BIN)
	@for t in $(TEST_BIN); do echo "Running $$t..."; ./$$t || exit 1; done
	@echo "All tests passed."

examples: libsse.a $(EXAMPLE_BIN)
	@for e in $(EXAMPLE_BIN); do echo "Running $$e..."; ./$$e; done

tests/%.exe: tests/%.c libsse.a
	$(CC) $(CFLAGS) $< -L. -lsse $(LDFLAGS) -o $@

examples/%.exe: examples/%.c libsse.a
	$(CC) $(CFLAGS) $< -L. -lsse $(LDFLAGS) -o $@

clean:
	rm -f src/*.o tests/*.o tests/*.exe examples/*.o examples/*.exe libsse.a

distclean: clean
	rm -f *~ include/*~ src/*~ tests/*~ examples/*~
