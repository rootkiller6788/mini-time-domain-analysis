# Makefile for mini-stability-routh-hurwitz
# Routh-Hurwitz Stability Criterion — Complete Engineering Implementation
#
# Usage:
#   make          - Build all object files
#   make test     - Build and run test suite
#   make examples - Build all example programs
#   make all-exe  - Build and run all examples
#   make clean    - Remove build artifacts
#   make count    - Count lines in include/ and src/

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

# Source files
SRC = src/routh_array.c \
      src/hurwitz_determinant.c \
      src/special_cases.c \
      src/root_distribution.c \
      src/stability_check.c \
      src/relative_stability.c \
      src/jury_stability.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_routh_hurwitz.c

EX_SRC = examples/ex1_dc_motor_stability.c \
         examples/ex2_aircraft_pitch_control.c \
         examples/ex3_spacecraft_attitude.c

EX_BIN = $(EX_SRC:.c=.exe)

.PHONY: all clean test examples all-exe count

all: $(OBJ)
	@echo "mini-stability-routh-hurwitz: all objects built successfully"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex1_dc_motor_stability.c $(OBJ) -o examples/ex1_dc_motor_stability.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex2_aircraft_pitch_control.c $(OBJ) -o examples/ex2_aircraft_pitch_control.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex3_spacecraft_attitude.c $(OBJ) -o examples/ex3_spacecraft_attitude.exe $(LDFLAGS)
	@echo "All examples built. Run individually or use 'make all-exe'"

all-exe: examples
	@echo "=== Running all examples ==="
	@echo ""
	./examples/ex1_dc_motor_stability.exe
	@echo ""
	./examples/ex2_aircraft_pitch_control.exe
	@echo ""
	./examples/ex3_spacecraft_attitude.exe

count:
	@echo "include/ lines: $$(find include -name '*.h' -exec cat {} + | wc -l)"
	@echo "src/ lines:     $$(find src -name '*.c' -exec cat {} + | wc -l)"
	@echo "Total (include+src): $$(find . \( -path '*/include/*.h' -o -path '*/src/*.c' \) -exec cat {} + | wc -l)"

clean:
	rm -f src/*.o tests/test_runner tests/test_runner.exe examples/*.exe
