# Makefile for mini-second-order-systems
# Time-Domain Analysis: Second-Order Systems
# L1-L9: Complete engineering implementation
#
# Usage:
#   make          - Build all object files
#   make test     - Build and run test suite
#   make examples - Build all example programs
#   make all-exe  - Build and run all examples
#   make clean    - Remove build artifacts

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

# Source files
SRC = src/second_order_core.c \
      src/transient_specs.c \
      src/response_computation.c \
      src/system_identification.c \
      src/canonical_models.c \
      src/second_order_design.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_second_order.c

EX_SRC = examples/ex1_mass_spring_damper.c \
         examples/ex2_rlc_circuit.c \
         examples/ex3_vehicle_suspension.c \
         examples/ex4_dc_motor_position.c \
         examples/ex5_system_identification.c

EX_BIN = $(EX_SRC:.c=.exe)

.PHONY: all clean test examples all-exe count

all: $(OBJ)
	@echo "mini-second-order-systems: all objects built successfully"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex1_mass_spring_damper.c $(OBJ) -o examples/ex1_mass_spring_damper.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex2_rlc_circuit.c $(OBJ) -o examples/ex2_rlc_circuit.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex3_vehicle_suspension.c $(OBJ) -o examples/ex3_vehicle_suspension.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex4_dc_motor_position.c $(OBJ) -o examples/ex4_dc_motor_position.exe $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex5_system_identification.c $(OBJ) -o examples/ex5_system_identification.exe $(LDFLAGS)
	@echo "All examples built. Run them individually or use 'make all-exe'"

all-exe: examples
	@echo "=== Running all examples ==="
	@echo ""
	./examples/ex1_mass_spring_damper.exe
	@echo ""
	./examples/ex2_rlc_circuit.exe
	@echo ""
	./examples/ex3_vehicle_suspension.exe
	@echo ""
	./examples/ex4_dc_motor_position.exe
	@echo ""
	./examples/ex5_system_identification.exe

count:
	@echo "include/ lines: $$(find include -name '*.h' -exec cat {} + | wc -l)"
	@echo "src/ lines:     $$(find src -name '*.c' -exec cat {} + | wc -l)"
	@echo "Total (include+src): $$(find . \\( -path '*/include/*.h' -o -path '*/src/*.c' \\) -exec cat {} + | wc -l)"

clean:
	rm -f src/*.o tests/test_runner examples/*.exe
