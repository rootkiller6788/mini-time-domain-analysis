# Makefile for mini-higher-order-reduction
# Higher-Order System Model Reduction for Control Theory
# L1-L9: Complete engineering implementation

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/reduction_core.c \
      src/reduction_dominant_pole.c \
      src/reduction_modal.c \
      src/reduction_balanced.c \
      src/reduction_routh.c \
      src/reduction_moment.c \
      src/reduction_demo.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_reduction.c
EX_SRC = examples/example_flexible_shaft.c \
         examples/example_aircraft_dynamics.c \
         examples/example_thermal_system.c

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-higher-order-reduction: all objects built"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

TEST_OBJ = src/reduction_core.o src/reduction_dominant_pole.o src/reduction_modal.o src/reduction_balanced.o src/reduction_routh.o src/reduction_moment.o

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(TEST_OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/example_flexible_shaft.c $(OBJ) -o examples/example_flexible_shaft $(LDFLAGS)
	$(CC) $(CFLAGS) examples/example_aircraft_dynamics.c $(OBJ) -o examples/example_aircraft_dynamics $(LDFLAGS)
	$(CC) $(CFLAGS) examples/example_thermal_system.c $(OBJ) -o examples/example_thermal_system $(LDFLAGS)

clean:
	rm -f src/*.o tests/test_runner examples/example_flexible_shaft examples/example_aircraft_dynamics examples/example_thermal_system
