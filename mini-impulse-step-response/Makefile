# Makefile for mini-impulse-step-response
# Impulse and Step Response -- Time-Domain Analysis
# L1-L9: Complete engineering implementation

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = src/time_response_common.c \
      src/impulse_response.c \
      src/step_response.c \
      src/convolution_response.c \
      src/response_metrics.c \
      src/discrete_response.c \
      src/partial_fraction.c \
      src/response_identification.c

OBJ = $(SRC:.c=.o)

TEST_SRC = tests/test_impulse_step.c

EX_SRC = examples/ex_first_order.c \
         examples/ex_second_order.c \
         examples/ex_motor_response.c

EX_BIN = $(EX_SRC:.c=)

.PHONY: all clean test examples

all: $(OBJ)
	@echo "mini-impulse-step-response: all objects built"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o tests/test_runner $(LDFLAGS)
	./tests/test_runner

examples: $(OBJ)
	$(CC) $(CFLAGS) examples/ex_first_order.c $(OBJ) -o examples/ex_first_order $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_second_order.c $(OBJ) -o examples/ex_second_order $(LDFLAGS)
	$(CC) $(CFLAGS) examples/ex_motor_response.c $(OBJ) -o examples/ex_motor_response $(LDFLAGS)

clean:
	rm -f src/*.o tests/test_runner examples/ex_first_order examples/ex_second_order examples/ex_motor_response
