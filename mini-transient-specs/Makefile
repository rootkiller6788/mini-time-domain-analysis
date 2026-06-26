#===============================================================================
# Makefile — mini-transient-specs
#
# Builds the transient response specifications library, test suite, and examples.
#
# Targets:
#   make          — build library and tests
#   make test     — build and run tests
#   make examples — build all examples
#   make all      — build everything
#   make clean    — remove build artifacts
#===============================================================================

CC       := gcc
CFLAGS   := -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS  := -lm

# Source files
CORE_SRC  := src/transient_specs_core.c
FO_SRC    := src/transient_first_order.c
SO_SRC    := src/transient_second_order.c
HO_SRC    := src/transient_higher_order.c
MET_SRC   := src/transient_metrics.c
ALL_SRC   := $(CORE_SRC) $(FO_SRC) $(SO_SRC) $(HO_SRC) $(MET_SRC)

# Object files
CORE_OBJ  := build/transient_specs_core.o
FO_OBJ    := build/transient_first_order.o
SO_OBJ    := build/transient_second_order.o
HO_OBJ    := build/transient_higher_order.o
MET_OBJ   := build/transient_metrics.o
ALL_OBJ   := $(CORE_OBJ) $(FO_OBJ) $(SO_OBJ) $(HO_OBJ) $(MET_OBJ)

# Test
TEST_SRC  := tests/test_transient_specs.c
TEST_BIN  := build/test_transient_specs

# Examples
EX1_SRC   := examples/example_design_specs.c
EX2_SRC   := examples/example_tradeoff.c
EX3_SRC   := examples/example_higher_order.c
EX1_BIN   := build/example_design_specs
EX2_BIN   := build/example_tradeoff
EX3_BIN   := build/example_higher_order

.PHONY: all clean test examples

all: builddir test examples

builddir:
	@mkdir -p build

#===============================================================================
# Library objects
#===============================================================================

$(CORE_OBJ): $(CORE_SRC) include/transient_specs.h
	$(CC) $(CFLAGS) -c $< -o $@

$(FO_OBJ): $(FO_SRC) include/transient_specs.h
	$(CC) $(CFLAGS) -c $< -o $@

$(SO_OBJ): $(SO_SRC) include/transient_specs.h
	$(CC) $(CFLAGS) -c $< -o $@

$(HO_OBJ): $(HO_SRC) include/transient_specs.h
	$(CC) $(CFLAGS) -c $< -o $@

$(MET_OBJ): $(MET_SRC) include/transient_specs.h
	$(CC) $(CFLAGS) -c $< -o $@

#===============================================================================
# Test
#===============================================================================

$(TEST_BIN): $(ALL_OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(ALL_OBJ) -o $@ $(LDFLAGS)

test: builddir $(TEST_BIN)
	@echo "=== Running transient_specs tests ==="
	./$(TEST_BIN)
	@echo "=== Tests complete ==="

#===============================================================================
# Examples
#===============================================================================

$(EX1_BIN): $(ALL_OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(ALL_OBJ) -o $@ $(LDFLAGS)

$(EX2_BIN): $(ALL_OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(ALL_OBJ) -o $@ $(LDFLAGS)

$(EX3_BIN): $(ALL_OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(ALL_OBJ) -o $@ $(LDFLAGS)

examples: builddir $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	@echo "=== Running example: Design from Specs ==="
	./$(EX1_BIN)
	@echo ""
	@echo "=== Running example: Design Trade-Off ==="
	./$(EX2_BIN)
	@echo ""
	@echo "=== Running example: Dominant Pole Reduction ==="
	./$(EX3_BIN)
	@echo "=== Examples complete ==="

#===============================================================================
# Utility
#===============================================================================

clean:
	rm -rf build

lines:
	@echo "=== Line counts ==="
	@echo "Headers:"
	@wc -l include/*.h
	@echo ""
	@echo "Sources:"
	@wc -l src/*.c
	@echo ""
	@echo "Lean:"
	@wc -l src/*.lean
	@echo ""
	@echo -n "Total include+src: "
	@cat include/*.h src/*.c src/*.lean 2>/dev/null | wc -l
