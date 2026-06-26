# Makefile for mini-sensitivity-analysis
# Domain: Time-Domain Analysis → Sensitivity Analysis (Control Theory)

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -I./include
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

SRCDIR   = src
INCDIR   = include
OBJDIR   = build/obj
LIBDIR   = build/lib
BINDIR   = build/bin

# Source files
SRCS_C   = $(wildcard $(SRCDIR)/*.c)
OBJS_C   = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS_C))

# Test files
TEST_SRC = tests/test_sensitivity.c
TEST_OBJ = $(OBJDIR)/test_sensitivity.o
TEST_BIN = $(BINDIR)/test_sensitivity

# Example files
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c, $(BINDIR)/%, $(EXAMPLE_SRCS))

# Library
LIBRARY  = $(LIBDIR)/libsensitivity.a

.PHONY: all clean test examples check

all: $(LIBRARY) $(TEST_BIN) $(EXAMPLE_BINS)

# Build directory structure
$(OBJDIR) $(LIBDIR) $(BINDIR):
	mkdir -p $@

# Compile C sources
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Build static library (excluding test object)
LIB_OBJS = $(filter-out $(OBJDIR)/test_sensitivity.o, $(OBJS_C))
$(LIBRARY): $(LIB_OBJS) | $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

# Build test binary
$(TEST_OBJ): $(TEST_SRC) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(LIB_OBJS) $(TEST_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Build example binaries
$(BINDIR)/%: examples/%.c $(LIBRARY) | $(BINDIR)
	$(CC) $(CFLAGS) $< -L$(LIBDIR) -lsensitivity -o $@ $(LDFLAGS)

# Run tests
test: $(TEST_BIN)
	@echo "Running sensitivity analysis tests..."
	$(TEST_BIN)
	@echo "All tests passed."

# Build examples
examples: $(EXAMPLE_BINS)

# Clean
clean:
	rm -rf build/

# Quick compile check (no linking)
check:
	@echo "Checking compilation..."
	@for src in $(SRCS_C); do \
		$(CC) $(CFLAGS) -fsyntax-only $$src || exit 1; \
	done
	@for hdr in $(INCDIR)/*.h; do \
		echo "#include \"$$(basename $$hdr)\"" | $(CC) $(CFLAGS) -fsyntax-only -xc - || exit 1; \
	done
	@echo "All files compile cleanly."

# Install: copy library and headers (optional)
PREFIX ?= /usr/local
install: $(LIBRARY)
	cp $(LIBRARY) $(PREFIX)/lib/
	cp $(INCDIR)/*.h $(PREFIX)/include/

# Help
help:
	@echo "Targets:"
	@echo "  all      - Build library, test, and examples"
	@echo "  test     - Build and run tests"
	@echo "  examples - Build example programs"
	@echo "  clean    - Remove build artifacts"
	@echo "  check    - Quick syntax check on all sources"
	@echo "  install  - Install to PREFIX=$(PREFIX)"
