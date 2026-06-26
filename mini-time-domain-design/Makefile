# Makefile for mini-time-domain-design
# Time-Domain Design - state-space control system design library
#
# Targets:
#   make          - build static library
#   make test     - build and run tests
#   make examples - build all example programs
#   make clean    - remove build artifacts
#   make all      - build everything and run tests

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

SRCDIR   = src
BUILDDIR = build
LIBNAME  = libtdd.a

SRCS     = $(wildcard $(SRCDIR)/*.c)
OBJS     = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))

TEST_SRC = tests/test_core.c
TEST_BIN = $(BUILDDIR)/test_core

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,$(BUILDDIR)/%,$(EXAMPLE_SRCS))

.PHONY: all clean test examples

all: $(BUILDDIR)/$(LIBNAME) test

# Build build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile source files to objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Build static library
$(BUILDDIR)/$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Build and run tests
test: $(BUILDDIR)/$(LIBNAME)
	$(CC) $(CFLAGS) $(TEST_SRC) -L$(BUILDDIR) -ltdd $(LDFLAGS) -o $(TEST_BIN)
	./$(TEST_BIN)

# Build all examples
examples: $(BUILDDIR)/$(LIBNAME) $(EXAMPLE_BINS)
	@echo "=== Examples built ==="
	@for e in $(EXAMPLE_BINS); do \
		echo "Running $$e..."; \
		$$e; \
		echo ""; \
	done

# Build each example
$(BUILDDIR)/%: examples/%.c $(BUILDDIR)/$(LIBNAME)
	$(CC) $(CFLAGS) $< -L$(BUILDDIR) -ltdd $(LDFLAGS) -o $@

# Clean build artifacts
clean:
	rm -rf $(BUILDDIR)

# Help
help:
	@echo "mini-time-domain-design Makefile"
	@echo "  make          - build library and run tests"
	@echo "  make test     - run unit tests"
	@echo "  make examples - build and run all examples"
	@echo "  make clean    - remove build artifacts"
