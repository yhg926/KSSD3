CC = gcc
CFLAGS = -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -ggdb -mavx2 -mbmi2 -fopenmp -march=native

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRCDIR := $(ROOT)/src
OBJDIR := $(ROOT)/obj
BINDIR := $(ROOT)/bin
PRONAME := kssd3
TARGET := $(BINDIR)/$(PRONAME)
PREFIX := /usr/local

all: $(TARGET)
	@echo "Build completed."

SRCS := $(wildcard $(SRCDIR)/*.c)

# sources requiring klib headers
KLIB_SRCS := $(SRCDIR)/command_ani.c \
             $(SRCDIR)/command_composite.c \
             $(SRCDIR)/command_matrix.c \
             $(SRCDIR)/command_operate.c \
             $(SRCDIR)/command_sketch.c \
             $(SRCDIR)/kssdlib_sort.c

KLIB_OBJS := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(KLIB_SRCS))

# generic source files (non-klib)
GEN_SRCS := $(filter-out $(KLIB_SRCS),$(SRCS))
GEN_OBJS := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(GEN_SRCS))

OBJS := $(GEN_OBJS) $(KLIB_OBJS)

$(TARGET): $(OBJS)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $^ -o $@ -lz -lm

$(KLIB_OBJS): CFLAGS += -Iklib

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

install:
	sudo install -m 755 $(TARGET) $(PREFIX)/bin/$(PRONAME)
	@echo "Installed $(PRONAME) to $(PREFIX)/bin"

uninstall:
	sudo rm -f $(PREFIX)/bin/$(PRONAME)
	@echo "Removed $(PRONAME) from $(PREFIX)/bin"

