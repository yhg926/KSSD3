CC = gcc
CFLAGS = -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -mavx2 -mbmi2 -fopenmp -march=native

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRCDIR := $(ROOT)/src
OBJDIR := $(ROOT)/obj
BINDIR := $(ROOT)/bin
PRONAME := kssd3
TARGET := $(BINDIR)/$(PRONAME)

# XGBoost paths
XGB_INCLUDE := $(ROOT)/xgboost/include
XGB_LIB := $(ROOT)/xgboost/lib

PREFIX := /usr/local

all: $(TARGET); \
    echo "Build completed."; \
    echo "👉 If you haven’t already, run: make install_env to set LD_LIBRARY_PATH for libxgboost"

SRCS := $(wildcard $(SRCDIR)/*.c)
OBJS := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

$(TARGET): $(OBJS); \
    mkdir -p $(BINDIR); \
    $(CC) $(CFLAGS) $^ -o $@ -L$(XGB_LIB) -lxgboost -lz -lm

$(OBJDIR)/%.o: $(SRCDIR)/%.c; \
    mkdir -p $(OBJDIR); \
    $(CC) $(CFLAGS) -Iklib -I$(XGB_INCLUDE) -c $< -o $@

clean:; rm -f $(TARGET) $(OBJS)

install_env:; \
    echo 'export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT)/xgboost/lib' >> ~/.bashrc; \
    echo "Added to .bashrc. Run 'source ~/.bashrc' to apply."

install: all; \
    sudo install -m 755 $(TARGET) $(PREFIX)/bin/$(PRONAME); \
    echo "Installed $(PRONAME) to $(PREFIX)/bin"

uninstall:; \
    sudo rm -f $(PREFIX)/bin/$(PRONAME); \
    echo "Removed $(PRONAME) from $(PREFIX)/bin"

