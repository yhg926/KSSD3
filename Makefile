CC = gcc

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(ROOT)/src
BIN_DIR := $(ROOT)/bin
OBJ_DIR := $(ROOT)/obj
PRONAME := kssd3
TARGET := $(BIN_DIR)/$(PRONAME)

# XGBoost paths
XGB_INCLUDE := $(ROOT)/xgboost/include
XGB_LIB := $(ROOT)/xgboost/lib

PREFIX := /usr/local

# Base flags shared by all compilations
BASE_CFLAGS := -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -fopenmp
LDLIBS := -lz -lm

SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))

# Source specific flags
$(OBJ_DIR)/command_sketch.o: CFLAGS += -march=native -mavx2 -mbmi2 -Iklib
$(OBJ_DIR)/command_ani.o: CFLAGS += -Iklib -I$(XGB_INCLUDE)
$(OBJ_DIR)/command_composite.o $(OBJ_DIR)/command_matrix.o $(OBJ_DIR)/command_operate.o: CFLAGS += -Iklib

all: $(TARGET)
	@echo "Build completed."
	@echo "👉 If you haven’t already, run: make install_env to set LD_LIBRARY_PATH for libxgboost"

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(BASE_CFLAGS) $^ -o $@ -L$(XGB_LIB) -lxgboost $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(BASE_CFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET)

install_env:
	echo 'export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT)/xgboost/lib' >> ~/.bashrc
	@echo "Added to .bashrc. Run 'source ~/.bashrc' to apply."

install: all
	sudo install -m 755 $(TARGET) $(PREFIX)/bin/$(PRONAME)
	@echo "Installed $(PRONAME) to $(PREFIX)/bin"

uninstall:
	sudo rm -f $(PREFIX)/bin/$(PRONAME)
	@echo "Removed $(PRONAME) from $(PREFIX)/bin"

