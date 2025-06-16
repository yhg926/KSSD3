CC = gcc
CFLAGS = -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -mavx2 -mbmi2 -fopenmp -march=native

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SOURCE := $(ROOT)/src
BIN := $(ROOT)/bin
PRONAME := kssd3
TARGET := $(BIN)/$(PRONAME)

# XGBoost paths
XGB_INCLUDE := $(ROOT)/xgboost/include
XGB_LIB := $(ROOT)/xgboost/lib

PREFIX := /usr/local

all: $(TARGET)

$(TARGET):
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) $(SOURCE)/*.c -o $(TARGET) \
	-Iklib -I$(XGB_INCLUDE) -L$(XGB_LIB) -lxgboost -lz -lm

clean:
	rm -f $(TARGET)

install_env:
	echo 'export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT)/lib' >> ~/.bashrc
	@echo "Added to .bashrc. Run 'source ~/.bashrc' to apply."

install: all
	sudo install -m 755 $(TARGET) $(PREFIX)/bin/$(PRONAME)
	@echo "Installed $(PRONAME) to $(PREFIX)/bin"

uninstall:
	sudo rm -f $(PREFIX)/bin/$(PRONAME)
	@echo "Removed $(PRONAME) from $(PREFIX)/bin"

