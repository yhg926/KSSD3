CC=gcc
CFLAGS= -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -mavx2 -mbmi2 -fopenmp -march=native 
SOURCE="./src"
BIN="./bin"
PRONAME="kssd3"

all:
	$(CC) $(CFLAGS) $(SOURCE)/*.c -o $(BIN)/$(PRONAME) -Iklib -lz -lm 

