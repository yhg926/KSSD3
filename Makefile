CFLAGS= -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -mavx2 -mbmi2 -fopenmp -march=native
SOURCE=./src
BIN=./bin
PRONAME=kssd3

# Add XGBoost include and lib paths
XGB_INCLUDE=/home/ubt/work1/tools/xgboost/include
XGB_LIB=/home/ubt/work1/tools/xgboost/lib

all:
	$(CC) $(CFLAGS) $(SOURCE)/*.c -o $(BIN)/$(PRONAME) \
	-Iklib -I$(XGB_INCLUDE) -L$(XGB_LIB) -lxgboost -lz -lm
clean:
	rm -f $(BIN)/$(PRONAME)


#old
#CC=gcc
#CFLAGS= -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -ggdb -mavx2 -mbmi2 -fopenmp -march=native 
#SOURCE="./src"
#BIN="./bin"
#PRONAME="kssd3"

#all:
#	$(CC) $(CFLAGS) $(SOURCE)/*.c -o $(BIN)/$(PRONAME) -Iklib -lz -lm 

