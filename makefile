CC=g++
CC_FLAGS=-std=c++11
PROG=hungarian
SRC=hungarian.cpp test.cpp

all: $(PROG)

$(PROG) : $(SRC)
		$(CC) $(CC_FLAGS) -o $(PROG) $(SRC)

clean : $(PROG)
		rm $(PROG)
