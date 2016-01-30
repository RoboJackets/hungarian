
CC=g++
CC_FLAGS=
PROG=hungarian
SRC=hungarian.cpp test.cpp

all: $(PROG)

$(PROG) : $(SRC)
		$(CC) $(CC_FLAGS) -o $(PROG) $(SRC)

clean : $(PROG)
		rm $(PROG)
