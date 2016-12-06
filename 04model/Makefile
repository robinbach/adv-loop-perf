# Makefile

CC = g++
LD = g++

CFLAGS=-O3 -std=c++11 -g

OBJS = LR.o dataImport.o

LR = LR.bin

default: $(LR)

debug: CFLAGS += -g
debug: $(LR)

$(LR): $(OBJS)
	$(LD) $(OBJS) -o $(LR)

LR.o: LR.cpp
	$(CC) $(CFLAGS) -I ./eigen/ LR.cpp

dataImport.o: dataImport.cpp dataImport.h
	$(CC) $(CFLAGS) dataImport.cpp

clean:
	rm -rf *.o $(LR)

binary:
	g++ LR.o dataImport.o -o LR.bin

run:
	./LR.bin covtype_blank.data
