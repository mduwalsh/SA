CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = pun allplot fplot sajob hist xsum

.PHONY: all clean

all: $(EXECUTABLES)

pun: pun.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
sajob: sajob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

fplot: fplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

allplot: allplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

hist: hist.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

xsum: xsum.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
