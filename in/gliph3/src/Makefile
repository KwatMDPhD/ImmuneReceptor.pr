CC=gcc -std=gnu99
JANSSONL= /usr/local/lib/libjansson.a
CFLAGS= -g -Wall -O2 -I/usr/local/include
LDFLAGS= -lm -lz -lpthread

ICTOOLS = ic_tools.o ic_clust.o ic_all.o ic_io.o ic_util.o \
					ic_type.o ic_config.o klib/kstring.o heap.o queue.o

all: ictools

ictools: $(ICTOOLS)
	$(CC) $(DEBUG) $^ -o $@ $(LDFLAGS) $(JANSSONL)
	mv ictools ../bin/

%.o: %.c
	$(CC) $(DEBUG) -c -o $@ $< $(CFLAGS)

klib/kstring.o:klib/kstring.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o klib/kstring.o

debug: DEBUG = -DDEBUG

debug: ictools
