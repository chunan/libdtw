CC = g++
CFLAGS = -Wall -Werror -ansi -pedantic -I /home/ugoc/include/

MACHINE = $(shell uname -m)
LDFLAGS = -ldtw \
					-lfeature \
					-lsegtree \
					-lutility

LDPATH = -L /home/ugoc/lib

.PHONY: all clean


all: CFLAGS +=  -O2
debug: CFLAGS += -g

all: example example2
debug: example example2

example: example.cpp
	${CC} ${CFLAGS} $(LDPATH) $< ${LDFLAGS} -o $@

example2: example2.cpp
	${CC} ${CFLAGS} $(LDPATH) $< ${LDFLAGS} -o $@

clean:
	rm -f example example2
