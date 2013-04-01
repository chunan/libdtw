CC = g++
CFLAGS = -Wall -Werror -ansi -pedantic \
				 -I ../include/ -I ../../libutility/include \
				 -I ../../libsegtree/include -I ../../libfeature/include

MACHINE = $(shell uname -m)
LDFLAGS = -ldtw \
					-lfeature \
					-lsegtree \
					-lutility

LDPATH = -L ../lib/$(MACHINE) -L ../../libutility/lib/$(MACHINE) \
				 -L ../../libsegtree/lib/$(MACHINE) -L ../../libfeature/lib/$(MACHINE)

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
