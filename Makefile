CC = g++
CXX = g++

CXXFLAGS = -Wall -Werror -ansi -pedantic -I include -I /home/ugoc/include

MACHINE = $(shell uname -m)

SRC = dtw_parm.cpp dtw_util.cpp
OBJ = $(addprefix obj/$(MACHINE)/,$(SRC:.cpp=.o))


vpath %.h   include
vpath %.cpp src
vpath %.o obj/$(MACHINE)
vpath %.a lib/$(MACHINE)

.PHONY: mk_machine_dir all clean allclean

all: CXXFLAGS += -O2

debug: CXXFLAGS += -g

all: mk_machine_dir lib/$(MACHINE)/libdtw.a

debug: mk_machine_dir lib/$(MACHINE)/libdtw.a


%.d: %.cpp
	$(CC) -M $(CXXFLAGS) $< > $@

lib/$(MACHINE)/libdtw.a: obj/$(MACHINE)/dtw_util.o obj/$(MACHINE)/dtw_parm.o
	$(AR) rucs $@ $^

obj/$(MACHINE)/%.o: src/%.cpp
	$(CC) -c $(CXXFLAGS) -o $@ $^

mk_machine_dir:
	@mkdir -p obj/$(MACHINE)
	@mkdir -p lib/$(MACHINE)

allclean: clean
	$(RM) lib/$(MACHINE)/libdtw.a

clean:
	$(RM) $(OBJ)
