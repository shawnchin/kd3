# Makefile to build test code.

SOURCES   = run_test.c kd3/kdtree.c 
HEADERS   = kd3/kdtree.h

GCC_CFLAGS_LVL1 = -Wall -pedantic 
GCC_CFLAGS_LVL2 = -Wextra -Wstrict-prototypes -Wmissing-prototypes -Wpointer-arith 
GCC_CFLAGS_LVL3 = -Wreturn-type -Wswitch -Wshadow -Wcast-align -Wunused 
GCC_CFLAGS_LVL4 = -Wwrite-strings -Wcast-qual 

CFLAGS    = -g -std=c99 
EXECUTABLE = run_test

CFLAGS += $(GCC_CFLAGS_LVL1)
CFLAGS += $(GCC_CFLAGS_LVL2)
CFLAGS += $(GCC_CFLAGS_LVL3)
CFLAGS += $(GCC_CFLAGS_LVL4)

DEPS      = $(HEADERS) Makefile 
OBJECTS   = $(SOURCES:.c=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

$(OBJECTS): $(DEPS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS) *.gcno *.gcda

