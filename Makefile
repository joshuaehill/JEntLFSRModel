#for static analysis, use 
#scan-build-21 --use-cc=/usr/bin/clang-21 make
#and setup the compiler as clang

#Clang 21
#CC=clang-21
#MUTE=-Wno-unsafe-buffer-usage
#DEBUGCFLAGS=-O1
#CFLAGS=-Wall -Wextra -Weverything -march=native -g -fno-omit-frame-pointer -fno-optimize-sibling-calls -Werror $(MUTE) $(DEBUGCFLAGS)

#GCC (tested with v11)
CC=gcc
MUTE=
DEBUGCFLAGS=-O1
CFLAGS=-Wall -Wextra -march=native -g -fno-omit-frame-pointer -fno-optimize-sibling-calls -Werror $(MUTE) $(DEBUGCFLAGS)

LDFLAGS=$(DEBUGLDFLAGS)

src = $(wildcard *.c)
obj = $(src:.c=.o)
dep = $(obj:.o=.d)  # one dependency file for each sourc

BINARIES=lfsr-matrix-evolution test-lfsr-matrix

all:	$(BINARIES)

-include $(dep)   # include all dep files in the makefile

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

sources:	$(src)

clean:
	rm -f *.o *~ $(BINARIES) a.out *.d

lfsr-matrix-evolution: lfsr-matrix-evolution.o
	$(CC) -o $@ $^ $(LDFLAGS)

test-lfsr-matrix: test-lfsr-matrix.o
	$(CC) -o $@ $^ $(LDFLAGS)
