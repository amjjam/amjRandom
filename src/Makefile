ifndef prefix
	prefix=/opt/amj
endif

FLAGS=-Wall -g -fPIC
CPPFLAGS=$(FLAGS) 
CPP=g++

AMJRANDOM_SRC=amjRandom.C
AMJRANDOM_OBJ=$(AMJRANDOM_SRC:.C=.o)

all: libamjRandom.so

libamjRandom.so: $(AMJRANDOM_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: libamjRandom.so ../include/amjRandom.H
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	cp libamjRandom.so $(prefix)/lib
	cp ../include/amjRandom.H $(prefix)/include

uninstall:
	- rm $(prefix)/lib/libamjRandom.so
	- rm $(prefix)/include/amjRandom.H
	- rmdir $(prefix)/lib
	- rmdir $(prefix)/include

clean:
	- rm libamjRandom.so $(AMJRANDOM_OBJ)
