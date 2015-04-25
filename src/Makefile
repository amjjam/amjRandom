libdir=$(prefix)/lib

CPPFLAGS=-Wall -g -fPIC
CPP=g++

RANDOM_SRC=Random.C
RANDOM_OBJ=$(RANDOM_SRC:.C=.o)

build: libRandom.a libRandom.so

libRandom.a: $(RANDOM_OBJ)
	ar rc $@ $^

libRandom.so: $(RANDOM_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: $(libdir)/libRandom.a $(libdir)/libRandom.so

$(libdir)/libRandom.a $(libdir)/libRandom.so: libRandom.a libRandom.so
	cp libRandom.a libRandom.so $(libdir)

uninstall:
	- rm -f $(libdir)/libRandom.a $(libdir)/libRandom.so

clean:
	- rm -f libRandom.a libRandom.so $(RANDOM_OBJ)

