INSTALLDIR=~/lib/

RANDOM_SRC=Random.C
RANDOM_OBJ=$(RANDOM_SRC:.C=.o)

install: libRandom.a
	cp libRandom.a $(INSTALLDIR)

uninstall:
	rm $(INSTALLDIR)/libRandom.a

clean:
	rm $(RANDOM_OBJ)

libRandom.a: $(RANDOM_OBJ)
	ar rc $@ $^
